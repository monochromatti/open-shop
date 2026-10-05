"""Lift a physical incumbent into all algebraic variables and audit feasibility."""
function _lift_start!(b, c, warm; tolerance = 1e-7)
    assigned=Dict{VariableRef,Float64}()
    variables=all_variables(b.m)
    for v in variables
        set_start_value(v, nothing)
    end
    byname=Dict(name(v)=>v for v in variables)
    function put(n, x)
        haskey(byname, n) || return
        y=Float64(x)
        isfinite(y) || error("nonfinite lift $n")
        assigned[byname[n]]=y
        set_start_value(byname[n], y)
    end
    function table(n, curve, x, lo, hi)
        haskey(byname, n) || return
        lo,hi=get(get(b.m.ext,:global_tables,Dict()),n,(lo,hi))
        cells=OpenSHOP._global_intervals(curve.x, lo, hi)
        k=findfirst(pair->pair[1]<=x<=pair[2], cells)
        k===nothing && error("table lift outside domain: $n")
        for (j, (l, r)) in enumerate(cells)
            λ=l==r ? 0.0 : (x-l)/(r-l)
            put("$(n)_cell[$j]", j==k ? 1.0 : 0.0)
            put("$(n)_left[$j]", j==k ? 1-λ : 0.0)
            put("$(n)_right[$j]", j==k ? λ : 0.0)
        end
        put(n, OpenSHOP.table_value(curve, x; extrapolation = :linear))
    end
    function turbine(n, curve, q, h, hlo, hhi, qmax)
        qlo,qmax,hlo,hhi=get(get(b.m.ext,:global_turbines,Dict()),n,(0.0,qmax,hlo,hhi))
        qc=OpenSHOP._global_intervals(curve.discharge, qlo, qmax)
        hc=OpenSHOP._global_intervals(curve.heads, hlo, hhi)
        cells=[(qr, hr) for qr in qc for hr in hc]
        cells=get(get(b.m.ext,:global_turbine_cells,Dict()),n,cells)
        normalized=get(get(b.m.ext,:global_turbine_normalized,Dict()),n,false)
        k=findfirst(pair->pair[1][1]<=q<=pair[1][2] && pair[2][1]<=h<=pair[2][2], cells)
        k===nothing && error("turbine lift outside domain: $n")
        for (j, ((ql, qr), (hl, hr))) in enumerate(cells)
            i=OpenSHOP._curve_segment(curve.discharge, ql+(qr-ql)/2, :linear)
            z=OpenSHOP._curve_segment(curve.heads, hl+(hr-hl)/2, :linear)
            put("$(n)_cell[$j]", j==k ? 1.0 : 0.0)
            put(
                "$(n)_q[$j]",
                j==k ? (normalized ? (ql==qr ? 0.0 : (q-ql)/(qr-ql)) : (q-curve.discharge[i])/(curve.discharge[i + 1]-curve.discharge[i])) :
                0.0,
            )
            put(
                "$(n)_h[$j]",
                j==k ? (normalized ? (hl==hr ? 0.0 : (h-hl)/(hr-hl)) : (h-curve.heads[z])/(curve.heads[z + 1]-curve.heads[z])) : 0.0,
            )
        end
        put(n, OpenSHOP.turbine_efficiency(curve, q, h; extrapolation = :linear))
    end
    try
        OpenSHOP.validate(c, warm)["valid"] || error("seed fails original physical audit")
        s=c.system
        T=length(c.prices)
        ix=OpenSHOP.nodeindex(s)
        for v in variables
            is_fixed(v) && put(name(v), fix_value(v))
        end
        for (i, r) in enumerate(s.reservoirs), t in 1:(T + 1)
            put("v[$i,$t]", warm["V"][i, t]/r.vmax)
        end
        for i in axes(warm["H"], 1), t in 1:T
            put("h[$i,$t]", warm["H"][i, t]/250)
        end
        for (i, e) in enumerate(s.tunnels), t in 1:T
            q=warm["tunnel_q"][i, t]
            put("q[$i,$t]", q/50)
            put("tunnel_positive_$(i)_$(t)", max(q, 0.0))
            put("tunnel_negative_$(i)_$(t)", max(-q, 0.0))
            put("tunnel_direction_$(i)_$(t)", q>=0 ? 1.0 : 0.0)
        end
        heads=Dict{Tuple{Symbol,Int},String}()
        for (i, g) in enumerate(s.generators), t in 1:T
            q=warm["generator_q"][i, t]
            p=warm["power"][i, t]
            u=warm["u"][i, t]
            put("gq[$i,$t]", q/50)
            put("p[$i,$t]", p/40)
            put("u[$i,$t]", u)
            prior=t==1 ? g.initial_on : warm["u"][i, t - 1]
            put("su[$i,$t]", max(0, u-prior))
            put("sd[$i,$t]", max(0, prior-u))
            pl=OpenSHOP.plantof(s, g)
            total=sum(
                warm["generator_q"][j, t] for
                (j, z) in enumerate(s.generators) if z.plant==pl.name
            )
            n="net_head_$(i)_$(t)"
            haskey(byname, n) && (heads[(pl.name, t)]=n)
            actual=heads[(pl.name, t)]
            hv=byname[actual]
            receiver=warm["H"][ix[pl.target], t]
            h=warm["H"][ix[pl.source], t]-OpenSHOP.outlet_head(pl, receiver)-OpenSHOP.tailwater(
                pl,
                total,
            )
            put(n, h)
            outlet="outlet_head_$(i)_$(t)"
            if haskey(byname, outlet)
                j=ix[pl.target]
                lo, hi=if j<=length(s.reservoirs)
                    r=s.reservoirs[j]
                    b.domains===nothing ? OpenSHOP._global_level_range(r, r.vmin, r.vmax) :
                    (b.domains.hlo[j, t], b.domains.hhi[j, t])
                elseif j<=length(s.reservoirs)+length(s.junctions)
                    r=s.junctions[j - length(s.reservoirs)]
                    (r.hmin, r.hmax)
                else
                    (receiver, receiver)
                end
                floor=pl.outlet_head_floor
                table(
                    outlet,
                    OpenSHOP.TableCurve([lo, floor, hi], [floor, floor, hi]),
                    receiver,
                    lo,
                    hi,
                )
            end
            pl.tailwater_curve!==nothing && table(
                "tailwater_$(i)_$(t)",
                pl.tailwater_curve,
                total,
                0.0,
                sum(z.qmax for z in s.generators if z.plant==pl.name),
            )
            if g.turbine_table===nothing
                put("eta_$(i)_$(t)", OpenSHOP.efficiency(g, q, h))
            else
                turbine(
                    "turbine_$(i)_$(t)",
                    g.turbine_table,
                    q,
                    h,
                    lower_bound(hv),
                    upper_bound(hv),
                    g.qmax,
                )
                table(
                    "qlo_$(i)_$(t)",
                    OpenSHOP.TableCurve(g.turbine_table.heads, g.turbine_table.qmin),
                    h,
                    lower_bound(hv),
                    upper_bound(hv),
                )
                table(
                    "qhi_$(i)_$(t)",
                    OpenSHOP.TableCurve(g.turbine_table.heads, g.turbine_table.qmax),
                    h,
                    lower_bound(hv),
                    upper_bound(hv),
                )
            end
            g.generator_efficiency_curve!==nothing &&
                table("electrical_$(i)_$(t)", g.generator_efficiency_curve, p, 0.0, g.pmax)
        end
        for (i, r) in enumerate(s.reservoirs), t in 1:T
            r.level_curve===nothing && continue
            lo=b.domains===nothing ? r.vmin :
               (b.domains.lower[i, t]+b.domains.lower[i, t + 1])/2
            hi=b.domains===nothing ? r.vmax :
               (b.domains.upper[i, t]+b.domains.upper[i, t + 1])/2
            table(
                "level_$(i)_$(t)",
                r.level_curve,
                (warm["V"][i, t]+warm["V"][i, t + 1])/2,
                lo,
                hi,
            )
        end
        shortfall=OpenSHOP.release_requirements(c, warm["river_release"]).shortfall
        for (i, r) in enumerate(s.rivers), t in 1:T
            put("rq[$i,$t]", warm["river_release"][i, t]/100)
            put("a[$i,$t]", warm["gate"][i, t])
            put("shortfall_release[$i,$t]", shortfall[i, t])
            r.law==:junction && continue
            h=warm["H"][ix[r.source], t]
            r.discharge_curve!==nothing && table(
                "river_law_$(i)_$(t)",
                r.discharge_curve,
                h,
                first(r.discharge_curve.x),
                last(r.discharge_curve.x),
            )
            wet=max(h-r.crest, 0.0)
            put("wet_$(i)_$(t)", wet)
            put("wet_root_$(i)_$(t)", sqrt(wet))
            put("wet_branch_$(i)_$(t)", h>=r.crest ? 1.0 : 0.0)
        end
        missing=[name(v) for v in variables if !haskey(assigned, v)]
        isempty(missing) ||
            return Dict("valid"=>false, "missing"=>missing, "assigned"=>length(assigned))
        residual=0.0
        worst=""
        count=0
        violation(x, set) =
            set isa MOI.EqualTo ? abs(x-set.value) :
            set isa MOI.LessThan ? max(0.0, x-set.upper) :
            set isa MOI.GreaterThan ? max(0.0, set.lower-x) :
            set isa MOI.Interval ? max(0.0, set.lower-x, x-set.upper) :
            set isa MOI.ZeroOne ? max(abs(x-round(x)), max(0.0, -x, x-1)) :
            set isa MOI.Integer ? abs(x-round(x)) :
            error("unsupported audit set $(typeof(set))")
        for ref in all_constraints(b.m; include_variable_in_set_constraints = true)
            obj=constraint_object(ref)
            x=obj.func isa Number ? obj.func : JuMP.value(v->assigned[v], obj.func)
            e=violation(x, obj.set)
            isfinite(e) || error("nonfinite constraint residual")
            count+=1
            if e>residual
                residual=e
                worst=string(ref)
            end
        end
        Dict(
            "valid"=>residual<=tolerance,
            "assigned"=>length(assigned),
            "variables"=>length(variables),
            "constraints"=>count,
            "max_residual"=>residual,
            "worst_constraint"=>worst,
            "tolerance"=>tolerance,
            "objective"=>JuMP.value(v->assigned[v], b.obj),
        )
    catch error
        Dict(
            "valid"=>false,
            "assigned"=>length(assigned),
            "error"=>sprint(showerror, error),
        )
    end
end

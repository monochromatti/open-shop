"""Read the versioned named-object profile; unsupported attributes are errors."""
function _input_keys(x, allowed, label)
    unknown=setdiff(Set(string.(keys(x))), Set(string.(allowed)))
    isempty(unknown) || throw(
        ArgumentError("Unknown $label attributes: $(join(sort!(collect(unknown)), ", "))"),
    )
end
function _read_table(x)
    x===nothing && return nothing
    _input_keys(x, ("x", "y"), "table")
    TableCurve(Float64.(x["x"]), Float64.(x["y"]))
end
function _read_turbine(x)
    x===nothing && return nothing
    _input_keys(
        x,
        (
            "heads",
            "discharge",
            "efficiency",
            "qmin",
            "qmax",
            "interpolation",
            "head_extrapolation",
        ),
        "turbine table",
    )
    vals=x["efficiency"]
    matrix=vals isa AbstractMatrix ? Float64.(vals) :
           permutedims(hcat([Float64.(row) for row in vals]...))
    TurbineTable(
        Float64.(x["heads"]),
        Float64.(x["discharge"]),
        matrix,
        Float64.(x["qmin"]),
        Float64.(x["qmax"]);
        interpolation = Symbol(get(x, "interpolation", "bilinear")),
        head_extrapolation = Symbol(get(x, "head_extrapolation", "error")),
    )
end
function _read_object(T, x)
    _input_keys(x, fieldnames(T), string(nameof(T)))
    vals=Dict{Symbol,Any}()
    for f in fieldnames(T)
        haskey(x, string(f)) || continue
        v=x[string(f)]
        if fieldtype(T, f)==Symbol
            v=Symbol(v)
        elseif f==:curves
            v=[
                begin
                    _input_keys(cv, ("reference_flow", "edges", "weights"), "delay curve")
                    RiverRouting.DelayCurve(
                        cv["reference_flow"],
                        Float64.(cv["edges"]),
                        Float64.(cv["weights"]),
                    )
                end for cv in v
            ]
        elseif f in (
            :level_curve,
            :tailwater_curve,
            :discharge_curve,
            :generator_efficiency_curve,
        )
            v=_read_table(v)
        elseif f==:turbine_table
            v=_read_turbine(v)
        elseif fieldtype(T, f)==Vector{Float64}
            v=Float64.(v)
        end
        vals[f]=v
    end
    T(; vals...)
end
function case_from_dict(d)
    get(d, "schema_version", 1) in (1, 2) ||
        throw(ArgumentError("Unsupported schema version"))
    objects=(
        reservoirs = Reservoir,
        junctions = Junction,
        boundaries = Boundary,
        tunnels = Tunnel,
        plants = Plant,
        generators = Generator,
        river_junctions = RiverJunction,
        rivers = River,
    )
    _input_keys(
        d,
        vcat(
            string.(collect(keys(objects))),
            ["schema_version", "name", "grid", "prices", "operations"],
        ),
        "case",
    )
    s=HydroSystem(;
        Dict(
            k=>[_read_object(T, v) for v in get(d, string(k), [])] for
            (k, T) in pairs(objects)
        )...,
    )
    s=normalize_river_connections(s)
    c=ScheduleCase(
        name = String(d["name"]),
        system = s,
        grid = Float64.(d["grid"]),
        prices = Float64.(d["prices"]),
        operations = [_read_object(OperationalSeries, x) for x in get(d, "operations", [])],
    )
    validate_inputs(c)
    c
end
readcase(path) = case_from_dict(JSON3.read(read(path, String)))
_table_dict(t::TableCurve) = Dict("x"=>t.x, "y"=>t.y)
_table_dict(t::TurbineTable) = Dict(
    string(f)=>getfield(t, f) for f in
    (:heads, :discharge, :efficiency, :qmin, :qmax, :interpolation, :head_extrapolation)
)
function component_dict(x)
    Dict(
        string(f)=>(
            getfield(x, f) isa Union{TableCurve,TurbineTable} ?
            _table_dict(getfield(x, f)) : getfield(x, f)
        ) for f in fieldnames(typeof(x))
    )
end
function case_dict(c)
    d=Dict{String,Any}(
        "schema_version"=>2,
        "name"=>c.name,
        "grid"=>c.grid,
        "prices"=>c.prices,
        "operations"=>component_dict.(c.operations),
    )
    for f in fieldnames(HydroSystem)
        d[string(f)]=component_dict.(getfield(c.system, f))
    end
    public=public_river_connections(c.system)
    d["river_junctions"]=component_dict.(public.river_junctions)
    d["rivers"]=[
        merge(component_dict(r), Dict("curves"=>component_dict.(r.curves))) for
        r in public.rivers
    ]
    d
end
function curve_samples(c)
    s=c.system
    Dict(
        "convention"=>"Explicit synthetic data; turbine and generator efficiencies are separate when tables are supplied; no measured calibration or SHOP file import",
        "reservoirs"=>[
            Dict(
                "name"=>r.name,
                "volume_Mm3"=>collect(range(r.vmin, r.vmax; length = 11)),
                "head_m"=>[head(r, v) for v in range(r.vmin, r.vmax; length = 11)],
            ) for r in s.reservoirs
        ],
        "generators"=>[
            Dict(
                "name"=>g.name,
                "plant"=>g.plant,
                "curves"=>[
                    Dict(
                        "reference_head_m"=>h,
                        "flow_m3s"=>[f*g.qmax for f in [0.35, 0.5, 0.65, 0.8, 0.85, 1.0]],
                        "turbine_efficiency_fraction"=>[
                            efficiency(g, f*g.qmax, h) for
                            f in [0.35, 0.5, 0.65, 0.8, 0.85, 1.0]
                        ],
                    ) for h in [g.hmin, clamp(g.hbest, g.hmin, g.hmax), g.hmax]
                ],
            ) for g in s.generators
        ],
    )
end

"""Identify cases requiring operational and tabulated-physics audit handling."""
function has_operational_data(c)
    !isempty(c.operations) ||
        any(r.level_curve!==nothing for r in c.system.reservoirs) ||
        any(
            g.turbine_table!==nothing ||
                g.generator_efficiency_curve!==nothing ||
                g.shutdown!=0.0 ||
                g.min_efficiency!=0.0 for g in c.system.generators
        ) ||
        any(
            p.tailwater_curve!==nothing || p.outlet_head_floor!==nothing for
            p in c.system.plants
        ) ||
        any(r.discharge_curve!==nothing || r.allow_dry for r in c.system.rivers)
end

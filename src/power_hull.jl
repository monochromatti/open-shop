_power_hull_value(x,assigned)=x isa Number ? Float64(x) : JuMP.value(v->assigned[v],x)
_power_hull_bounds(v)=is_fixed(v) ? (fix_value(v),fix_value(v)) : (lower_bound(v),upper_bound(v))
_power_hull_fixed_off(u)=u isa Number ? u==0 : is_fixed(u)&&fix_value(u)==0

function _power_hull_outward_bounds(lo,hi)
    slack=1e-10*max(1.0,abs(lo),abs(hi))
    (lo-slack,hi+slack)
end

# Conditional ranges change only the on-state box; the source graph and its
# full off-inclusive bounds remain intact. Between head knots the graph is
# linear, so discharge-polynomial extrema at the head knots bound the box.
function _power_hull_eta_bounds(g,qbox,hbox)
    if g.turbine_table===nothing
        return _power_hull_outward_bounds(_global_analytic_eta_bounds(g,qbox...,hbox...)...)
    end
    table=g.turbine_table
    qs=_global_tensor_nodes(table.discharge,qbox...)
    hs=_global_tensor_nodes(table.heads,hbox...)
    lo=Inf;hi=-Inf
    for h in hs
        if length(qs)==1
            value=turbine_efficiency(table,only(qs),h;extrapolation=:linear)
            lo=min(lo,value);hi=max(hi,value)
        else
            for i in 1:(length(qs)-1)
                polynomial=_global_tensor_polynomial(table,qs[i],qs[i+1],h)
                a,b=_global_polynomial_range(polynomial,0.0,1.0)
                lo=min(lo,a);hi=max(hi,b)
            end
        end
    end
    _power_hull_outward_bounds(lo,hi)
end

function _power_hull_on_boxes(g,qbox,hbox,etabox,pmin,electrical_max)
    function restrict_eta(headbox,previous)
        lo,hi=_power_hull_eta_bounds(g,qbox,headbox)
        all(isfinite,(lo,hi)) || return previous
        (max(previous[1],lo),min(previous[2],hi))
    end
    etabox=restrict_eta(hbox,etabox)
    etabox[1]<=etabox[2] || return hbox,etabox
    electrical_upper=_power_hull_outward_bounds(electrical_max,electrical_max)[2]
    denominator=0.00981*qbox[2]*etabox[2]*electrical_upper
    if pmin>0 && electrical_max>0 && isfinite(denominator) && denominator>0
        headlo=pmin/denominator
        if isfinite(headlo)
            # Guard downward: rounding must not exclude the limiting on point.
            headlo-=1e-10*max(1.0,abs(headlo))
            hbox=(max(hbox[1],headlo),hbox[2])
            hbox[1]<=hbox[2] && (etabox=restrict_eta(hbox,etabox))
        end
    end
    hbox,etabox
end

function _power_hull_binary_product_rows!(m,z,u,x,lo,hi)
    @constraint(m,z>=lo*u)
    @constraint(m,z<=hi*u)
    @constraint(m,z>=x-hi*(1-u))
    @constraint(m,z<=x-lo*(1-u))
end

function _power_hull_unit!(m,q,head,eta,power,u,qbox,hbox,etabox,fullhead,fulleta,
    electrical_min,electrical_max;name)
    axes=map(box->box[1]==box[2] ? [box[1]] : [box[1],box[2]],(qbox,hbox,etabox))
    corners=[(qv,hv,ev) for qv in axes[1] for hv in axes[2] for ev in axes[3]]
    weights=@variable(m,[1:length(corners)],lower_bound=0,upper_bound=1,base_name="$(name)_corner")
    onhead=@variable(m,lower_bound=0,upper_bound=hbox[2],base_name="$(name)_onhead")
    oneta=@variable(m,lower_bound=0,upper_bound=etabox[2],base_name="$(name)_oneta")
    @constraint(m,sum(weights)==u)
    @constraint(m,q==sum(corners[k][1]*weights[k] for k in eachindex(corners)))
    @constraint(m,onhead==sum(corners[k][2]*weights[k] for k in eachindex(corners)))
    @constraint(m,oneta==sum(corners[k][3]*weights[k] for k in eachindex(corners)))
    _power_hull_binary_product_rows!(m,onhead,u,head,fullhead...)
    _power_hull_binary_product_rows!(m,oneta,u,eta,fulleta...)
    product=sum(prod(corners[k])*weights[k] for k in eachindex(corners))
    @constraint(m,power<=0.00981*electrical_max*product)
    @constraint(m,power>=0.00981*electrical_min*product)
    (q=q,head=head,eta=eta,power=power,u=u,axes=axes,corners=corners,
     weights=weights,onhead=onhead,oneta=oneta)
end

"""Add the joint trilinear box hull as a redundant envelope of exact power.

On-state head and efficiency moments are linked to the original variables using
their full off-inclusive bounds. Corner weights add no new integer decisions.
The original turbine graphs and nonlinear electrical power equality remain.
Exact conditional efficiency ranges and a guarded minimum-power head bound
shrink only the box used by the added corner rows.
"""
function _add_power_hull!(b,c)
    m=b.m
    haskey(m.ext,:global_power_hulls) && throw(ArgumentError("power hull has already been added"))
    records=Dict{String,Any}()
    m.ext[:global_power_hulls]=records
    before_variables=num_variables(m)
    before_constraints=num_constraints(m;count_variable_in_set_constraints=false)
    skipped=Dict("fixed_off"=>0,"empty_on_box"=>0,"electrical"=>0)
    for (i,g) in enumerate(c.system.generators),t in eachindex(c.prices)
        n="power_hull_$(i)_$(t)"
        u=b.u[i,t]
        if _power_hull_fixed_off(u)
            skipped["fixed_off"]+=1;continue
        end
        hd=b.shared_heads[(plantof(c.system,g).name,t)]
        tensor=get(get(m.ext,:global_tensor_turbines,Dict()),"turbine_$(i)_$(t)",nothing)
        eta=tensor===nothing ? variable_by_name(m,"eta_$(i)_$(t)") : tensor.eta
        fullhead=_power_hull_bounds(hd);fulleta=_power_hull_bounds(eta)
        qbox=(max(0.0,opinterval(c,g.name,:qmin,t,g.qmin)),
              min(g.qmax,opinterval(c,g.name,:qmax,t,g.qmax)))
        # Positive on-state power and nonnegative efficiencies force positive head.
        hbox=(max(0.0,g.hmin,fullhead[1]),min(g.hmax,fullhead[2]))
        etabox=(max(0.0,g.min_efficiency,fulleta[1]),min(1.0,fulleta[2]))
        if !all(box->all(isfinite,box)&&box[1]<=box[2],(qbox,hbox,etabox,fullhead,fulleta))
            skipped["empty_on_box"]+=1;continue
        end
        plo=max(0.0,opinterval(c,g.name,:pmin,t,g.pmin))
        phi=min(g.pmax,opinterval(c,g.name,:pmax,t,g.pmax))
        electrical=if g.generator_efficiency_curve===nothing
            [1.0]
        elseif isfinite(plo)&&isfinite(phi)&&plo<=phi
            curve=g.generator_efficiency_curve
            [table_value(curve,p;extrapolation=:linear) for p in
                _global_tensor_nodes(curve.x,plo,phi)]
        else
            Float64[]
        end
        if isempty(electrical)||!all(isfinite,electrical)||minimum(electrical)<0
            skipped["electrical"]+=1;continue
        end
        hbox,etabox=_power_hull_on_boxes(g,qbox,hbox,etabox,plo,maximum(electrical))
        if hbox[1]>hbox[2] || etabox[1]>etabox[2]
            skipped["empty_on_box"]+=1;continue
        end
        records[n]=_power_hull_unit!(m,b.GQ[i,t],hd,eta,b.P[i,t],u,qbox,hbox,etabox,
            fullhead,fulleta,minimum(electrical),maximum(electrical);name=n)
    end
    profile=Dict("units_added"=>length(records),
        "variables_added"=>num_variables(m)-before_variables,
        "constraints_added"=>num_constraints(m;count_variable_in_set_constraints=false)-before_constraints,
        "corner_weights"=>sum(length(r.weights) for r in values(records);init=0),
        "skipped"=>skipped)
    m.ext[:global_power_hull_profile]=profile
    profile
end

function _power_hull_endpoint_weights(axis,value)
    if length(axis)==1
        abs(value-only(axis))<=1e-7 || throw(DomainError(value,"power-hull singleton start mismatch"))
        return [1.0]
    end
    t=(value-axis[1])/(axis[2]-axis[1])
    # Do not clip weights: inconsistent source values must fail the final audit.
    [1-t,t]
end

"""Lift physical points into rank-one corner weights, including zero off mass."""
function _lift_power_hulls!(m,assigned)
    for record in values(get(m.ext,:global_power_hulls,Dict()))
        u=_power_hull_value(record.u,assigned)
        head=_power_hull_value(record.head,assigned);eta=_power_hull_value(record.eta,assigned)
        assigned[record.onhead]=u*head
        assigned[record.oneta]=u*eta
        if u==0
            for weight in record.weights
                assigned[weight]=0.0
            end
            continue
        end
        coordinate=map(_power_hull_endpoint_weights,record.axes,(_power_hull_value(record.q,assigned),head,eta))
        k=1
        for qi in eachindex(coordinate[1]),hi in eachindex(coordinate[2]),ei in eachindex(coordinate[3])
            assigned[record.weights[k]]=u*coordinate[1][qi]*coordinate[2][hi]*coordinate[3][ei]
            k+=1
        end
    end
    assigned
end

module ProductHull
using JuMP
import OpenSHOP
export add_product_hull!, lift_product_hull!

_value(x,assigned)=x isa Number ? Float64(x) : JuMP.value(v->assigned[v],x)
_bounds(v)=is_fixed(v) ? (fix_value(v),fix_value(v)) : (lower_bound(v),upper_bound(v))
_fixed_off(u)=u isa Number ? u==0 : is_fixed(u)&&fix_value(u)==0

function _binary_product_rows!(m,z,u,x,lo,hi)
    @constraint(m,z>=lo*u)
    @constraint(m,z<=hi*u)
    @constraint(m,z>=x-hi*(1-u))
    @constraint(m,z<=x-lo*(1-u))
end

function _add_unit!(m,q,head,eta,power,u,qbox,hbox,etabox,fullhead,fulleta,
    electrical_min,electrical_max;name,lower_power=false)
    axes=map(box->box[1]==box[2] ? [box[1]] : [box[1],box[2]],(qbox,hbox,etabox))
    corners=[(qv,hv,ev) for qv in axes[1] for hv in axes[2] for ev in axes[3]]
    weights=@variable(m,[1:length(corners)],lower_bound=0,upper_bound=1,base_name="$(name)_corner")
    onhead=@variable(m,lower_bound=0,upper_bound=hbox[2],base_name="$(name)_onhead")
    oneta=@variable(m,lower_bound=0,upper_bound=etabox[2],base_name="$(name)_oneta")
    @constraint(m,sum(weights)==u)
    @constraint(m,q==sum(corners[k][1]*weights[k] for k in eachindex(corners)))
    @constraint(m,onhead==sum(corners[k][2]*weights[k] for k in eachindex(corners)))
    @constraint(m,oneta==sum(corners[k][3]*weights[k] for k in eachindex(corners)))
    _binary_product_rows!(m,onhead,u,head,fullhead...)
    _binary_product_rows!(m,oneta,u,eta,fulleta...)
    product=sum(prod(corners[k])*weights[k] for k in eachindex(corners))
    @constraint(m,power<=0.00981*electrical_max*product)
    lower_power && @constraint(m,power>=0.00981*electrical_min*product)
    (q=q,head=head,eta=eta,power=power,u=u,axes=axes,corners=corners,
     weights=weights,onhead=onhead,oneta=oneta)
end

"""Add the joint trilinear box hull as a redundant envelope of exact power.

On-state head and efficiency moments are linked to the original variables using
their full off-inclusive bounds. Corner weights add no new integer decisions.
The original turbine graphs and nonlinear electrical power equality remain.
"""
function add_product_hull!(b,c;lower_power=false,selected=nothing,max_units=typemax(Int))
    max_units>=0 || throw(ArgumentError("product-hull size limit must be nonnegative"))
    m=b.m
    haskey(m.ext,:product_hull) && throw(ArgumentError("product hull has already been added"))
    records=Dict{String,Any}()
    m.ext[:product_hull]=records
    before_variables=num_variables(m)
    before_constraints=num_constraints(m;count_variable_in_set_constraints=false)
    wanted=selected===nothing ? nothing : Set(string.(selected))
    skipped=Dict("fixed_off"=>0,"empty_on_box"=>0,"electrical"=>0,"size"=>0)
    for (i,g) in enumerate(c.system.generators),t in eachindex(c.prices)
        n="power_hull_$(i)_$(t)"
        wanted!==nothing && !(n in wanted) && continue
        u=b.u[i,t]
        if _fixed_off(u)
            skipped["fixed_off"]+=1;continue
        end
        hd=b.shared_heads[(OpenSHOP.plantof(c.system,g).name,t)]
        tensor=get(get(m.ext,:global_tensor_turbines,Dict()),"turbine_$(i)_$(t)",nothing)
        eta=tensor===nothing ? variable_by_name(m,"eta_$(i)_$(t)") : tensor.eta
        fullhead=_bounds(hd);fulleta=_bounds(eta)
        qbox=(max(0.0,OpenSHOP.opinterval(c,g.name,:qmin,t,g.qmin)),
              min(g.qmax,OpenSHOP.opinterval(c,g.name,:qmax,t,g.qmax)))
        hbox=(max(0.0,g.hmin,fullhead[1]),min(g.hmax,fullhead[2]))
        etabox=(max(0.0,g.min_efficiency,fulleta[1]),min(1.0,fulleta[2]))
        if !all(box->all(isfinite,box)&&box[1]<=box[2],(qbox,hbox,etabox,fullhead,fulleta))
            skipped["empty_on_box"]+=1;continue
        end
        plo=max(0.0,OpenSHOP.opinterval(c,g.name,:pmin,t,g.pmin))
        phi=min(g.pmax,OpenSHOP.opinterval(c,g.name,:pmax,t,g.pmax))
        electrical=if g.generator_efficiency_curve===nothing
            [1.0]
        elseif isfinite(plo)&&isfinite(phi)&&plo<=phi
            curve=g.generator_efficiency_curve
            [OpenSHOP.table_value(curve,p;extrapolation=:linear) for p in
                OpenSHOP._global_tensor_nodes(curve.x,plo,phi)]
        else
            Float64[]
        end
        if isempty(electrical)||!all(isfinite,electrical)||minimum(electrical)<0
            skipped["electrical"]+=1;continue
        end
        if length(records)>=max_units
            skipped["size"]+=1;continue
        end
        records[n]=_add_unit!(m,b.GQ[i,t],hd,eta,b.P[i,t],u,qbox,hbox,etabox,
            fullhead,fulleta,minimum(electrical),maximum(electrical);name=n,lower_power)
    end
    profile=Dict("units_added"=>length(records),"lower_power"=>lower_power,
        "variables_added"=>num_variables(m)-before_variables,
        "constraints_added"=>num_constraints(m;count_variable_in_set_constraints=false)-before_constraints,
        "corner_weights"=>sum(length(r.weights) for r in values(records);init=0),
        "unit_names"=>sort!(collect(keys(records))),"skipped"=>skipped)
    m.ext[:product_hull_profile]=profile
    !isempty(records) && push!(get!(m.ext,:experiment_start_lifters,Any[]),
        assigned->lift_product_hull!(m,assigned))
    profile
end

function _endpoint_weights(axis,value)
    if length(axis)==1
        abs(value-only(axis))<=1e-7 || throw(DomainError(value,"product-hull singleton start mismatch"))
        return [1.0]
    end
    t=(value-axis[1])/(axis[2]-axis[1])
    # Do not clip weights: inconsistent source values must fail the final audit.
    [1-t,t]
end

"""Lift physical points into rank-one corner weights, including zero off mass."""
function lift_product_hull!(m,assigned)
    for record in values(get(m.ext,:product_hull,Dict()))
        u=_value(record.u,assigned)
        head=_value(record.head,assigned);eta=_value(record.eta,assigned)
        assigned[record.onhead]=u*head
        assigned[record.oneta]=u*eta
        if u==0
            for weight in record.weights
                assigned[weight]=0.0
            end
            continue
        end
        coordinate=map(_endpoint_weights,record.axes,(_value(record.q,assigned),head,eta))
        k=1
        for qi in eachindex(coordinate[1]),hi in eachindex(coordinate[2]),ei in eachindex(coordinate[3])
            assigned[record.weights[k]]=u*coordinate[1][qi]*coordinate[2][hi]*coordinate[3][ei]
            k+=1
        end
    end
    assigned
end
end

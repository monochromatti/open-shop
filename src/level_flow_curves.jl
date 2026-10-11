"""Discharge curves indexed by a reference reservoir level (masl).

Each curve maps upstream minus downstream level (m) to signed discharge
(m³/s). Delta grids may differ. Curves must be passive, nondecreasing, and
contain (0,0); forward and reverse discharge need not be symmetric. Values
are interpolated linearly between reference levels, and references outside
the supplied range use the nearest curve. A single curve is independent of
reference level. Delta extrapolation requires `extrapolation=:linear`.
"""
struct LevelFlowCurve
    reference_levels::Vector{Float64}
    curves::Vector{TableCurve}
    function LevelFlowCurve(reference_levels, curves)
        refs=Float64.(reference_levels)
        cc=TableCurve[curves...]
        !isempty(refs) && length(refs)==length(cc) ||
            throw(ArgumentError("level-flow tables need one curve per reference level"))
        all(isfinite,refs) && all(diff(refs).>0) ||
            throw(ArgumentError("reference levels must be finite and strictly increasing"))
        for c in cc
            all(diff(c.y).>=0) ||
                throw(ArgumentError("delta-level discharge must be nondecreasing"))
            origin=findfirst(iszero,c.x)
            origin!==nothing && iszero(c.y[origin]) ||
                throw(ArgumentError("every delta-level discharge table must contain (0,0)"))
            all((x>=0 && y>=0) || (x<=0 && y<=0) for (x,y) in zip(c.x,c.y)) ||
                throw(ArgumentError("delta-level discharge must have the sign of the level difference"))
        end
        new(refs,cc)
    end
end
LevelFlowCurve(;reference_levels,curves)=LevelFlowCurve(reference_levels,curves)

# Knot derivatives follow TableCurve's right-segment convention, except at
# the final reference. Beyond the reference range all reference derivatives
# are zero because the nearest curve is held constant.
function _level_flow_reference(c::LevelFlowCurve,reference)
    isfinite(reference) || throw(DomainError(reference,"reference level must be finite"))
    n=length(c.reference_levels)
    (n==1 || reference<first(c.reference_levels)) && return (1,1,0.0,0.0)
    reference>last(c.reference_levels) && return (n,n,0.0,0.0)
    j=clamp(searchsortedlast(c.reference_levels,reference),1,n-1)
    width=c.reference_levels[j+1]-c.reference_levels[j]
    (j,j+1,(reference-c.reference_levels[j])/width,1/width)
end

function level_flow_value(c::LevelFlowCurve,delta,reference;extrapolation=:error)
    j,k,t,_=_level_flow_reference(c,reference)
    t==1 && return table_value(c.curves[k],delta;extrapolation)
    left=table_value(c.curves[j],delta;extrapolation)
    (j==k || t==0) && return left
    right=table_value(c.curves[k],delta;extrapolation)
    left+t*(right-left)
end

"""Analytic derivatives within the active piecewise bilinear cell.

Returns `(delta, reference, delta_reference)`; both diagonal second
derivatives are zero. Derivatives at knots use the selected segment.
"""
function level_flow_derivatives(c::LevelFlowCurve,delta,reference;extrapolation=:error)
    j,k,t,scale=_level_flow_reference(c,reference)
    left=table_value(c.curves[j],delta;extrapolation)
    dl=table_slope(c.curves[j],delta;extrapolation)
    j==k && return (delta=dl,reference=zero(dl),delta_reference=zero(dl))
    right=table_value(c.curves[k],delta;extrapolation)
    dr=table_slope(c.curves[k],delta;extrapolation)
    (delta=dl+t*(dr-dl),reference=(right-left)*scale,delta_reference=(dr-dl)*scale)
end

"""JuMP operator with explicit derivatives and linear delta extension.
The physical model must constrain or audit its intended delta domain.
"""
function level_flow_operator(m,c::LevelFlowCurve;name=gensym(:level_flow))
    function gradient(g,delta,reference)
        d=level_flow_derivatives(c,delta,reference;extrapolation=:linear)
        g[1]=d.delta
        g[2]=d.reference
        nothing
    end
    function hessian(h,delta,reference)
        d=level_flow_derivatives(c,delta,reference;extrapolation=:linear)
        h[1,1]=0
        h[2,1]=d.delta_reference
        h[2,2]=0
        nothing
    end
    JuMP.add_nonlinear_operator(m,2,
        (delta,reference)->level_flow_value(c,delta,reference;extrapolation=:linear),
        gradient,hessian;name)
end

"""Exact discharge extrema on a rectangle, including linear delta extension."""
function level_flow_bounds(c::LevelFlowCurve,dlo,dhi,rlo,rhi)
    deltas=_global_tensor_nodes(sort!(unique(vcat((t.x for t in c.curves)...))),dlo,dhi)
    refs=_global_tensor_nodes(c.reference_levels,rlo,rhi)
    extrema(level_flow_value(c,d,r;extrapolation=:linear) for d in deltas,r in refs)
end

# The union of delta knots makes each reference column linear on every delta
# cell, even when the supplied curves have different grids. Shared independent
# SOS2 coordinates then describe the exact bilinear graph with one auxiliary
# per reference column, avoiding a variable for every pair of weights.
function _global_level_flow!(m,c::LevelFlowCurve,delta,reference,dlo,dhi,rlo,rhi;
    name=gensym(:level_flow))
    all(isfinite,(rlo,rhi)) && rlo<=rhi ||
        throw(ArgumentError("level-flow reference bounds must be finite and ordered"))
    if length(c.curves)==1
        flow=_global_tensor_table!(m,only(c.curves),delta,dlo,dhi;name)
        table=m.ext[:global_tensor_tables][string(name)]
        get!(m.ext,:global_level_flow_tables,Dict{String,Any}())[string(name)]=
            (dcoordinate=table.coordinate,rcoordinate=nothing,
             nodal=reshape(table.values,:,1),columns=[flow],flow=flow)
        return flow
    end
    knots=sort!(unique(vcat((t.x for t in c.curves)...)))
    dc=_global_tensor_coordinate!(m,delta,_global_tensor_nodes(knots,dlo,dhi);name="$(name)_delta")
    rc=_global_tensor_coordinate!(m,reference,_global_tensor_nodes(c.reference_levels,rlo,rhi);name="$(name)_reference")
    nd,nr=length(dc.nodes),length(rc.nodes)
    nodal=[level_flow_value(c,d,r;extrapolation=:linear) for d in dc.nodes,r in rc.nodes]
    all(isfinite,nodal) || throw(ArgumentError("level-flow extension bounds overflow"))
    columns=@variable(m,[1:nr],base_name="$(name)_column")
    for j in 1:nr
        set_lower_bound(columns[j],minimum(view(nodal,:,j)))
        set_upper_bound(columns[j],maximum(view(nodal,:,j)))
        @constraint(m,columns[j]==sum(nodal[i,j]*dc.weights[i] for i in 1:nd))
    end
    flow=@variable(m,lower_bound=minimum(nodal),upper_bound=maximum(nodal),base_name=string(name))
    @constraint(m,flow==sum(columns[j]*rc.weights[j] for j in 1:nr))
    get!(m.ext,:global_level_flow_tables,Dict{String,Any}())[string(name)]=
        (dcoordinate=dc,rcoordinate=rc,nodal=nodal,columns=columns,flow=flow)
    flow
end

"""Values of the lifted graph variables for a physical delta/reference point."""
function _global_level_flow_values(data,delta,reference)
    dw=_global_tensor_coordinate_weights(data.dcoordinate.nodes,delta)
    isfinite(reference) || throw(DomainError(reference,"reference level must be finite"))
    rw=data.rcoordinate===nothing ? [1.0] :
        _global_tensor_coordinate_weights(data.rcoordinate.nodes,reference)
    values=Dict{VariableRef,Float64}()
    for (coord,weights) in ((data.dcoordinate,dw),(data.rcoordinate,rw))
        coord===nothing && continue
        for (v,x) in zip(coord.weights,weights)
            v isa VariableRef && (values[v]=x)
        end
    end
    flow=0.0
    for j in eachindex(data.columns)
        q=sum(data.nodal[i,j]*dw[i] for i in eachindex(dw))
        values[data.columns[j]]=q
        flow+=q*rw[j]
    end
    values[data.flow]=flow
    values
end

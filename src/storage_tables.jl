# Exact integral of a piecewise-linear level curve in the endpoint's existing
# SOS2 coordinate. A quadratic's chord exceeds it by slope*width²*t*(1-t)/2.
struct _StoragePotentialGraph
    coordinate::_TensorCoordinate
    y::VariableRef
    reference_volume::Float64
    reference_head::Float64
    scale::Float64
end

function _global_storage_potential!(m,c,level,reference_volume,scale;name)
    coordinate=level.coordinate
    reference_head=table_value(c,reference_volume;extrapolation=:linear)
    values=[_level_potential(c,reference_volume,reference_head,x)/scale
        for x in coordinate.nodes]
    weights=coordinate.weights
    nodes=coordinate.nodes
    correction=[table_slope(c,nodes[k];extrapolation=:linear)*
        (nodes[k+1]-nodes[k])^2/(2scale) for k in 1:(length(nodes)-1)]
    y=@variable(m,lower_bound=minimum(values),upper_bound=maximum(values),base_name=string(name))
    @constraint(m,y==sum(values[k]*weights[k] for k in eachindex(values))-
        sum(correction[k]*weights[k]*weights[k+1] for k in eachindex(correction);init=0.0))
    _StoragePotentialGraph(coordinate,y,reference_volume,reference_head,scale)
end

function _global_storage_mean_slope!(m,r,a,b,h,level,potential,lo,hi;name)
    c=r.level_curve
    slopes=[table_slope(c,x;extrapolation=:linear) for x in
        _global_tensor_nodes(c.x,lo,hi)]
    lower,upper=extrema(slopes)
    slope=@variable(m,lower_bound=lower*r.vmax/500,upper_bound=upper*r.vmax/500,
        base_name=string(name))
    change=(b-a)/r.vmax
    # H = Z(a) + D*(b-a), D in [min Z'/2,max Z'/2]. This bounded
    # divided difference fixes H=Z(a) at b=a, where the energy row alone is 0=0.
    @constraint(m,h/250==level.y/250+slope*change)
    @constraint(m,change*(h-potential[1].reference_head)/250==potential[2].y-potential[1].y)
    slope
end

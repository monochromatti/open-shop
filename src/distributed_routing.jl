# Shared reference interpolation for distributed river transport. The nonlinear
# dispatch evaluates the original neighboring-reference polynomial; the global
# model selects the same reference pair through one shared SOS2 coordinate.

"""A shared exact reference coordinate for all transfers of one release.

One or two curves need no coordinate: their released transfers are linear or
quadratic. Identical delay distributions also need no coordinate. Otherwise the
caller keeps the physical release bounds on `q` and reuses this coordinate for
arrivals, terminal inventory and observation rows.
"""
function distributed_reference_coordinate!(m, curves, q; name=gensym(:delay_reference))
    RiverRouting.check_curves(curves)
    length(curves)<=2 && return nothing
    all(c->c.edges==first(curves).edges && c.weights==first(curves).weights,curves) &&
        return nothing
    references=[c.reference_flow for c in curves]
    coordinates=get!(m.ext,:distributed_reference_coordinates,Dict{Symbol,Any}())
    key=Symbol(name)
    if haskey(coordinates,key)
        existing=coordinates[key]
        isequal(existing.x,q) && existing.nodes==references ||
            throw(ArgumentError("delay reference name is already used by another release"))
        return existing
    end
    coordinate=_global_tensor_coordinate!(m,q,references;name)
    coordinates[key]=coordinate
    coordinate
end

"""Populate exact SOS2 reference weights from an assigned physical release."""
function _lift_distributed_references!(m, assigned)
    for coordinate in values(get(m.ext,:distributed_reference_coordinates,Dict()))
        q=JuMP.value(v->assigned[v],coordinate.x)
        weights=_global_tensor_coordinate_weights(coordinate.nodes,q)
        for (variable,value) in zip(coordinate.weights,weights)
            variable isa VariableRef && (assigned[variable]=value)
        end
    end
    assigned
end

"""Exact released-transfer expression, in the same flow units as `q`.

With a shared `coordinate`, this is `q*sum(coefficient*weight)` on its selected
SOS2 pair. Otherwise a scalar operator evaluates the neighboring-reference
piecewise quadratic with analytic derivatives. Numerical trial points use the
endpoint polynomial continuation; physical bounds retain the supplied domain.
"""
function distributed_transfer_expression(m, curves, q, coefficients; coordinate=nothing)
    RiverRouting.check_transfer_coefficients(curves,coefficients)
    values=Float64.(collect(coefficients))
    if length(curves)==1 || all(==(first(values)),values)
        return first(values)*q
    end
    if coordinate!==nothing
        length(coordinate.weights)==length(curves) ||
            throw(ArgumentError("delay coordinate does not match its reference curves"))
        return q*sum(values[j]*coordinate.weights[j] for j in eachindex(values))
    end
    references=[c.reference_flow for c in curves]
    slopes=diff(values)./diff(references)
    if length(curves)==2
        return q*(values[1]+slopes[1]*(q-references[1]))
    end
    cache=get!(m.ext,:distributed_transfer_operators,Dict{Any,Any}())
    key=(Tuple(references),Tuple(values))
    operator=get!(cache,key) do
        value(x)=begin
            j=RiverRouting.reference_segment(references,x)
            x*(values[j]+slopes[j]*(x-references[j]))
        end
        derivative(x)=begin
            j=RiverRouting.reference_segment(references,x)
            values[j]+slopes[j]*(x+(x-references[j]))
        end
        second_derivative(x)=2slopes[RiverRouting.reference_segment(references,x)]
        JuMP.add_nonlinear_operator(m,1,value,derivative,second_derivative;
            name=gensym(:delay_transfer))
    end
    operator(q)
end

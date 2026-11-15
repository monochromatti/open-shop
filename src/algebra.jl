# Convex hull of z = state*x for a binary state and lo <= x <= hi.
# These rows also give its standard continuous relaxation.
function _binary_product_rows!(m,z,state,x,lo,hi)
    @constraint(m,z>=lo*state)
    @constraint(m,z<=hi*state)
    @constraint(m,z>=x-hi*(1-state))
    @constraint(m,z<=x-lo*(1-state))
    nothing
end

function _bounded_binary_product!(m,x,state,lo,hi;name)
    state isa Number && return state*x
    z=@variable(m,lower_bound=min(0.0,lo),upper_bound=max(0.0,hi),base_name=name)
    _binary_product_rows!(m,z,state,x,lo,hi)
    z
end

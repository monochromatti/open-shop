# Convex hull of z = state*x for a binary state and lo <= x <= hi.
# These rows also give its standard continuous relaxation.
function _binary_product_rows!(m,z,state,x,lo,hi)
    @constraint(m,z>=lo*state)
    @constraint(m,z<=hi*state)
    @constraint(m,z>=x-hi*(1-state))
    @constraint(m,z<=x-lo*(1-state))
    nothing
end

"""Register first and second derivatives of a small scalar nonlinear operator.

Missing multivariate Hessians make Ipopt approximate curvature for the entire
model. Supply both callbacks and reuse their ForwardDiff workspaces. Each
operator owns its buffers; the supplied gradient, when present, is retained.
"""
function autodiff_operator(m,n,f,provided_gradient=nothing;name)
    value=x->f(x...)
    point=zeros(n)
    hessian=zeros(n,n)
    gradient_config=ForwardDiff.GradientConfig(value,point)
    hessian_config=ForwardDiff.HessianConfig(value,point)
    gradient=provided_gradient===nothing ? (g,args...)->begin
        copyto!(point,args)
        ForwardDiff.gradient!(g,value,point,gradient_config)
        nothing
    end : provided_gradient
    second_derivative=(H,args...)->begin
        copyto!(point,args)
        ForwardDiff.hessian!(hessian,value,point,hessian_config)
        for i in 1:n,j in 1:i
            H[i,j]=hessian[i,j]
        end
        nothing
    end
    JuMP.add_nonlinear_operator(m,n,f,gradient,second_derivative;name)
end

function _bounded_binary_product!(m,x,state,lo,hi;name)
    state isa Number && return state*x
    z=@variable(m,lower_bound=min(0.0,lo),upper_bound=max(0.0,hi),base_name=name)
    _binary_product_rows!(m,z,state,x,lo,hi)
    z
end

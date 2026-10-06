"""Lift turbine discharge times efficiency into one bounded physical variable.

The cuts require nonnegative discharge, zero discharge while off, and operating
efficiency in `[min_efficiency,1]` while on. The caller supplies those constraints.
Efficiency itself retains its numerical off-state continuation, including
negative values; neither its graph nor its bounds are changed here.
"""
function _global_effective_flow!(m,q,eta::VariableRef,qmax,min_efficiency;
    name=gensym(:effective_flow))
    isfinite(qmax) && qmax>=0 || throw(ArgumentError("effective flow needs a finite nonnegative discharge limit"))
    isfinite(min_efficiency) && 0<=min_efficiency<=1 ||
        throw(ArgumentError("operating minimum efficiency must lie in [0,1]"))
    has_upper_bound(eta) && isfinite(upper_bound(eta)) ||
        throw(ArgumentError("effective flow needs a finite efficiency upper bound"))
    upper_efficiency=min(1.0,upper_bound(eta))
    flow=@variable(m,lower_bound=0,upper_bound=qmax*max(0.0,upper_efficiency),base_name=string(name))
    @constraint(m,flow==q*eta)
    @constraint(m,flow>=min_efficiency*q)
    @constraint(m,flow<=upper_efficiency*q)
    get!(m.ext,:global_effective_flows,Dict{String,Any}())[string(name)]=
        (q=q,eta=eta,flow=flow,qmax=Float64(qmax),
         min_efficiency=Float64(min_efficiency),max_efficiency=upper_efficiency)
    flow
end

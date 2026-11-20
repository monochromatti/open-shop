using OpenSHOP, JuMP, Ipopt, LinearAlgebra
BLAS.set_num_threads(1)
source=read(joinpath(@__DIR__,"../test/head_flow_rules_tests.jl"),String)
a=findfirst("function flow_rule_case",source).start
b=findnext("@testset",source,a).start
include_string(Main,source[a:prevind(source,b)])

for opening in (0.5,0.85), method in (:none,:adaptive,:recalc,:relaxed)
    c=flow_rule_case(river_kwargs=(min_allowed_flow=1.0,flow_change_cost=2.0,
        flow_start_cost=7.0,flow_stop_cost=3.0,))
    c=OpenSHOP._river_replace(c;grid=[0.0,0.5,2.0],prices=zeros(2))
    u=zeros(Int,0,2)
    warm=dispatch_from_controls(c,u,zeros(0,2),reshape([opening,0.0],1,2))
    model=OpenSHOP._build_dispatch(c;u,warm,river_on=reshape([1,0],1,2))
    OpenSHOP._local_optimizer!(model.m,"Ipopt",120.0)
    if method==:adaptive
        set_optimizer_attribute(model.m,"mu_strategy","adaptive")
    elseif method==:recalc
        set_optimizer_attribute(model.m,"recalc_y","yes")
    elseif method==:relaxed
        set_optimizer_attribute(model.m,"bound_relax_factor",1e-8)
    end
    unset_silent(model.m)
    set_optimizer_attribute(model.m,"print_level",4)
    println("BEGIN_RIVER ",opening," ",method)
    if method in (:physical,:adaptive,:recalc,:relaxed)
        OpenSHOP.retain_feasible_dispatch!(model)
    elseif method!=:none
        point=Float64[];lower=Float64[];upper=Float64[];violation=Float64[]
        callback=function(mode,args...)
            mode==0 || return true
            native=unsafe_backend(model.m).inner
            if method in (:violations,:both)
                OpenSHOP.ipopt_primal_residual!(lower,upper,violation,native)
            end
            if method in (:iterate,:both)
                resize!(point,native.n)
                Ipopt.GetIpoptCurrentIterate(native,false,native.n,point,
                    C_NULL,C_NULL,native.m,C_NULL,C_NULL)
            end
            true
        end
        MOI.set(model.m,Ipopt.CallbackFunction(),callback)
    end
    elapsed=@elapsed optimize!(model.m)
    println("RIVER_DIAGNOSTIC ",repr((;opening,method,
        status=termination_status(model.m),objective=value(model.obj),
        release=value.(model.RQ),gate=value.(model.variables.a),elapsed,
        iterations=MOI.get(model.m,MOI.BarrierIterations()));context=:limit=>false))
    raw=OpenSHOP._dispatch_values(model)
    raw.diagnostics["validation"]=validate(c,raw)
    report=Dict{String,Any}()
    reconstructed=OpenSHOP._reconstruct_dispatch!(c,raw,report)
    println("RECONSTRUCTED ",repr((;opening,method,objective=reconstructed.objective,
        valid=get(report,"forward_reconstructed",false),report);context=:limit=>false))
    flush(stdout)
end

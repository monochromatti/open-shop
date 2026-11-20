using OpenSHOP, JuMP, Test
source=read(joinpath(@__DIR__,"../test/head_flow_rules_tests.jl"),String)
a=findfirst("function flow_rule_case",source).start
b=findnext("@testset",source,a).start
include_string(Main,source[a:prevind(source,b)])
@testset "Linux river optimum from feasible starting releases" begin
    for opening in (0.25,0.28,0.37,0.5,0.58,0.85,1.0)
        c=flow_rule_case(river_kwargs=(min_allowed_flow=1.0,flow_change_cost=2.0,
            flow_start_cost=7.0,flow_stop_cost=3.0,))
        c=OpenSHOP._river_replace(c;grid=[0.0,0.5,2.0],prices=zeros(2))
        u=zeros(Int,0,2)
        warm=dispatch_from_controls(c,u,zeros(0,2),reshape([opening,0.0],1,2))
        result=solve_case(c;u,warm,river_on=reshape([1,0],1,2),time_limit=120.0)
        println("RIVER_FIXED ",(;opening,status=result.status,objective=result.objective))
        flush(stdout)
        @test result.accepted
        @test result.objective≈-14.0 atol=1e-3
        @test validate(c,result.solution)["valid"]
        @test validate_replay(c,result.solution)["valid"]
    end
end

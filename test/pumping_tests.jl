include("../examples/pumped_storage.jl")
using .PumpedStorageExample: pumped_storage_case

pump_replace(x;kwargs...)=OpenSHOP._river_replace(x;kwargs...)

@testset "Pump work, water conservation and signed electrical exchange" begin
    c=pumped_storage_case()
    u=[0 0 1 1;1 1 0 0]
    q=20.0 .* u
    gate=zeros(0,4)
    x=dispatch_from_controls(c,u,q,gate)
    @test x.diagnostics["validation"]["valid"]
    @test x.storage[1,:]≈[5.0,5.072,5.144,5.072,5.0] atol=1e-10
    @test x.storage[2,:]≈[15.0,14.928,14.856,14.928,15.0] atol=1e-10
    @test vec(sum(x.storage;dims=1))≈fill(20.0,5) atol=1e-10
    expected=zeros(2,4)
    for t in 1:4
        upper=299.5+0.1*(x.storage[1,t]+x.storage[1,t+1])/2
        lower=48.5+0.1*(x.storage[2,t]+x.storage[2,t+1])/2
        h=upper-lower
        expected[1,t]=0.00981*q[1,t]*h*0.91
        expected[2,t]=0.00981*q[2,t]*h/(0.88*0.96)
    end
    @test x.power≈expected atol=1e-10
    @test net_power(c,x)≈vec(expected[1,:]-expected[2,:])
    @test x.objective≈sum(c.prices .* net_power(c,x))-20.0 atol=1e-8
    @test sum(expected[1,:])<sum(expected[2,:])
    @test power(c.system.generators[2],0.0,0.0)==0.0
    idle=dispatch_from_controls(c,zeros(Int,2,4),zeros(2,4),gate)
    @test idle.diagnostics["validation"]["valid"]
    @test idle.power==zeros(2,4)
    @test net_power(c,idle)==zeros(4)
    @test replay_audit(c,x)["valid"]
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    @test case_from_dict(case_dict(c)).system.generators[2].mode==:pump
    @test schedule_dict(x)["power"]≈expected atol=1e-10
    @test_throws DimensionMismatch net_power(c,zeros(1,4))

    continued=restart_case(c,x,1.0)
    rest=dispatch_from_controls(continued,u[:,2:end],q[:,2:end],zeros(0,3))
    @test rest.diagnostics["validation"]["valid"]
    @test rest.storage≈x.storage[:,2:end] atol=1e-10
    @test rest.power≈x.power[:,2:end] atol=1e-10
    @test continued.system.generators[2].initial_on==1
    @test continued.system.generators[2].initial_power≈x.power[2,1]
    @test continued.system.commitment_exclusions==c.system.commitment_exclusions
    @test replay_audit(continued,rest)["valid"]
end

@testset "Hydraulic pump tables and connected pressure network" begin
    c=pumped_storage_case(;table=true,connected=true)
    u=[0 0 1 1;1 1 0 0]
    q=20.0 .* u
    x=dispatch_from_controls(c,u,q,zeros(0,4))
    @test x.diagnostics["validation"]["valid"]
    @test vec(x.tunnel_q)≈[20.0,20.0,0.0,0.0] atol=1e-9
    for t in 1:2
        lift=x.head[3,t]-x.head[2,t]
        @test x.head[3,t]-x.head[1,t]≈0.002*20.0^2 atol=1e-8
        eta=0.86+0.04*(lift-200.0)/100.0
        @test x.power[2,t]≈0.00981*20.0*lift/(eta*0.96) atol=1e-8
    end
    @test replay_audit(c,x)["valid"]
    local_result=solve_case(c;u,warm=x,feasibility_only=true,time_limit=20.0)
    @test local_result.accepted
    @test local_result.solution.storage≈x.storage atol=1e-8
    @test local_result.solution.power≈x.power atol=1e-7
    # The polynomial global model keeps the reciprocal-efficiency equality;
    # it must accept pump consumption above every generation power support.
    global_result=solve(c;initial=x,fixed_u=u,fixed_tunnel_opening=ones(1,4),
        time_limit=40.0,relative_gap=1e-4)
    @test global_result.accepted
    @test global_result.solution.power≈x.power atol=1e-5
end

@testset "Reversible mode selection and independently rejected schedules" begin
    c=pumped_storage_case()
    u=[0 0 1 1;1 1 0 0]
    x=dispatch_from_controls(c,u,20.0 .* u,zeros(0,4))
    result=solve(c;initial=x,time_limit=60.0,relative_gap=1e-4)
    @test result.accepted
    @test result.global_certificate
    @test result.solution.commitment==u
    @test result.solution.objective≈x.objective atol=1e-3
    @test all(vec(sum(result.solution.commitment;dims=1)).<=1)
    proposal=OpenSHOP.propose_commitment(c;time_limit=10.0)
    @test haskey(proposal,"u")
    @test all(vec(sum(proposal["u"];dims=1)).<=1)
    @test proposal["u"]==u
    @test_throws ArgumentError dispatch_from_controls(c,ones(Int,2,4),fill(20.0,2,4),zeros(0,4))
    @test_throws ArgumentError simulate(c;generator_q=fill(20.0,2,4),gate=zeros(0,4))
    bad=schedule_dict(x);bad["power"][2,1]*=(0.88*0.96)^2
    @test !validate(c,bad)["valid"]
    bad=schedule_dict(x);bad["u"][1,1]=1
    @test !validate(c,bad)["valid"]
    bad=schedule_dict(x);bad["objective"]+=2c.prices[1]*x.power[2,1]
    @test !validate(c,bad)["valid"]

    units=copy(c.system.generators)
    units[2]=pump_replace(units[2];initial_on=1,initial_age=0.25,minup=1.5)
    history=pump_replace(c;system=pump_replace(c.system;generators=units))
    @test !OpenSHOP.admissible(history,zeros(Int,2,4))
    hu=[0 0 1 1;1 1 0 0]
    hx=dispatch_from_controls(history,hu,20.0 .* hu,zeros(0,4))
    @test hx.diagnostics["validation"]["valid"]
    continued=restart_case(history,hx,1.0)
    @test continued.system.generators[2].initial_age≈1.25
    @test !OpenSHOP.admissible(continued,[1 1 1;0 0 0])
    @test OpenSHOP.admissible(continued,hu[:,2:end])
end

@testset "Pump scope and exclusion input admission" begin
    c=pumped_storage_case()
    pump=c.system.generators[2]
    for wrong in (
        pump_replace(pump;mode=:unknown),
        pump_replace(pump;motor_efficiency=0.0),
        pump_replace(pump;motor_efficiency=1.01),
        pump_replace(pump;hmin=0.0),
        pump_replace(pump;generator_efficiency_curve=TableCurve([0.0,65.0],[0.9,0.95])),
        pump_replace(pump;plant=:Generating),
    )
        invalid=pump_replace(c;system=pump_replace(c.system;generators=[c.system.generators[1],wrong]))
        @test_throws ArgumentError OpenSHOP.validate_inputs(invalid)
    end
    for rule in (
        CommitmentExclusion(name=:Bad,groups=[[:Turbine]]),
        CommitmentExclusion(name=:Bad,groups=[[:Turbine],Symbol[]]),
        CommitmentExclusion(name=:Bad,groups=[[:Turbine],[:Turbine]]),
        CommitmentExclusion(name=:Bad,groups=[[:Turbine],[:Missing]]),
    )
        invalid=pump_replace(c;system=pump_replace(c.system;commitment_exclusions=[rule]))
        @test_throws ArgumentError OpenSHOP.validate_inputs(invalid)
    end
    units=[pump_replace(g;initial_on=1) for g in c.system.generators]
    @test_throws ArgumentError OpenSHOP.validate_inputs(pump_replace(c;
        system=pump_replace(c.system;generators=units)))

    # Exclusion chooses a mode group, not one unit globally: parallel pumps
    # in the same group can run together.
    extra=pump_replace(pump;name=:ParallelPump)
    grouped=pump_replace(c;system=pump_replace(c.system;
        generators=[c.system.generators;extra],
        plants=[c.system.plants[1],pump_replace(c.system.plants[2];pmax=130.0)],
        commitment_exclusions=[CommitmentExclusion(name=:Reversible,
            groups=[[:Turbine],[:Pump,:ParallelPump]])]))
    parallel_u=[0 0 0 0;1 1 0 0;1 1 0 0]
    @test OpenSHOP.admissible(grouped,parallel_u)
    parallel=dispatch_from_controls(grouped,parallel_u,20.0 .* parallel_u,zeros(0,4))
    @test parallel.diagnostics["validation"]["valid"]
    @test parallel.storage[1,3]≈5.288 atol=1e-10
    @test !OpenSHOP.admissible(grouped,[1 0 0 0;1 1 0 0;0 0 0 0])
end

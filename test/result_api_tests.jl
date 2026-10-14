@testset "Typed schedules and one solve-result interface" begin
    fixture=MixedSourceNetworkExample.mixed_source_small_fixture()
    c=fixture.case
    x=dispatch_from_controls(c; commitment=fixture.u, generator_q=fixture.q, gate=fixture.gates)
    @test x isa Schedule
    @test x.storage==x["V"]
    @test x.commitment==x["u"]
    @test_throws MethodError setindex!(x,0.0,"objective")
    profile=schedule_dict(x)
    profile["power"][1]+=1.0
    @test profile["power"]!=x.power
    @test validate(c,x)["valid"]
    profile=schedule_dict(x)
    @test OpenSHOP._schedule_input(OpenSHOP.jsonready(x);case=c).storage==x.storage
    @test validate(c,OpenSHOP.jsonready(x))["valid"]
    for key in ("generator_q","V","u","objective")
        missing=deepcopy(profile)
        delete!(missing,key)
        @test !validate(c,missing)["valid"]
    end
    invalid=deepcopy(profile)
    invalid["power"][1]+=1.0
    bad=OpenSHOP._schedule_input(invalid;case=c)
    bad.diagnostics["validation"]=validate(c,bad)
    @test !bad.diagnostics["validation"]["valid"]
    result=SolveResult(Dict("status"=>"LOCALLY_SOLVED","accepted"=>false);solution=bad)
    @test result.objective==bad.objective
    @test result.feasible_lower_bound===nothing
    @test result["feasible_lower_bound"]===nothing
    @test OpenSHOP.jsonready(result)["feasible_lower_bound"]===nothing
    empty=SolveResult(Dict("status"=>"TIME_LIMIT"))
    @test empty.solution===nothing
    @test !haskey(empty,"objective")
    @test empty.feasible_lower_bound===nothing
    @test occursin("TIME_LIMIT",sprint(show,MIME"text/plain"(),empty))
    invalid=deepcopy(profile)
    invalid["u"]=Float64.(invalid["u"])
    invalid["u"][1]=0.5
    @test OpenSHOP._schedule_input(invalid;case=c).commitment[1]==0.5
    @test !validate(c,invalid)["valid"]
    invalid=deepcopy(profile)
    invalid["u"]=BigFloat.(invalid["u"])
    invalid["u"][1]=big"1e-400"
    @test OpenSHOP._schedule_input(invalid;case=c).commitment[1]==big"1e-400"
    @test !validate(c,invalid)["valid"]
    invalid["u"][1]=prevfloat(BigFloat(1))
    @test OpenSHOP._schedule_input(invalid;case=c).commitment[1]!=1
    @test !validate(c,invalid)["valid"]
    invalid=deepcopy(profile)
    invalid["u"]=vec(invalid["u"])
    @test_throws ArgumentError OpenSHOP._schedule_input(invalid;case=c)
    @test !validate(c,invalid)["valid"]
    invalid=deepcopy(profile)
    invalid["generator_q"]=BigFloat.(invalid["generator_q"])
    invalid["generator_q"][1]=big"1e-400"
    @test_throws ArgumentError OpenSHOP._schedule_input(invalid;case=c)
    @test !validate(c,invalid)["valid"]
end

@testset "Native cases resolve direct river connections" begin
    f=MixedSourceNetworkExample.mixed_source_fixture()
    c=f.case
    public=public_river_connections(c.system)
    input=OpenSHOP._river_replace(c.system;
        rivers=public.rivers,river_junctions=public.river_junctions)
    before=[(r.source,r.target) for r in input.rivers]
    native=ScheduleCase(name=c.name,system=input,grid=c.grid,prices=c.prices,operations=c.operations)
    positional=ScheduleCase(c.name,input,c.grid,c.prices,c.operations,c.flow_requirements)
    repeated=ScheduleCase(name=c.name,system=native.system,grid=c.grid,prices=c.prices,operations=c.operations)
    imported=case_from_dict(case_dict(native))
    @test [(r.source,r.target) for r in input.rivers]==before
    @test case_dict(native)==case_dict(positional)==case_dict(repeated)==case_dict(imported)
    @test OpenSHOP.validate_inputs(native)
    @test simulate(native;generator_q=f.q,gate=f.gates)==simulate(imported;generator_q=f.q,gate=f.gates)
end

@testset "Missing optional accounting and explicit selectable gates" begin
    c=tg_flow_case()
    tunnel=only(c.system.tunnels)
    tunnel=OpenSHOP._river_replace(tunnel;gate_adjustment_cost=10.0,initial_opening=0.0)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;tunnels=[tunnel]))
    x=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3);
        tunnel_opening=reshape([0.0,0.5,1.0],1,3))
    @test x.gate_adjustment_cost>0
    profile=schedule_dict(x)
    delete!(profile,"gate_adjustment_cost")
    @test validate(c,profile)["valid"]
    @test !haskey(OpenSHOP._schedule_input(profile;case=c),"gate_adjustment_cost")
    delete!(profile,"tunnel_opening")
    @test !validate(c,profile)["valid"]
end

@testset "Soft release accounting remains required at zero shortfall" begin
    lake=Reservoir(name=:Lake,z0=100.0,slope=0.0,v0=5.0,vmin=0.0,vmax=10.0,water_value=0.0)
    river=River(name=:Outlet,source=:Lake,target=:Sea,law=:controlled,
        capacity=10.0,water_value=0.0,curves=RiverRouting.DelayCurve[],
        deterministic_delay=0.0,history_grid=[-1.0,0.0],history_release=[0.0])
    system=HydroSystem(reservoirs=[lake],rivers=[river],boundaries=[Boundary(name=:Sea,head=0.0)])
    c=ScheduleCase(name="release accounting",system=system,grid=[0.0,1.0],prices=[0.0],operations=[
        OperationalSeries(object=:Outlet,attribute=:min_release,times=[0.0],values=[5.0]),
        OperationalSeries(object=:Outlet,attribute=:release_penalty,times=[0.0],values=[10.0]),
    ])
    x=dispatch_from_controls(c,zeros(Int,0,1),zeros(0,1),ones(1,1))
    @test x.diagnostics["validation"]["valid"]
    @test x.release_penalty_cost==0.0
    for key in ("shortfall_release","release_penalty_cost")
        profile=schedule_dict(x)
        delete!(profile,key)
        audit=validate(c,profile)
        @test !audit["valid"]
        @test "missing soft release accounting" in audit["errors"]
    end
end

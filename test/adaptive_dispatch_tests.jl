@testset "Adaptive integration resolves small passive reservoirs" begin
    for scale in (0.01,1.0,100.0)
        c=passive_replay_case(scale)
        initial=dispatch_from_controls(c,zeros(Int,0,1),zeros(0,1),zeros(0,1))
        result=solve_verified(c;u=zeros(Int,0,1),initial,max_refinements=24,
            feasibility_only=true,time_limit=10.0)
        @test result["accepted"]
        if result["accepted"]
            x=result["solution"]
            @test validate(result["problem"],x)["valid"]
            @test replay_audit(result["problem"],x)["valid"]
            @test x.storage[1,end]/scale ≈ 0.0002 atol=0.2e-6
            @test all(t in x.grid for t in c.grid)
            @test length(x.grid)>length(c.grid)
        end
    end
end

@testset "Coarse overshoot does not reject an admissible passive trajectory" begin
    base=passive_replay_case(1.0)
    reservoir=OpenSHOP._river_replace(only(base.system.reservoirs);vmin=0.000199)
    c=OpenSHOP._river_replace(base;system=OpenSHOP._river_replace(base.system;reservoirs=[reservoir]))
    coarse=dispatch_from_controls(c,zeros(Int,0,1),zeros(0,1),zeros(0,1))
    @test !validate(c,coarse)["valid"]
    result=solve_verified(c;u=zeros(Int,0,1),warm=coarse,feasibility_only=true,
        max_refinements=24,time_limit=1e-9)
    @test result["accepted"]
    if result["accepted"]
        @test validate(result["problem"],result["solution"])["valid"]
        @test minimum(result["solution"].storage)>=reservoir.vmin-1e-8
        @test result["solution"].storage[1,end]≈0.0002 atol=0.1e-6
    end
end

@testset "Local refinement respects an executed gate schedule" begin
    base=passive_replay_case(1.0)
    tunnel=OpenSHOP._river_replace(only(base.system.tunnels);gate_positions=[0.5,1.0],
        initial_opening=0.5,gate_adjustment_cost=7.0)
    flows=[OpenSHOP.tunnel_flow(tunnel,2.0,a) for a in (0.5,1.0)]
    inflow=OperationalSeries(object=:Forebay,attribute=:inflow,times=[0.0,1.0],values=flows)
    c=OpenSHOP._river_replace(base;grid=[0.0,1.0,2.0],prices=[0.0,0.0],operations=[inflow],
        system=OpenSHOP._river_replace(base.system;tunnels=[tunnel]))
    opening=reshape([0.5,1.0],1,2)
    initial=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),zeros(0,2);tunnel_opening=opening)
    result=solve_verified(c;u=zeros(Int,0,2),initial,tunnel_opening=opening,
        feasibility_only=true,max_refinements=24,time_limit=10.0)
    @test result["accepted"]
    if result["accepted"]
        x=result["solution"]
        @test validate(result["problem"],x)["valid"]
        mesh=OpenSHOP.TimeDiscretization(c.grid,x.grid)
        @test x.tunnel_opening==opening[:,mesh.interval]
        @test objective_components(c,x).gate_cost ≈ 3.5
        @test replay_audit(result["problem"],x)["valid"]
    end
    # An unchecked selected position does not replace the case's prescribed
    # default when dispatch did not request that position.
    default=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;
        tunnels=[OpenSHOP._river_replace(tunnel;opening=1.0)]))
    if result["accepted"]
        resumed=solve_verified(default;u=zeros(Int,0,2),warm=result["solution"],
            time_limit=1e-9,max_refinements=0)
        @test !resumed["accepted"] || all(resumed["solution"].tunnel_opening.==1.0)
    end
end

@testset "A numerical start survives an exhausted optimization budget only after reconstruction" begin
    c=warm_dispatch_fixture(;rivers=false)
    u=ones(Int,2,2)
    warm=dispatch_from_controls(c,u,[3.0 3.0;2.8 1.8],zeros(0,2))
    expected=warm.objective
    result=solve_verified(c;u,warm,time_limit=1e-9,max_refinements=0)
    @test result["accepted"]
    if result["accepted"]
        @test result["solution"].objective>=expected-1e-8
        @test validate(result["problem"],result["solution"])["valid"]
    end
    @test warm.objective==expected
    worse=dispatch_from_controls(c,u,[1.0 1.0;2.5 1.5],zeros(0,2))
    @test validate(c,worse)["valid"]
    @test worse.objective<warm.objective
    interrupted=(case;kwargs...)->SolveResult(Dict("status"=>"TIME_LIMIT");solution=worse)
    retained=OpenSHOP._verified_dispatch(c,interrupted;u,warm,max_refinements=0)
    @test retained.accepted
    @test retained.solution.objective>=expected-1e-8
    @test validate(retained.case,retained.solution)["valid"]
    lake=OpenSHOP._river_replace(only(c.system.reservoirs);inflow=12.0)
    changed=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;reservoirs=[lake]),
        prices=[100.0,20.0])
    @test !validate(changed,warm)["valid"]
    updated=solve_verified(changed;u,warm,time_limit=1e-9,max_refinements=0)
    @test updated["accepted"]
    if updated["accepted"]
        @test validate(changed,updated["solution"])["valid"]
        @test updated["solution"].storage != warm.storage
        @test objective_components(changed,updated["solution"]).total ≈ updated["solution"].objective
    end
    expired=solve_verified(changed;u,warm,deadline=time()-1)
    @test !expired["accepted"]

    base=warm_dispatch_fixture()
    reach=OpenSHOP._river_replace(only(base.system.rivers);law=:controlled,min_allowed_flow=1.0)
    fixed=OpenSHOP._river_replace(base;system=OpenSHOP._river_replace(base.system;rivers=[reach]))
    off=dispatch_from_controls(fixed,zeros(Int,2,2),zeros(2,2),zeros(1,2))
    @test validate(fixed,off)["valid"]
    changed_state=solve_verified(fixed;u=zeros(Int,2,2),warm=off,river_on=ones(Int,1,2),
        time_limit=1e-9,max_refinements=0)
    @test !changed_state["accepted"]
end

@testset "Adaptive dispatch retains original power schedules and reserve decisions" begin
    c=resolved_operating_case(reserve_groups=[ReserveGroup(name=:Reserve,generators=[:Unit],
        up_requirement=0.1,up_max=0.1,response_hours=0.25)])
    grid=[0.0,0.2,1.0,1.25,2.8,3.0]
    u=ones(Int,1,2)
    allocation=ReserveAllocation(up=fill(0.1,1,1,2),down=zeros(1,1,2),generator_q=[1.0 2.0])
    warm=dispatch_from_controls(c,u,[1.0 2.0],zeros(0,2);grid,reserves=allocation)
    averages=OpenSHOP.operating_window_values(c,warm)
    schedule=OperationalSeries(object=:Unit,attribute=:power,times=[0.0,1.0],values=vec(averages.power))
    c=OpenSHOP._river_replace(c;operations=[schedule])
    result=solve_verified(c;u,warm,time_limit=10.0,max_refinements=4)
    @test result["accepted"]
    if result["accepted"]
        x=result["solution"]
        @test validate(result["problem"],x)["valid"]
        @test OpenSHOP.operating_window_values(c,x).power ≈ averages.power atol=1e-4
        @test objective_components(c,x).total ≈ x.objective
        @test all(isapprox.(x.reserves.up,0.1;atol=1e-6))
    end
end

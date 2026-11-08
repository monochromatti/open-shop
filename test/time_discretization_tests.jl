@testset "Operating observations preserve durations and water volumes" begin
    mesh=OpenSHOP.TimeDiscretization([0.0,1.0,3.0],[0.0,0.1,0.6,1.0,2.5,3.0])
    rates=[2.0 4.0 8.0 3.0 7.0]
    means=OpenSHOP.interval_means(mesh,rates)
    @test means ≈ [5.4 4.0]
    water=0.0036 .* rates .* permutedims(mesh.dt)
    @test OpenSHOP.interval_sums(mesh,water) ≈ 0.0036 .* means .* [1.0 2.0]
    @test sum(OpenSHOP.interval_sums(mesh,water)) ≈ sum(water)
    same=OpenSHOP.TimeDiscretization(mesh.grid)
    @test OpenSHOP.interval_means(same,rates) ≈ rates
    model=Model()
    @variable(model,q[1:1,1:5])
    set_start_value.(q,rates)
    @test value.(start_value,OpenSHOP.interval_means(mesh,q)) ≈ means
    @test_throws ArgumentError OpenSHOP.TimeDiscretization([0.0,1.0,3.0],[0.0,0.5,3.0])
    @test_throws ArgumentError OpenSHOP.TimeDiscretization([0.0,1.0],[0.0,2.0])
end

@testset "Generated meshes share vertices without merging input events" begin
    operating=[0.0,0.4,1.2,2.0,3.3,4.2,6.0]
    transport=OpenSHOP.refined_grid(operating;factor=4)
    physical=OpenSHOP.refined_grid(operating;existing=transport)
    for level in 1:6
        transport=OpenSHOP.refined_grid(transport;existing=physical)
        physical=OpenSHOP.bulk_refined_grid(physical,ones(length(physical)-1);existing=transport)
        @test all(t in transport for t in physical)
        @test all(t in physical for t in operating)
        @test minimum(diff(transport))≈0.4/(4*2^level) rtol=1e-10
    end
    input=[0.0,0.6,nextfloat(0.6),1.0]
    refined=OpenSHOP.refined_grid(input)
    @test all(t in refined for t in input)
    @test all(diff(refined).>0)
    shared=OpenSHOP.refined_grid([0.4,0.8];existing=[0.4,0.6,0.8])
    @test minimum(diff(sort!(unique!(vcat(shared,[0.6])))))≈0.2
end

@testset "Reconstruction retains operating controls on unequal hydraulic cells" begin
    c=fixed_pump_case()
    grid=[0.0,0.1,0.5,0.8,1.25,1.3,2.0]
    x=dispatch_from_controls(c,[1 0 1],zeros(1,3),zeros(0,3);grid)
    @test validate(c,x)["valid"]
    @test x.generator_q[1,1] != x.generator_q[1,2]
    @test x.commitment == [1 1 0 0 1 1]
    @test x.objective ≈ -sum(c.prices[OpenSHOP.TimeDiscretization(c.grid,grid).interval] .*
        diff(grid).*vec(x.power)) atol=1e-8
    @test validate(c,schedule_dict(x))["valid"]
    mktemp() do path,io
        close(io)
        writejson(path,x)
        restored=OpenSHOP.JSON3.read(read(path,String),Dict{String,Any})
        @test validate(c,restored)["valid"]
    end
    @test objective_components(c,x).total ≈ x.objective atol=1e-8
    restarted=restart_case(c,x,0.5)
    average=sum(x.power[1,1:2].*diff(grid)[1:2])/0.5
    @test restarted.system.generators[1].initial_power ≈ average
    @test restarted.system.generators[1].initial_discharge ≈
        sum(x.generator_q[1,1:2].*diff(grid)[1:2])/0.5
    continued=dispatch_from_controls(restarted,[0 1],zeros(1,2),zeros(0,2);grid=grid[3:end])
    @test validate(restarted,continued)["valid"]
    @test continued.storage ≈ x.storage[:,3:end] atol=1e-9
end

function resolved_operating_case(;operations=OperationalSeries[],soft_constraints=SoftConstraint[],
    flow_requirements=FlowRequirement[],level_windows=LevelWindow[],reserve_groups=ReserveGroup[])
    lake=Reservoir(;name=:Lake,z0=100.0,slope=1.0,v0=5.0,vmin=1.0,vmax=10.0,
        water_value=0.0,level_windows)
    unit=Generator(name=:Unit,plant=:Station,qmin=0.0,qmax=10.0,pmin=0.0,pmax=15.0,
        qbest=5.0,hbest=100.0,hmin=90.0,hmax=120.0,minup=0.0,mindown=0.0,
        efficiency=1.0,qcurvature=0.0,hcurvature=0.0,startup=7.0,shutdown=3.0,initial_on=0)
    system=HydroSystem(reservoirs=[lake],boundaries=[Boundary(name=:Sea,head=0.0)],
        plants=[Plant(name=:Station,source=:Lake,target=:Sea,pmax=15.0)],generators=[unit])
    ScheduleCase(;name="resolved operating windows",system,grid=[0.0,1.0,3.0],prices=[10.0,100.0],
        operations,soft_constraints,flow_requirements,reserve_groups)
end

@testset "Resolved optimization preserves schedules, penalties and assessment windows" begin
    grid=[0.0,0.2,1.0,1.25,2.8,3.0]
    u=ones(Int,1,2)
    base=resolved_operating_case()
    reference=dispatch_from_controls(base,u,[1.0 2.0],zeros(0,2);grid)
    averages=OpenSHOP.operating_window_values(base,reference)
    schedule=OperationalSeries(object=:Unit,attribute=:power,times=[0.0,1.0],values=vec(averages.power))
    soft=SoftConstraint(name=:Target,object=:Lake,quantity=:volume,time=1.0,
        lower=5.0,max_lower_slack=0.01,lower_cost=100.0)
    gauge=FlowRequirement(name=:Gauge,generators=[:Unit],min_flow=0.0,max_flow=2.0,
        volume_windows=[VolumeWindow(start_time=0.5,end_time=2.5,min_volume=0.01,max_volume=0.02)])
    c=resolved_operating_case(operations=[schedule],soft_constraints=[soft],flow_requirements=[gauge],
        level_windows=[LevelWindow(start_time=0.0,end_time=3.0,max_change=0.1)])
    warm=dispatch_from_controls(c,u,[1.0 2.0],zeros(0,2);grid)
    @test validate(c,warm)["valid"]
    result=solve_case(c;u,grid,warm,time_limit=60.0,feasibility_only=true)
    @test result.accepted
    if result.solution!==nothing
        x=result.solution
        @test validate(result.problem,x)["valid"]
        @test OpenSHOP.operating_window_values(c,x).power ≈ averages.power atol=1e-5
        @test x.soft_slacks[1].lower ≈ 0.0036 atol=1e-6
        @test objective_components(c,x).transition_cost ≈ 7.0 atol=1e-8
        @test objective_components(c,x).total ≈ x.objective atol=1e-7
        @test x.generator_q[1,1] ≈ x.generator_q[1,2] atol=1e-9
        @test replay_audit(c,x)["valid"]
    end
end

@testset "Refined primary pumping optimization" begin
    c=fixed_pump_case();grid=[0.0,0.1,0.5,0.8,1.25,1.3,2.0]
    warm=dispatch_from_controls(c,[1 0 1],zeros(1,3),zeros(0,3);grid)
    result=solve_case(c;u=[1 0 1],grid,warm,feasibility_only=true,time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        @test validate(c,result.solution)["valid"]
        @test result.solution.generator_q[1,1] != result.solution.generator_q[1,2]
    end
end

@testset "Resolved reserve capacity is held while its physical witness changes" begin
    group=ReserveGroup(name=:Spinning,generators=[:Unit],up_requirement=0.1,up_max=2.0,
        up_price=5.0,response_hours=0.25)
    c=resolved_operating_case(reserve_groups=[group])
    grid=[0.0,0.2,1.0,1.25,2.8,3.0]
    allocation=ReserveAllocation(up=fill(0.1,1,1,2),down=zeros(1,1,2),generator_q=[1.0 2.0])
    warm=dispatch_from_controls(c,ones(Int,1,2),[1.0 2.0],zeros(0,2);grid,reserves=allocation)
    @test validate(c,warm)["valid"]
    @test objective_components(c,warm).reserve_revenue ≈ 1.5 atol=1e-10
    result=solve_case(c;u=ones(Int,1,2),grid,warm,feasibility_only=true,time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        x=result.solution
        @test validate(c,x)["valid"]
        @test x.reserves.up[1,1,1] ≈ x.reserves.up[1,1,2] atol=1e-9
        @test objective_components(c,x).total ≈ x.objective atol=1e-7
    end
    response=OperationalSeries(object=:Unit,attribute=:ramp_up,times=[0.0,1.0],values=[4.0,0.4])
    response_group=OpenSHOP._river_replace(group;response_hours=0.5)
    restricted=OpenSHOP._river_replace(c;operations=[response],reserve_groups=[response_group])
    allocation=ReserveAllocation(up=reshape([1.8,0.1],1,1,2),down=zeros(1,1,2),generator_q=[1.0 2.0])
    coarse=dispatch_from_controls(restricted,ones(Int,1,2),[1.0 2.0],zeros(0,2);reserves=allocation)
    fine=dispatch_from_controls(restricted,ones(Int,1,2),[1.0 2.0],zeros(0,2);grid,reserves=allocation)
    @test validate(restricted,coarse)["valid"]
    @test !validate(restricted,fine)["valid"]
end

@testset "Resolved river continuation retains executed cohorts" begin
    reservoirs=[Reservoir(name=name,z0=100.0,slope=1.0,v0=5.0,vmin=1.0,vmax=10.0,
        water_value=0.0) for name in (:Upper,:Lower)]
    river=River(name=:Reach,source=:Upper,target=:Lower,law=:orifice,coefficient=1.0,
        crest=0.0,capacity=20.0,water_value=0.0,deterministic_delay=0.3,
        curves=RiverRouting.DelayCurve[],history_grid=[-1.0,0.0],history_release=[0.0])
    c=ScheduleCase(name="resolved river restart",system=HydroSystem(;reservoirs,rivers=[river]),
        grid=[0.0,1.0,2.0],prices=[0.0,0.0])
    grid=[0.0,0.1,0.9,1.0,1.05,1.4,2.0]
    x=dispatch_from_controls(c,zeros(Int,0,2),zeros(0,2),ones(1,2);grid)
    @test validate(c,x)["valid"]
    continued=restart_case(c,x,1.0)
    y=dispatch_from_controls(continued,zeros(Int,0,1),zeros(0,1),ones(1,1);grid=grid[4:end])
    @test validate(continued,y)["valid"]
    @test y.storage ≈ x.storage[:,4:end] atol=1e-9
    @test y.arrival_volume ≈ x.arrival_volume[:,4:end] atol=1e-9
    @test y.terminal_transit ≈ x.terminal_transit atol=1e-9
    @test continued.system.rivers[1].initial_release ≈
        sum(x.river_release[1,1:3].*diff(grid)[1:3]) atol=1e-9
end

@testset "Physical release deficits cannot be hidden by an interval average" begin
    reservoir=Reservoir(name=:Lake,z0=-400.0,slope=100.0,v0=5.0,vmin=4.5,vmax=6.0,water_value=0.0)
    river=River(name=:Outlet,source=:Lake,target=:Sea,law=:orifice,coefficient=1.0,
        capacity=10.0,gate_min=0.5,water_value=0.0,deterministic_delay=0.1,
        curves=RiverRouting.DelayCurve[],history_grid=[-1.0,0.0],history_release=[5.0])
    system=HydroSystem(reservoirs=[reservoir],boundaries=[Boundary(name=:Sea,head=0.0)],rivers=[river])
    base=ScheduleCase(name="release assessment window",system=system,grid=[0.0,2.0],prices=[0.0])
    grid=[0.0,0.5,1.0,1.25,2.0]
    x=dispatch_from_controls(base,zeros(Int,0,1),zeros(0,1),fill(0.5,1,1);grid)
    requirement=OpenSHOP.operating_window_values(base,x).river_release[1,1]
    operations=[OperationalSeries(object=:Outlet,attribute=attribute,times=[0.0],values=[value])
        for (attribute,value) in ((:gate_max,0.5),(:min_release,requirement),(:release_penalty,1000.0))]
    c=OpenSHOP._river_replace(base;operations)
    x=dispatch_from_controls(c,zeros(Int,0,1),zeros(0,1),fill(0.5,1,1);grid)
    @test validate(c,x)["valid"]
    @test minimum(x.river_release) < requirement
    @test OpenSHOP.operating_window_values(c,x).river_release[1,1] ≈ requirement atol=1e-12
    expected=1000*0.0036*sum(diff(grid).*max.(0.0,requirement.-vec(x.river_release)))
    @test expected>0.0
    @test x.release_penalty_cost ≈ expected atol=1e-10
    @test objective_components(c,x).release_penalty ≈ expected atol=1e-10
    result=solve_case(c;u=zeros(Int,0,1),grid,warm=x,feasibility_only=true,time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        @test result.solution.release_penalty_cost ≈ expected atol=1e-6
    end
end

@testset "A plant capacity limits each physical cell" begin
    c=resolved_operating_case()
    grid=[0.0,0.2,1.0,1.25,2.8,3.0]
    u=ones(Int,1,2)
    x=dispatch_from_controls(c,u,fill(10.0,1,2),zeros(0,2);grid)
    average=maximum(OpenSHOP.operating_window_values(c,x).power)
    limit=(average+maximum(x.power))/2
    plant=OpenSHOP._river_replace(only(c.system.plants);pmax=limit)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;plants=[plant]))
    @test maximum(OpenSHOP.operating_window_values(c,x).power) < limit
    @test !validate(c,x)["valid"]
    result=solve_case(c;u,grid,warm=x,time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        @test maximum(result.solution.power)<=limit+1e-4
    end
end

@testset "Level windows include an interior physical storage maximum" begin
    window=LevelWindow(start_time=0.0,end_time=1.0,max_change=0.1)
    base=resolved_operating_case(level_windows=[window])
    lake=OpenSHOP._river_replace(only(base.system.reservoirs);z0=-400.0,slope=100.0)
    river=River(name=:Return,source=:Upstream,target=:Lake,law=:controlled,capacity=2.0,
        water_value=0.0,deterministic_delay=1.0,curves=RiverRouting.DelayCurve[],
        history_grid=[-2.0,-1.0,-0.5,0.0],history_release=[0.0,2.0,0.0])
    upstream=Reservoir(name=:Upstream,z0=100.0,slope=1.0,v0=5.0,vmin=1.0,vmax=10.0,water_value=0.0)
    system=OpenSHOP._river_replace(base.system;reservoirs=[lake,upstream],rivers=[river])
    c=OpenSHOP._river_replace(base;system,grid=[0.0,1.0],prices=[10.0])
    x=dispatch_from_controls(c,ones(Int,1,1),ones(1,1),zeros(1,1);grid=[0.0,0.5,1.0])
    @test x.storage[1,end] ≈ x.storage[1,1] atol=1e-10
    @test head(lake,x.storage[1,2])-head(lake,x.storage[1,1]) ≈ 0.18 atol=1e-8
    @test !validate(c,x)["valid"]
    @test validate(c,x)["residuals"]["level_window_Lake_1"] ≈ 0.08 atol=1e-8
end

@testset "Running gross head and startup eligibility retain operating windows" begin
    c=resolved_operating_case()
    unit=OpenSHOP._river_replace(only(c.system.generators);gross_head_min=104.9,
        gross_head_max=106.0,startup_gross_head_min=104.99)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;generators=[unit]),
        operations=[OperationalSeries(object=:Unit,attribute=:discharge,times=[0.0,1.0],values=[10.0,10.0])])
    grid=[0.0,0.2,1.0,1.25,2.8,3.0]
    warm=dispatch_from_controls(c,ones(Int,1,2),fill(10.0,1,2),zeros(0,2);grid)
    @test validate(c,warm)["valid"]
    result=solve_case(c;u=ones(Int,1,2),grid,warm,feasibility_only=true,time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        @test validate(c,result.solution)["valid"]
        @test OpenSHOP.operating_window_values(c,result.solution).generator_q ≈ fill(10.0,1,2) atol=1e-6
    end
end

@testset "Physical refinement creates no extra control opportunities" begin
    c=resolved_operating_case()
    grid=[0.0,0.2,1.0,1.25,2.8,3.0]
    finer=with_grid(c,grid)
    changing=dispatch_from_controls(finer,ones(Int,1,5),[1.0 2.0 2.0 3.0 3.0],zeros(0,5))
    @test validate(finer,changing)["valid"]
    @test !validate(c,changing)["valid"]
    held=dispatch_from_controls(c,ones(Int,1,2),[1.0 2.0],zeros(0,2);grid)
    @test validate(c,held)["valid"]
    @test solve_verified(c;u=ones(Int,1,2),initial=held,max_refinements=0)["accepted"]
    @test_throws ArgumentError solve(c;initial=held,preparation_time_limit=0.0,time_limit=5.0)
    result=solve_case(c;u=ones(Int,1,2),grid,warm=held,feasibility_only=true,
        replay_grid=OpenSHOP.refined_grid(grid),time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        @test validate(c,result.solution)["valid"]
        @test replay_audit(c,result.solution)["valid"]
    end
end

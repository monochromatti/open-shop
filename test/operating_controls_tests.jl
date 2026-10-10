# These fixtures use explicit water balances and constant efficiency. Expected
# powers, transition costs, and ramp integrals are computed independently here.
const OPERATING_TEST_ALPHA = 0.00981 * 100.0 * 0.9

operating_series(object, attribute, times, values) = OperationalSeries(
    object=object, attribute=attribute, times=Float64.(times), values=Float64.(values))

function operating_case(; grid=[0.0, 0.5, 2.0, 3.0], units=1, plants=1,
    operations=OperationalSeries[], qmax=10.0, initial_on=1,
    reservoir_kwargs=(;), generator_kwargs=(;), plant_kwargs=(;), river_kwargs=(;))
    lake=Reservoir(; name=:Lake, z0=99.5, slope=0.1, v0=5.0, vmin=4.0,
        vmax=6.0, water_value=0.0, reservoir_kwargs...)
    ps=[Plant(; name=Symbol("Plant", i), source=:Lake, target=:Sea,
        pmax=20.0, ramp=100.0, plant_kwargs...) for i in 1:plants]
    gs=[Generator(; name=Symbol("Unit", i), plant=ps[1+mod(i-1,plants)].name,
        qmin=1.0, qmax=qmax, pmin=0.1, pmax=10.0, efficiency=0.9,
        qbest=3.0, hbest=100.0, hmin=80.0, hmax=120.0,
        qcurvature=0.0, hcurvature=0.0, initial_on=initial_on,
        initial_age=10.0, minup=0.0, mindown=0.0, startup=20.0,
        shutdown=5.0, generator_kwargs...) for i in 1:units]
    reach=River(; name=:River, source=:Lake, target=:Sea, law=:controlled,
        capacity=4.0, curves=RiverRouting.DelayCurve[], deterministic_delay=0.0,
        water_value=0.0, history_grid=[-2.0,-1.0,0.0], history_release=[1.0,1.0],
        river_kwargs...)
    system=HydroSystem(; reservoirs=[lake], junctions=Junction[],
        boundaries=[Boundary(name=:Sea,head=0.0)], tunnels=Tunnel[],
        plants=ps, generators=gs, rivers=[reach])
    ScheduleCase(; name="independent operating controls", system,
        grid=Float64.(grid), prices=fill(50.0,length(grid)-1), operations)
end

function operating_replace(c; reservoirs=c.system.reservoirs,
    generators=c.system.generators, plants=c.system.plants, rivers=c.system.rivers,
    operations=c.operations, prices=c.prices)
    system=OpenSHOP._river_replace(c.system; reservoirs, generators, plants, rivers)
    OpenSHOP._river_replace(c; system, operations, prices)
end

function operating_balanced(c, q; releases=fill(1.0,length(c.prices)))
    inflow=operating_series(:Lake,:inflow,c.grid[1:end-1],
        vec(sum(q;dims=1))+releases)
    operating_replace(c; operations=vcat(
        [z for z in c.operations if !(z.object==:Lake && z.attribute==:inflow)], [inflow]))
end

# To test commitment tampering, reconstruct an actual physical trajectory with
# only state controls removed, then independently audit the original restrictions.
function operating_physical_without_states(c,u,q,gate)
    cc=operating_replace(c;operations=[z for z in c.operations
        if z.attribute ∉ (:maintenance,:forced_on,:power,:discharge)])
    dispatch_from_controls(cc,u,q,gate)
end

function operating_two_plant_case()
    alpha=OPERATING_TEST_ALPHA
    c=operating_case(;units=4,plants=2,qmax=6.0,initial_on=0)
    # Group consecutive units in a plant; their combined controls are prescribed.
    generators=[OpenSHOP._river_replace(g;plant=Symbol("Plant",i<=2 ? 1 : 2))
        for (i,g) in enumerate(c.system.generators)]
    plants=[OpenSHOP._river_replace(p;pmin=0.5,qmin=1.0,qmax=6.0)
        for p in c.system.plants]
    knot=[0.0,0.5,2.0]
    operations=[
        operating_series(:Unit1,:maintenance,knot,[0,1,0]),
        operating_series(:Unit1,:forced_on,knot,[1,-1,1]),
        operating_series(:Unit2,:forced_on,knot,[0,1,1]),
        operating_series(:Plant2,:maintenance,knot,[0,1,0]),
        operating_series(:Plant2,:forced_on,knot,[1,-1,1]),
        operating_series(:Unit4,:maintenance,[0.0],[1]),
        operating_series(:Plant1,:discharge,knot,[2,2,3]),
        # Before .5h this optional hard schedule must be absent, not zero.
        operating_series(:Unit1,:discharge,[0.5,2.0],[0,2]),
        operating_series(:Plant2,:power,knot,[2alpha,0,2alpha]),
        operating_series(:Unit1,:startup,[0.0,2.0],[50,80]),
        operating_series(:Unit1,:shutdown,knot,[5,15,5]),
        operating_series(:River,:gate_min,[0.0],[0.25]),
        operating_series(:River,:gate_max,[0.0],[0.25]),
    ]
    u=[1 0 1; 0 1 1; 1 0 1; 0 0 0]
    q=[2.0 0.0 2.0; 0.0 2.0 1.0; 2.0 0.0 2.0; 0.0 0.0 0.0]
    c=operating_replace(c;generators,plants,operations,prices=[30.0,90.0,60.0])
    operating_balanced(c,q),u,q,fill(0.25,1,3)
end

@testset "Operating schedules, maintenance and operational transition costs" begin
    c,u,q,gate=operating_two_plant_case()
    @test OpenSHOP.validate_inputs(c)
    @test admissible(c,u)
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    @test OpenSHOP.optional_operation(c,:Unit1,:discharge,1)===nothing
    x=dispatch_from_controls(c,u,q,gate)
    @test x["validation"]["valid"]
    @test x["V"]≈fill(5.0,1,4) atol=1e-10
    @test x["power"]≈OPERATING_TEST_ALPHA*q atol=1e-9
    expected_revenue=OPERATING_TEST_ALPHA*(30*0.5*4+90*1.5*2+60*1*5)
    # Unit1 50+15+80; Unit2 20; Unit3 20+5+20; Unit4 no transitions.
    @test x["objective"]≈expected_revenue-210.0 atol=1e-8
    @test OpenSHOP.transition_costs(c,u)==210.0
    @test replay_audit(c,x)["valid"]
    continued=restart_case(c,x,0.5)
    remainder=dispatch_from_controls(continued,u[:,2:end],q[:,2:end],gate[:,2:end])
    @test remainder["validation"]["valid"]
    @test replay_audit(continued,remainder)["valid"]
    @test remainder["V"]≈x["V"][:,2:end] atol=1e-10
    @test remainder["power"]≈x["power"][:,2:end] atol=1e-9
    # Historical output windows use the executed .5h interval at this edge.
    @test all(g.initial_interval_hours==0.5 for g in continued.system.generators)
    @test all(p.initial_interval_hours==0.5 for p in continued.system.plants)

    # Changing a unit schedule while preserving the plant sum is still rejected.
    changed=copy(q); changed[1,3]=1.5; changed[2,3]=1.5
    bad=dispatch_from_controls(c,u,changed,gate)
    @test !bad["validation"]["valid"]
    @test bad["validation"]["residuals"]["schedule_Unit1"]≈0.5 atol=1e-9
    changed=copy(q); changed[1,1]=3.0
    bad=dispatch_from_controls(c,u,changed,gate)
    @test !bad["validation"]["valid"]
    @test bad["validation"]["residuals"]["schedule_Plant1"]≈1.0 atol=1e-9
    changed=copy(q); changed[3,3]=2.5
    bad=dispatch_from_controls(c,u,changed,gate)
    @test !bad["validation"]["valid"]
    @test bad["validation"]["residuals"]["schedule_Plant2"]>0.4
    changed_u=copy(u); changed_u[1,2]=1; changed_u[2,2]=0
    changed=copy(q); changed[1,2]=2.0; changed[2,2]=0.0
    @test !admissible(c,changed_u)
    bad=operating_physical_without_states(c,changed_u,changed,gate)
    @test !validate(c,bad)["valid"]
    changed_u=copy(u); changed_u[3,2]=1
    changed=copy(q); changed[3,2]=2.0
    @test !admissible(c,changed_u)
    bad=operating_physical_without_states(c,changed_u,changed,gate)
    @test !validate(c,bad)["valid"]

    # Both numerical optimizers receive the same physical schedules, and the
    # independent audits check their outputs instead of trusting solver status.
    local_result=solve_case(c;u,warm=x,time_limit=45.0)
    @test haskey(local_result,"validation")
    if haskey(local_result,"validation")
        @test local_result["validation"]["valid"]
        @test local_result["objective"]≈x["objective"] atol=1e-3
        @test replay_audit(c,local_result)["valid"]
    end
    global_result=solve(c;initial=x,time_limit=30.0,relative_gap=1e-3)
    @test global_result["accepted"]
    @test global_result["global_certificate"]
    if global_result["solution"]!==nothing
        @test global_result["solution"]["objective"]≈x["objective"] atol=1e-3
    end
    proposed=OpenSHOP.propose_commitment(c;reference=x,time_limit=10.0)
    @test haskey(proposed,"u")
    if haskey(proposed,"u")
        @test admissible(c,proposed["u"])
        @test proposed["u"]==u
    end
end

@testset "Aggregate plant operating limits apply only to a running plant" begin
    c=operating_case(;units=2,initial_on=0,qmax=6.0)
    plant=only(c.system.plants)
    plant=OpenSHOP._river_replace(plant;pmin=2.0,qmin=3.0,qmax=6.0)
    c=operating_replace(c;plants=[plant])
    u=[1 0 1;1 0 0]
    q=[2.0 0.0 3.0;2.0 0.0 0.0]
    c=operating_balanced(c,q)
    @test dispatch_from_controls(c,u,q,fill(0.25,1,3))["validation"]["valid"]
    for (replacement, expected) in (([2.0,5.0],1.0),([1.0,1.0],1.0))
        changed=copy(q); changed[:,1]=replacement
        bad=dispatch_from_controls(c,u,changed,fill(0.25,1,3))
        @test !bad["validation"]["valid"]
        @test bad["validation"]["residuals"]["plant_capacity_Plant1"]>=expected-1e-8
    end
    lower_power=operating_replace(c;plants=[OpenSHOP._river_replace(plant;pmax=3.0)])
    bad=dispatch_from_controls(lower_power,u,q,fill(0.25,1,3))
    @test !bad["validation"]["valid"]
    @test bad["validation"]["residuals"]["plant_capacity_Plant1"]>0.5
end

@testset "Nonuniform midpoint ramps integrate piecewise rates in physical hours" begin
    for (attribute, values) in ((:discharge_ramp_up,[2.0,4.4,6.8]),
        (:discharge_ramp_down,[6.8,4.4,2.0]))
        ops=[operating_series(object,attribute,[0.0,0.5,2.0],[1.0,3.0,0.5])
            for object in (:Unit1,:Plant1)]
        c=operating_case(;operations=ops)
        q=reshape(values,1,3)
        c=operating_balanced(c,q)
        u=ones(Int,1,3); gate=fill(0.25,1,3)
        x=dispatch_from_controls(c,u,q,gate)
        @test x["validation"]["valid"]
        @test replay_audit(c,x)["valid"]
        for object in (only(c.system.generators),only(c.system.plants))
            @test OpenSHOP.ramp_allowance(c,object,attribute,0.25,1.25)≈2.5 atol=1e-12
            @test OpenSHOP.ramp_allowance(c,object,attribute,1.25,2.5)≈2.5 atol=1e-12
        end
        changed=copy(q); changed[1,2]+=attribute==:discharge_ramp_up ? 0.2 : -0.2
        bad=dispatch_from_controls(c,u,changed,gate)
        @test !bad["validation"]["valid"]
        @test bad["validation"]["residuals"]["discharge_ramp_Unit1"]≈0.1 atol=1e-9
        @test bad["validation"]["residuals"]["discharge_ramp_Plant1"]≈0.1 atol=1e-9
    end
    # Power rate integration uses exactly the same midpoint windows, in MW/h.
    ops=[operating_series(object,:ramp_up,[0.0,0.5,2.0],
        OPERATING_TEST_ALPHA*[1.0,3.0,0.5]) for object in (:Unit1,:Plant1)]
    c=operating_case(;operations=ops)
    q=reshape([2.0,4.4,6.8],1,3); c=operating_balanced(c,q)
    x=dispatch_from_controls(c,ones(Int,1,3),q,fill(0.25,1,3))
    @test x["validation"]["valid"]
    @test replay_audit(c,x)["valid"]
    changed=copy(q); changed[1,2]+=0.2
    bad=dispatch_from_controls(c,ones(Int,1,3),changed,fill(0.25,1,3))
    @test !bad["validation"]["valid"]
    @test bad["validation"]["residuals"]["unit_ramp_Unit1"]>0.08
end

@testset "Historical unit, plant and river midpoint windows" begin
    c=operating_case(;generator_kwargs=(initial_discharge=1.0,initial_interval_hours=2.0,
        discharge_ramp_up=1.0), plant_kwargs=(initial_discharge=1.0,
        initial_interval_hours=2.0,discharge_ramp_up=1.0))
    q=fill(2.2,1,3); c=operating_balanced(c,q)
    x=dispatch_from_controls(c,ones(Int,1,3),q,fill(0.25,1,3))
    @test x["validation"]["valid"] # -1h to .25h: allowance1.25, not .5.
    @test replay_audit(c,x)["valid"]
    changed=copy(q); changed[1,1]=2.3
    bad=dispatch_from_controls(c,ones(Int,1,3),changed,fill(0.25,1,3))
    @test !bad["validation"]["valid"]
    @test bad["validation"]["residuals"]["discharge_ramp_Unit1"]≈0.05 atol=1e-9
    @test bad["validation"]["residuals"]["discharge_ramp_Plant1"]≈0.05 atol=1e-9
    c=operating_case(;river_kwargs=(history_grid=[-2.0,0.0],history_release=[1.0],ramp_up=1.0))
    q=fill(2.0,1,3); c=operating_balanced(c,q;releases=fill(2.2,3))
    x=dispatch_from_controls(c,ones(Int,1,3),q,fill(2.2/4,1,3))
    @test x["validation"]["valid"]
    @test replay_audit(c,x)["valid"]
    bad=dispatch_from_controls(c,ones(Int,1,3),q,fill(2.3/4,1,3))
    @test !bad["validation"]["valid"]
    @test bad["validation"]["residuals"]["river_ramp_River"]≈0.05 atol=1e-9
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
end

@testset "Discharge ramps bind in both numerical dispatch models" begin
    operations=[operating_series(object,:discharge_ramp_up,[0.0,0.5,2.0],
        [1.0,3.0,0.5]) for object in (:Unit1,:Plant1)]
    append!(operations,[operating_series(:River,:gate_min,[0.0],[0.25]),
        operating_series(:River,:gate_max,[0.0],[0.25])])
    c=operating_case(;operations,generator_kwargs=(initial_discharge=1.25,
        discharge_ramp_up=1.0))
    q=reshape([2.0,4.5,7.0],1,3);c=operating_balanced(c,q)
    u=ones(Int,1,3);x=dispatch_from_controls(c,u,q,fill(0.25,1,3))
    @test x["validation"]["valid"]
    local_result=solve_case(c;u,warm=x,time_limit=20.0)
    @test get(get(local_result,"validation",Dict()),"valid",false)
    if haskey(local_result,"generator_q")
        @test local_result["generator_q"]≈q atol=2e-4
        @test replay_audit(c,local_result)["valid"]
    end
    global_result=solve(c;initial=x,fixed_u=u,time_limit=20.0)
    @test global_result["accepted"]
    @test global_result["global_certificate"]
    if global_result["solution"]!==nothing
        @test global_result["solution"]["generator_q"]≈q atol=2e-4
    end
end

@testset "Tabulated and nonlinear reservoir level ramps enter both models" begin
    for kind in (:tabulated,:quadratic)
        c=operating_case(;grid=[0.0,0.5,2.0],
            reservoir_kwargs=(level_ramp_down=0.02,),operations=[
                operating_series(:River,:gate_min,[0.0],[0.25]),
                operating_series(:River,:gate_max,[0.0],[0.25])])
        r=only(c.system.reservoirs)
        r=kind==:tabulated ? OpenSHOP._river_replace(r;
            level_curve=TableCurve([4.0,5.0,6.0],[99.0,100.0,102.0])) :
            OpenSHOP._river_replace(r;z0=97.5,slope=0.0,curvature=0.1)
        # The independent witness lowers the level by .02m/hour, at vertices.
        levels=[100.0,99.99,99.96]
        volumes=kind==:tabulated ? 5.0 .+ (levels .-100.0) :
            sqrt.((levels .-97.5)./0.1)
        q=reshape(1.0 .- diff(volumes)./(0.0036 .* diff(c.grid)),1,2)
        c=operating_replace(c;reservoirs=[r],operations=vcat(c.operations,
            [operating_series(:Lake,:inflow,[0.0],[2.0])]))
        u=ones(Int,1,2);x=dispatch_from_controls(c,u,q,fill(0.25,1,2))
        @test x["validation"]["valid"]
        @test head.(Ref(r),vec(x["V"]))≈levels atol=1e-9
        @test replay_audit(c,x)["valid"]
        local_result=solve_case(c;u,warm=x,time_limit=20.0)
        @test get(get(local_result,"validation",Dict()),"valid",false)
        if haskey(local_result,"V")
            @test head.(Ref(r),vec(local_result["V"]))≈levels atol=2e-4
            @test replay_audit(c,local_result)["valid"]
        end
        global_result=solve(c;initial=x,fixed_u=u,time_limit=20.0)
        @test global_result["accepted"]
        @test global_result["global_certificate"]
        if global_result["solution"]!==nothing
            @test head.(Ref(r),vec(global_result["solution"]["V"]))≈levels atol=2e-4
        end
    end
end

@testset "Storage and level rates use adjacent vertices" begin
    # Positive and negative storage changes remain well inside physical bounds.
    for (attribute, delta) in ((:volume_ramp_up,0.006),(:volume_ramp_down,-0.006),
        (:level_ramp_up,0.006),(:level_ramp_down,-0.006))
        level=attribute in (:level_ramp_up,:level_ramp_down)
        rate=level ? 0.001 : 0.01
        rkwargs=NamedTuple{(attribute,)}((rate,))
        c=operating_case(;reservoir_kwargs=rkwargs)
        q=fill(2.0,1,3); dt=diff(c.grid)
        desired=delta .* dt
        inflows=3.0 .+ desired ./ (0.0036 .* dt)
        c=operating_replace(c;operations=[operating_series(:Lake,:inflow,c.grid[1:end-1],inflows)])
        x=dispatch_from_controls(c,ones(Int,1,3),q,fill(0.25,1,3))
        @test x["validation"]["valid"]
        @test vec(diff(x["V"];dims=2))≈desired atol=1e-10
        @test replay_audit(c,x)["valid"]
        bad_rate=level ? 0.0005 : 0.005
        r=OpenSHOP._river_replace(only(c.system.reservoirs);NamedTuple{(attribute,)}((bad_rate,))...)
        cc=operating_replace(c;reservoirs=[r])
        audit=validate(cc,x)
        @test !audit["valid"]
        label=level ? "level_ramp_Lake" : "volume_ramp_Lake"
        @test audit["residuals"][label]≈(abs(delta)-(level ? bad_rate/0.1 : bad_rate))*maximum(dt)*(level ? 0.1 : 1.0) atol=1e-9
        # Restrictions add linear rows, without storage/head state auxiliaries.
        plain=operating_case()
        @test JuMP.num_variables(OpenSHOP._build_dispatch(c).m)==
            JuMP.num_variables(OpenSHOP._build_dispatch(plain).m)
    end
end

@testset "Changing minimum power permits start/stop jumps but never discharge jumps" begin
    alpha=OPERATING_TEST_ALPHA
    c=operating_case(;initial_on=0,generator_kwargs=(ramp_up=0.0,ramp_down=0.0),
        plant_kwargs=(ramp_up=0.0,ramp_down=0.0),operations=[
            operating_series(:Unit1,:pmin,[0.0,0.5,2.0],[0.1,2alpha,0.1])])
    u=reshape([0,1,0],1,3);q=reshape([0.0,2.0,0.0],1,3)
    c=operating_balanced(c,q)
    x=dispatch_from_controls(c,u,q,fill(0.25,1,3))
    @test x["validation"]["valid"]
    @test replay_audit(c,x)["valid"]
    # Current pmin at startup, previous pmin at shutdown: both jumps2alpha.
    @test x["validation"]["residuals"]["unit_ramp_Unit1"]<1e-9
    @test x["validation"]["residuals"]["plant_ramp_Plant1"]<1e-9
    continued=restart_case(c,x,2.0)
    @test only(continued.system.generators).initial_power≈2alpha atol=1e-9
    @test only(continued.system.generators).initial_interval_hours==1.5
    remainder=dispatch_from_controls(continued,u[:,3:end],q[:,3:end],fill(0.25,1,1))
    @test remainder["validation"]["valid"]
    @test remainder["objective"]≈-5.0 atol=1e-9
    @test replay_audit(continued,remainder)["valid"]
    strict=operating_replace(c;generators=[OpenSHOP._river_replace(only(c.system.generators);
        discharge_ramp_up=0.0,discharge_ramp_down=0.0)])
    @test !validate(strict,x)["valid"]
    @test validate(strict,x)["residuals"]["discharge_ramp_Unit1"]≈2.0 atol=1e-9
    local_result=solve_case(c;u,warm=x,time_limit=20.0)
    @test get(get(local_result,"validation",Dict()),"valid",false)
    global_result=solve(c;initial=x,fixed_u=u,time_limit=20.0)
    @test global_result["accepted"]
    @test global_result["global_certificate"]
    # A historical shutdown uses the static preceding minimum, not today's one.
    initial=operating_case(;generator_kwargs=(pmin=alpha,initial_power=alpha,
        ramp_down=0.0),plant_kwargs=(initial_power=alpha,ramp_down=0.0),
        operations=[operating_series(:Unit1,:pmin,[0.0],[2alpha])])
    q=zeros(1,3);initial=operating_balanced(initial,q)
    x=dispatch_from_controls(initial,zeros(Int,1,3),q,fill(0.25,1,3))
    @test x["validation"]["valid"]
    too_large=operating_replace(initial;generators=[OpenSHOP._river_replace(
        only(initial.system.generators);initial_power=2alpha)])
    @test !validate(too_large,x)["valid"]
    @test validate(too_large,x)["residuals"]["unit_ramp_Unit1"]≈alpha atol=1e-9
end

@testset "Optional ramps reject invalid history and avoid unused state variables" begin
    c=operating_case()
    off=operating_case(;initial_on=0)
    for field in (:initial_power,:initial_discharge)
        history=OpenSHOP._river_replace(only(off.system.plants);NamedTuple{(field,)}((1.0,))...)
        @test_throws ArgumentError OpenSHOP.validate_inputs(operating_replace(off;plants=[history]))
    end
    g=only(c.system.generators)
    @test_throws ArgumentError OpenSHOP.validate_inputs(operating_replace(c;
        generators=[OpenSHOP._river_replace(g;ramp_up=-1.0)]))
    @test_throws ErrorException OpenSHOP.validate_inputs(operating_replace(c;
        generators=[OpenSHOP._river_replace(g;initial_interval_hours=0.0)]))
    two=operating_case(;units=2)
    generators=[OpenSHOP._river_replace(z;initial_power=1.0,initial_discharge=2.0,
        initial_interval_hours=Float64(i)) for (i,z) in enumerate(two.system.generators)]
    two=operating_replace(two;generators)
    @test_throws ArgumentError OpenSHOP._build_dispatch(two)
    # Explicit aggregate histories permit different unit observation durations.
    plants=[OpenSHOP._river_replace(only(two.system.plants);initial_power=2.0,
        initial_discharge=4.0,initial_interval_hours=1.5)]
    explicit=operating_replace(two;plants)
    @test OpenSHOP._build_dispatch(explicit).m isa JuMP.Model
    # Unit discharge histories need not share a window when no aggregate
    # discharge ramp consumes them.
    only_discharge=[OpenSHOP._river_replace(z;initial_power=nothing) for z in generators]
    @test OpenSHOP._build_dispatch(operating_replace(two;generators=only_discharge)).m isa JuMP.Model
    no_power_ramp=operating_replace(c;plants=[OpenSHOP._river_replace(
        only(c.system.plants);ramp=nothing)])
    @test case_dict(case_from_dict(case_dict(no_power_ramp)))==case_dict(no_power_ramp)
    @test JuMP.num_variables(OpenSHOP._build_dispatch(no_power_ramp).m)==
        JuMP.num_variables(OpenSHOP._build_dispatch(c).m)
    @test JuMP.num_variables(OpenSHOP._build_global_dispatch(no_power_ramp;joint=true).m)==
        JuMP.num_variables(OpenSHOP._build_global_dispatch(c;joint=true).m)
end

@testset "Confluence restart retains interval release history and exact pulse cohorts" begin
    lakes=[Reservoir(name=name,z0=99.5,slope=0.1,v0=5.0,vmin=4.0,vmax=6.0,
        water_value=0.0) for name in (:TributaryLakeA,:TributaryLakeB)]
    reaches=[
        River(name=:DelayedTributary,source=:TributaryLakeA,target=:Merge,
            law=:controlled,capacity=10.0,curves=RiverRouting.DelayCurve[],
            deterministic_delay=0.5,water_value=0.0,
            history_grid=[-2.0,-1.5,-1.0,0.0],history_release=[0.0,10.0,0.0]),
        River(name=:QuietTributary,source=:TributaryLakeB,target=:Merge,
            law=:controlled,capacity=10.0,curves=RiverRouting.DelayCurve[],
            deterministic_delay=0.0,water_value=0.0,
            history_grid=[-2.0,0.0],history_release=[0.0]),
        # Previous [-1,0]h actual pulses are10 then0: their average is5.
        # A zero ramp applies to control-window averages, not individual pulses.
        River(name=:MergedRiver,source=:Merge,target=:Sea,law=:junction,
            capacity=10.0,curves=RiverRouting.DelayCurve[],
            deterministic_delay=0.0,water_value=0.0,
            history_grid=[-2.0,-1.0,-0.5,0.0],history_release=[0.0,10.0,0.0],
            initial_release=5.0,initial_interval_hours=1.0,
            ramp_up=0.0,ramp_down=0.0),
    ]
    operations=[
        OperationalSeries(object=:TributaryLakeA,attribute=:inflow,
            times=[0.0,1.0],values=[10.0,0.0]),
        OperationalSeries(object=:DelayedTributary,attribute=:gate_min,
            times=[0.0,1.0],values=[1.0,0.0]),
        OperationalSeries(object=:DelayedTributary,attribute=:gate_max,
            times=[0.0,1.0],values=[1.0,0.0]),
        OperationalSeries(object=:QuietTributary,attribute=:gate_min,
            times=[0.0],values=[0.0]),
        OperationalSeries(object=:QuietTributary,attribute=:gate_max,
            times=[0.0],values=[0.0]),
    ]
    system=HydroSystem(reservoirs=lakes,junctions=Junction[],
        boundaries=[Boundary(name=:Sea,head=0.0)],tunnels=Tunnel[],
        plants=Plant[],generators=Generator[],
        river_junctions=[RiverJunction(name=:Merge)],rivers=reaches)
    c=ScheduleCase(name="river average versus pulse restart",system=system,
        grid=[0.0,1.0,2.0],prices=zeros(2),operations=operations)
    u=zeros(Int,0,2);q=zeros(0,2);gate=[1.0 0.0;0.0 0.0;0.0 0.0]
    @test OpenSHOP.validate_inputs(c)
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    x=dispatch_from_controls(c,u,q,gate)
    @test x["validation"]["valid"]
    @test x["V"]≈fill(5.0,2,3) atol=1e-10
    @test x["river_release"]≈[10.0 0.0;0.0 0.0;5.0 5.0] atol=1e-10
    @test replay_audit(c,x)["valid"]

    fine=OpenSHOP.without_interval_controls(with_grid(c,[0.0,0.5,1.0,1.5,2.0]))
    fine_x=dispatch_from_controls(fine,zeros(Int,0,4),zeros(0,4),gate[:,[1,1,2,2]])
    @test fine_x["validation"]["valid"]
    @test vec(fine_x["river_release"][3,:])≈[0.0,10.0,10.0,0.0] atol=1e-10

    continued=restart_case(c,x,1.0)
    merged=continued.system.rivers[3]
    @test merged.initial_release≈5.0 atol=1e-10
    @test merged.initial_interval_hours==1.0
    @test last(merged.history_release)≈10.0 atol=1e-10
    @test merged.history_grid[end]-merged.history_grid[end-1]==0.5
    remainder=dispatch_from_controls(continued,zeros(Int,0,1),zeros(0,1),gate[:,2:end])
    @test remainder["validation"]["valid"]
    @test remainder["river_release"]≈x["river_release"][:,2:end] atol=1e-10
    @test remainder["terminal_transit"]≈x["terminal_transit"] atol=1e-10
    @test replay_audit(continued,remainder)["valid"]
    @test case_dict(case_from_dict(case_dict(continued)))==case_dict(continued)

    # Reproduce the former restart behavior while leaving routing cohorts intact:
    # it incorrectly compares the next average5 to the final pulse10.
    old_rivers=copy(continued.system.rivers)
    old_rivers[3]=OpenSHOP._river_replace(merged;initial_release=nothing)
    old_system=OpenSHOP._river_replace(continued.system;rivers=old_rivers)
    old_case=OpenSHOP._river_replace(continued;system=old_system)
    rejected=validate(old_case,remainder)
    @test !rejected["valid"]
    @test rejected["residuals"]["river_ramp_MergedRiver"]≈5.0 atol=1e-10

    local_result=solve_case(continued;u=zeros(Int,0,1),warm=remainder,time_limit=45.0)
    @test get(get(local_result,"validation",Dict()),"valid",false)
    if haskey(local_result,"river_release")
        @test local_result["river_release"]≈remainder["river_release"] atol=1e-7
        @test replay_audit(continued,local_result)["valid"]
    end
    global_result=solve(continued;initial=remainder,time_limit=30.0)
    @test global_result["accepted"]
    @test global_result["global_certificate"]
end

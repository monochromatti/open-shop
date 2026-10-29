function warm_dispatch_fixture(; rivers = true, minimum_power = 0.1)
    lake=Reservoir(name = :WarmLake, z0 = 100.0, slope = 1.0, curvature = 0.05,
        v0 = 2.0, vmin = 0.1, vmax = 4.0, inflow = 8.0, water_value = 10.0)
    generators=[Generator(name = name, plant = :WarmPlant, qmin = 1.0, qmax = 5.0,
        pmin = minimum_power, pmax = 5.0, efficiency = 0.9, qbest = 3.0,
        qcurvature = 0.0, hbest = 100.0, hcurvature = 0.0, hmin = 70.0,
        hmax = 130.0, minup = 0.0, mindown = 0.0, startup = 0.0) for name in (:A, :B)]
    reaches=rivers ? [River(name = :WarmReach, source = :WarmLake, target = :WarmSea,
        curves = RiverRouting.DelayCurve[], deterministic_delay = 0.25,
        capacity = 5.0, law = :orifice, coefficient = 0.1, water_value = 10.0)] : River[]
    system=HydroSystem(reservoirs = [lake], junctions = [Junction(name = :WarmIntake)],
        boundaries = [Boundary(name = :WarmSea, head = 0.0)],
        tunnels = [Tunnel(name = :WarmTunnel, source = :WarmLake, target = :WarmIntake,
            resistance = 0.2, capacity = 20.0)],
        plants = [Plant(name = :WarmPlant, source = :WarmIntake, target = :WarmSea,
            pmax = 10.0, tailwater_curve = TableCurve([0.0, 5.0, 10.0], [0.0, 0.5, 2.0]))],
        generators = generators, rivers = reaches)
    operations=[OperationalSeries(object = :B, attribute = :qmin, times = [0.0, 1.0],
        values = [2.5, 1.5]), OperationalSeries(object = :B, attribute = :qmax,
        times = [0.0, 1.0], values = [3.0, 2.0])]
    ScheduleCase(name = "changed commitment hydraulics", system = system, grid = [0.0, 1.0, 2.0],
        prices = [50.0, 90.0], operations = operations)
end

@testset "Preparation improves an accepted schedule and preserves its feasibility" begin
    c=warm_dispatch_fixture(;rivers=false)
    off=dispatch_from_controls(c,zeros(Int,2,2),zeros(2,2),zeros(0,2))
    @test replay_audit(c,off)["valid"]
    prepared=schedule_case(c;initial=off,time_limit=15.0,proposal_time_limit=2.0,
        nlp_time_limit=5.0,max_refinements=0,refine_grid=false)
    @test prepared["accepted"]
    if prepared["accepted"]
        @test prepared["solution"]["objective"]>off.objective+1.0
        @test validate(c,prepared["solution"])["valid"]
        @test replay_audit(c,prepared["solution"])["valid"]
    end
    damaged=OpenSHOP.schedule_dict(off)
    damaged["power"][1]=1.0
    @test_throws ArgumentError schedule_case(c;initial=damaged,time_limit=15.0)
end

@testset "Shared preparation deadline" begin
    c=warm_dispatch_fixture()
    u=zeros(Int,2,2)
    expired=solve_verified(c;u,deadline=time()-1,max_refinements=2)
    @test !expired["accepted"]
    @test expired["solution"]===nothing
    @test expired["case"].grid==c.grid

    # Auditing a supplied schedule remains necessary even when no optimization
    # time remains. A deadline must never turn stale acceptance into evidence.
    x=dispatch_from_controls(c,u,zeros(2,2),ones(1,2))
    audited=solve_verified(c;u,initial=x,deadline=time()-1,max_refinements=0)
    @test audited["accepted"]
    damaged=OpenSHOP.schedule_dict(x)
    damaged["power"][1]=1.0
    rejected=solve_verified(c;u,initial=damaged,deadline=time()-1,max_refinements=0)
    @test !rejected["accepted"]

    prepared=schedule_case(c;time_limit=1e-9)
    @test !prepared["accepted"]
    @test prepared["solution"]===nothing
    for invalid in (0.0,-1.0,NaN)
        @test_throws ArgumentError schedule_case(c;time_limit=invalid)
    end
    @test_throws ArgumentError solve_verified(c;deadline=NaN)
end

@testset "Forecast updates use numerical starts without retaining stale feasibility" begin
    c=warm_dispatch_fixture(;rivers=false)
    x=dispatch_from_controls(c,ones(Int,2,2),[3.0 3.0;2.8 1.8],zeros(0,2))
    lake=OpenSHOP._river_replace(only(c.system.reservoirs);inflow=12.0)
    updated=OpenSHOP._river_replace(c;prices=[100.0,20.0],
        system=OpenSHOP._river_replace(c.system;reservoirs=[lake]),
        operations=vcat(c.operations,[OperationalSeries(object=:A,attribute=:maintenance,
            times=[0.0],values=[1.0])]))
    @test !validate(updated,x)["valid"]
    prepared=schedule_case(updated;warm=x,time_limit=15.0,
        proposal_time_limit=2.0,nlp_time_limit=5.0,max_refinements=0)
    @test prepared["accepted"]
    if prepared["accepted"]
        y=prepared["solution"]
        @test validate(updated,y)["valid"]
        @test replay_audit(updated,y)["valid"]
        @test objective_components(updated,y).total≈y.objective
        @test all(iszero,y.commitment[1,:])
    end
    expired=schedule_case(updated;warm=x,time_limit=1e-9)
    @test !expired["accepted"]
    @test expired["solution"]===nothing
    wrong_grid=schedule_dict(x)
    wrong_grid["grid"]=[0.0,0.5,2.0]
    @test_throws ArgumentError schedule_case(updated;warm=wrong_grid)
end

@testset "Finer hydraulic constraints preserve the requested decision grid" begin
    base=warm_dispatch_fixture(;rivers=false)
    lake=OpenSHOP._river_replace(only(base.system.reservoirs);inflow=80.0)
    units=[OpenSHOP._river_replace(g;pmax=2.0) for g in base.system.generators]
    system=OpenSHOP._river_replace(base.system;reservoirs=[lake],generators=units)
    c=OpenSHOP._river_replace(base;system)
    u=[1 1;0 0]
    verified=solve_verified(c;u,max_refinements=1,refine_grid=false,time_limit=10.0)
    @test verified["accepted"]
    @test verified["case"].grid==c.grid
    @test verified["solution"]["u"]==u
    @test replay_audit(c,verified["solution"])["valid"]
end

@testset "Finer replay measures the power shortfall in physical units" begin
    lake=Reservoir(name=:Lake,z0=100.0,slope=1.0,v0=2.0,vmin=1.0,vmax=4.0,
        inflow=80.0,water_value=0.0)
    # One Mm³ raises head by one metre. Storage and power rise linearly, so
    # the last half-hour's mean power exceeds the hour's mean by gain*rise/4.
    # No optimizer is needed to know that discrepancy.
    q=2.0
    gain=0.00981*0.9*q
    rise=0.0036*(lake.inflow-q)
    capacity=gain*(lake.z0+lake.v0+rise/2)
    unit=Generator(name=:Unit,plant=:Plant,qmin=0.0,qmax=5.0,pmin=0.0,pmax=capacity,
        efficiency=0.9,qbest=3.0,qcurvature=0.0,hbest=100.0,hcurvature=0.0,
        hmin=90.0,hmax=120.0,minup=0.0,mindown=0.0,startup=0.0)
    system=HydroSystem(reservoirs=[lake],generators=[unit],
        plants=[Plant(name=:Plant,source=:Lake,target=:Sea,pmax=10.0)],
        boundaries=[Boundary(name=:Sea,head=0.0)])
    c=ScheduleCase(name="known replay shortfall",system=system,grid=[0.0,1.0],prices=[50.0])
    x=dispatch_from_controls(c,ones(Int,1,1),fill(q,1,1),zeros(0,1))
    @test validate(c,x)["valid"]
    check=replay_audit(c,x)
    @test !check["valid"]
    @test check["shortfalls"].power≈gain*rise/4 atol=1e-10
    @test check["shortfalls"].arrival==0.0
    @test !replay_audit(c,x;grid=[0.0,0.5,1.5])["valid"]
    @test_throws MethodError replay_audit(c,x;grid=Any[0.0,"midpoint",1.0])
    verified=solve_verified(c;u=x.commitment,warm=x,time_limit=10.0,
        max_refinements=1,refine_grid=false)
    @test verified["accepted"]
    @test verified["case"].grid==c.grid
    @test replay_audit(c,verified["solution"])["valid"]
end

@testset "Local reconstruction corrects controls and recomputes physical values" begin
    c=warm_dispatch_fixture()
    u=zeros(Int, 2, 2)
    physical=dispatch_from_controls(c, u, zeros(2, 2), ones(1, 2))
    raw=OpenSHOP.schedule_dict(physical)
    raw["status"]="TIME_LIMIT"
    raw["generator_q"].=1e-7
    raw["gate"].+=1e-14
    raw["objective"]=123456.0
    raw["power"].=100.0
    raw["validation"]=validate(c, raw)
    @test !raw["validation"]["valid"]
    original=deepcopy(raw)
    rebuilt, correction=OpenSHOP._reconstruct_candidate(c, raw)
    @test raw==original
    @test correction.flow==1e-7
    @test correction.gate>0
    @test rebuilt["validation"]["valid"]
    report=Dict{String,Any}("status"=>raw["status"])
    raw=OpenSHOP._repair_dispatch!(c, OpenSHOP._schedule_input(raw;case=c), report)
    @test report["forward_reconstructed"]
    @test report["status"]=="TIME_LIMIT"
    @test report["raw_solver_validation"]==original["validation"]
    @test report["raw_solver_objective"]==123456.0
    @test raw["generator_q"]==zeros(2, 2)
    @test raw["gate"]==ones(1, 2)
    @test raw["power"]==zeros(2, 2)
    @test raw["V"]≈physical["V"] atol=1e-10
    @test raw["objective"]≈physical["objective"] atol=1e-8
    @test replay_audit(c, raw)["valid"]
    for (key, value) in (("generator_q", 1e-3), ("generator_q", -1e-3),
            ("gate", 1.0001), ("gate", NaN))
        invalid=deepcopy(original)
        invalid[key][1]=value
        @test_throws ArgumentError OpenSHOP._reconstruct_candidate(c, invalid)
        report=Dict{String,Any}()
        invalid=OpenSHOP._repair_dispatch!(c, OpenSHOP._schedule_input(invalid;case=c), report)
        @test !get(report, "forward_reconstructed", false)
        @test haskey(report, "forward_reconstruction_error")
        @test invalid["objective"]==original["objective"]
    end
end

@testset "Warm starts require compatible grids and finite arrays" begin
    for rivers in (false, true)
        c=warm_dispatch_fixture(; rivers)
        u=[1 1; 0 0]
        initial=dispatch_from_controls(c, u, [4.0 4.0; 0.0 0.0],
            fill(0.8, length(c.system.rivers), 2))
        compatible, reason=OpenSHOP._dispatch_start(c, u, initial)
        @test compatible.generator_q==initial.generator_q
        @test compatible.objective==initial.objective
        @test reason["used"]
        for grid in ([0.0, 0.5, 2.0], [0.0, 0.5, 1.0, 1.5, 2.0])
            wrong=OpenSHOP.schedule_dict(initial)
            wrong["grid"]=grid
            discarded, reason=OpenSHOP._dispatch_start(c, u, wrong)
            @test discarded===nothing
            @test occursin("grid", reason["ignored_reason"])
        end
        legacy=OpenSHOP.schedule_dict(initial)
        delete!(legacy, "grid")
        @test first(OpenSHOP._dispatch_start(c,u,legacy)) isa Schedule
        for key in ("generator_q", "power", "V", "gate")
            for damage in (:shape, :missing)
                wrong=deepcopy(legacy)
                damage==:shape ?
                    (wrong[key]=zeros(size(wrong[key], 1), size(wrong[key], 2)+1)) :
                    delete!(wrong, key)
                discarded, reason=OpenSHOP._dispatch_start(c, u, wrong)
                @test discarded===nothing
                @test !reason["used"]
            end
        end
        wrong=OpenSHOP.schedule_dict(initial)
        wrong["H"][1]=NaN
        @test first(OpenSHOP._dispatch_start(c, u, wrong))===nothing
        @test first(OpenSHOP._dispatch_start(c, u, nothing))===nothing
    end
end

@testset "Passive river semantics survive guarded reconstruction" begin
    c=warm_dispatch_fixture()
    weir=OpenSHOP._river_replace(only(c.system.rivers); law = :weir, coefficient = 0.001)
    c=OpenSHOP._river_replace(c; system = OpenSHOP._river_replace(c.system; rivers = [weir]))
    x=dispatch_from_controls(c, zeros(Int, 2, 2), zeros(2, 2), zeros(1, 2))
    @test !x["validation"]["valid"]
    report=Dict{String,Any}()
    repaired=OpenSHOP._repair_dispatch!(c, x, report)
    @test repaired.generator_q==x.generator_q
    @test repaired.objective==x.objective
    @test !get(report, "forward_reconstructed", false)
    @test !report["forward_reconstruction_audit"]["valid"]
    @test x["gate"]==zeros(1, 2)
end

@testset "Verified dispatch uses a fresh original audit" begin
    c=warm_dispatch_fixture()
    u=zeros(Int, 2, 2)
    physical=dispatch_from_controls(c, u, zeros(2, 2), ones(1, 2))
    @test physical["validation"]["valid"]
    for key in ("power", "V")
        bad=OpenSHOP.schedule_dict(physical)
        bad[key][end]+=1.0
        @test bad["validation"]["valid"]
        result=solve_verified(c; u, initial=bad, max_refinements=0)
        @test !result["accepted"]
        @test !result["attempts"][1]["model_valid"]
        @test !isempty(result["attempts"][1]["model_errors"])
        @test !haskey(result["audit"], "replay_grid")
    end
    for key in ("power", "gate")
        bad=OpenSHOP.schedule_dict(physical)
        bad[key]=zeros(1, 3)
        result=solve_verified(c; u, initial=bad, max_refinements=0)
        @test !result["accepted"]
        @test "missing or malformed $key" in result["attempts"][1]["model_errors"]
    end
    valid=OpenSHOP.schedule_dict(physical)
    valid["validation"]=Dict("valid"=>false, "errors"=>["stale diagnostic"])
    valid["raw_solver_validation"]=Dict("valid"=>false)
    result=solve_verified(c; u, initial=valid, max_refinements=0)
    @test result["accepted"]
    @test result["audit"]["original_audit"]["valid"]
    @test isempty(result["attempts"][1]["model_errors"])
    @test haskey(result["audit"], "replay_grid")
    wrong_grid=OpenSHOP.schedule_dict(physical)
    wrong_grid["grid"]=[0.0, 0.5, 2.0]
    @test_throws ArgumentError solve_verified(c; u, initial=wrong_grid, max_refinements=0)
end

@testset "Preparation finds and improves a cold physical schedule" begin
    c=warm_dispatch_fixture()
    prepared=schedule_case(c;time_limit=20.0,proposal_time_limit=3.0,
        nlp_time_limit=5.0,max_refinements=0,refine_grid=false)
    @test prepared["accepted"]
    if prepared["accepted"]
        @test validate(c,prepared["solution"])["valid"]
        @test replay_audit(c,prepared["solution"])["valid"]
        @test prepared["case"]["grid"]==c.grid
        @test prepared["solution"]["objective"]>0.0
    end
end

@testset "Physical proposal domains enclose curved heads and signed flows" begin
    lakes=[Reservoir(name=name,z0=z,slope=0.5,curvature=0.02,v0=2.0,
        vmin=1.0,vmax=3.0,water_value=0.0) for (name,z) in ((:A,100.0),(:B,99.0))]
    link=Tunnel(name=:Link,source=:A,target=:B,resistance=0.5,capacity=10.0,
        gate_positions=[0.0,0.5,1.0],
        gate_loss_factor_curve=TableCurve([0.0,0.5,1.0],[2.0,1.0,0.0]))
    delta=LevelFlowCurve([99.0,103.0],[TableCurve([-4.0,0.0,4.0],[-2.0,0.0,3.0]),
        TableCurve([-4.0,0.0,4.0],[-1.0,0.0,4.0])])
    channel=River(name=:Channel,source=:A,target=:B,law=:delta,capacity=10.0,
        deterministic_delay=0.0,curves=RiverRouting.DelayCurve[],water_value=0.0,
        delta_level_ref_up_flow_curve=delta)
    c=ScheduleCase(name="physical domains",system=HydroSystem(reservoirs=lakes,
        tunnels=[link],rivers=[channel]),grid=[0.0,1.0],prices=[0.0])
    domains=OpenSHOP._proposal_flow_bounds(c)
    for a in range(1.0,3.0;length=9),b in range(1.0,3.0;length=9),p in link.gate_positions
        ha,hb=head(lakes[1],a),head(lakes[2],b)
        q=OpenSHOP.tunnel_flow(link,ha-hb,p)
        r=OpenSHOP.channel_flow(channel,ha,hb)
        @test domains.tunnel_lower[1]<=q<=domains.tunnel_upper[1]
        @test domains.river_lower[1]<=r<=domains.river_upper[1]
    end
end

@testset "A joint neighborhood counts switches and categorical gates once" begin
    base=warm_dispatch_fixture()
    tunnel=OpenSHOP._river_replace(only(base.system.tunnels);gate_positions=[0.0,0.5,1.0])
    river=OpenSHOP._river_replace(only(base.system.rivers);law=:controlled,min_allowed_flow=1.0)
    c=OpenSHOP._river_replace(base;operations=OperationalSeries[],
        system=OpenSHOP._river_replace(base.system;tunnels=[tunnel],rivers=[river]))
    x=dispatch_from_controls(c,ones(Int,2,2),fill(2.0,2,2),fill(0.5,1,2);
        tunnel_opening=ones(1,2))
    @test replay_audit(c,x)["valid"]
    for radius in (0,1,3)
        p=OpenSHOP.propose_commitment(c;reference=x,incumbent=x,radius,time_limit=5.0)
        @test haskey(p,"u")
        if haskey(p,"u")
            changes=sum(p["u"].!=x.commitment)+sum(p["tunnel_opening"].!=x.tunnel_opening)+
                sum(p["river_on"].!=OpenSHOP.river_on_schedule(c,x.river_release))
            @test changes<=radius
            @test admissible(c,p["u"])
            @test validate_tunnel_openings(c,p["tunnel_opening"],p["u"])
        end
    end
    @test_throws ArgumentError OpenSHOP.propose_commitment(c;radius=1)
    @test_throws ArgumentError OpenSHOP.propose_commitment(c;incumbent=x,radius=-1)
end

function fixed_power_replay_case()
    lake=Reservoir(name=:Lake,z0=98.0,slope=1.0,v0=2.0,vmin=1.0,vmax=4.0,
        inflow=80.0,water_value=10.0)
    output=0.00981*0.9*200
    unit=Generator(name=:Unit,plant=:Plant,qmin=0.0,qmax=5.0,pmin=output,pmax=output,
        efficiency=0.9,qbest=3.0,qcurvature=0.0,hbest=100.0,hcurvature=0.0,
        hmin=90.0,hmax=120.0,minup=0.0,mindown=0.0,startup=0.0)
    outlet=River(name=:Outlet,source=:Lake,target=:Sea,law=:controlled,capacity=80.0,
        deterministic_delay=0.0,curves=RiverRouting.DelayCurve[],water_value=0.0)
    system=HydroSystem(reservoirs=[lake],generators=[unit],rivers=[outlet],
        plants=[Plant(name=:Plant,source=:Lake,target=:Sea,pmax=10.0)],
        boundaries=[Boundary(name=:Sea,head=0.0)])
    ScheduleCase(name="fixed power with changing storage",system=system,
        grid=[0.0,1.0],prices=[50.0])
end

@testset "Finer equations recover fixed power without contradictory margins" begin
    c=fixed_power_replay_case()
    u=ones(Int,1,1)
    witness=dispatch_from_controls(c,u,fill(2.0,1,1),fill(78/80,1,1))
    @test replay_audit(c,witness)["valid"]
    b=100+0.0018*80
    q=400/(b+sqrt(b^2-4*0.0018*200))
    coarse=dispatch_from_controls(c,u,fill(q,1,1),zeros(1,1))
    @test validate(c,coarse)["valid"]
    @test !replay_audit(c,coarse)["valid"]
    for refine_grid in (false,true)
        result=solve_verified(c;u,initial=coarse,max_refinements=1,refine_grid,time_limit=15.0)
        @test result["accepted"]
        @test result["case"].grid==c.grid
        @test result["solution"].grid==c.grid
        @test result["solution"].generator_q≈fill(2.0,1,1) atol=1e-5
        @test result["solution"].gate≈fill(78/80,1,1) atol=1e-4
        @test replay_audit(c,result["solution"])["valid"]
        # A nonzero explicitly requested inset also preserves a fixed value.
        inset=solve_case(c;u,warm=witness,operational_margin=0.1,time_limit=10.0)
        @test inset.accepted
        @test inset.solution.power≈witness.power atol=1e-5
    end
end

@testset "Finer replay keeps discharge and power errors in their own units" begin
    lower=Reservoir(name=:Lower,z0=0.0,slope=1.0,v0=2.0,vmin=1.0,vmax=4.0,water_value=0.0)
    upper=Reservoir(name=:Upper,z0=10.0,slope=1.0,v0=2.0,vmin=1.0,vmax=4.0,water_value=0.0)
    characteristic=PumpCharacteristic([5.0,15.0],[5.0,1.0],[0.9,0.9])
    unit=Pump(name=:Pump,plant=:Lift,characteristic=characteristic,qmin=0.0,qmax=10.0,
        pmin=0.0,pmax=10.0,hmin=5.0,hmax=15.0,qbest=3.0,hbest=10.0,
        minup=0.0,mindown=0.0,startup=0.0)
    q=3/(1+0.4*0.0036)
    system=HydroSystem(reservoirs=[lower,upper],generators=[unit],
        plants=[Plant(name=:Lift,source=:Lower,target=:Upper,pmax=10.0,qmax=q)])
    c=ScheduleCase(name="discharge error only",system=system,grid=[0.0,1.0],prices=[0.0])
    x=dispatch_from_controls(c,ones(Int,1,1),fill(q,1,1),zeros(0,1))
    @test validate(c,x)["valid"]
    check=replay_audit(c,x)
    @test !check["valid"]
    @test check["shortfalls"].power==0.0
    @test check["shortfalls"].discharge≈0.0021553422261328947 atol=1e-9
end

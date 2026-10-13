function delta_channel_case(; direction=1.0, reference=:up, gain=4.0)
    upstream=Reservoir(name=:UpperChannel,z0=100.0,slope=1.0,v0=3.0+direction,
        vmin=1.0,vmax=6.0,water_value=10.0)
    downstream=Reservoir(name=:LowerChannel,z0=100.0,slope=1.0,v0=3.0,
        vmin=1.0,vmax=6.0,water_value=10.0)
    table=LevelFlowCurve([100.0],[TableCurve([-10.0,0.0,10.0],[-10gain,0.0,10gain])])
    kwargs=reference==:up ? (;delta_level_ref_up_flow_curve=table) :
                          (;delta_level_ref_down_flow_curve=table)
    river=River(name=:LevelChannel,source=:UpperChannel,target=:LowerChannel,
        law=:delta,curves=RiverRouting.DelayCurve[],deterministic_delay=0.0,
        capacity=50.0,water_value=0.0;kwargs...)
    system=HydroSystem(reservoirs=[upstream,downstream],junctions=Junction[],
        boundaries=Boundary[],tunnels=Tunnel[],plants=Plant[],generators=Generator[],rivers=[river])
    ScheduleCase(name="signed level channel",system=system,grid=[0.0,0.5,1.0,2.0],prices=zeros(3))
end

# Independent solution of the two-reservoir linear conductance midpoint equation.
function channel_linear_trajectory(c,gain)
    volume=[r.v0 for r in c.system.reservoirs]
    storage=reshape(copy(volume),:,1)
    release=Float64[]
    for dt in diff(c.grid)
        q=gain*(volume[1]-volume[2])/(1+0.0036dt*gain)
        volume+=[-0.0036dt*q,0.0036dt*q]
        storage=hcat(storage,volume)
        push!(release,q)
    end
    storage,reshape(release,1,:)
end

@testset "Signed delta channels conserve and restart in both directions" begin
    for direction in (-1.0,0.0,1.0), reference in (:up,:down)
        c=delta_channel_case(;direction,reference)
        @test OpenSHOP.validate_inputs(c)
        @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
        expected,flow=channel_linear_trajectory(c,4.0)
        x=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(1,3))
        @test x["validation"]["valid"]
        @test x["V"]≈expected atol=1e-10
        @test x["river_release"]≈flow atol=1e-9
        @test x["arrival_volume"]≈flow .* reshape(0.0036 .* diff(c.grid),1,:) atol=1e-12
        @test maximum(abs,x["terminal_transit"])<1e-12
        @test vec(sum(x["V"];dims=1))≈fill(6.0+direction,4) atol=1e-10
        @test replay_audit(c,x)["valid"]
        @test transport_audit(c,x)["accepted"]
        continued=restart_case(c,x,1.0)
        rest=dispatch_from_controls(continued,zeros(Int,0,1),zeros(0,1),zeros(1,1))
        @test rest["validation"]["valid"]
        @test rest["V"]≈x["V"][:,3:end] atol=1e-10
        @test rest["river_release"]≈flow[:,3:end] atol=1e-9
        if direction!=0
            broken=OpenSHOP.schedule_dict(x)
            broken["river_release"].=-broken["river_release"]
            @test !validate(c,broken)["valid"]
        end
        # Public routing/query APIs also preserve signed instantaneous transfer.
        routed=route_network(c,x["river_release"])
        @test routed["arrival_volume"]≈x["arrival_volume"] atol=1e-12
        @test point_arrival(only(c.system.rivers),c.grid,vec(flow),0.25)≈flow[1,1]
    end
end

@testset "Delta channels have no downstream transport or implicit gate" begin
    c=delta_channel_case()
    r=only(c.system.rivers)
    bad_case(river;operations=c.operations)=OpenSHOP._river_replace(c;
        system=OpenSHOP._river_replace(c.system;rivers=[river]),operations)
    for river in (
        OpenSHOP._river_replace(r;law=:controlled),
        OpenSHOP._river_replace(r;deterministic_delay=0.1),
        OpenSHOP._river_replace(r;inflow=1.0),
        OpenSHOP._river_replace(r;min_arrival=1.0),
        OpenSHOP._river_replace(r;gate_min=0.5),
        OpenSHOP._river_replace(r;target=:UpperChannel),
        OpenSHOP._river_replace(r;discharge_curve=TableCurve([100.0,110.0],[0.0,10.0])),
        OpenSHOP._river_replace(r;delta_level_ref_down_flow_curve=r.delta_level_ref_up_flow_curve),
        OpenSHOP._river_replace(r;delta_level_ref_up_flow_curve=nothing))
        @test_throws ArgumentError OpenSHOP.validate_inputs(bad_case(river))
    end
    op=OperationalSeries(object=r.name,attribute=:inflow,times=[0.0],values=[1.0])
    @test_throws ArgumentError OpenSHOP.validate_inputs(bad_case(r;operations=[op]))
end

@testset "Local and global delta equations retain reverse schedules" begin
    for direction in (-1.0,1.0), reference in (:up,:down)
        c=delta_channel_case(;direction,reference)
        expected,flow=channel_linear_trajectory(c,4.0)
        initial=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(1,3))
        local_result=solve_case(c;warm=initial,time_limit=10.0)
        @test local_result["validation"]["valid"]
        @test local_result["river_release"]≈flow atol=1e-6
        b=OpenSHOP._build_global_dispatch(c;warm=initial)
        @test OpenSHOP._lift_start!(b,c,initial)["valid"]
        JuMP.set_optimizer(b.m,OpenSHOP.SCIP.Optimizer)
        JuMP.set_silent(b.m)
        JuMP.set_optimizer_attribute(b.m,"limits/time",15.0)
        JuMP.set_optimizer_attribute(b.m,"numerics/feastol",1e-8)
        JuMP.optimize!(b.m)
        @test JuMP.has_values(b.m)
        if JuMP.has_values(b.m)
            raw=OpenSHOP._dispatch_values(b)
            @test raw.commitment==zeros(Int,0,3)
            @test validate(c,raw)["valid"]
            @test raw["river_release"]≈flow atol=1e-5
            @test raw["V"]≈expected atol=1e-7
        end
        domains=b.domains
        @test all(domains.lower .<= initial["V"] .<= domains.upper)
        proposal=OpenSHOP.propose_commitment(c;reference=initial,time_limit=5.0)
        @test proposal["status"] in ("OPTIMAL","TIME_LIMIT")
    end
end

@testset "Absolute-level references change the physical channel flow" begin
    for reference in (:up,:down), direction in (-1.0,1.0)
        c=delta_channel_case(;direction,reference)
        original=only(c.system.rivers)
        family=LevelFlowCurve([100.0,105.0,110.0],[
            TableCurve([-10.0,0.0,10.0],[-10.0,0.0,20.0]),
            TableCurve([-10.0,0.0,10.0],[-30.0,0.0,40.0]),
            TableCurve([-10.0,0.0,10.0],[-50.0,0.0,60.0])])
        changes=reference==:up ? (;delta_level_ref_up_flow_curve=family) :
                                (;delta_level_ref_down_flow_curve=family)
        river=OpenSHOP._river_replace(original;changes...)
        c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;rivers=[river]))
        x=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(1,3))
        @test x["validation"]["valid"]
        volumes=[r.v0 for r in c.system.reservoirs]
        for (t,dt) in enumerate(diff(c.grid))
            alpha=0.0036dt
            residual(q)=begin
                ref=100.0+(reference==:up ? volumes[1]-alpha*q/2 : volumes[2]+alpha*q/2)
                conductance=(direction>0 ? 2.0 : 1.0)+0.4*(ref-100.0)
                q-conductance*(volumes[1]-volumes[2]-alpha*q)
            end
            lo,hi=-50.0,50.0
            for _ in 1:80
                mid=(lo+hi)/2
                residual(mid)>0 ? (hi=mid) : (lo=mid)
            end
            q=(lo+hi)/2
            @test x["river_release"][1,t]≈q atol=1e-9
            volumes+=[-alpha*q,alpha*q]
            @test x["V"][:,t+1]≈volumes atol=1e-10
        end
        b=OpenSHOP._build_global_dispatch(c;warm=x)
        @test OpenSHOP._lift_start!(b,c,x)["valid"]
        local_result=solve_case(c;warm=x,time_limit=10.0)
        @test local_result["validation"]["valid"]
        @test local_result["river_release"]≈x["river_release"] atol=1e-6
    end
end

function channel_network_case(; distributed=false)
    lake(name,z,slope,v)=Reservoir(;name,z0=z,slope,v0=v,vmin=1.0,vmax=6.0,water_value=0.0)
    reservoirs=[lake(:LakeA,100.0,10.0,2.0),lake(:LakeB,100.0,10.0,2.05),
        lake(:LakeC,100.0,10.0,2.02),lake(:LakeD,100.0,5.0,4.0),
        lake(:LowerLake,60.0,2.0,3.0),lake(:TributaryLake,130.0,2.0,2.0)]
    response=LevelFlowCurve([100.0],[TableCurve([-10.0,0.0,10.0],[-30.0,0.0,40.0])])
    channel(name,source,target)=River(;name,source,target,law=:delta,
        delta_level_ref_up_flow_curve=response,curves=RiverRouting.DelayCurve[],
        deterministic_delay=0.0,capacity=50.0,water_value=0.0)
    reaches=[channel(:AB,:LakeA,:LakeB),channel(:BC,:LakeB,:LakeC),channel(:CA,:LakeC,:LakeA),
        River(name=:Tributary,source=:TributaryLake,target=:CombinedRiver,law=:controlled,
            curves=RiverRouting.DelayCurve[],deterministic_delay=0.25,capacity=8.0,water_value=0.0),
        River(name=:CombinedRiver,target=:LowerLake,law=:junction,inflow=2.0,
            curves=distributed ? [RiverRouting.DelayCurve(0.0,[0.1,0.5,1.0],[0.6,0.4])] : RiverRouting.DelayCurve[],
            deterministic_delay=distributed ? nothing : 0.5,capacity=25.0,water_value=0.0)]
    plants=[Plant(name=:StationB,source=:LakeB,target=:LowerLake,pmax=10.0,
        discharge_river=:CombinedRiver),Plant(name=:StationD,source=:LakeD,target=:LowerLake,
        discharge_river=:CombinedRiver,pmax=10.0)]
    generators=[Generator(name=Symbol("Unit",k),plant=p.name,qmin=0.1,qmax=5.0,
        pmin=0.01,pmax=10.0,hmin=20.0,hmax=80.0,qbest=3.0,hbest=50.0,
        efficiency=0.9,qcurvature=0.0,hcurvature=0.0,minup=0.0,mindown=0.0,
        startup=0.0,initial_on=1) for (k,p) in enumerate(plants)]
    tunnel=Tunnel(name=:ACPressureLink,source=:LakeA,target=:LakeC,resistance=0.5,
        capacity=50.0,opening=0.5,gate_loss_factor_curve=TableCurve([0.0,1.0],[2.0,0.0]))
    system=normalize_river_connections(HydroSystem(;reservoirs,junctions=Junction[],
        boundaries=Boundary[],tunnels=[tunnel],plants,generators,rivers=reaches))
    grid=[0.0,0.25,0.5,1.0,1.5,2.0]
    operations=[OperationalSeries(object=:LakeA,attribute=:inflow,times=[0.5],values=[30.0]),
        OperationalSeries(object=:ACPressureLink,attribute=:opening,times=[1.0],values=[0.0]),
        OperationalSeries(object=:AB,attribute=:capacity,times=[1.0],values=[40.0])]
    ScheduleCase(name="lake channels and delayed confluence",system=system,grid=grid,
        prices=[40.0,50.0,100.0,80.0,30.0],operations=operations)
end

@testset "Hydraulic channel loops with pressure links and delayed confluences" begin
    for distributed in (false,true)
        c=channel_network_case(;distributed)
        @test OpenSHOP.validate_inputs(c)
        x=dispatch_from_controls(c,ones(Int,2,5),fill(1.0,2,5),
            [zeros(3,5);fill(0.5,1,5);zeros(1,5)])
        @test x["validation"]["valid"]
        @test x["river_release"][1,1]<0
        @test x["river_release"][1,end]>0
        @test maximum(abs,x["terminal_transit"][1:3])<1e-12
        @test x["tunnel_q"][1,4:5]≈zeros(2) atol=1e-9
        @test x["terminal_transit"][5]>0
        continued=restart_case(c,x,1.0)
        rest=dispatch_from_controls(continued,ones(Int,2,2),fill(1.0,2,2),x["gate"][:,4:5])
        @test rest["validation"]["valid"]
        @test rest["V"]≈x["V"][:,4:end] atol=1e-9
        @test rest["terminal_transit"]≈x["terminal_transit"] atol=1e-9
        b=OpenSHOP._build_global_dispatch(c;warm=x)
        @test OpenSHOP._lift_start!(b,c,x)["valid"]
        result=solve_case(c;warm=x,u=ones(Int,2,5),time_limit=15.0)
        @test result["validation"]["valid"]
        @test result["objective"]>=x["objective"]-1e-3
    end
end

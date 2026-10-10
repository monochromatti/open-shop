module MixedSourceNetworkExample
using OpenSHOP
export mixed_source_fixture, mixed_source_small_fixture

series(object,attribute,times,values)=OperationalSeries(
    object=object,attribute=attribute,times=Float64.(times),values=Float64.(values))

function delay_curves(capacity,delay)
    # Three distinct flow references: only the two neighbors of a release blend.
    [RiverRouting.DelayCurve(reference,
        [offset,offset+0.25,offset+0.8],[0.65,0.35])
        for (reference,offset) in ((0.0,delay+0.8),
            (capacity/2,delay+0.35),(capacity,delay))]
end

function reach(name,source,target,capacity,delay,history;
    distributed=false,law=:junction,inflow=0.0,coefficient=0.0,crest=0.0)
    River(;name,source,target,capacity,law,inflow,coefficient,crest,
        deterministic_delay=distributed ? nothing : delay,
        curves=distributed ? delay_curves(capacity,delay) : RiverRouting.DelayCurve[],
        history_grid=[-8.0,-4.0,-2.0,0.0],history_release=fill(history,3),
        water_value=150.0)
end

function unit(name,plant;lowhead=false,table=false)
    heads=lowhead ? [50.0,75.0,100.0] : [100.0,150.0,200.0]
    discharge=[1.0,5.0,10.0,15.0]
    efficiencies=[0.82 0.84 0.83;0.89 0.91 0.90;0.92 0.94 0.93;0.90 0.92 0.91]
    turbine=table ? TurbineTable(heads,discharge,efficiencies,
        fill(1.0,3),fill(lowhead ? 8.0 : 15.0,3);interpolation=:pchip_discharge) : nothing
    Generator(;name,plant,qmin=1.0,qmax=lowhead ? 8.0 : 15.0,
        pmin=0.1,pmax=40.0,efficiency=0.91,qbest=lowhead ? 4.0 : 9.0,
        hbest=lowhead ? 75.0 : 150.0,hmin=first(heads),hmax=last(heads),
        qcurvature=0.0,hcurvature=0.0,turbine_table=turbine,
        initial_on=1,initial_age=8.0,minup=0.3,mindown=0.3,startup=40.0,shutdown=10.0)
end

"""A connected six-reservoir basin, with routed plant and tunnel outfalls.

The four passive supply tunnels form a pressurized loop. The canal's target
and every plant target remain their hydraulic head references; discharged water
enters its named river once and reaches the pond/lower reservoirs after routing.
Natural inflows vary independently of the operating schedule.
"""
function mixed_source_fixture(;distributed=false)
    grid=[0.0,0.4,1.2,2.0,3.3,4.2,6.0]
    reservoirs=[
        Reservoir(name=:UpperA,z0=297.0,slope=0.2,v0=15.0,vmin=10.0,vmax=20.0,
            inflow=4.0,water_value=350.0),
        Reservoir(name=:UpperB,z0=297.6,slope=0.2,v0=12.0,vmin=8.0,vmax=16.0,
            inflow=3.0,water_value=350.0),
        Reservoir(name=:Pond,z0=148.4,slope=0.8,v0=2.0,vmin=1.0,vmax=4.0,
            water_value=220.0),
        Reservoir(name=:Lower,z0=78.8,slope=0.2,v0=6.0,vmin=4.0,vmax=9.0,
            inflow=0.3,water_value=150.0),
        Reservoir(name=:Estuary,z0=9.7,slope=0.1,v0=3.0,vmin=2.0,vmax=5.0,
            water_value=50.0),
        Reservoir(name=:SideLake,z0=169.0,slope=0.2,v0=5.0,vmin=3.0,vmax=7.0,
            inflow=1.5,water_value=220.0),
    ]
    tunnels=[
        Tunnel(name=:SupplyA,source=:UpperA,target=:IntakeA,resistance=0.03,capacity=30.0),
        Tunnel(name=:SupplyB,source=:UpperB,target=:IntakeB,resistance=0.04,capacity=30.0),
        Tunnel(name=:IntakeCross,source=:IntakeA,target=:IntakeB,resistance=0.1,capacity=20.0),
        Tunnel(name=:UpperCross,source=:UpperA,target=:UpperB,resistance=0.05,capacity=20.0),
        Tunnel(name=:RoutedCanal,source=:UpperA,target=:Pond,resistance=1.5,
            capacity=15.0,discharge_river=:CanalRiver),
    ]
    plants=[
        Plant(name=:HighA,source=:IntakeA,target=:Pond,pmax=40.0,discharge_river=:TailA),
        Plant(name=:HighB,source=:IntakeB,target=:Pond,pmax=40.0,discharge_river=:TailB),
        Plant(name=:LowStation,source=:Pond,target=:Lower,pmax=15.0,discharge_river=:LowTail),
    ]
    generators=[unit(:UnitA,:HighA;table=true),unit(:UnitB,:HighB),
        unit(:LowUnit,:LowStation;lowhead=true,table=true)]
    rivers=[
        reach(:TailA,:auto,:MountainMerge,25.0,0.15,10.4;distributed,inflow=0.4),
        reach(:TailB,:auto,:MountainMerge,25.0,1.75,6.0;distributed),
        reach(:CanalRiver,:auto,:UpperMerge,20.0,0.25,10.0;distributed),
        reach(:SideRelease,:SideLake,:UpperMerge,6.0,0.1,2.0;distributed,law=:controlled),
        reach(:MountainMerge,:auto,:UpperMerge,40.0,0.35,17.1;distributed,inflow=0.7),
        reach(:UpperMerge,:auto,:Pond,60.0,0.2,29.4;distributed,inflow=0.3),
        reach(:LowTail,:auto,:Lower,15.0,1.25,4.2;distributed,inflow=0.2),
        reach(:PondOutlet,:Pond,:Lower,10.0,0.4,2.0;distributed,law=:controlled),
        reach(:PondSpill,:Pond,:Lower,10.0,0.15,0.6;distributed,
            law=:weir,coefficient=0.6,crest=149.0),
        reach(:LowerRelease,:Lower,:Estuary,15.0,1.1,6.0;distributed,law=:controlled),
        reach(:EstuaryOutlet,:Estuary,:Sea,10.0,0.6,sqrt(10.0);distributed,
            law=:orifice,coefficient=1.0),
    ]
    operations=[series(:TailA,:inflow,[0.0,2.0,4.2],[0.4,0.8,0.3]),
        series(:MountainMerge,:inflow,[0.0,1.2,3.3],[0.7,1.1,0.6]),
        series(:RoutedCanal,:opening,[0.0,2.0,4.2],[1.0,0.8,1.0])]
    system=HydroSystem(;reservoirs,junctions=[Junction(name=:IntakeA,hmin=285.0,hmax=305.0),
        Junction(name=:IntakeB,hmin=285.0,hmax=305.0)],
        boundaries=[Boundary(name=:Sea,head=0.0)],tunnels,plants,generators,rivers)
    raw=ScheduleCase(;name=distributed ? "mixed_sources_distributed" : "mixed_sources_exact",
        system,grid,prices=[35.0,80.0,45.0,100.0,30.0,70.0],operations)
    # Round trip through the public representation compiles direct connections.
    case=case_from_dict(case_dict(raw))
    u=[1 1 1 0 1 1;1 1 1 1 1 1;1 0 1 1 0 1]
    q=[10.0 9.0 11.0 0.0 10.0 8.0;6.0 6.5 6.0 6.0 5.5 6.0;
        4.0 0.0 3.0 3.5 0.0 4.0]
    gates=zeros(length(rivers),6)
    gate_for=Dict(:SideRelease=>2.0/6.0,:PondOutlet=>2.0/10.0,
        :PondSpill=>1.0,:LowerRelease=>6.0/15.0,:EstuaryOutlet=>1.0)
    for (i,r) in enumerate(rivers)
        gates[i,:].=get(gate_for,r.name,0.0)
    end
    (;case,u,q,gates)
end

"""Small exact-routing case with two scheduled units and independent inflows.
The discharge schedules fix commitment; SCIP's certificate still covers the
declared complete case, since the fixed states are operational requirements.
"""
function mixed_source_small_fixture(;mixed=false)
    grid=[0.0,0.5,1.5,2.5]
    reservoirs=[
        Reservoir(name=:UpperA,z0=99.0,slope=0.2,v0=5.0,vmin=4.0,vmax=6.0,water_value=0.0),
        Reservoir(name=:UpperB,z0=119.0,slope=0.2,v0=5.0,vmin=4.0,vmax=6.0,
            inflow=3.0,water_value=0.0),
        Reservoir(name=:Receiver,z0=49.0,slope=0.2,v0=5.0,vmin=4.0,vmax=6.0,water_value=0.0),
    ]
    plants=[Plant(name=Symbol("Station",i),source=Symbol("Upper",i==1 ? "A" : "B"),
        target=:Receiver,pmax=10.0,discharge_river=Symbol("Tail",i)) for i in 1:2]
    generators=[Generator(name=Symbol("Unit",i),plant=plants[i].name,qmin=1.0,qmax=5.0,
        pmin=0.1,pmax=10.0,efficiency=0.9,qbest=3.0,hbest=i==1 ? 50.0 : 70.0,
        hmin=40.0,hmax=80.0,qcurvature=0.0,hcurvature=0.0,
        initial_on=1,initial_age=8.0,minup=0.0,mindown=0.0,startup=0.0,shutdown=0.0)
        for i in 1:2]
    rivers=[
        reach(:Tail1,:auto,:MainRiver,8.0,0.25,2.25;inflow=0.25,distributed=mixed),
        reach(:Tail2,:auto,:MainRiver,8.0,0.75,3.0),
        reach(:MainRiver,:auto,:Receiver,20.0,0.3,5.75;inflow=0.5),
        reach(:ReceiverOutlet,:Receiver,:Sea,10.0,0.4,3.0;law=:controlled),
    ]
    # This fixture tests electricity revenue alone; changing river inventories
    # are still physically conserved and audited, with zero economic water value.
    rivers=[River(;Dict(f=>getfield(r,f) for f in fieldnames(River) if f!=:water_value)...,
        water_value=0.0) for r in rivers]
    operations=[series(:UpperA,:inflow,grid[1:end-1],[2.0,0.0,2.0]),
        series(:Unit1,:discharge,grid[1:end-1],[2.0,0.0,2.0]),
        series(:Unit2,:discharge,[0.0],[3.0]),
        series(:ReceiverOutlet,:gate_min,[0.0],[0.3]),
        series(:ReceiverOutlet,:gate_max,[0.0],[0.3])]
    system=HydroSystem(;reservoirs,junctions=Junction[],
        boundaries=[Boundary(name=:Sea,head=0.0)],tunnels=Tunnel[],plants,generators,rivers)
    raw=ScheduleCase(;name=mixed ? "mixed_sources_small_heterogeneous" : "mixed_sources_small_exact",system,grid,
        prices=[40.0,90.0,55.0],operations)
    case=case_from_dict(case_dict(raw))
    u=[1 0 1;1 1 1];q=[2.0 0.0 2.0;3.0 3.0 3.0]
    gates=zeros(4,3);gates[4,:].=0.3
    (;case,u,q,gates)
end

function main()
    large=mixed_source_fixture()
    witness=dispatch_from_controls(large.case,large.u,large.q,large.gates)
    witness["validation"]["valid"] || error(join(witness["validation"]["errors"],"; "))
    replay_audit(large.case,witness)["valid"] || error("network finer replay failed")
    println("Six-reservoir mixed-source network: original and finer audits passed.")
    small=mixed_source_small_fixture()
    initial=dispatch_from_controls(small.case,small.u,small.q,small.gates)
    result=solve(small.case;initial,time_limit=60.0,relative_gap=1e-4)
    println("Small exact-routing case accepted: ",result["accepted"],
        "; global certificate: ",result["global_certificate"])
    println("Objective interval: [",result["feasible_lower_bound"],", ",result["global_bound"],"]")
    result["accepted"] && result["global_certificate"] || error("small case failed")
end
end

if abspath(PROGRAM_FILE)==@__FILE__
    MixedSourceNetworkExample.main()
end

module RiverSplitsExample
using OpenSHOP
export river_split_case

"""One generating unit feeds a 40/60 river fork that rejoins downstream.

The left branch also receives natural inflow and a pressure-link outfall.
Those local sources are not shared with the right branch. Different travel
times preserve each branch's arrival times before the downstream merge.
"""
function river_split_case()
    reservoirs=[Reservoir(name=name,z0=head,slope=0.2,v0=5.0,
        vmin=4.0,vmax=6.0,water_value=0.0)
        for (name,head) in ((:Upper,100.0),(:Receiver,50.0))]
    plant=Plant(name=:Station,source=:Upper,target=:Receiver,pmax=5.0,
        discharge_river=:Incoming)
    unit=Generator(name=:Unit,plant=:Station,qmin=1.0,qmax=5.0,
        pmin=0.1,pmax=5.0,hmin=40.0,hmax=60.0,qbest=3.0,hbest=50.0,
        efficiency=0.9,qcurvature=0.0,hcurvature=0.0,
        minup=0.0,mindown=0.0,startup=0.0,shutdown=0.0)
    # A separate supply enters Left through a fixed-head pressure link.
    tunnel=Tunnel(name=:Canal,source=:Supply,target=:CanalOutlet,
        resistance=50.0,capacity=2.0,discharge_river=:Left)
    rivers=[River(;name,source,target,allocation,deterministic_delay=delay,
        curves=RiverRouting.DelayCurve[],capacity,inflow,water_value=0.0,
        history_grid=[-4.0,0.0],history_release=[history])
        for (name,source,target,allocation,delay,capacity,inflow,history) in (
            (:Incoming,:auto,:Fork,1.0,0.25,10.0,0.0,3.5),
            (:Left,:Fork,:Merge,0.4,0.35,3.6,1.0,3.4),
            (:Right,:Fork,:Merge,0.6,0.8,10.0,0.0,2.1),
            (:Downstream,:Merge,:Receiver,1.0,0.2,15.0,0.5,6.0))]
    system=HydroSystem(;reservoirs,
        boundaries=[Boundary(name=:Supply,head=100.0),Boundary(name=:CanalOutlet,head=50.0)],
        plants=[plant],generators=[unit],
        tunnels=[tunnel],rivers)
    ScheduleCase(name="prescribed river fork",system=system,
        grid=[0.0,0.5,1.5,2.5,4.0],prices=[40.0,90.0,55.0,70.0])
end

function main()
    case=river_split_case()
    result=solve(case;time_limit=180.0,relative_gap=1e-4)
    result.accepted || error("No schedule passed validation and finer replay: $(result.status)")
    println("Unit discharge (m³/s): ",vec(result.solution.generator_q))
    println("Left branch discharge (m³/s): ",vec(result.solution.river_release[2,:]))
    println("Right branch discharge (m³/s): ",vec(result.solution.river_release[3,:]))
    println("Global certificate: ",result.global_certificate)
end
end

if abspath(PROGRAM_FILE)==@__FILE__
    RiverSplitsExample.main()
end

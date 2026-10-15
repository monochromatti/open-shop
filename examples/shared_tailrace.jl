module SharedTailraceExample
using OpenSHOP
export shared_tailrace_case

"""Two independent penstocks discharging through one pressurized tailrace.

Each Plant is a hydraulic supply group. Their outlet node is shared, so the
tailrace tunnel carries the sum of their flows. Station-wide operating rules
are not implicitly pooled across the two Plants.
"""
function shared_tailrace_case()
    system=HydroSystem(
        reservoirs=[Reservoir(name=:Lake,z0=149.0,slope=0.01,v0=100.0,
            vmin=90.0,vmax=110.0,water_value=0.0)],
        junctions=[Junction(name=:PenstockA,hmin=130.0,hmax=151.0),
            Junction(name=:PenstockB,hmin=120.0,hmax=151.0),
            Junction(name=:Tail,hmin=0.0,hmax=15.0)],
        boundaries=[Boundary(name=:Sea,head=0.0)],
        tunnels=[
            Tunnel(name=:SupplyA,source=:Lake,target=:PenstockA,
                resistance=0.08,capacity=10.0),
            Tunnel(name=:SupplyB,source=:Lake,target=:PenstockB,
                resistance=0.20,capacity=10.0),
            Tunnel(name=:Tailrace,source=:Tail,target=:Sea,
                resistance=0.03,capacity=10.0)],
        plants=[Plant(name=:BranchA,source=:PenstockA,target=:Tail,pmax=20.0,ramp=nothing),
            Plant(name=:BranchB,source=:PenstockB,target=:Tail,pmax=20.0,ramp=nothing)],
        generators=[Generator(name=name,plant=plant,qmin=1.0,qmax=10.0,
            pmin=0.1,pmax=20.0,efficiency=eta,qbest=5.0,hbest=140.0,
            hmin=100.0,hmax=160.0,qcurvature=0.0,hcurvature=0.0,
            initial_on=0,minup=0.0,mindown=0.0,startup=cost,shutdown=0.0)
            for (name,plant,eta,cost) in ((:UnitA,:BranchA,0.90,2.0),
                (:UnitB,:BranchB,0.92,3.0))],
        rivers=River[])
    ScheduleCase(name="Shared tailrace",system=system,
        grid=[0.0,0.5,1.5],prices=[60.0,100.0])
end

function main()
    case=shared_tailrace_case()
    initial=dispatch_from_controls(case,ones(Int,2,2),[5.0 5.0;3.0 3.0],zeros(0,2))
    initial["validation"]["valid"] || error("Initial hydraulic schedule failed")
    result=solve(case;initial,time_limit=60.0,relative_gap=1e-4)
    result.accepted || error("Shared-tailrace optimization failed")
    result.solution["objective"]>initial["objective"]+100.0 ||
        error("Optimization did not improve the initial schedule")
    println("Unit discharges (m³/s): ",result.solution["generator_q"])
    println("Shared tailrace discharge (m³/s): ",result.solution["tunnel_q"][3,:])
    println("Objective: ",result.solution["objective"])
end
end

if abspath(PROGRAM_FILE)==@__FILE__
    SharedTailraceExample.main()
end

module PumpedStorageExample
using OpenSHOP
export pumped_storage_case

"""Synthetic 250 m lift, two 20 m³/s operating modes and four hourly prices.

The efficiency and reservoir data are illustrative, not measured calibration.
Equal discharge bounds prescribe a fixed-flow operating point; this example
does not claim a fixed rotational-speed pump characteristic.
"""
function pumped_storage_case(; table=false, connected=false)
    reservoirs=[
        Reservoir(name=:Upper,z0=299.5,slope=0.1,v0=5.0,vmin=1.0,vmax=20.0,
            water_value=18000.0),
        Reservoir(name=:Lower,z0=48.5,slope=0.1,v0=15.0,vmin=1.0,vmax=20.0,
            water_value=1000.0),
    ]
    pump_target=connected ? :Delivery : :Upper
    plants=[
        Plant(name=:Generating,source=:Upper,target=:Lower,pmax=60.0,ramp=nothing),
        Plant(name=:Pumping,source=:Lower,target=pump_target,pmax=65.0,ramp=nothing),
    ]
    pump_table=table ? TurbineTable([200.0,300.0],[0.0,20.0],
        [0.86 0.90;0.86 0.90],[20.0,20.0],[20.0,20.0]) : nothing
    units=[
        Generator(name=:Turbine,plant=:Generating,qmin=20.0,qmax=20.0,pmin=40.0,pmax=60.0,
            hmin=200.0,hmax=300.0,efficiency=0.91,qbest=20.0,hbest=250.0,
            qcurvature=0.0,hcurvature=0.0,initial_on=0,minup=0.0,mindown=0.0,startup=10.0),
        Pump(name=:Pump,plant=:Pumping,qmin=20.0,qmax=20.0,pmin=50.0,pmax=65.0,
            hmin=200.0,hmax=300.0,efficiency=0.88,motor_efficiency=0.96,
            qbest=20.0,hbest=250.0,qcurvature=0.0,hcurvature=0.0,minup=0.0,mindown=0.0,
            startup=10.0,turbine_table=pump_table),
    ]
    system=HydroSystem(reservoirs=reservoirs,plants=plants,generators=units,
        junctions=connected ? [Junction(name=:Delivery,hmin=290.0,hmax=330.0)] : Junction[],
        tunnels=connected ? [Tunnel(name=:DeliveryPipe,source=:Delivery,target=:Upper,
            resistance=0.002,capacity=40.0)] : Tunnel[],
        commitment_exclusions=[CommitmentExclusion(name=:Reversible,
            groups=[[:Turbine],[:Pump]])])
    ScheduleCase(name="synthetic_pumped_storage",system=system,
        grid=collect(0.0:4.0),prices=[10.0,15.0,110.0,120.0])
end

function run_example()
    case=pumped_storage_case()
    result=solve(case;time_limit=180.0,relative_gap=1e-4)
    result.accepted || error("pumped-storage schedule was not accepted: $(result.status)")
    println("Unit order: ", [g.name for g in case.system.generators])
    println("Net electrical exchange (MW): ",net_power(case,result.solution))
    println("Objective: ",result.solution.objective)
    result
end
end

if abspath(PROGRAM_FILE)==@__FILE__
    PumpedStorageExample.run_example()
end

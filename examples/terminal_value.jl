module TerminalValueExample
using OpenSHOP
export terminal_value_case

"""A two-station cascade with explicit reference-based future water values.

Upper water can generate at both stations. Water already between stations can
generate only below. Historical water in the half-hour reach is valued too.
All prices, costs and values use the same caller-chosen currency.
"""
function terminal_value_case()
    future_price=30.0
    lower_value=volume_water_value(future_price,energy_conversion_factor(80.0,0.92))
    upper_value=volume_water_value(future_price,energy_conversion_factor(120.0,0.90);
        downstream_value=lower_value)
    lakes=[Reservoir(name=:UpperLake,z0=198.0,slope=1.0,v0=2.0,vmin=1.0,vmax=3.0,
            water_value=upper_value),
        Reservoir(name=:LowerLake,z0=78.0,slope=1.0,v0=2.0,vmin=1.0,vmax=3.0,
            water_value=lower_value)]
    units=[Generator(name=:UpperUnit,plant=:UpperStation,qmin=1.0,qmax=4.0,
            pmin=0.1,pmax=8.0,hmin=110.0,hmax=130.0,efficiency=0.90,qbest=4.0,hbest=120.0,
            qcurvature=0.0,hcurvature=0.0,initial_on=0,minup=0.0,mindown=0.0,startup=3.0,shutdown=4.0),
        Generator(name=:LowerUnit,plant=:LowerStation,qmin=1.0,qmax=2.0,
            pmin=0.1,pmax=4.0,hmin=70.0,hmax=90.0,efficiency=0.92,qbest=2.0,hbest=80.0,
            qcurvature=0.0,hcurvature=0.0,initial_on=0,minup=0.0,mindown=0.0,startup=5.0,shutdown=7.0)]
    reach=River(name=:BetweenStations,source=:Outfall,target=:LowerLake,law=:junction,
        capacity=10.0,deterministic_delay=0.5,curves=RiverRouting.DelayCurve[],
        history_grid=[-1.0,0.0],history_release=[2.0],water_value=lower_value)
    system=HydroSystem(reservoirs=lakes,generators=units,
        plants=[Plant(name=:UpperStation,source=:UpperLake,target=:LowerLake,pmax=8.0,
                discharge_river=:BetweenStations,ramp=nothing),
            Plant(name=:LowerStation,source=:LowerLake,target=:Sea,pmax=4.0,ramp=nothing)],
        river_junctions=[RiverJunction(name=:Outfall)],rivers=[reach],
        boundaries=[Boundary(name=:Sea,head=0.0)])
    ScheduleCase(name="energy-based terminal water value",system=system,
        grid=[0.0,1.0,2.0],prices=[100.0,100.0])
end

function main()
    c=terminal_value_case()
    result=solve(c;time_limit=60.0)
    result.accepted || error("No schedule passed validation and finer replay")
    println("Reservoir values (currency/Mm³): ",[r.water_value for r in c.system.reservoirs])
    println("Objective components: ",objective_components(c,result.solution))
end
end

if abspath(PROGRAM_FILE)==@__FILE__
    TerminalValueExample.main()
end

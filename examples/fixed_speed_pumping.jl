module FixedSpeedPumpingExample
using OpenSHOP
export fixed_speed_case

"""Synthetic fixed-speed lift with one operating point at each supplied head.

The characteristic is illustrative, not equipment calibration. The two basin
levels change sufficiently to make the operating flow change during pumping.
"""
function fixed_speed_case()
    lower=Reservoir(name=:Lower,z0=-10.0,v0=10.0,vmin=8.0,vmax=12.0,water_value=0.0)
    upper=Reservoir(name=:Upper,z0=90.0,v0=10.0,vmin=8.0,vmax=12.0,water_value=0.0)
    characteristic=PumpCharacteristic([90.0,100.0,110.0],
        [10.2,10.0,9.8],[0.845,0.85,0.855])
    pump=Pump(name=:Pump,plant=:Lift,characteristic=characteristic,qmin=9.0,qmax=11.0,
        pmin=1.0,pmax=20.0,hmin=90.0,hmax=110.0,qbest=10.0,hbest=100.0,
        motor_efficiency=0.95,minup=0.0,mindown=0.0,startup=0.0)
    system=HydroSystem(reservoirs=[lower,upper],generators=[pump],plants=[
        Plant(name=:Lift,source=:Lower,target=:Upper,pmax=20.0,ramp=nothing)])
    ScheduleCase(name="synthetic_fixed_speed_pumping",system=system,
        grid=collect(0.0:4.0),prices=[-30.0,-20.0,50.0,60.0])
end

function run_example()
    case=fixed_speed_case()
    # Commitment controls the fixed-speed pump. Zero flow guesses are allowed.
    initial=dispatch_from_controls(case,[1 1 0 0],zeros(1,4),zeros(0,4))
    result=solve(case;initial,time_limit=60.0,relative_gap=1e-4)
    result.accepted || error("fixed-speed pumping schedule was not accepted")
    x=result.solution
    println("On/off schedule: ",x.commitment)
    println("Characteristic discharge (m³/s): ",x.generator_q)
    println("Net electrical exchange (MW): ",net_power(case,x))
    for t in eachindex(case.prices)
        x.commitment[1,t]==0 && continue
        h=x.head[2,t]-x.head[1,t]
        eta=0.8+0.0005h
        expected=0.00981*x.generator_q[1,t]*h/(eta*0.95)
        isapprox(x.power[1,t],expected;atol=1e-8) || error("consumption mismatch")
    end
    result
end
end

if abspath(PROGRAM_FILE)==@__FILE__
    FixedSpeedPumpingExample.run_example()
end

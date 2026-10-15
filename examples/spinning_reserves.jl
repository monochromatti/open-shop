using OpenSHOP

# A small synthetic station. Qualified activation energy is external input;
# reserve capacity is not an additional flow in the scheduled water balance.
lake=Reservoir(name=:Lake,z0=100.0,slope=0.01,v0=100.0,vmin=10.0,vmax=1000.0,
    water_value=0.0)
unit=Generator(name=:Unit,plant=:Station,qmin=2.0,qmax=10.0,pmin=1.7658,pmax=8.829,
    efficiency=0.9,qcurvature=0.0,hcurvature=0.0,qbest=5.0,hbest=100.0,
    hmin=50.0,hmax=150.0,initial_on=0,initial_age=8.0,minup=0.0,mindown=0.0,
    startup=5.0,ramp_up=20.0,ramp_down=20.0,
    reserve_up_energy=1.0,reserve_down_energy=0.8)
station=Plant(name=:Station,source=:Lake,target=:Sea,pmax=8.829,ramp=nothing)
system=HydroSystem(reservoirs=[lake],generators=[unit],plants=[station],
    boundaries=[Boundary(name=:Sea,head=0.0)])
reserve=ReserveGroup(name=:Spinning,generators=[:Unit],up_requirement=3.0,
    down_requirement=1.0,up_price=100.0,down_price=25.0,up_max=4.0,down_max=4.0,
    response_hours=0.25,activation_hours=0.25)
case=ScheduleCase(name="energy and spinning capacity",system=system,
    grid=[0.0,1.0,2.0],prices=[50.0,50.0],reserve_groups=[reserve])

result=solve(case;time_limit=180.0)
@assert result.accepted
schedule=result.solution
@assert validate(case,schedule)["valid"]
println("Status: ",result.status,"; accepted: ",result.accepted,
    "; global certificate: ",result.global_certificate)
println("Scheduled power, MW: ",vec(schedule.power))
println("Upward capacity, MW: ",vec(schedule.reserves.up))
println("Downward capacity, MW: ",vec(schedule.reserves.down))
println("Objective: ",schedule.objective)

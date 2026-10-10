using OpenSHOP

# Synthetic two-unit station with a maintenance window and a required release
# schedule. Storage and head evolve under the original nonlinear equations.
lake=Reservoir(name=:Lake,z0=98.0,slope=0.2,v0=10.0,vmin=8.0,vmax=12.0,
    inflow=10.0,water_value=2000.0,level_ramp_down=0.025)
station=Plant(name=:Station,source=:Lake,target=:Tailwater,pmax=40.0,
    pmin=2.0,qmax=35.0,ramp_up=12.0,ramp_down=20.0,
    discharge_ramp_up=10.0,discharge_ramp_down=20.0,
    initial_power=0.0,initial_discharge=0.0,
    initial_on=0,initial_age=8.0,minup=2.0,mindown=1.0)
units=[Generator(name=name,plant=:Station,qmin=4.0,qmax=20.0,
    pmin=2.0,pmax=20.0,efficiency=0.94,qbest=12.0,qcurvature=0.1,
    hbest=100.0,hcurvature=0.01,hmin=90.0,hmax=110.0,
    initial_on=0,initial_age=8.0,minup=1.0,mindown=1.0,
    startup=30.0,shutdown=0.0) for name in (:UnitA,:UnitB)]
system=HydroSystem(reservoirs=[lake],junctions=Junction[],
    boundaries=[Boundary(name=:Tailwater,head=0.0)],tunnels=Tunnel[],
    plants=[station],generators=units,rivers=River[])
case=ScheduleCase(name="station operating rules",system=system,
    grid=collect(0.0:1.0:6.0),prices=[40.0,30.0,80.0,100.0,60.0,20.0],
    operations=[
        OperationalSeries(object=:UnitB,attribute=:maintenance,
            times=[0.0,3.0],values=[1.0,0.0]),
        OperationalSeries(object=:Station,attribute=:discharge,
            times=collect(0.0:1.0:5.0),values=[0.0,10.0,20.0,30.0,20.0,0.0]),
    ])

result=solve(case;time_limit=60.0,relative_gap=1e-3)
result["accepted"] || error("no validated schedule: $(result["status"])")
schedule=result["solution"]
println("Unit states: ",schedule["u"])
println("Plant discharge: ",vec(sum(schedule["generator_q"];dims=1)))
println("Objective: ",schedule["objective"])
println("Global certificate: ",result["global_certificate"])

using OpenSHOP

# Two connected lakes feed one station. The supplied response describes an
# asymmetric open-channel connection; a parallel pressure link has a scheduled gate.
response=LevelFlowCurve([100.0,110.0],[
    TableCurve([-5.0,0.0,5.0],[-15.0,0.0,20.0]),
    TableCurve([-5.0,0.0,5.0],[-25.0,0.0,30.0]),
])
lakes=[
    Reservoir(name=:WestLake,z0=100.0,slope=1.0,v0=3.0,vmin=1.0,vmax=6.0,water_value=0.0),
    Reservoir(name=:EastLake,z0=100.0,slope=1.0,v0=4.0,vmin=1.0,vmax=6.0,water_value=0.0),
]
channel=River(name=:LakeConnection,source=:WestLake,target=:EastLake,
    law=:delta,delta_level_ref_up_flow_curve=response,
    curves=RiverRouting.DelayCurve[],deterministic_delay=0.0,capacity=40.0,water_value=0.0)
tunnel=Tunnel(name=:PressureLink,source=:WestLake,target=:EastLake,
    resistance=1.0,capacity=10.0,opening=0.5,
    gate_loss_factor_curve=TableCurve([0.0,0.5,1.0],[9.0,3.0,0.0]))
plant=Plant(name=:Station,source=:WestLake,target=:Tail,pmax=5.0,ramp=nothing)
unit=Generator(name=:Unit,plant=:Station,qmin=0.5,qmax=2.0,pmin=0.1,pmax=5.0,
    hmin=40.0,hmax=60.0,qbest=1.0,hbest=53.0,efficiency=0.9,
    qcurvature=0.0,hcurvature=0.0,minup=0.0,mindown=0.0,
    initial_on=0,startup=0.0)
system=HydroSystem(reservoirs=lakes,junctions=Junction[],
    boundaries=[Boundary(name=:Tail,head=50.0)],tunnels=[tunnel],
    plants=[plant],generators=[unit],rivers=[channel])
case=ScheduleCase(name="supplied lake responses",system=system,grid=collect(0.0:4.0),
    prices=[-20.0,80.0,120.0,-10.0],operations=[
        OperationalSeries(object=:PressureLink,attribute=:opening,times=[2.0],values=[0.0])])
result=OpenSHOP.solve(case;time_limit=60.0,relative_gap=1e-3)
result["accepted"] || error("No schedule passed physical validation and finer replay")
solution=result["solution"]
println("Objective: ",solution["objective"])
println("Unit commitment: ",vec(solution["u"]))
println("Signed channel discharge (m³/s): ",vec(solution["river_release"]))
println("Pressure-link discharge (m³/s): ",vec(solution["tunnel_q"]))
println("Global certificate: ",result["global_certificate"])

module DirectionalInterlocksExample
using OpenSHOP
export directional_interlock_case

"""Two reservoir intakes, with different plant shutdown times on a handover.

The next intake may open immediately. The station must stay stopped for 1.5 h
after HighIntake closes and .5 h after LowIntake closes. Gate schedules describe
the handovers; unit commitment and nonlinear discharge remain decisions.
"""
function directional_interlock_case()
    reservoirs=[Reservoir(name=name,z0=head-0.5,slope=0.1,v0=5.0,
        vmin=4.0,vmax=6.0,water_value=0.0)
        for (name,head) in ((:HighLake,100.0),(:LowLake,80.0))]
    tunnels=[Tunnel(name=name,source=source,target=:Intake,resistance=1.0,
        capacity=6.0,opening=initial,initial_opening=initial,
        gate_positions=[0.0,0.5,1.0],gate_adjustment_cost=5.0)
        for (name,source,initial) in ((:HighIntake,:HighLake,1.0),(:LowIntake,:LowLake,0.0))]
    unit=Generator(name=:Unit,plant=:Station,qmin=1.0,qmax=5.0,pmin=0.1,pmax=5.0,
        hmin=40.0,hmax=110.0,efficiency=0.9,qbest=5.0,hbest=75.0,
        qcurvature=0.0,hcurvature=0.0,initial_on=1,initial_age=8.0,
        minup=0.0,mindown=0.0,startup=0.0,shutdown=0.0)
    interlock=TunnelInterlock(name=:SupplyHandover,tunnels=[:HighIntake,:LowIntake],
        min_open=1,max_open=1,forward_tunnels=[:HighIntake],backward_tunnels=[:LowIntake],
        plants=[:Station],forward_switch_time=1.5,backward_switch_time=0.5)
    system=HydroSystem(;reservoirs,junctions=[Junction(name=:Intake,hmin=40.0,hmax=101.0)],
        boundaries=[Boundary(name=:Tail,head=0.0)],tunnels,
        plants=[Plant(name=:Station,source=:Intake,target=:Tail,pmax=5.0,ramp=nothing)],
        generators=[unit],tunnel_interlocks=[interlock])
    operations=[OperationalSeries(object=name,attribute=:opening,
        times=[0.0,1.0,3.0],values=values)
        for (name,values) in ((:HighIntake,[1.0,0.0,1.0]),(:LowIntake,[0.0,1.0,0.0]))]
    ScheduleCase(;name="directional supply handovers",system,
        grid=[0.0,0.5,1.0,2.0,2.5,3.0,3.5,4.5,5.0],
        prices=[40.0,50.0,90.0,100.0,60.0,110.0,75.0,35.0],operations)
end

function main()
    case=directional_interlock_case()
    result=solve(case;time_limit=90.0,relative_gap=1e-4)
    result.accepted || error("No schedule passed validation and finer replay")
    println("Grid (hours): ",case.grid)
    println("Unit commitment: ",vec(result.solution.commitment))
    println("Unit discharge (m³/s): ",vec(result.solution.generator_q))
    println("Gate positions: ",result.solution.tunnel_opening)
    println("Global certificate: ",result.global_certificate)
end
end

if abspath(PROGRAM_FILE)==@__FILE__
    DirectionalInterlocksExample.main()
end

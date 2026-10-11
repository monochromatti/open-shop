module TunnelGatesExample
using OpenSHOP

export tunnel_gate_case

"""Two reservoirs supply one station, with one intake gate open at a time."""
function tunnel_gate_case()
    reservoirs=[
        Reservoir(name=:HighLake,z0=99.5,slope=0.1,v0=5.0,vmin=4.0,vmax=6.0,water_value=0.0),
        Reservoir(name=:LowLake,z0=79.5,slope=0.1,v0=5.0,vmin=4.0,vmax=6.0,water_value=0.0),
    ]
    tunnels=[
        Tunnel(name=:HighIntake,source=:HighLake,target=:Intake,resistance=1.0,
            capacity=6.0,opening=1.0,gate_positions=[0.0,1.0],
            initial_opening=1.0,gate_adjustment_cost=5.0),
        Tunnel(name=:LowIntake,source=:LowLake,target=:Intake,resistance=1.0,
            capacity=6.0,opening=0.0,gate_positions=[0.0,1.0],
            initial_opening=0.0,gate_adjustment_cost=5.0),
    ]
    unit=Generator(name=:Unit,plant=:Station,qmin=1.0,qmax=5.0,pmin=0.1,pmax=5.0,
        hmin=40.0,hmax=110.0,efficiency=0.9,qbest=5.0,hbest=75.0,
        qcurvature=0.0,hcurvature=0.0,initial_on=1,initial_age=8.0,
        minup=1.0,mindown=1.0,startup=0.0,shutdown=0.0)
    system=HydroSystem(;reservoirs,junctions=[Junction(name=:Intake,hmin=40.0,hmax=101.0)],
        boundaries=[Boundary(name=:Tail,head=0.0)],tunnels,
        plants=[Plant(name=:Station,source=:Intake,target=:Tail,pmax=5.0,ramp=nothing)],
        generators=[unit],rivers=River[],
        tunnel_interlocks=[TunnelInterlock(name=:IntakeSelection,
            tunnels=[:HighIntake,:LowIntake],min_open=1,max_open=1)])
    ScheduleCase(;name="optimized reservoir supply gates",system,
        grid=[0.0,1.0,2.0,3.0],prices=[40.0,-40.0,40.0],
        operations=[OperationalSeries(object=:Unit,attribute=:discharge,
            times=[0.0],values=[5.0])])
end

function main()
    case=tunnel_gate_case()
    initial=dispatch_from_controls(case,ones(Int,1,3),fill(5.0,1,3),zeros(0,3);
        tunnel_opening=[1.0 1.0 1.0;0.0 0.0 0.0])
    result=OpenSHOP.solve(case;initial,time_limit=60.0,relative_gap=1e-4)
    result["accepted"] || error("No gate schedule passed physical validation and finer replay")
    schedule=result["solution"]
    # At q=5, each open intake loses R*q²=25 m. The expected switch is
    # HighLake → LowLake → HighLake; reservoir depletion lowers these heads.
    head=[74.9991,54.9991,74.9973]
    expected=sum(case.prices .* (0.00981*0.9*5.0 .* head))-20.0
    isapprox(schedule["objective"],expected;atol=1e-3) ||
        error("Returned schedule did not reach the independently calculated optimum")
    println("Tunnel gate positions (rows: HighIntake, LowIntake): ",schedule["tunnel_opening"])
    println("Adjustment cost: ",schedule["gate_adjustment_cost"])
    println("Objective: ",schedule["objective"])
    println("Accepted: ",result["accepted"])
    println("Global certificate: ",result["global_certificate"])
end
end

if abspath(PROGRAM_FILE)==@__FILE__
    TunnelGatesExample.main()
end

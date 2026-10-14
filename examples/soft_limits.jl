module SoftLimitsExample
using OpenSHOP
export soft_limit_case

"""Dispatch may spend up to .018 Mm³ of a terminal storage target, at a price."""
function soft_limit_case(;penalty=1000.0,allowance=0.018)
    reservoir=Reservoir(name=:Lake,z0=99.0,slope=1.0,v0=1.0,vmin=0.95,vmax=2.0,water_value=0.0)
    generator=Generator(name=:Unit,plant=:Station,qmin=1.0,qmax=5.0,pmin=0.1,pmax=5.0,
        hmin=90.0,hmax=110.0,efficiency=0.9,qbest=3.0,hbest=100.0,
        qcurvature=0.0,hcurvature=0.0,initial_on=0,initial_age=8.0,
        minup=0.0,mindown=0.0,startup=0.0,shutdown=0.0)
    system=HydroSystem(reservoirs=[reservoir],boundaries=[Boundary(name=:Tail,head=0.0)],
        plants=[Plant(name=:Station,source=:Lake,target=:Tail,pmax=5.0)],generators=[generator])
    target=SoftConstraint(name=:TerminalStorage,object=:Lake,quantity=:volume,time=1.0,
        lower=1.0,max_lower_slack=allowance,lower_cost=penalty)
    ScheduleCase(;name="priced terminal storage target",system,grid=[0.0,1.0],prices=[40.0],
        soft_constraints=[target])
end

function main()
    for penalty in (1000.0,20000.0)
        c=soft_limit_case(;penalty)
        initial=dispatch_from_controls(c,zeros(Int,1,1),zeros(1,1),zeros(0,1))
        result=solve(c;initial,time_limit=60.0,relative_gap=1e-4)
        result.accepted || error("No schedule passed physical and penalty checks")
        println("Penalty price: ",penalty," currency/Mm³; discharge: ",result.solution.generator_q[1,1])
        println("Terminal shortfall: ",result.solution.soft_slacks[1].lower," Mm³; cost: ",result.solution.soft_penalty_cost)
        println("Net objective: ",result.solution.objective)
    end
end
end
if abspath(PROGRAM_FILE)==@__FILE__
    SoftLimitsExample.main()
end

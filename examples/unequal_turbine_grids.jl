module UnequalTurbineGridsExample
using OpenSHOP

function unequal_grid_case(;interpolation=:pchip_discharge)
    # Original discharge/efficiency samples at 100 m and 200 m. Each head keeps
    # its own interpolation; these curves are not refitted on a common grid.
    table=TurbineTable(heads=[100.0,200.0],curves=[
        TableCurve([1.0,3.0,5.0],[0.80,0.90,0.88]),
        TableCurve([1.0,2.0,5.0],[0.82,0.87,0.94])],
        qmin=[1.0,1.0],qmax=[5.0,5.0];interpolation)
    lake=Reservoir(name=:Lake,z0=148.0,slope=1.0,v0=2.0,vmin=1.0,vmax=3.0,
        inflow=2.5,water_value=0.0)
    unit=Generator(name=:Unit,plant=:Station,qmin=1.0,qmax=5.0,pmin=0.1,pmax=10.0,
        hmin=100.0,hmax=200.0,qbest=3.0,hbest=150.0,minup=0.0,mindown=0.0,
        startup=0.0,turbine_table=table)
    system=HydroSystem(reservoirs=[lake],generators=[unit],
        boundaries=[Boundary(name=:Sea,head=0.0)],
        plants=[Plant(name=:Station,source=:Lake,target=:Sea,pmax=10.0,ramp=nothing)])
    ScheduleCase(name="unequal turbine discharge grids",system=system,
        grid=[0.0,1.0,2.0],prices=[40.0,80.0])
end

function main()
    case=unequal_grid_case()
    table=only(case.system.generators).turbine_table
    println("Efficiency at q=2.5 m³/s, h=150 m: ",OpenSHOP.turbine_efficiency(table,2.5,150.0))
    # q=5 is a supplied knot at both heads. It gives the hand reference below;
    # the reservoir loses .0036*(5-2.5)=.009 Mm³ in each scheduling hour.
    reference=dispatch_from_controls(case,ones(Int,1,2),fill(5.0,1,2),zeros(0,2))
    heads=[149.9955,149.9865]
    efficiency=0.88 .+ (heads .- 100.0).*(0.94-0.88)/100.0
    expected=0.00981 .* 5.0 .* heads .* efficiency
    @assert validate(case,reference)["valid"]
    @assert isapprox(vec(reference.power),expected;atol=1e-10)
    result=solve(case;initial=reference,time_limit=60.0)
    @assert result.accepted
    println("Accepted objective: ",result.solution.objective)
    println("Global certificate: ",result.global_certificate)
end

end

if abspath(PROGRAM_FILE)==@__FILE__
    UnequalTurbineGridsExample.main()
end

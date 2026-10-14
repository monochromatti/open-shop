using OpenSHOP

# A synthetic table with one reference head. The constant policy assumes that
# efficiency and discharge limits do not change with head. Generation still
# depends on the actual reservoir level. Efficiencies are fractions, not percent.
table = TurbineTable(
    heads = [100.0],
    discharge = [1.0, 3.0],
    efficiency = [0.88, 0.92],
    qmin = [0.5],
    qmax = [4.0],
    head_extrapolation = :constant,
    discharge_extrapolation = :linear,
)

lake = Reservoir(name=:Lake, z0=98.0, slope=1.0, v0=2.0, vmin=1.0,
    vmax=3.0, inflow=2.0, water_value=0.0)
unit = Generator(name=:Unit, plant=:Plant, qmin=0.5, qmax=4.0,
    pmin=0.1, pmax=5.0, hmin=90.0, hmax=110.0, qbest=2.0, hbest=100.0,
    minup=0.0, mindown=0.0, startup=0.0, turbine_table=table,
    generator_efficiency_curve=TableCurve([0.0, 5.0], [0.95, 0.95]))
system = HydroSystem(reservoirs=[lake], generators=[unit],
    boundaries=[Boundary(name=:Sea, head=0.0)],
    plants=[Plant(name=:Plant, source=:Lake, target=:Sea, pmax=5.0)])
case = ScheduleCase(name="single_head_turbine", system=system,
    grid=[0.0, 1.0, 2.0], prices=[50.0, 100.0])

# At q=4, the declared secant gives eta=.94, beyond the final efficiency knot.
# With inflow=2, each one-hour interval loses .0036*(4-2)=.0072 Mm³.
# Thus the midpoint heads are 99.9964 and 99.9892 m.
reference = dispatch_from_controls(case, ones(Int, 1, 2), fill(4.0, 1, 2), zeros(0, 2))
expected_power = 0.00981 .* 4.0 .* [99.9964, 99.9892] .* 0.94 .* 0.95
@assert reference.diagnostics["validation"]["valid"]
@assert isapprox(vec(reference.power), expected_power; atol=1e-9)

result = OpenSHOP.solve(case; initial=reference, time_limit=30.0)
@assert result.accepted
println("Accepted schedule value: ", result.solution.objective)
println("Global certificate: ", result.global_certificate)

using OpenSHOP

lake = Reservoir(
    name = :Lake,
    z0 = 40.0,
    slope = 20.0,
    v0 = 3.0,
    vmin = 1.0,
    vmax = 4.0,
    water_value = 1000.0,
)
plant = Plant(name = :Station, source = :Lake, target = :Tailwater, pmax = 25.0)
unit = Generator(
    name = :Unit,
    plant = :Station,
    qmin = 5.0,
    qmax = 20.0,
    pmin = 1.0,
    pmax = 25.0,
    efficiency = 0.9,
    qbest = 10.0,
    qcurvature = 0.0,
    hbest = 80.0,
    hcurvature = 0.0,
    hmin = 35.0,
    hmax = 105.0,
    initial_on = 0,
    initial_age = 8.0,
    minup = 0.0,
    mindown = 0.0,
    startup = 0.0,
)
system = HydroSystem(
    reservoirs = [lake],
    plants = [plant],
    generators = [unit],
    junctions = Junction[],
    boundaries = [Boundary(name = :Tailwater, head = 20.0)],
    tunnels = Tunnel[],
    rivers = River[],
)
case = ScheduleCase(
    name = "single_reservoir",
    system = system,
    grid = [0.0, 1.0],
    prices = [100.0],
)

result = solve(case; time_limit = 60.0, relative_gap = 1e-4)
println("Status: ", result["status"])
println("Accepted: ", result["accepted"])
println("Certified: ", result["global_certificate"])
println(
    "Objective interval: [",
    result["feasible_lower_bound"],
    ", ",
    result["global_bound"],
    "]",
)
result["accepted"] && result["global_certificate"] || error("example did not pass")

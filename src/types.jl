"""
    Reservoir(; name, v0, vmin, vmax, water_value, ...)

Storage in Mm³ with a level–storage curve in metres. Supply `level_curve`, or
use `z0 + slope*v + curvature*v^2`. `inflow` is m³/s; `water_value` is currency
per Mm³. Initial storage is `v0`, and `vmin`/`vmax` bound the scheduling domain.
"""
Base.@kwdef struct Reservoir
    name::Symbol
    z0::Float64 = 0.0
    slope::Float64 = 1.0
    curvature::Float64 = 0.0
    v0::Float64
    vmin::Float64
    vmax::Float64
    inflow::Float64 = 0.0
    water_value::Float64
    level_curve::Union{Nothing, TableCurve} = nothing
    volume_ramp_up::Union{Nothing, Float64} = nothing # Mm³/hour, between storage vertices
    volume_ramp_down::Union{Nothing, Float64} = nothing
    level_ramp_up::Union{Nothing, Float64} = nothing # metres/hour, between level vertices
    level_ramp_down::Union{Nothing, Float64} = nothing
end
"""
    Junction(; name, hmin=0.0, hmax=700.0)

A hydraulic node with no storage. Its pressure head (metres) is determined by
connected tunnel losses and flow conservation within `[hmin, hmax]`.
"""
Base.@kwdef struct Junction
    name::Symbol
    hmin::Float64 = 0.0
    hmax::Float64 = 700.0
end
"""
    Boundary(; name, head)

An external hydraulic node at fixed head in metres. Water may enter or leave
the modeled system here; it has no modeled storage.
"""
Base.@kwdef struct Boundary
    name::Symbol
    head::Float64
end
"""
    Tunnel(; name, source, target, resistance, capacity, ...)

A signed hydraulic connection with flow in m³/s and quadratic head loss in
metres. Positive flow is from `source` to `target`. `opening` prescribes a gate
position in [0,1]; nonempty `gate_positions` lets the global solve choose among
those positions. Zero is closed. `gate_loss_factor_curve` adds resistance as a
function of position. `discharge_river` changes the water destination while
`target` remains the hydraulic head reference.
"""
Base.@kwdef struct Tunnel
    name::Symbol
    source::Symbol
    target::Symbol
    resistance::Float64
    capacity::Float64
    opening::Float64 = 1.0
    gate_positions::Vector{Float64} = Float64[] # empty keeps opening prescribed
    initial_opening::Union{Nothing, Float64} = nothing
    gate_adjustment_cost::Float64 = 0.0 # cost per absolute position change
    gate_loss_factor_curve::Union{Nothing, TableCurve} = nothing # added q|q| resistance by gate position
    discharge_river::Union{Nothing, Symbol} = nothing # forward outfall; target remains hydraulic head reference
end
"""
    TunnelInterlock(; name, tunnels, min_open=0, max_open=length(tunnels))

Bounds on the number of listed tunnels with positive gate position in each
interval. A positive position counts as open even if its water flow is zero.
Closing a tunnel in `forward_tunnels` or `backward_tunnels` requires the listed
`plants` to stay stopped for the corresponding switch time, in hours. Initial
waits are the remaining parts of pre-horizon closure obligations, also in hours.
"""
Base.@kwdef struct TunnelInterlock
    name::Symbol
    tunnels::Vector{Symbol}
    min_open::Int = 0
    max_open::Int = length(tunnels)
    forward_tunnels::Vector{Symbol} = Symbol[]
    backward_tunnels::Vector{Symbol} = Symbol[]
    plants::Vector{Symbol} = Symbol[]
    forward_switch_time::Float64 = 0.0
    backward_switch_time::Float64 = 0.0
    initial_forward_wait::Float64 = 0.0
    initial_backward_wait::Float64 = 0.0
end
"""
    Plant(; name, source, target, pmax, ...)

A station sharing intake head, outlet head, losses and operating limits
across its units. All members have the same generating or pumping mode;
endpoints follow the actual water direction. Power is MW, flow is m³/s, and dwell times are hours.
`discharge_river` optionally routes its water into a delayed river.
"""
Base.@kwdef struct Plant
    name::Symbol
    source::Symbol
    target::Symbol
    pmax::Float64
    discharge_river::Union{Nothing, Symbol} = nothing # water destination, independent of outlet head
    ramp::Union{Nothing, Float64} = 100.0 # MW/hour; nothing disables the symmetric default
    initial_power::Union{Nothing, Float64} = nothing
    initial_interval_hours::Float64 = 1.0
    tailwater_curve::Union{Nothing, TableCurve} = nothing
    outlet_head_floor::Union{Nothing, Float64} = nothing
    pmin::Float64 = 0.0
    qmin::Float64 = 0.0
    qmax::Union{Nothing, Float64} = nothing
    ramp_up::Union{Nothing, Float64} = nothing # nothing inherits the symmetric power ramp
    ramp_down::Union{Nothing, Float64} = nothing
    discharge_ramp_up::Union{Nothing, Float64} = nothing # (m³/s)/hour
    discharge_ramp_down::Union{Nothing, Float64} = nothing
    initial_discharge::Union{Nothing, Float64} = nothing
    minup::Float64 = 0.0
    mindown::Float64 = 0.0
    initial_on::Union{Nothing, Int} = nothing # otherwise any initially running member unit
    initial_age::Union{Nothing, Float64} = nothing
end
"""
    Generator(; name, plant, qmin, qmax, pmin, pmax, hmin, hmax, ...)

A unit that is either off or within its on-state operating range.
`mode=:generate` produces power; [`Pump`](@ref) constructs `mode=:pump`.
Discharge and power are positive magnitudes in m³/s and MW, effective head is
metres, and efficiencies are fractions. Consumption is subtracted in [`net_power`](@ref).
Supply `turbine_table` for measured efficiency data. The analytic efficiency
parameters otherwise define a quadratic response. Startup/shutdown costs use
the same currency as the objective; minimum dwell times are in hours.
"""
Base.@kwdef struct Generator
    name::Symbol
    plant::Symbol
    mode::Symbol = :generate
    qmin::Float64
    qmax::Float64
    pmin::Float64
    pmax::Float64
    efficiency::Float64 = 0.94
    min_efficiency::Float64 = 0.0
    qbest::Float64 = 1.0
    qcurvature::Float64 = 0.28
    hbest::Float64 = 1.0
    hcurvature::Float64 = 0.015
    hmin::Float64
    hmax::Float64
    initial_on::Int = 1
    initial_age::Float64 = 8.0
    minup::Float64 = 3.0
    mindown::Float64 = 2.0
    startup::Float64 = 120.0
    shutdown::Float64 = 0.0
    turbine_table::Union{Nothing, TurbineTable} = nothing
    generator_efficiency_curve::Union{Nothing, TableCurve} = nothing
    motor_efficiency::Float64 = 1.0
    ramp_up::Union{Nothing, Float64} = nothing
    ramp_down::Union{Nothing, Float64} = nothing
    discharge_ramp_up::Union{Nothing, Float64} = nothing
    discharge_ramp_down::Union{Nothing, Float64} = nothing
    initial_power::Union{Nothing, Float64} = nothing
    initial_discharge::Union{Nothing, Float64} = nothing
    initial_interval_hours::Float64 = 1.0
    reserve_available::Bool = true
    reserve_up_energy::Union{Nothing,Float64} = nothing # qualified activation budget, MWh
    reserve_down_energy::Union{Nothing,Float64} = nothing
end

"""
    Pump(; name, plant, qmin, qmax, pmin, pmax, hmin, hmax, ...)

A pumping mode of the common unit record. `plant.source` is the water intake
and `plant.target` is the pumped-water destination. Discharge and electrical
consumption are positive magnitudes. Hydraulic efficiency uses `efficiency`
or `turbine_table`; `motor_efficiency` is a constant fraction in (0,1].
Pumps default to initially off. Electrical efficiency curves are not supported
for pumps. Use [`CommitmentExclusion`](@ref) for alternative reversible modes.
"""
Pump(; initial_on=0, kwargs...) = Generator(; mode=:pump, initial_on, kwargs...)

"""
    CommitmentExclusion(; name, groups)

At most one group of named units may run in any interval. Several units in
the same group may run together. For a reversible machine, supply
`groups=[[:Turbine], [:Pump]]`. This restricts commitment, not signed flow.
"""
Base.@kwdef struct CommitmentExclusion
    name::Symbol
    groups::Vector{Vector{Symbol}}
end
Base.@kwdef struct RiverJunction
    name::Symbol
end
"""
    River(; name, target, curves, capacity, water_value, ...)

A routed reach with flow in m³/s, travel time in hours and inventory in Mm³.
`source` names a reservoir or explicit river confluence. A direct river
connection uses the upstream `target` to name the downstream reach, whose
`source` is `:auto`; `ScheduleCase` resolves these connections during construction. Plants and tunnels can feed a reach through `discharge_river`.
A shared source/target label also defines an internal zero-storage river node.
At a split, each outgoing reach supplies its `allocation` fraction of upstream
arrivals; fractions sum to one. Local inflow and outfalls enter after allocation.
Supply `deterministic_delay` and `curves=[]` for a fixed travel time, or
`curves` for distributed responses. Reservoir outlets use the selected
`law`; a `:delta` channel instead exchanges water instantaneously according to
one supplied delta-level response curve.
"""
Base.@kwdef struct River
    name::Symbol
    source::Symbol = :auto
    target::Symbol
    allocation::Float64 = 1.0 # prescribed share of arrivals at a zero-storage river node
    curves::Vector{RiverRouting.DelayCurve}
    capacity::Float64
    inflow::Float64 = 0.0 # natural water entering the top of this reach
    law::Symbol = :junction # :orifice, :weir, or :junction
    coefficient::Float64 = 0.0 # q=C*a*sqrt(H-crest), or C*(H-crest)^1.5
    crest::Float64 = 0.0
    min_arrival::Float64 = 0.0
    arrival_policy::Symbol = :interval_average # or :pointwise for exact cohort-wave minima
    arrival_window_grid::Vector{Float64} = Float64[] # empty selects scheduling grid; preserved on refinement
    deterministic_delay::Union{Nothing, Float64} = nothing
    gate_min::Float64 = 0.0 # optional sufficient environmental-release guard
    water_value::Float64 # inventory value per Mm³, explicit even at river junctions
    history_grid::Vector{Float64} = [-8.0, -6.0, -4.0, -2.0, 0.0]
    history_release::Vector{Float64} = zeros(4)
    discharge_curve::Union{Nothing, TableCurve} = nothing
    delta_level_ref_up_flow_curve::Union{Nothing, LevelFlowCurve} = nothing
    delta_level_ref_down_flow_curve::Union{Nothing, LevelFlowCurve} = nothing
    allow_dry::Bool = false
    ramp_up::Union{Nothing, Float64} = nothing # release (m³/s)/hour
    ramp_down::Union{Nothing, Float64} = nothing
    initial_release::Union{Nothing, Float64} = nothing # previous control-window average, separate from cohorts
    initial_interval_hours::Float64 = 1.0
end
"""
    HydroSystem(; reservoirs=[], plants=[], generators=[], ...)

Named physical objects in their array order. Omitted collections are empty;
include only the object types present in the watercourse. Names identify
connections. Solve and simulation arrays follow the order of these collections.
`generators` is the common, ordered unit collection for both `Generator` and
`Pump` constructors. The numerical algorithms do not mutate the supplied objects.
"""
Base.@kwdef struct HydroSystem
    reservoirs::Vector{Reservoir} = Reservoir[]
    junctions::Vector{Junction} = Junction[]
    boundaries::Vector{Boundary} = Boundary[]
    tunnels::Vector{Tunnel} = Tunnel[]
    tunnel_interlocks::Vector{TunnelInterlock} = TunnelInterlock[]
    plants::Vector{Plant} = Plant[]
    generators::Vector{Generator} = Generator[]
    commitment_exclusions::Vector{CommitmentExclusion} = CommitmentExclusion[]
    river_junctions::Vector{RiverJunction} = RiverJunction[]
    rivers::Vector{River} = River[]
end
"""Piecewise-constant operational input at absolute hour knots; final value is held."""
Base.@kwdef struct OperationalSeries
    object::Symbol
    attribute::Symbol
    times::Vector{Float64}
    values::Vector{Float64}
end
"""Hard minimum of a named flow observation, in m³/s.

Contributions are generator discharge and river release at the observation;
this records an operating rule and adds no water or hydraulic node. Delayed
arrivals must be constrained on the corresponding River instead.
"""
Base.@kwdef struct FlowRequirement
    name::Symbol
    generators::Vector{Symbol} = Symbol[]
    rivers::Vector{Symbol} = Symbol[]
    inflow::Float64 = 0.0
    min_flow::Float64 = 0.0
end
"""A bounded, explicitly penalized operating restriction.

Volume/level are assessed at `time`; power/discharge are averaged over
`[time,end_time]`. Times are absolute and survive refinement and restart.
Slacks use Mm³, metres, MW or m³/s respectively. Costs use currency/Mm³,
currency/metre, currency/MWh or currency/Mm³. Existing operating inputs and
physical domains are never softened by this record.
"""
Base.@kwdef struct SoftConstraint
    name::Symbol
    object::Symbol
    quantity::Symbol
    time::Float64
    end_time::Union{Nothing,Float64} = nothing
    lower::Union{Nothing,Float64} = nothing
    upper::Union{Nothing,Float64} = nothing
    max_lower_slack::Float64 = 0.0
    max_upper_slack::Float64 = 0.0
    lower_cost::Float64 = 0.0
    upper_cost::Float64 = 0.0
end

"""Spinning generation reserve capacity, MW, on a named group of units.

Prices and optional shortfall costs are currency/(MW·hour), paid on reserved
capacity. `response_hours` limits response through existing unit/plant ramps.
Positive `activation_hours` additionally requires qualified unit energy budgets
in MWh; it does not simulate activation or certify network delivery.
"""
Base.@kwdef struct ReserveGroup
    name::Symbol
    generators::Vector{Symbol}
    up_requirement::Float64 = 0.0
    down_requirement::Float64 = 0.0
    up_price::Float64 = 0.0
    down_price::Float64 = 0.0
    up_max::Union{Nothing,Float64} = nothing
    down_max::Union{Nothing,Float64} = nothing
    shortfall_cost::Union{Nothing,Float64} = nothing # nothing is a hard obligation
    response_hours::Float64 = 0.25
    activation_hours::Float64 = 0.0
end
"""
    ScheduleCase(; name, system, grid, prices, operations=[], flow_requirements=[], soft_constraints=[], reserve_groups=[])

A watercourse and its scheduling horizon. `grid` contains increasing absolute
hours; `prices[t]` is currency/MWh on `[grid[t], grid[t+1]]`. Prices and controls
have one column/value per interval. Storage has one column per grid vertex.
`operations` supplies time-dependent restrictions on named physical objects.
"""
Base.@kwdef struct ScheduleCase
    name::String
    system::HydroSystem
    grid::Vector{Float64}
    prices::Vector{Float64}
    operations::Vector{OperationalSeries} = OperationalSeries[]
    flow_requirements::Vector{FlowRequirement} = FlowRequirement[]
    soft_constraints::Vector{SoftConstraint} = SoftConstraint[]
    reserve_groups::Vector{ReserveGroup} = ReserveGroup[]
    function ScheduleCase(name::String, system::HydroSystem, grid::Vector{Float64},
        prices::Vector{Float64}, operations::Vector{OperationalSeries},
        flow_requirements::Vector{FlowRequirement}, soft_constraints::Vector{SoftConstraint},
        reserve_groups::Vector{ReserveGroup})
        new(name, normalize_river_connections(system), grid, prices, operations, flow_requirements,
            soft_constraints, reserve_groups)
    end
end
ScheduleCase(name, system::HydroSystem, grid, prices, operations, flow_requirements,
    soft_constraints=SoftConstraint[], reserve_groups=ReserveGroup[]) =
    ScheduleCase(String(name), system, Vector{Float64}(grid), Vector{Float64}(prices),
        Vector{OperationalSeries}(operations), Vector{FlowRequirement}(flow_requirements),
        Vector{SoftConstraint}(soft_constraints), Vector{ReserveGroup}(reserve_groups))

nodes(s) = vcat(
    [r.name for r in s.reservoirs],
    [r.name for r in s.junctions],
    [r.name for r in s.boundaries],
)
nodeindex(s) = Dict(n=>i for (i, n) in enumerate(nodes(s)))
plantof(s, g) = only(p for p in s.plants if p.name==g.plant)

function Base.show(io::IO, ::MIME"text/plain", c::ScheduleCase)
    s = c.system
    print(
        io,
        "ScheduleCase(",
        repr(c.name),
        ")\n",
        "  ",
        length(c.prices),
        " intervals, ",
        length(s.reservoirs),
        " reservoirs, ",
        length(s.generators),
        " generators\n",
        "  ",
        length(s.tunnels),
        " tunnels, ",
        length(s.rivers),
        " rivers",
    )
end

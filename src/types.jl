Base.@kwdef struct Reservoir
    name::Symbol
    z0::Float64=0.0
    slope::Float64=1.0
    curvature::Float64=0.0
    v0::Float64
    vmin::Float64
    vmax::Float64
    inflow::Float64=0.0
    water_value::Float64
    level_curve::Union{Nothing,TableCurve}=nothing
    volume_ramp_up::Union{Nothing,Float64}=nothing # Mm³/hour, between storage vertices
    volume_ramp_down::Union{Nothing,Float64}=nothing
    level_ramp_up::Union{Nothing,Float64}=nothing # metres/hour, between level vertices
    level_ramp_down::Union{Nothing,Float64}=nothing
end
Base.@kwdef struct Junction
    name::Symbol
    hmin::Float64=0.0
    hmax::Float64=700.0
end
Base.@kwdef struct Boundary
    name::Symbol
    head::Float64
end
Base.@kwdef struct Tunnel
    name::Symbol
    source::Symbol
    target::Symbol
    resistance::Float64
    capacity::Float64
    opening::Float64=1.0
    discharge_river::Union{Nothing,Symbol}=nothing # forward outfall; target remains hydraulic head reference
end
Base.@kwdef struct Plant
    name::Symbol
    source::Symbol
    target::Symbol
    pmax::Float64
    discharge_river::Union{Nothing,Symbol}=nothing # water destination, independent of outlet head
    ramp::Union{Nothing,Float64}=100.0 # MW/hour; nothing disables the symmetric default
    initial_power::Union{Nothing,Float64}=nothing
    initial_interval_hours::Float64=1.0
    tailwater_curve::Union{Nothing,TableCurve}=nothing
    outlet_head_floor::Union{Nothing,Float64}=nothing
    pmin::Float64=0.0
    qmin::Float64=0.0
    qmax::Union{Nothing,Float64}=nothing
    ramp_up::Union{Nothing,Float64}=nothing # nothing inherits the symmetric power ramp
    ramp_down::Union{Nothing,Float64}=nothing
    discharge_ramp_up::Union{Nothing,Float64}=nothing # (m³/s)/hour
    discharge_ramp_down::Union{Nothing,Float64}=nothing
    initial_discharge::Union{Nothing,Float64}=nothing
    minup::Float64=0.0
    mindown::Float64=0.0
    initial_on::Union{Nothing,Int}=nothing # otherwise any initially running member unit
    initial_age::Union{Nothing,Float64}=nothing
end
Base.@kwdef struct Generator
    name::Symbol
    plant::Symbol
    qmin::Float64
    qmax::Float64
    pmin::Float64
    pmax::Float64
    efficiency::Float64=0.94
    min_efficiency::Float64=0.0
    qbest::Float64=1.0
    qcurvature::Float64=0.28
    hbest::Float64=1.0
    hcurvature::Float64=0.015
    hmin::Float64
    hmax::Float64
    initial_on::Int=1
    initial_age::Float64=8.0
    minup::Float64=3.0
    mindown::Float64=2.0
    startup::Float64=120.0
    shutdown::Float64=0.0
    turbine_table::Union{Nothing,TurbineTable}=nothing
    generator_efficiency_curve::Union{Nothing,TableCurve}=nothing
    ramp_up::Union{Nothing,Float64}=nothing
    ramp_down::Union{Nothing,Float64}=nothing
    discharge_ramp_up::Union{Nothing,Float64}=nothing
    discharge_ramp_down::Union{Nothing,Float64}=nothing
    initial_power::Union{Nothing,Float64}=nothing
    initial_discharge::Union{Nothing,Float64}=nothing
    initial_interval_hours::Float64=1.0
end
Base.@kwdef struct RiverJunction
    name::Symbol
end
Base.@kwdef struct River
    name::Symbol
    source::Symbol=:auto
    target::Symbol
    curves::Vector{RiverRouting.DelayCurve}
    capacity::Float64
    inflow::Float64=0.0 # natural water entering the top of this reach
    law::Symbol=:junction # :orifice, :weir, or :junction
    coefficient::Float64=0.0 # q=C*a*sqrt(H-crest), or C*(H-crest)^1.5
    crest::Float64=0.0
    min_arrival::Float64=0.0
    arrival_policy::Symbol=:interval_average # or :pointwise for exact cohort-wave minima
    arrival_window_grid::Vector{Float64}=Float64[] # empty selects scheduling grid; preserved on refinement
    deterministic_delay::Union{Nothing,Float64}=nothing
    gate_min::Float64=0.0 # optional sufficient environmental-release guard
    water_value::Float64 # inventory value per Mm³, explicit even at river junctions
    history_grid::Vector{Float64}=[-8.0, -6.0, -4.0, -2.0, 0.0]
    history_release::Vector{Float64}=zeros(4)
    discharge_curve::Union{Nothing,TableCurve}=nothing
    allow_dry::Bool=false
    ramp_up::Union{Nothing,Float64}=nothing # release (m³/s)/hour
    ramp_down::Union{Nothing,Float64}=nothing
    initial_release::Union{Nothing,Float64}=nothing # previous control-window average, separate from cohorts
    initial_interval_hours::Float64=1.0
end
Base.@kwdef struct HydroSystem
    reservoirs::Vector{Reservoir}
    junctions::Vector{Junction}
    boundaries::Vector{Boundary}
    tunnels::Vector{Tunnel}
    plants::Vector{Plant}
    generators::Vector{Generator}
    river_junctions::Vector{RiverJunction}=RiverJunction[]
    rivers::Vector{River}
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
    generators::Vector{Symbol}=Symbol[]
    rivers::Vector{Symbol}=Symbol[]
    inflow::Float64=0.0
    min_flow::Float64=0.0
end
Base.@kwdef struct ScheduleCase
    name::String
    system::HydroSystem
    grid::Vector{Float64}
    prices::Vector{Float64}
    operations::Vector{OperationalSeries}=OperationalSeries[]
    flow_requirements::Vector{FlowRequirement}=FlowRequirement[]
end
nodes(s) = vcat(
    [r.name for r in s.reservoirs],
    [r.name for r in s.junctions],
    [r.name for r in s.boundaries],
)
nodeindex(s) = Dict(n=>i for (i, n) in enumerate(nodes(s)))
plantof(s, g) = only(p for p in s.plants if p.name==g.plant)

"""Operating-bound exceedances from a physical hydraulic step, in Mm³ and m³/s."""
struct HydraulicBoundViolations
    storage::Float64
    tunnel_flow::Float64
    river_flow::Float64
end

"""One independently integrated midpoint step.

`storage` is endpoint storage in Mm³, `head` is midpoint head in metres,
`power` is MW, and all flow arrays are m³/s. Each step owns its output arrays;
later integration steps do not overwrite them. Operating-bound exceedances
are reported separately from convergence of the physical equations.
"""
struct HydraulicStep
    storage::Vector{Float64}
    head::Vector{Float64}
    tunnel_q::Vector{Float64}
    tunnel_opening::Vector{Float64}
    river_release::Vector{Float64}
    river_injection::Vector{Float64}
    river_arrival::Vector{Float64}
    generator_head::Vector{Float64}
    generator_q::Vector{Float64}
    power::Vector{Float64}
    residual::Float64
    converged::Bool
    iterations::Int
    boundary_net_inflow::Vector{Float64}
    routing_domain_valid::Bool
    bound_violations::HydraulicBoundViolations
    law_domain_valid::Bool
end

"""Chronological physical trajectory on a specified time grid.

Columns of `storage` and `transit` are grid vertices; all other matrices use
interval columns. Storage and transported volumes are Mm³, heads are metres,
power is nonnegative production/consumption in MW, and releases, injections and
boundary outflows are m³/s. Generating and pumping units share the unit axis.
"""
struct HydraulicTrajectory
    grid::Vector{Float64}
    storage::Matrix{Float64}
    head::Matrix{Float64}
    tunnel_q::Matrix{Float64}
    tunnel_opening::Matrix{Float64}
    generator_q::Matrix{Float64}
    gate::Matrix{Float64}
    power::Matrix{Float64}
    river_release::Matrix{Float64}
    river_injection::Matrix{Float64}
    arrival_volume::Matrix{Float64}
    transit::Matrix{Float64}
    boundary_outflow::Vector{Float64}
    converged::Bool
    residual::Float64
    steps::Vector{HydraulicStep}
end

# Read-only access for the original public string keys. Numerical code uses
# fields directly; these records are not dictionaries or mutable key stores.
function Base.getindex(x::HydraulicStep, key::AbstractString)
    field = key == "Vnew" ? :storage : key == "H" ? :head : Symbol(key)
    hasfield(HydraulicStep, field) || throw(KeyError(key))
    getfield(x, field)
end
function Base.getindex(x::HydraulicTrajectory, key::AbstractString)
    field = key == "V" ? :storage : key == "H" ? :head : Symbol(key)
    hasfield(HydraulicTrajectory, field) || throw(KeyError(key))
    getfield(x, field)
end
function Base.getindex(x::HydraulicBoundViolations, key::AbstractString)
    field = Symbol(key)
    hasfield(HydraulicBoundViolations, field) || throw(KeyError(key))
    getfield(x, field)
end

function Base.:(==)(a::T, b::T) where {
    T<:Union{HydraulicBoundViolations, HydraulicStep, HydraulicTrajectory},
}
    all(getfield(a, i) == getfield(b, i) for i in 1:fieldcount(T))
end

# Serialization retains the original public spellings and row-oriented
# matrices. Converting to dictionaries belongs only at this I/O boundary.
function jsonready(x::HydraulicBoundViolations)
    jsonready(Dict(
        "storage"=>x.storage,
        "tunnel_flow"=>x.tunnel_flow,
        "river_flow"=>x.river_flow,
    ))
end

function jsonready(x::HydraulicStep)
    jsonready(Dict{String,Any}(
        "Vnew"=>x.storage,
        "H"=>x.head,
        "tunnel_q"=>x.tunnel_q,
        "tunnel_opening"=>x.tunnel_opening,
        "river_release"=>x.river_release,
        "river_injection"=>x.river_injection,
        "river_arrival"=>x.river_arrival,
        "generator_head"=>x.generator_head,
        "generator_q"=>x.generator_q,
        "power"=>x.power,
        "residual"=>x.residual,
        "converged"=>x.converged,
        "iterations"=>x.iterations,
        "boundary_net_inflow"=>x.boundary_net_inflow,
        "routing_domain_valid"=>x.routing_domain_valid,
        "bound_violations"=>x.bound_violations,
        "law_domain_valid"=>x.law_domain_valid,
    ))
end

function jsonready(x::HydraulicTrajectory)
    jsonready(Dict{String,Any}(
        "grid"=>x.grid,
        "V"=>x.storage,
        "H"=>x.head,
        "tunnel_q"=>x.tunnel_q,
        "tunnel_opening"=>x.tunnel_opening,
        "generator_q"=>x.generator_q,
        "gate"=>x.gate,
        "power"=>x.power,
        "river_release"=>x.river_release,
        "arrival_volume"=>x.arrival_volume,
        "transit"=>x.transit,
        "boundary_outflow"=>x.boundary_outflow,
        "converged"=>x.converged,
        "residual"=>x.residual,
        "steps"=>x.steps,
    ))
end

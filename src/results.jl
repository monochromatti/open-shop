"""
    Schedule

A physical schedule in the case's object order. Storage has one column per
grid vertex; rates, heads, power and commitment have one column per interval.
Storage and river inventory are Mm³, discharge is m³/s, head is metres and
power is MW. `diagnostics` contains validation and reporting data, not equations.
"""
struct Schedule{U<:AbstractMatrix}
    grid::Vector{Float64}
    commitment::U
    storage::Matrix{Float64}
    head::Matrix{Float64}
    tunnel_q::Matrix{Float64}
    tunnel_opening::Matrix{Float64}
    generator_q::Matrix{Float64}
    power::Matrix{Float64}
    gate::Matrix{Float64}
    river_release::Matrix{Float64}
    arrival_volume::Matrix{Float64}
    terminal_transit::Vector{Float64}
    shortfall_release::Union{Nothing,Matrix{Float64}}
    objective::Float64
    gate_adjustment_cost::Union{Nothing,Float64}
    release_penalty_cost::Union{Nothing,Float64}
    diagnostics::Dict{String,Any}
end

# These aliases are an input/output compatibility boundary. Numerical routines
# access the fields above directly.
const _SCHEDULE_KEYS = (
    "grid"=>:grid, "u"=>:commitment, "V"=>:storage, "H"=>:head,
    "tunnel_q"=>:tunnel_q, "tunnel_opening"=>:tunnel_opening,
    "generator_q"=>:generator_q, "power"=>:power, "gate"=>:gate,
    "river_release"=>:river_release, "arrival_volume"=>:arrival_volume,
    "terminal_transit"=>:terminal_transit, "shortfall_release"=>:shortfall_release,
    "objective"=>:objective, "gate_adjustment_cost"=>:gate_adjustment_cost,
    "release_penalty_cost"=>:release_penalty_cost,
)
function Base.getindex(x::Schedule, key::AbstractString)
    for (name, field) in _SCHEDULE_KEYS
        key == name && return getfield(x, field)
    end
    x.diagnostics[key]
end
function Base.haskey(x::Schedule,key::AbstractString)
    for (name,field) in _SCHEDULE_KEYS
        key==name && return getfield(x,field)!==nothing
    end
    haskey(x.diagnostics,key)
end
Base.get(x::Schedule, key::AbstractString, default) = haskey(x, key) ? x[key] : default

"""Return the serializable schedule profile, with independent copies of arrays."""
function schedule_dict(x::Schedule)
    result = deepcopy(x.diagnostics)
    for (key, field) in _SCHEDULE_KEYS
        value = getfield(x, field)
        value===nothing && continue
        result[key] = value isa AbstractArray ? copy(value) : value
    end
    result
end

# JSON encodes a matrix as rows. Empty matrices need the interval count because
# an empty JSON array cannot express its number of columns.
function _schedule_matrix(x, columns; control=false)
    A=if x isa AbstractMatrix
        Matrix{Float64}(x)
    else
        x isa AbstractVector || throw(ArgumentError("schedule arrays must be matrices or JSON rows"))
        isempty(x) && return zeros(0,columns)
        all(row->row isa AbstractVector,x) || throw(ArgumentError("schedule matrix needs rows"))
        permutedims(reduce(hcat,Vector{Float64}.(x)))
    end
    if control
        unchanged=x isa AbstractMatrix ?
            all(i->iszero(x[i])==iszero(A[i]),eachindex(A)) :
            all(iszero(value)==iszero(A[i,j]) for (i,row) in enumerate(x) for (j,value) in enumerate(row))
        unchanged || throw(ArgumentError("nonzero controls must remain nonzero in Float64"))
    end
    A
end

_schedule_input(x::Schedule; case=nothing) = x
function _schedule_input(x::AbstractDict; case=nothing)
    grid = Float64.(get(x, "grid", case === nothing ? Float64[] : case.grid))
    T = length(grid)-1
    T >= 0 || throw(ArgumentError("schedule needs a grid"))
    matrix(key) = _schedule_matrix(x[key], key == "V" ? T+1 : T;
        control=key in ("u","generator_q","gate","tunnel_opening"))
    # Preserve noninteger commitment input for validation. This conversion does
    # not round unit states or repair malformed schedules.
    commitment=if x["u"] isa AbstractMatrix
        Matrix(x["u"])
    else
        rows=x["u"]
        rows isa AbstractVector && all(row -> row isa AbstractVector, rows) ||
            throw(ArgumentError("commitment matrix needs rows"))
        isempty(rows) ? zeros(Int,0,T) : permutedims(reduce(hcat,rows))
    end
    release = matrix("river_release")
    opening = haskey(x, "tunnel_opening") ? matrix("tunnel_opening") :
        case === nothing ? throw(ArgumentError("schedule needs tunnel_opening")) :
        _tunnel_opening_input(case,nothing)
    diagnostics = Dict{String,Any}(string(k)=>v for (k,v) in x
        if !any(p -> first(p) == string(k), _SCHEDULE_KEYS))
    Schedule(grid, commitment, matrix("V"), matrix("H"), matrix("tunnel_q"),
        opening, matrix("generator_q"), matrix("power"), matrix("gate"), release,
        matrix("arrival_volume"), Vector{Float64}(x["terminal_transit"]),
        haskey(x,"shortfall_release") ? matrix("shortfall_release") : nothing,
        Float64(x["objective"]),
        haskey(x,"gate_adjustment_cost") ? Float64(x["gate_adjustment_cost"]) : nothing,
        haskey(x,"release_penalty_cost") ? Float64(x["release_penalty_cost"]) : nothing, diagnostics)
end

"""
    SolveResult

Optimization status, a physical `solution` (or `nothing`), and optional global
bounds. `accepted` describes the checks performed by the called solver;
`global_certificate` additionally requires an enclosing global bound.
Detailed solver work counters and validation reports are in `diagnostics`.
"""
struct SolveResult{S}
    status::String
    accepted::Bool
    global_certificate::Bool
    solution::S
    global_bound::Union{Nothing,Float64}
    relative_gap::Union{Nothing,Float64}
    absolute_gap::Union{Nothing,Float64}
    total_seconds::Float64
    diagnostics::Dict{String,Any}
end

const _SOLVE_KEYS = ("status", "accepted", "global_certificate", "solution",
    "global_bound", "relative_gap", "absolute_gap", "total_seconds")
function SolveResult(report::AbstractDict; solution=get(report,"solution",nothing))
    schedule = solution === nothing ? nothing : _schedule_input(solution)
    diagnostics = Dict{String,Any}(string(k)=>v for (k,v) in report
        if !(string(k) in _SOLVE_KEYS) && !(string(k) in ("objective","feasible_lower_bound")) &&
            !any(p -> first(p)==string(k),_SCHEDULE_KEYS))
    SolveResult(string(report["status"]), get(report,"accepted",false),
        get(report,"global_certificate",false), schedule, get(report,"global_bound",nothing),
        get(report,"relative_gap",nothing), get(report,"absolute_gap",nothing),
        get(report,"total_seconds",0.0), diagnostics)
end

function Base.getproperty(x::SolveResult, field::Symbol)
    if field in (:objective, :feasible_lower_bound)
        solution = getfield(x,:solution)
        solution===nothing && return nothing
        field==:feasible_lower_bound &&
            !get(get(solution.diagnostics,"validation",Dict()),"valid",false) && return nothing
        return solution.objective
    end
    getfield(x,field)
end
Base.propertynames(x::SolveResult) = (fieldnames(typeof(x))..., :objective, :feasible_lower_bound)
function Base.getindex(x::SolveResult, key::AbstractString)
    key in _SOLVE_KEYS && return getfield(x,Symbol(key))
    key in ("objective","feasible_lower_bound") && return getproperty(x,Symbol(key))
    haskey(x.diagnostics,key) && return x.diagnostics[key]
    x.solution === nothing && throw(KeyError(key))
    x.solution[key]
end
function Base.haskey(x::SolveResult, key::AbstractString)
    key == "objective" && return x.solution !== nothing
    key == "feasible_lower_bound" && return true
    key in _SOLVE_KEYS || haskey(x.diagnostics,key) ||
        (x.solution !== nothing && haskey(x.solution,key))
end
Base.get(x::SolveResult, key::AbstractString, default) = haskey(x,key) ? x[key] : default
_schedule_input(x::SolveResult; case=nothing) =
    x.solution === nothing ? throw(ArgumentError("result has no physical solution")) : x.solution

function jsonready(x::Schedule)
    jsonready(schedule_dict(x))
end
function jsonready(x::SolveResult)
    report = copy(x.diagnostics)
    for key in _SOLVE_KEYS
        report[key] = x[key]
    end
    report["feasible_lower_bound"]=x.feasible_lower_bound
    x.solution===nothing || (report["objective"]=x.objective)
    jsonready(report)
end

function Base.show(io::IO, ::MIME"text/plain", x::Schedule)
    print(io,"Schedule(",length(x.grid)-1," intervals, ",size(x.commitment,1),
        " units; objective = ",x.objective,")")
end
function Base.show(io::IO, ::MIME"text/plain", x::SolveResult)
    print(io,"SolveResult(",x.status,"; accepted = ",x.accepted,
        ", global_certificate = ",x.global_certificate)
    x.solution === nothing || print(io,", objective = ",x.objective)
    print(io,")")
end

Base.:(==)(a::Schedule,b::Schedule) =
    all(field -> getfield(a,field)==getfield(b,field),fieldnames(typeof(a)))

schedule_dict(x::SolveResult) = schedule_dict(_schedule_input(x))

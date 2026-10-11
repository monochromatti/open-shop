# One interpretation of prescribed and selectable tunnel gates for every model
# and independent physical audit. Opening zero is closed; every positive listed
# position counts as open, regardless of its hydraulic resistance.

tunnel_initial_opening(e::Tunnel) =
    e.initial_opening===nothing ? e.opening : e.initial_opening

"""Index of a declared gate position, or `nothing` when no position matches."""
function gate_position_index(e::Tunnel,x; tolerance=0)
    isfinite(tolerance) && tolerance>=0 ||
        throw(ArgumentError("gate-position tolerance must be finite and nonnegative"))
    x isa Real && isfinite(x) || return nothing
    positions=e.gate_positions
    isempty(positions) && return nothing
    j=searchsortedfirst(positions,x)
    j<=length(positions) && positions[j]==x && return j
    best=nothing
    distance=Inf
    for k in (j-1,j)
        1<=k<=length(positions) || continue
        d=abs(positions[k]-x)
        if d<=tolerance && d<distance
            best,distance=k,d
        end
    end
    best
end

function _validate_tunnel_positions(e)
    all(isfinite,e.gate_positions) && all(x->0<=x<=1,e.gate_positions) &&
        all(diff(e.gate_positions).>0) ||
        throw(ArgumentError("$(e.name): gate positions must be finite, unique and increasing in [0,1]"))
    isfinite(tunnel_initial_opening(e)) && 0<=tunnel_initial_opening(e)<=1 ||
        throw(ArgumentError("$(e.name): initial opening must lie in [0,1]"))
    isfinite(e.gate_adjustment_cost) && e.gate_adjustment_cost>=0 ||
        throw(ArgumentError("$(e.name): gate adjustment cost must be nonnegative"))
    if !isempty(e.gate_positions)
        gate_position_index(e,e.opening)!==nothing &&
            gate_position_index(e,tunnel_initial_opening(e))!==nothing ||
            throw(ArgumentError("$(e.name): static and initial openings must be declared gate positions"))
    end
    nothing
end

"""Permitted positions on an interval; an active opening series is a hard schedule.

An optimized gate remains free before the first series knot. Prescribed gates
retain the existing interval opening, including its static fallback.
"""
function tunnel_allowed_positions(c,e::Tunnel,t)
    isempty(e.gate_positions) &&
        return [opinterval(c,e.name,:opening,t,e.opening)]
    scheduled=optional_operation(c,e.name,:opening,t)
    scheduled===nothing ? e.gate_positions : [scheduled]
end

"""Complete default opening trajectory in physical gate-position units.

For an optimized gate this is an initial witness, not a restriction before an
opening schedule begins. Explicit optimized controls may select any listed
position on those free intervals.
"""
function tunnel_opening_schedule(c)
    [isempty(e.gate_positions) ? opinterval(c,e.name,:opening,t,e.opening) :
        opvalue(c,e.name,:opening,c.grid[t],e.opening)
        for e in c.system.tunnels,t in eachindex(c.prices)]
end

function _tunnel_interlock_indices(c)
    index=Dict(e.name=>i for (i,e) in enumerate(c.system.tunnels))
    [begin
        length(unique(rule.tunnels))==length(rule.tunnels) ||
            throw(ArgumentError("duplicate tunnel in interlock $(rule.name)"))
        all(n->haskey(index,n),rule.tunnels) ||
            throw(ArgumentError("unknown tunnel in interlock $(rule.name)"))
        [index[n] for n in rule.tunnels]
    end for rule in c.system.tunnel_interlocks]
end

function _tunnel_interlock_bounds(c,rule,t)
    lo=opvalue(c,rule.name,:min_open,c.grid[t],rule.min_open)
    hi=opvalue(c,rule.name,:max_open,c.grid[t],rule.max_open)
    isinteger(lo) && isinteger(hi) && 0<=lo<=hi<=length(rule.tunnels) ||
        throw(ArgumentError("invalid open-count bounds for interlock $(rule.name)"))
    Int(lo),Int(hi)
end

function _validate_tunnel_operating_controls(c)
    indices=_tunnel_interlock_indices(c)
    for (j,rule) in enumerate(c.system.tunnel_interlocks),t in eachindex(c.prices)
        lo,hi=_tunnel_interlock_bounds(c,rule,t)
        positions=[tunnel_allowed_positions(c,c.system.tunnels[i],t) for i in indices[j]]
        forced_open=count(p->all(>(0),p),positions)
        possible_open=count(p->any(>(0),p),positions)
        forced_open<=hi && possible_open>=lo ||
            throw(ArgumentError("fixed tunnel controls conflict with interlock $(rule.name) at $(c.grid[t])"))
    end
    nothing
end

"""Check complete physical openings against schedules, positions and interlocks.

Tolerance applies to declared-position and prescribed-schedule comparisons.
The physical domain remains [0,1], and only exact zero counts as closed.
"""
function validate_tunnel_openings(c,A;tolerance=1e-8)
    isfinite(tolerance) && tolerance>=0 ||
        throw(ArgumentError("gate-position tolerance must be finite and nonnegative"))
    A isa AbstractMatrix && size(A)==(length(c.system.tunnels),length(c.prices)) ||
        throw(DimensionMismatch("tunnel openings need one value per tunnel and control interval"))
    for (i,e) in enumerate(c.system.tunnels),t in eachindex(c.prices)
        x=A[i,t]
        x isa Real && isfinite(x) && 0<=x<=1 ||
            throw(ArgumentError("$(e.name): opening must be finite and within [0,1]"))
        if isempty(e.gate_positions)
            prescribed=opinterval(c,e.name,:opening,t,e.opening)
            abs(x-prescribed)<=tolerance ||
                throw(ArgumentError("$(e.name): opening differs from its prescribed schedule"))
            (x==0)==(prescribed==0) ||
                throw(ArgumentError("$(e.name): tolerance cannot change whether a gate is closed"))
        else
            position=gate_position_index(e,x;tolerance)
            position!==nothing ||
                throw(ArgumentError("$(e.name): opening is not a declared gate position"))
            (x==0)==(e.gate_positions[position]==0) ||
                throw(ArgumentError("$(e.name): tolerance cannot change whether a gate is closed"))
            scheduled=optional_operation(c,e.name,:opening,t)
            scheduled===nothing || abs(x-scheduled)<=tolerance ||
                throw(ArgumentError("$(e.name): opening differs from its hard schedule"))
        end
    end
    indices=_tunnel_interlock_indices(c)
    for (j,rule) in enumerate(c.system.tunnel_interlocks),t in eachindex(c.prices)
        lo,hi=_tunnel_interlock_bounds(c,rule,t)
        count_open=count(i->A[i,t]>0,indices[j])
        lo<=count_open<=hi ||
            throw(ArgumentError("tunnel openings violate interlock $(rule.name) at $(c.grid[t])"))
    end
    true
end

"""Cost of gate changes at interval starts; no duration multiplier is applied."""
function tunnel_adjustment_cost(c,A)
    A isa AbstractMatrix && size(A)==(length(c.system.tunnels),length(c.prices)) ||
        throw(DimensionMismatch("tunnel openings need one value per tunnel and control interval"))
    all(x->x isa Real && isfinite(x) && 0<=x<=1,A) ||
        throw(ArgumentError("tunnel openings must be finite and within [0,1]"))
    cost=0.0
    for (i,e) in enumerate(c.system.tunnels)
        previous=tunnel_initial_opening(e)
        for t in eachindex(c.prices)
            coefficient=opvalue(c,e.name,:gate_adjustment_cost,c.grid[t],e.gate_adjustment_cost)
            cost+=coefficient*abs(A[i,t]-previous)
            previous=A[i,t]
        end
    end
    cost
end

"""Compile selectable positions with K−1 state variables and shared interlocks.

Position one is implicit. Fixed, prescribed, scheduled and singleton controls
stay numeric; a fixed trajectory introduces no gate variables. Positive costs
use an absolute-change epigraph only when the change is an affine expression.
"""
function _tunnel_gate_states!(m,c;fixed=nothing,relaxed=false)
    fixed===nothing || (fixed=_physical_tunnel_openings(c,fixed))
    E,T=length(c.system.tunnels),length(c.prices)
    opening=Matrix{Any}(undef,E,T)
    open=similar(opening)
    weights=Dict{Tuple{Int,Int},Vector{Any}}()
    changes=Dict{Tuple{Int,Int},VariableRef}()
    cost=AffExpr(0.0)
    for (i,e) in enumerate(c.system.tunnels)
        previous=tunnel_initial_opening(e)
        for t in 1:T
            positions=fixed===nothing ? tunnel_allowed_positions(c,e,t) : [Float64(fixed[i,t])]
            if length(positions)==1
                opening[i,t]=only(positions)
                open[i,t]=Float64(only(positions)>0)
            else
                selectors=[@variable(m,lower_bound=0,upper_bound=1,binary=!relaxed,
                    base_name="gate_select_$(i)_$(t)_$(k)") for k in 2:length(positions)]
                implicit=1-sum(selectors)
                @constraint(m,implicit>=0)
                weights[(i,t)]=Any[implicit;selectors]
                opening[i,t]=positions[1]+sum((positions[k]-positions[1])*selectors[k-1]
                    for k in 2:length(positions))
                open[i,t]=positions[1]>0 ? 1.0 : sum(selectors)
            end
            coefficient=opvalue(c,e.name,:gate_adjustment_cost,c.grid[t],e.gate_adjustment_cost)
            if coefficient>0
                difference=opening[i,t]-previous
                if difference isa Real
                    add_to_expression!(cost,coefficient*abs(difference))
                else
                    change=@variable(m,lower_bound=0,upper_bound=1,
                        base_name="gate_change_$(i)_$(t)")
                    @constraint(m,change>=difference)
                    @constraint(m,change>=-difference)
                    changes[(i,t)]=change
                    add_to_expression!(cost,coefficient,change)
                end
            end
            previous=opening[i,t]
        end
    end
    indices=_tunnel_interlock_indices(c)
    for (j,rule) in enumerate(c.system.tunnel_interlocks),t in 1:T
        lo,hi=_tunnel_interlock_bounds(c,rule,t)
        count_open=sum(open[i,t] for i in indices[j];init=0.0)
        if count_open isa Real
            lo<=count_open<=hi || throw(ArgumentError(
                "fixed tunnel controls conflict with interlock $(rule.name) at $(c.grid[t])"))
        else
            @constraint(m,lo<=count_open<=hi)
        end
    end
    (;opening,open,weights,cost,changes,positions=[e.gate_positions for e in c.system.tunnels])
end

"""Extract exact positions from integral selectors, independent of position spacing.

The implicit first selector is selected only when all explicit selectors round
to zero. Position-valued expressions are never used to infer a discrete state.
"""
function _tunnel_gate_solution(data,value_function=JuMP.value;tolerance=1e-6)
    isfinite(tolerance) && 0<=tolerance<0.5 ||
        throw(ArgumentError("invalid gate-selector integrality tolerance"))
    opening=Matrix{Float64}(value_function.(data.opening))
    for ((i,t),weights) in data.weights
        selected=1
        for k in 2:length(weights)
            x=value_function(weights[k])
            isfinite(x) && abs(x-round(x))<=tolerance && round(x) in (0,1) ||
                throw(ArgumentError("fractional tunnel gate selector"))
            if round(x)==1
                selected==1 || throw(ArgumentError("multiple selected tunnel gate positions"))
                selected=k
            end
        end
        implicit=value_function(first(weights))
        isfinite(implicit) && abs(implicit-Float64(selected==1))<=tolerance ||
            throw(ArgumentError("inconsistent implicit tunnel gate selector"))
        opening[i,t]=data.positions[i][selected]
    end
    opening
end

"""Lift a physical opening trajectory into selector and change variables."""
function _tunnel_gate_values(data,c,A)
    validate_tunnel_openings(c,A;tolerance=0)
    for i in axes(A,1),t in axes(A,2)
        fixed=data.opening[i,t]
        fixed isa Real && A[i,t]!=fixed &&
            throw(ArgumentError("opening trajectory differs from the compiled fixed controls"))
    end
    assigned=Dict{VariableRef,Float64}()
    for ((i,t),weights) in data.weights
        position=gate_position_index(c.system.tunnels[i],A[i,t])
        for k in 2:length(weights)
            assigned[weights[k]]=Float64(k==position)
        end
    end
    for ((i,t),variable) in data.changes
        previous=t==1 ? tunnel_initial_opening(c.system.tunnels[i]) : A[i,t-1]
        assigned[variable]=abs(A[i,t]-previous)
    end
    assigned
end

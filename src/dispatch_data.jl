# Numerical data and routing expressions for local dispatch construction.

"""Prepare capacity and gate bounds on the case's control intervals."""
function _dispatch_outlet_bounds(c)
    D, T = length(c.system.rivers), length(c.prices)
    capacity = Matrix{Float64}(undef, D, T)
    gate_lower = similar(capacity)
    gate_upper = similar(capacity)
    for (i, r) in enumerate(c.system.rivers), t in 1:T
        capacity[i, t] = opinterval(c, r.name, :capacity, t, r.capacity)
        gate_lower[i, t] = opinterval(c, r.name, :gate_min, t, r.gate_min)
        gate_upper[i, t] = opinterval(c, r.name, :gate_max, t, 1.0)
    end
    (; capacity, gate_lower, gate_upper)
end

function _dispatch_transport_expressions(m,c,release;transport=nothing,
    injections=river_injections(c,release),coordinate=(r,q,i,k)->nothing,distributed=nothing)
    exact=all(r->r.deterministic_delay!==nothing,c.system.rivers)
    !exact && transport!==nothing && throw(ArgumentError(
        "compiled deterministic transport cannot serve distributed dispatch"))
    nd=exact ? _transport_data(c,transport) : nothing
    rd=exact ? nothing : distributed===nothing ? DistributedRouting(c) : distributed
    D,T=size(release)
    arrivals=Matrix{Union{AffExpr,QuadExpr,NonlinearExpr}}(undef,D,T)
    terminal=Vector{Union{AffExpr,QuadExpr,NonlinearExpr}}(undef,D)
    rates=Vector{Vector{Any}}(undef,D)
    coordinates=Vector{Vector{Any}}(undef,D)
    starts=Pair{VariableRef,Union{AffExpr,QuadExpr,NonlinearExpr}}[]
    if exact
        for i in 1:D,t in 1:T
            arrivals[i,t]=AffExpr(nd.arrival_history[i,t])+sum(
                v*injections[d,k] for (d,k,v) in nd.arrival_terms[i,t];init=0.0)
        end
        for i in 1:D
            terminal[i]=AffExpr(nd.terminal_history[i])+sum(
                v*injections[d,k] for (d,k,v) in nd.terminal_terms[i];init=0.0)
        end
    else
        K=length(rd.grid)-1
        fine_arrival=Matrix{Union{AffExpr,QuadExpr,NonlinearExpr}}(undef,D,K)
        reservoirs=river_source_indices(c.system).reservoirs
        for i in river_order(c)
            r=c.system.rivers[i]
            if rd.source_cohorts[i]
                rates[i]=Any[reservoirs[i] ? release[i,t] : injections[i,t] for t in 1:T]
            else
                q=@variable(m,[1:K],lower_bound=0.0,base_name="routing_release_$i")
                rates[i]=Any[q[k] for k in 1:K]
                upstream=findall(e->e.target==r.source,c.system.rivers)
                for k in 1:K
                    t=rd.control[k]
                    set_upper_bound(q[k],opinterval(c,r.name,:capacity,t,r.capacity))
                    duration=rd.grid[k+1]-rd.grid[k]
                    incoming=r.allocation*sum(fine_arrival[d,k] for d in upstream;init=0.0)
                    @constraint(m,(0.0036duration*(q[k]-injections[i,t])-incoming)/0.3==0)
                    push!(starts,q[k]=>injections[i,t]+incoming/(0.0036duration))
                end
            end
            coordinates[i]=Any[coordinate(r,q,i,k) for (k,q) in enumerate(rates[i])]
            cohort_dt=diff(rd.cohort_grid[i])
            B=rd.coefficients[i]
            for k in 1:K
                fine_arrival[i,k]=AffExpr(rd.history_arrival[i,k])+sum(
                    0.0036cohort_dt[j]*river_transfer_expression(m,r,rates[i][j],
                        collect(view(B,k,j,1:river_curve_count(r)));coordinate=coordinates[i][j])
                    for j in eachindex(rates[i]) if any(!iszero,view(B,k,j,:));init=0.0)
            end
            for t in 1:T
                arrivals[i,t]=sum(fine_arrival[i,k]
                    for k in rd.vertices[t]:(rd.vertices[t+1]-1);init=AffExpr(0.0))
            end
            terminal[i]=AffExpr(rd.history_terminal[i])+sum(
                0.0036cohort_dt[j]*river_transfer_expression(m,r,rates[i][j],
                    [1-sum(view(B,:,j,l)) for l in 1:river_curve_count(r)];coordinate=coordinates[i][j])
                for j in eachindex(rates[i]);init=0.0)
        end
        m.ext[:distributed_routing]=(;mesh=rd,rates)
    end
    routing=(;exact,nd,rd,arrivals,terminal,rates,coordinates,injections,starts)
    constrain_river_stream_capacities!(m,c,routing)
    routing
end

function routing_volume_expression(m,c,routing,i,a,b)
    if routing.exact
        expression=_pulse_expression(routing.nd.arrival_pulses[i],a,b)
        return expression.history+sum(v*routing.injections[d,k]
            for (d,k,v) in expression.terms;init=0.0)
    end
    r=c.system.rivers[i]
    grid=routing.rd.cohort_grid[i]
    B=transfer_coefficients(r,[a,b],grid)
    history=only(route_volumes(r,[a,b],r.history_grid,r.history_release))
    history+sum(0.0036*(grid[k+1]-grid[k])*river_transfer_expression(m,r,routing.rates[i][k],
        collect(view(B,1,k,1:river_curve_count(r)));coordinate=routing.coordinates[i][k])
        for k in eachindex(routing.rates[i]) if any(!iszero,view(B,1,k,:));init=0.0)
end

function routing_point_expression(m,c,routing,i,time;side)
    if routing.exact
        expression=deterministic_point_data(routing.nd,i,time;side)
        return expression.history+sum(v*routing.injections[d,k]
            for (d,k,v) in expression.terms;init=0.0)
    end
    r=c.system.rivers[i]
    grid=routing.rd.cohort_grid[i]
    B=point_coefficients(r,grid,time;side)
    history=cohort_point_arrival(r,r.history_grid,r.history_release,time;side)
    history+sum(river_transfer_expression(m,r,routing.rates[i][k],
        collect(view(B,k,1:river_curve_count(r)));coordinate=routing.coordinates[i][k])
        for k in eachindex(routing.rates[i]) if any(!iszero,view(B,k,:));init=0.0)
end

function routing_arrival_knots(c,routing,i)
    routing.exact ? sort!(unique!(vcat(deterministic_knots(routing.nd,i),c.grid))) :
        operational_arrival_knots(c,c.system.rivers[i],routing.rd.cohort_grid[i])
end

"""Bound actual stream rates, including confluence waves before cell averaging.
Each wave is affine between its cohort/response and operating knots. Both sides
of jumps are retained. Existing release coordinates serve these rows too.
"""
function constrain_river_stream_capacities!(m,c,routing)
    for (i,r) in enumerate(c.system.rivers)
        if routing.exact
            knots=sort!(unique!(vcat(deterministic_knots(routing.nd,i;kind=:release),c.grid)))
            for time in knots
                sides=time==first(c.grid) ? (:right,) :
                    time==last(c.grid) ? (:left,) : (:left,:right)
                for side in sides
                    expression=deterministic_point_data(routing.nd,i,time;side,kind=:release)
                    rate=expression.history+sum(v*routing.injections[d,k]
                        for (d,k,v) in expression.terms;init=0.0)
                    @constraint(m,rate<=opvalue(c,r.name,:capacity,time,r.capacity;side))
                end
            end
            continue
        end
        routing.rd.source_cohorts[i] && continue
        upstream=findall(e->e.target==r.source,c.system.rivers)
        for time in confluence_knots(c,routing.rd.cohort_grid,upstream)
            sides=time==first(c.grid) ? (:right,) :
                time==last(c.grid) ? (:left,) : (:left,:right)
            for side in sides
                t=river_control_interval(c,time,side)
                rate=routing.injections[i,t]+r.allocation*sum(
                    routing_point_expression(m,c,routing,d,time;side) for d in upstream;init=0.0)
                @constraint(m,rate<=opvalue(c,r.name,:capacity,time,r.capacity;side))
            end
        end
    end
    nothing
end

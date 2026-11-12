"""River equations and their transport discretization.

`arrivals`, `terminal` and `initial_transit` are volumes in Mm³. Window and
point observations, auxiliary starts and global domains use the same transport
coefficients for deterministic and distributed reaches.
`compiled` is the reusable deterministic transport cache, or `nothing`.
"""
struct RiverTransport{D,I}
    _data::D
    arrivals::Matrix{Union{AffExpr,QuadExpr,NonlinearExpr}}
    terminal::Vector{Union{AffExpr,QuadExpr,NonlinearExpr}}
    initial_transit::Vector{Float64}
    compiled::Union{Nothing,CompiledTransport}
    _rates::Vector{Vector{Any}}
    _coordinates::Vector{Vector{Any}}
    _primitives::Vector{Union{Nothing,ReferenceIntegrals}}
    _injections::I
    _starts::Vector{Pair{VariableRef,Union{AffExpr,QuadExpr,NonlinearExpr}}}
    _volume_rows::Dict{Tuple{Int,Float64,Float64},RiverTransferRow}
    _point_rows::Dict{Tuple{Int,Float64,Symbol},RiverTransferRow}
end

function transport_discretization(c,transport=nothing)
    if all(r->r.deterministic_delay!==nothing,c.system.rivers)
        return _transport_data(c,transport)
    end
    transport===nothing || throw(ArgumentError(
        "compiled deterministic transport cannot serve distributed dispatch"))
    DistributedRouting(c)
end

function RiverTransport(m,c,release;transport=nothing,
    injections=river_injections(c,release),coordinate=(r,q,i,k,window)->nothing,
    discretization=transport_discretization(c,transport),integral_grid=c.grid)
    exact=discretization isa CompiledTransport
    nd=exact ? discretization : nothing
    rd=exact ? nothing : discretization
    D,T=size(release)
    arrivals=Matrix{Union{AffExpr,QuadExpr,NonlinearExpr}}(undef,D,T)
    terminal=Vector{Union{AffExpr,QuadExpr,NonlinearExpr}}(undef,D)
    rates=Vector{Vector{Any}}(undef,D)
    coordinates=Vector{Vector{Any}}(undef,D)
    primitives=Vector{Union{Nothing,ReferenceIntegrals}}(nothing,D)
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
            coordinates[i]=Any[coordinate(r,q,i,k,(rd.cohort_grid[i][k],rd.cohort_grid[i][k+1]))
                for (k,q) in enumerate(rates[i])]
            # Fixed mixtures are linear coefficient maps. Variable mixtures
            # use anchored integrals to keep fine observation rows short.
            if !rd.source_cohorts[i] && r.deterministic_delay===nothing
                varying=any(curve->curve.edges!=first(r.curves).edges ||
                    curve.weights!=first(r.curves).weights,r.curves)
                fixed=all(z->z!==nothing && all(w->w isa Real,z.weights),coordinates[i])
                if varying && !fixed
                    primitives[i]=reference_integrals!(m,integral_grid,r,rates[i],coordinates[i],rd.cohort_grid[i],starts)
                end
            end
            primitive=primitives[i]
            cohort_dt=diff(rd.cohort_grid[i])
            B=rd.coefficients[i]
            for k in 1:K
                fine_arrival[i,k]=AffExpr(rd.history_arrival[i,k])+(primitive!==nothing ?
                    reference_arrival_volume(primitive,r,rd.grid[k],rd.grid[k+1]) : sum(
                    0.0036cohort_dt[j]*river_transfer_expression(m,r,rates[i][j],
                        fractions;coordinate=coordinates[i][j])
                    for (j,fractions) in arrival_terms(B,k);init=0.0))
            end
            intervals=rd.interval_coefficients[i]
            for t in 1:T
                arrivals[i,t]=AffExpr(sum(view(rd.history_arrival,i,
                    rd.vertices[t]:(rd.vertices[t+1]-1))))+(primitive!==nothing ?
                    reference_arrival_volume(primitive,r,c.grid[t],c.grid[t+1]) : sum(
                    0.0036cohort_dt[j]*river_transfer_expression(m,r,rates[i][j],
                        fractions;coordinate=coordinates[i][j])
                    for (j,fractions) in arrival_terms(intervals,t);init=0.0))
            end
            terminal[i]=AffExpr(rd.history_terminal[i])+(primitive!==nothing ?
                reference_transit(primitive,r,last(c.grid)) : sum(
                0.0036cohort_dt[j]*river_transfer_expression(m,r,rates[i][j],
                    view(B.terminal,:,j);coordinate=coordinates[i][j])
                for j in eachindex(rates[i]);init=0.0))
        end
    end
    volume_rows=Dict{Tuple{Int,Float64,Float64},RiverTransferRow}()
    point_rows=Dict{Tuple{Int,Float64,Symbol},RiverTransferRow}()
    routing=RiverTransport(discretization,arrivals,terminal,
        exact ? nd.initial_transit : rd.history_initial,nd,
        rates,coordinates,primitives,injections,starts,volume_rows,point_rows)
    constrain_river_stream_capacities!(m,c,routing)
    routing
end

function routing_volume_expression(m,c,routing,i,a,b)
    if routing._data isa CompiledTransport
        expression=_pulse_expression(routing._data.arrival_pulses[i],a,b)
        return expression.history+sum(v*routing._injections[d,k]
            for (d,k,v) in expression.terms;init=0.0)
    end
    r=c.system.rivers[i]
    primitive=routing._primitives[i]
    if primitive!==nothing
        history=only(route_volumes(r,[a,b],r.history_grid,r.history_release))
        return history+reference_arrival_volume(primitive,r,a,b)
    end
    grid=routing._data.cohort_grid[i]
    row=get!(routing._volume_rows,(i,Float64(a),Float64(b))) do
        volume_transfer_row(r,grid,a,b)
    end
    history=only(route_volumes(r,[a,b],r.history_grid,r.history_release))
    history+sum(0.0036*(grid[k+1]-grid[k])*river_transfer_expression(m,r,routing._rates[i][k],
        view(row.fractions,:,j);coordinate=routing._coordinates[i][k])
        for (j,k) in enumerate(row.cohorts);init=0.0)
end

function routing_point_expression(m,c,routing,i,time;side)
    if routing._data isa CompiledTransport
        expression=deterministic_point_data(routing._data,i,time;side)
        return expression.history+sum(v*routing._injections[d,k]
            for (d,k,v) in expression.terms;init=0.0)
    end
    r=c.system.rivers[i]
    primitive=routing._primitives[i]
    if primitive!==nothing
        history=cohort_point_arrival(r,r.history_grid,r.history_release,time;side)
        return history+reference_arrival_rate(primitive,r,time)
    end
    grid=routing._data.cohort_grid[i]
    row=get!(routing._point_rows,(i,Float64(time),side)) do
        point_transfer_row(r,grid,time;side)
    end
    history=cohort_point_arrival(r,r.history_grid,r.history_release,time;side)
    history+sum(river_transfer_expression(m,r,routing._rates[i][k],
        view(row.fractions,:,j);coordinate=routing._coordinates[i][k])
        for (j,k) in enumerate(row.cohorts);init=0.0)
end

function routing_arrival_knots(c,routing,i)
    (routing._data isa CompiledTransport) ? sort!(unique!(vcat(deterministic_knots(routing._data,i),c.grid))) :
        operational_arrival_knots(c,c.system.rivers[i],routing._data.cohort_grid[i])
end

"""Bound actual stream rates, including confluence waves before cell averaging.
Each wave is affine between its cohort/response and operating knots. Both sides
of jumps are retained. Existing release coordinates serve these rows too.
"""
function constrain_river_stream_capacities!(m,c,routing)
    for (i,r) in enumerate(c.system.rivers)
        if routing._data isa CompiledTransport
            knots=sort!(unique!(vcat(deterministic_knots(routing._data,i;kind=:release),c.grid)))
            for time in knots
                sides=time==first(c.grid) ? (:right,) :
                    time==last(c.grid) ? (:left,) : (:left,:right)
                for side in sides
                    expression=deterministic_point_data(routing._data,i,time;side,kind=:release)
                    rate=expression.history+sum(v*routing._injections[d,k]
                        for (d,k,v) in expression.terms;init=0.0)
                    @constraint(m,rate<=opvalue(c,r.name,:capacity,time,r.capacity;side))
                end
            end
            continue
        end
        routing._data.source_cohorts[i] && continue
        upstream=findall(e->e.target==r.source,c.system.rivers)
        for time in confluence_knots(c,routing._data.cohort_grid,upstream)
            sides=time==first(c.grid) ? (:right,) :
                time==last(c.grid) ? (:left,) : (:left,:right)
            for side in sides
                t=river_control_interval(c,time,side)
                rate=routing._injections[i,t]+r.allocation*sum(
                    routing_point_expression(m,c,routing,d,time;side) for d in upstream;init=0.0)
                @constraint(m,rate<=opvalue(c,r.name,:capacity,time,r.capacity;side))
            end
        end
    end
    nothing
end

function storage_observation_times(c,routing::RiverTransport,incoming)
    data=routing._data
    data isa CompiledTransport ? sort!(unique!(vcat(c.grid,
        [routing_arrival_knots(c,routing,j) for j in incoming]...))) : data.grid
end

"""Initialize transport auxiliaries in upstream-to-downstream equation order."""
function start_transport!(routing::RiverTransport)
    for (variable,expression) in routing._starts
        set_start_value(variable,JuMP.value(v->something(start_value(v),0.0),expression))
    end
    nothing
end

"""Lift prescribed physical cohorts into this model's auxiliary variables."""
function start_transport!(routing::RiverTransport,c,physical,put)
    routing._data isa CompiledTransport && return
    routed=routed_schedule(c,physical.river_release;
        generator_q=physical.generator_q,tunnel_q=physical.tunnel_q)
    for i in eachindex(routing._rates),k in eachindex(routing._rates[i])
        variable=routing._rates[i][k]
        variable isa VariableRef && put(variable,routed.cohorts[i].release[k])
    end
    for (i,data) in enumerate(routing._primitives)
        data===nothing && continue
        start_reference_integrals!(data,c.system.rivers[i],routed.cohorts[i].release,put)
    end
    nothing
end

# Global search needs finite domains. These are the integrals of the maximum
# admissible reference flow over each control window, including its time moment.
function bound_transport!(c,routing::RiverTransport)
    routing._data isa CompiledTransport && return
    for (i,data) in enumerate(routing._primitives)
        data===nothing && continue
        r=c.system.rivers[i]
        J=size(data.variables,1)
        for t in 1:length(data.control)-1
            left,right=data.control[t],data.control[t+1]
            cap=opvalue(c,r.name,:capacity,left,r.capacity)
            for series in c.operations
                series.object==r.name && series.attribute==:capacity || continue
                for (time,capacity) in zip(series.times,series.values)
                    left<time<right && (cap=max(cap,capacity))
                end
            end
            H=right-left
            for l in 1:J
                unit=zeros(J);unit[l]=1.0
                _,maximum_flux=river_transfer_range(r,unit,0.0,cap)
                for k in data.edges[t]:data.edges[t+1]-1
                    elapsed=data.grid[k+1]-data.control[t]
                    set_lower_bound(data.variables[l,1,k],0.0)
                    set_lower_bound(data.variables[l,2,k],0.0)
                    set_upper_bound(data.variables[l,1,k],0.0036elapsed*maximum_flux/0.3)
                    set_upper_bound(data.variables[l,2,k],0.0036elapsed^2/(2H)*maximum_flux/0.3)
                end
            end
        end
    end
end

"""Arrival-volume enclosures over all admissible releases, in Mm³."""
function routing_arrival_bounds(c,data,release_cap)
    s=c.system
    T=length(c.prices)
    al=zeros(length(s.rivers), T)
    au=similar(al)
    for (i, r) in enumerate(s.rivers), t in 1:T
        if data isa CompiledTransport
            lo=hi=data.arrival_history[i, t]
            for (d, k, v) in data.arrival_terms[i, t]
                z=v*release_cap[d, k]
                lo+=is_channel(s.rivers[d]) ? -abs(z) : min(0.0,z)
                hi+=max(0.0, z)
            end
        else
            rows=data.vertices[t]:(data.vertices[t+1]-1)
            lo=hi=sum(view(data.history_arrival,i,rows))
            B=data.interval_coefficients[i]
            cohort=data.cohort_grid[i]
            source=data.source_cohorts[i]
            for (k,coefficients) in arrival_terms(B,t)
                period=source ? k : data.control[k]
                cap=release_cap[i,period]
                l,h=river_transfer_range(r,coefficients,river_flow_lower(r,cap),cap)
                l*=0.0036*(cohort[k+1]-cohort[k])
                h*=0.0036*(cohort[k+1]-cohort[k])
                lo+=l
                hi+=h
            end
        end
        al[i, t]=lo
        au[i, t]=hi
    end
    al,au
end

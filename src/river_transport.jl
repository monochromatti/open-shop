"""River transport mesh, independent of control and hydraulic decision times.

A supplied mesh is retained. Otherwise distributed confluences use four transport
cells per control interval; source cohorts retain their control grid.
Exact deterministic networks do not need a mixing mesh. Finer independent routing
and hydraulic replay still determine whether this resolution is sufficient.
"""
function river_routing_grid(c;physical_grid=c.grid)
    grid=if !isempty(c.routing_grid)
        c.routing_grid
    elseif any(!,river_cohort_sources(c.system)) && any(r->r.deterministic_delay===nothing,c.system.rivers)
        refined_grid(c.grid;factor=4,existing=physical_grid)
    else
        c.grid
    end
    sort!(unique!(vcat(grid,physical_grid)))
end

"""Finite-volume river transport data, without control or state decisions.

Source reaches retain the hydraulic/control cohort grid. Zero-storage relay
reaches use `grid`; their rates follow subinterval conservation. Coefficients
map each reach's own cohort grid into transport-cell arrival volumes.
"""
struct DistributedRouting
    grid::Vector{Float64}
    control::Vector{Int}
    vertices::Vector{Int}
    source_cohorts::BitVector
    cohort_grid::Vector{Vector{Float64}}
    coefficients::Vector{RiverTransfer}
    interval_coefficients::Vector{RiverTransfer}
    history_arrival::Matrix{Float64}
    history_initial::Vector{Float64}
    history_terminal::Vector{Float64}
end

function river_cohort_sources(s)
    reservoirs=Set(r.name for r in s.reservoirs)
    BitVector(r.source in reservoirs || !any(e->e.target==r.source,s.rivers)
        for r in s.rivers)
end

function DistributedRouting(c;grid=river_routing_grid(c))
    RiverRouting.check_grid(grid)
    first(grid)==first(c.grid) && last(grid)==last(c.grid) && all(t in grid for t in c.grid) ||
        throw(ArgumentError("transport mesh must retain the hydraulic grid"))
    source=river_cohort_sources(c.system)
    control=[searchsortedlast(c.grid,(grid[k]+grid[k+1])/2) for k in 1:(length(grid)-1)]
    vertices=[searchsortedfirst(grid,t) for t in c.grid]
    cohorts=[source[d] ? copy(c.grid) : copy(grid) for d in eachindex(source)]
    B=[RiverTransfer(r,grid,cohorts[d]) for (d,r) in enumerate(c.system.rivers)]
    intervals=grid==c.grid ? B :
        [RiverTransfer(r,c.grid,cohorts[d]) for (d,r) in enumerate(c.system.rivers)]
    history=zeros(length(source),length(grid)-1)
    initial=zeros(length(source));terminal=similar(initial)
    for (d,r) in enumerate(c.system.rivers)
        history[d,:].=route_volumes(r,grid,r.history_grid,r.history_release)
        initial[d]=remaining_volume(r,first(grid),r.history_grid,r.history_release)
        terminal[d]=remaining_volume(r,last(grid),r.history_grid,r.history_release)
    end
    DistributedRouting(copy(grid),control,vertices,source,cohorts,B,intervals,history,initial,terminal)
end

"""Known arrivals and current-cohort transfer rows for one hydraulic step."""
struct DistributedRoutingStep
    known_arrival::Matrix{Float64}
    transfer::Vector{RiverTransfer}
    block::UnitRange{Int}
    interval::Int
    duration::Vector{Float64}
    source_cohorts::BitVector
end

function DistributedRoutingStep(mesh::DistributedRouting,s,t,known)
    block=mesh.vertices[t]:(mesh.vertices[t+1]-1)
    source=mesh.source_cohorts
    DistributedRoutingStep(copy(known[:,block]),mesh.coefficients,block,t,
        diff(view(mesh.grid,first(block):(last(block)+1))),source)
end

# Trial-point continuation agrees with the existing release-time reference law.
# Only the blend coordinate is held at an endpoint outside its physical domain.
function released_fraction(r,q,coefficients)
    (r.deterministic_delay!==nothing || length(r.curves)==1) && return q*first(coefficients)
    curves=r.curves
    a=clamp(q,first(curves).reference_flow,last(curves).reference_flow)
    j=RiverRouting.reference_segment(curves,a;reference=c->c.reference_flow)
    weight=(a-curves[j].reference_flow)/(curves[j+1].reference_flow-curves[j].reference_flow)
    q*((1-weight)*coefficients[j]+weight*coefficients[j+1])
end

"""Independent routing, observed on the original control intervals."""
function routed_schedule(c,release;transport=nothing,injections=nothing,generator_q=nothing,tunnel_q=nothing)
    injections=injections===nothing ? river_injections(c,release,generator_q,tunnel_q) : injections
    exact=all(r->r.deterministic_delay!==nothing,c.system.rivers)
    if exact
        routed=route_network_exact(c,release;transport,injections,generator_q,tunnel_q)
        return (;release=routed["release"],arrival_volume=routed["arrival_volume"],
            transit=routed["transit"],cohorts=routed["release_cohorts"],injections)
    end
    grid=river_routing_grid(c)
    routed=_route_network(c,release;grid,injections,generator_q,tunnel_q)
    edges=[searchsortedfirst(grid,t) for t in c.grid]
    R,A=routed.release,routed.arrival_volume
    mean_release=[sum(R[d,k]*(grid[k+1]-grid[k]) for k in edges[t]:(edges[t+1]-1))/
        (c.grid[t+1]-c.grid[t]) for d in axes(R,1),t in eachindex(c.prices)]
    arrivals=[sum(view(A,d,edges[t]:(edges[t+1]-1))) for d in axes(A,1),t in eachindex(c.prices)]
    (;release=mean_release,arrival_volume=arrivals,transit=routed.transit[:,edges],
        routed.cohorts,injections)
end

function reference_flux!(values,r,q)
    fill!(values,0.0)
    if r.deterministic_delay!==nothing
        values[1]=q
    else
        for (l,w) in RiverRouting.blend(r.curves,q)
            values[l]=q*w
        end
    end
    values
end

"""Refine transport cohorts for prescribed controls before the next NLP.

Compare cumulative arrivals on nested meshes in receiving-reservoir volume
units. The returned mesh is the coarse one, so use the second-order 4/3
correction for its estimated error, as in coupled hydraulic replay.
Mark entry cells by variation of q*theta(q), the reference flux that is
projected onto constant cohorts. Bulk marking splits cells covering half of
that variation. Control times and all original transport vertices are retained.
This is a numerical error estimate; independent coupled replay remains required.
"""

function refine_transport_grid(c,x,storage_limits,absolute;deadline=Inf,max_refinements=16)
    grid=river_routing_grid(c)
    all(r->r.deterministic_delay!==nothing,c.system.rivers) && return grid
    error_factor=4/3
    s=c.system
    sources=river_cohort_sources(s)
    injections=river_injections(c,x.river_release,x.generator_q,x.tunnel_q)
    incidence=[Float64(reach.target==r.name) for r in s.reservoirs,reach in s.rivers]
    current=_route_network(c,injections;grid,injections)
    for iteration in 1:max_refinements
        time()<deadline || break
        probe=sort!(unique!(refined_grid(grid)))
        fine=_route_network(c,injections;grid=probe,injections)
        coarse=zeros(size(fine.arrival_volume))
        for (d,r) in enumerate(s.rivers)
            cohort=current.cohorts[d]
            coarse[d,:].=route_volumes(r,probe,cohort.grid,cohort.release)+
                route_volumes(r,probe,r.history_grid,r.history_release)
        end
        cumulative=cumsum(fine.arrival_volume-coarse;dims=2)
        errors=vec(maximum(abs.(incidence*cumulative);dims=2))
        all(error_factor .*errors .<=storage_limits) &&
            error_factor*maximum(abs,cumulative;init=0.0)<=absolute/10 &&
            error_factor*maximum(abs,current.transit[:,end]-fine.transit[:,end];init=0.0)<=absolute && break
        edges=[searchsortedfirst(probe,t) for t in grid]
        duration=diff(probe)
        indicator=zeros(length(grid)-1)
        for (d,r) in enumerate(s.rivers)
            sources[d] && continue
            J=river_curve_count(r)
            mean_flux=zeros(J);flux=similar(mean_flux)
            for k in eachindex(indicator)
                block=edges[k]:(edges[k+1]-1)
                rates=view(fine.release,d,block)
                mean=sum(fine.release[d,j]*duration[j] for j in block)/(grid[k+1]-grid[k])
                mean=clamp(mean,minimum(rates),maximum(rates))
                reference_flux!(mean_flux,r,mean)
                for j in block
                    reference_flux!(flux,r,fine.release[d,j])
                    indicator[k]+=0.0036duration[j]*sum(abs.(flux-mean_flux))
                end
            end
        end
        total=sum(indicator)
        total>0 || break
        marked=Int[];covered=0.0
        for k in sortperm(indicator;rev=true)
            push!(marked,k);covered+=indicator[k]
            covered>=total/2 && break
        end
        vertices=[refinement_vertex(grid[k],grid[k+1],0.5,probe) for k in marked]
        grid=sort!(unique!(vcat(grid,vertices)))
        current=_route_network(c,injections;grid,injections)
    end
    grid
end

"""Arrival and capacity shortfalls in m³/s, indexed by reach.

`mean_arrival` uses the declared observation windows; `pointwise_arrival` uses
the actual wave. `capacity` also checks confluence peaks before cell averaging.
Channels have zero entries: their signed hydraulic laws are checked separately.
"""
function river_arrival_residuals(c,routed)
    D=length(c.system.rivers)
    mean_arrival=zeros(D)
    capacity=zeros(D)
    pointwise_arrival=zeros(D)
    exact=all(r->r.deterministic_delay!==nothing,c.system.rivers)
    for (d,r) in enumerate(c.system.rivers)
        is_channel(r) && continue
        cohort=routed.cohorts[d]
        windows=isempty(r.arrival_window_grid) ? c.grid : r.arrival_window_grid
        volume=route_volumes(r,windows,cohort.grid,cohort.release)+
            route_volumes(r,windows,r.history_grid,r.history_release)
        mean_arrival[d]=maximum(
            opaverage(c,r.name,:min_arrival,windows[k],windows[k+1],r.min_arrival)-
            volume[k]/(0.0036*(windows[k+1]-windows[k]))
            for k in eachindex(volume);init=0.0)
        capacity_error=maximum(cohort.release[k]-opinterval(c,r.name,:capacity,t,r.capacity)
            for k in eachindex(cohort.release),t in eachindex(c.prices)
            if min(cohort.grid[k+1],c.grid[t+1])>max(cohort.grid[k],c.grid[t]);init=0.0)
        upstream=findall(e->e.target==r.source,c.system.rivers)
        if !exact && r.law==:junction && !isempty(upstream)
            grids=[cohort.grid for cohort in routed.cohorts]
            for time in confluence_knots(c,grids,upstream)
                sides=time==first(c.grid) ? (:right,) :
                    time==last(c.grid) ? (:left,) : (:left,:right)
                for side in sides
                    t=river_control_interval(c,time,side)
                    rate=routed.injections[d,t]+r.allocation*sum(
                        point_arrival(c.system.rivers[i],routed.cohorts[i].grid,
                            routed.cohorts[i].release,time;side) for i in upstream;init=0.0)
                    capacity_error=max(capacity_error,rate-opvalue(c,r.name,:capacity,time,r.capacity;side))
                end
            end
        end
        capacity[d]=capacity_error
        pointwise_arrival[d]=r.arrival_policy==:pointwise ?
            arrival_requirement_violation(c,r,cohort.grid,cohort.release) : 0.0
    end
    (;mean_arrival,capacity,pointwise_arrival)
end

river_control_interval(c,time,side)=clamp(side==:left ? searchsortedfirst(c.grid,time)-1 :
    searchsortedlast(c.grid,time),1,length(c.prices))

function confluence_knots(c,grids,upstream)
    sort!(unique!(vcat(c.grid,[arrival_knots(c.system.rivers[d],grids[d]) for d in upstream]...)))
end

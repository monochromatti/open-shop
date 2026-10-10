# Local injections and hydraulic incidence are compiled from explicit destinations.
# Directed outfalls transport water; their hydraulic targets supply head only.
water_incidence(x, node) = Int(x.discharge_river===nothing && x.target==node)-Int(x.source==node)

function river_source_indices(s)
    plants=Dict(p.name=>p.discharge_river for p in s.plants)
    (
        generators=[findall(g->plants[g.plant]==r.name,s.generators) for r in s.rivers],
        tunnels=[findall(e->e.discharge_river==r.name,s.tunnels) for r in s.rivers],
        reservoirs=BitVector(any(x.name==r.source for x in s.reservoirs) for r in s.rivers),
    )
end

function river_local_sources(c)
    idx=river_source_indices(c.system)
    BitVector(idx.reservoirs[d] || !isempty(idx.generators[d]) || !isempty(idx.tunnels[d]) ||
        r.inflow>0 || any(z.object==r.name && z.attribute==:inflow && any(>(0),z.values) for z in c.operations)
        for (d,r) in enumerate(c.system.rivers))
end

"""Local source rates for compiled transport, excluding incoming river arrivals.

Reservoir outlet rows contain total reach input (including natural inflow).
Other rows contain plant/tunnel discharges and natural inflow; confluence
arrivals are added by the transport operator. Each stream enters exactly once.
"""
function river_injections(c, release, generator_q=nothing, tunnel_q=nothing)
    s=c.system
    indices=river_source_indices(s)
    T=length(c.prices)
    size(release)==(length(s.rivers),T) || throw(DimensionMismatch("river releases"))
    any(!isempty,indices.generators) && generator_q===nothing &&
        throw(ArgumentError("routing plant outfalls requires generator_q"))
    any(!isempty,indices.tunnels) && tunnel_q===nothing &&
        throw(ArgumentError("routing tunnel outfalls requires tunnel_q"))
    generator_q===nothing || size(generator_q)==(length(s.generators),T) ||
        throw(DimensionMismatch("generator discharges"))
    tunnel_q===nothing || size(tunnel_q)==(length(s.tunnels),T) ||
        throw(DimensionMismatch("tunnel discharges"))
    [indices.reservoirs[d] ? release[d,t] :
        opinterval(c,r.name,:inflow,t,r.inflow)+
        sum(generator_q[i,t] for i in indices.generators[d];init=0.0)+
        sum(tunnel_q[i,t] for i in indices.tunnels[d];init=0.0)
        for (d,r) in enumerate(s.rivers),t in 1:T]
end
function _river_injections(s, release, generator_q, tunnel_q;indices=river_source_indices(s))
    [indices.reservoirs[d] ? release[d] :
        r.inflow+sum(generator_q[i] for i in indices.generators[d];init=0.0)+
        sum(tunnel_q[i] for i in indices.tunnels[d];init=0.0)
        for (d,r) in enumerate(s.rivers)]
end

# Ordinary conservation rows at zero-storage confluences. The same local source
# matrix also drives the exact routing operator and chronological reconstruction.
function constrain_river_sources!(m,c,release,arrivals,injections)
    for j in c.system.river_junctions
        outgoing=only(findall(r->r.source==j.name,c.system.rivers))
        incoming=findall(r->r.target==j.name,c.system.rivers)
        for t in eachindex(c.prices)
            @constraint(m,(0.0036*(c.grid[t+1]-c.grid[t])*(release[outgoing,t]-injections[outgoing,t])-
                sum(arrivals[d,t] for d in incoming;init=0.0))/0.3==0)
        end
    end
    nothing
end

# Fixed-delay reaches remain linear even inside a distributed network.
river_curve_count(r) = r.deterministic_delay===nothing ? length(r.curves) : 1
river_reference_coordinate!(m,r,q;kwargs...) = r.deterministic_delay===nothing ?
    distributed_reference_coordinate!(m,r.curves,q;kwargs...) : nothing
river_transfer_expression(m,r,q,b;kwargs...) = r.deterministic_delay===nothing ?
    distributed_transfer_expression(m,r.curves,q,b;kwargs...) : first(b)*q
river_transfer_range(r,b,lo,hi) = r.deterministic_delay===nothing ?
    RiverRouting.transfer_range(r.curves,b,lo,hi) : extrema((first(b)*lo,first(b)*hi))

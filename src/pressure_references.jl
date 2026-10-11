# Open passive tunnels determine head differences. A junction component with
# no reservoir/boundary has a constant-pressure nullspace. Active generation
# preserves that freedom only if both endpoints share the same component and
# there is no absolute outlet-head floor. No threshold substitutes for closure.
function _floating_pressure_components(s,opening,generator_q;ix=nodeindex(s))
    nr,nj=length(s.reservoirs),length(s.junctions)
    nj==0 && return Vector{Int}[]
    parent=collect(1:length(ix))
    function root(i)
        while parent[i]!=i
            parent[i]=parent[parent[i]]
            i=parent[i]
        end
        i
    end
    for (i,e) in enumerate(s.tunnels)
        opening[i]>0 || continue
        a,b=root(ix[e.source]),root(ix[e.target])
        parent[b]=a
    end
    component=root.(eachindex(parent))
    anchored=falses(length(parent))
    for i in eachindex(parent)
        nr<i<=nr+nj || (anchored[component[i]]=true)
    end
    for (i,g) in enumerate(s.generators)
        generator_q[i]==0 && continue
        p=plantof(s,g)
        a,b=component[ix[p.source]],component[ix[p.target]]
        if a!=b || p.outlet_head_floor!==nothing
            anchored[a]=anchored[b]=true
        end
    end
    groups=Vector{Int}[]
    group_index=zeros(Int,length(parent))
    for i in (nr+1):(nr+nj)
        r=component[i]
        anchored[r] && continue
        if group_index[r]==0
            push!(groups,Int[])
            group_index[r]=length(groups)
        end
        push!(groups[group_index[r]],i)
    end
    groups
end

# Choose the midpoint of the compatible reference-shift interval. This moves
# an entire nullspace component by one constant, never individual heads or
# physical flows. Incompatible bounds remain violated and fail the usual audit.
function _normalize_pressure_references!(heads,s,groups)
    nr=length(s.reservoirs)
    for group in groups
        lo=maximum(s.junctions[i-nr].hmin-heads[i] for i in group)
        hi=minimum(s.junctions[i-nr].hmax-heads[i] for i in group)
        lo<=hi || continue
        shift=isfinite(lo) && isfinite(hi) ? lo/2+hi/2 : clamp(0.0,lo,hi)
        for i in group
            heads[i]+=shift
        end
    end
    heads
end

# Compare replay heads modulo the same physical nullspace, while the ordinary
# equation/bound audits still test every supplied head and every head difference.
function _pressure_replay_difference(s,opening,generator_q,heads,replayed)
    difference=heads-replayed
    ix=nodeindex(s)
    for t in axes(heads,2)
        groups=_floating_pressure_components(s,view(opening,:,t),view(generator_q,:,t);ix)
        for group in groups
            shift=sum(difference[i,t] for i in group)/length(group)
            for i in group
                difference[i,t]-=shift
            end
        end
    end
    difference
end

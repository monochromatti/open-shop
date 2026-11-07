scaled_defect(error,limit) = iszero(error) ? 0.0 : error/limit

"""Bisect cells accounting for half the squared local defect indicators."""
function bulk_refined_grid(grid, indicator)
    length(indicator)==length(grid)-1 || throw(DimensionMismatch("cell indicators"))
    total=sum(abs2,indicator)
    total>0 || return grid
    marked=Int[]
    covered=0.0
    for k in sortperm(indicator;rev=true)
        push!(marked,k)
        covered+=indicator[k]^2
        covered>=total/2 && break
    end
    sort!(unique!(vcat(grid,[(grid[k]+grid[k+1])/2 for k in marked])))
end

"""Compare two half steps with each represented cell, retaining executed cohorts.

These local defects locate integration error without propagating an earlier
state error through every later cell. They choose a mesh; coupled whole-horizon
replay remains the acceptance test. Existing transport vertices are retained.
"""
function hydraulic_defects(c,x;storage_limits,level_limit,transport_limit,deadline=Inf)
    physical=hydraulic_case(c,x.grid)
    s=physical.system
    routed=routed_schedule(physical,x.river_release;generator_q=x.generator_q,tunnel_q=x.tunnel_q)
    transport_grid=river_routing_grid(physical)
    structure=(;topology=HydraulicTopology(s),sources=river_source_indices(s),
        active_rivers=findall(r->r.law!=:junction,s.rivers),river_order=_forward_river_order(s))
    indicator=zeros(length(x.grid)-1)
    for k in eachindex(indicator)
        time()<deadline || return nothing
        a,b=x.grid[k:k+1]
        middle=(a+b)/2
        reservoirs=[_river_replace(r;v0=x.storage[i,k]) for (i,r) in enumerate(s.reservoirs)]
        rivers=[begin
            if is_channel(r)
                _river_replace(r;history_grid=[a-1,a],history_release=[0.0])
            else
                cohort=routed.cohorts[d]
                past=release_history(cohort.grid,cohort.release,a)
                history_grid=isempty(r.history_grid) ? past.grid : vcat(r.history_grid,past.grid[2:end])
                _river_replace(r;history_grid,history_release=vcat(r.history_release,past.release))
            end
        end for (d,r) in enumerate(s.rivers)]
        system=_river_replace(s;reservoirs,rivers)
        routing_grid=transport_grid[searchsortedfirst(transport_grid,a):searchsortedlast(transport_grid,b)]
        cell=ScheduleCase(;name="local hydraulic step",system,grid=[a,b],prices=[0.0],
            operations=c.operations,routing_grid)
        q=x.generator_q[:,k:k]
        gate=x.gate[:,k:k]
        u=x.commitment[:,k:k]
        opening=x.tunnel_opening[:,k:k]
        reference=x.head[:,k:k]
        fine=_simulate_hydraulics(cell,q,gate,u,opening,reference,[a,middle,b];structure...)
        if !fine.converged
            indicator[k]=Inf
            continue
        end
        for (i,r) in enumerate(s.reservoirs)
            indicator[k]=max(indicator[k],
                scaled_defect(abs(fine.storage[i,end]-x.storage[i,k+1]),storage_limits[i]),
                scaled_defect(abs(head(r,fine.storage[i,end])-head(r,x.storage[i,k+1])),level_limit))
        end
        for (d,r) in enumerate(s.rivers)
            is_channel(r) && continue
            J=river_curve_count(r)
            flux=zeros(J)
            left=similar(flux)
            right=similar(flux)
            reference_flux!(flux,r,x.river_release[d,k])
            reference_flux!(left,r,fine.river_release[d,1])
            reference_flux!(right,r,fine.river_release[d,2])
            variation=0.0036*(b-a)/2*(sum(abs,left-flux)+sum(abs,right-flux))
            indicator[k]=max(indicator[k],scaled_defect(variation,transport_limit))
        end
    end
    indicator
end

"""Reintegrate held controls, with a new commitment or a new hydraulic mesh."""
function controlled_trajectory(c,x;u=nothing,opening=nothing,grid=x.grid)
    mesh=TimeDiscretization(c.grid,x.grid)
    cells=first.(mesh.cells)
    u=u===nothing ? x.commitment[:,cells] : u
    q=copy(x.generator_q[:,cells])
    q[u.==0].=0.0
    opening=opening===nothing ? x.tunnel_opening[:,cells] : opening
    dispatch_from_controls(c,u,q,x.gate[:,cells];grid,tunnel_opening=opening,
        reserves=_reserve_slice(x.reserves,cells),head_reference=interval_means(mesh,x.head))
end

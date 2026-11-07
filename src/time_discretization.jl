"""A hydraulic partition with unchanged operating intervals.

`interval[k]` locates hydraulic cell k on the operating grid; `cells[j]`
and `vertices[j]` locate an operating interval and its starting state.
"""
struct TimeDiscretization
    grid::Vector{Float64}
    interval::Vector{Int}
    vertices::Vector{Int}
    cells::Vector{UnitRange{Int}}
    dt::Vector{Float64}
    duration::Vector{Float64}
end

function TimeDiscretization(operating, hydraulic=operating)
    RiverRouting.check_grid(operating)
    grid=Float64.(hydraulic)
    RiverRouting.check_grid(grid)
    first(grid)==first(operating) && last(grid)==last(operating) ||
        throw(ArgumentError("hydraulic mesh must span the operating horizon"))
    vertices=searchsortedfirst.(Ref(grid),operating)
    all(grid[vertices].==operating) ||
        throw(ArgumentError("hydraulic mesh must retain every operating boundary"))
    cells=[vertices[j]:(vertices[j+1]-1) for j in 1:length(operating)-1]
    interval=Vector{Int}(undef,length(grid)-1)
    for j in eachindex(cells)
        interval[cells[j]].=j
    end
    TimeDiscretization(grid,interval,vertices,cells,diff(grid),Float64.(diff(operating)))
end

"""Physical interval data with the original river-transport resolution retained."""
function hydraulic_case(c,grid)
    physical=with_grid(c,grid)
    _river_replace(physical;routing_grid=river_routing_grid(c;physical_grid=grid))
end

"""Duration-weighted interval means, for numbers or optimization expressions."""
function interval_means(mesh::TimeDiscretization,values::AbstractMatrix)
    size(values,2)==length(mesh.dt) || throw(DimensionMismatch("hydraulic rate columns"))
    length(mesh.cells)==length(mesh.dt) && return values
    [sum(values[i,k]*mesh.dt[k] for k in mesh.cells[j])/mesh.duration[j]
        for i in axes(values,1), j in eachindex(mesh.cells)]
end

"""Add extensive cell quantities, such as transported water volumes."""
function interval_sums(mesh::TimeDiscretization,values::AbstractMatrix)
    size(values,2)==length(mesh.dt) || throw(DimensionMismatch("hydraulic volume columns"))
    length(mesh.cells)==length(mesh.dt) && return values
    [sum(values[i,k] for k in mesh.cells[j])
        for i in axes(values,1), j in eachindex(mesh.cells)]
end

function operating_window_values(mesh::TimeDiscretization,z,u)
    (;commitment=u,power=interval_means(mesh,z.power),
        generator_q=interval_means(mesh,z.generator_q),
        river_release=interval_means(mesh,z.river_release),
        storage=z.storage[:,mesh.vertices],tunnel_q=interval_means(mesh,z.tunnel_q),
        arrival_volume=interval_sums(mesh,z.arrival_volume))
end

function operating_window_values(c,z,grid,indices,u)
    mesh=TimeDiscretization(c.grid,grid)
    indices==mesh.interval || throw(ArgumentError("inconsistent operating interval indices"))
    operating_window_values(mesh,z,u)
end

function operating_window_values(c,x::Schedule)
    mesh=TimeDiscretization(c.grid,x.grid)
    firstcell=first.(mesh.cells)
    operating_window_values(mesh,x,x.commitment[:,firstcell])
end

"""Maximum departure from each operating interval's first control value."""
function held_control_difference(mesh::TimeDiscretization,values::AbstractMatrix)
    difference=0.0
    for cells in mesh.cells,k in cells,i in axes(values,1)
        difference=max(difference,abs(values[i,k]-values[i,first(cells)]))
    end
    difference
end

"""Deviation from controls held over the original operating intervals."""
function held_control_residuals(c,x::Schedule,mesh::TimeDiscretization)
    deviation(values)=held_control_difference(mesh,values)
    adjustable=findall(g->g.characteristic===nothing,c.system.generators)
    (;commitment=deviation(x.commitment),discharge=deviation(view(x.generator_q,adjustable,:)),
        gate=deviation(x.gate),tunnel_opening=deviation(x.tunnel_opening))
end

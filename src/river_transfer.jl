"""One observation row: release indices and reference × release fractions."""
struct RiverTransferRow
    release_indices::Vector{Int}
    fractions::Matrix{Float64}
end

"""Constant-release transfer fractions in compressed column and row order.

`fractions[:,p]` contains all reference curves for one nonzero arrival/release
pair. Columns are interval releases; the row index borrows the same coefficient
buffer. Thus dispatch and chronological replay traverse only overlapping pairs
without duplicating coefficients or introducing reference-flow variables.
`terminal[:,j]` is the already-released fraction still in transit at the horizon.
"""
struct RiverTransfer
    fractions::Matrix{Float64}
    release_offsets::Vector{Int}
    arrival_indices::Vector{Int}
    arrival_offsets::Vector{Int}
    release_indices::Vector{Int}
    fraction_indices::Vector{Int}
    terminal::Matrix{Float64}
end

arrival_terms(B::RiverTransfer,t) =
    ((B.release_indices[p],view(B.fractions,:,B.fraction_indices[p]))
        for p in B.arrival_offsets[t]:(B.arrival_offsets[t+1]-1))
release_terms(B::RiverTransfer,j) =
    ((B.arrival_indices[p],view(B.fractions,:,p))
        for p in B.release_offsets[j]:(B.release_offsets[j+1]-1))

function release_terms(B::RiverTransfer,j,rows::UnitRange{Int})
    start,stop=B.release_offsets[j],B.release_offsets[j+1]-1
    column=view(B.arrival_indices,start:stop)
    lo=start-1+searchsortedfirst(column,first(rows))
    hi=start-1+searchsortedlast(column,last(rows))
    ((B.arrival_indices[p],view(B.fractions,:,p)) for p in lo:hi)
end

function river_delay_support(r)
    if r.deterministic_delay!==nothing
        d=r.deterministic_delay
        isfinite(d) && d>=0 || throw(ArgumentError("deterministic delay must be finite and nonnegative"))
        return d,d
    end
    RiverRouting.check_curves(r.curves)
    minimum(first(c.edges) for c in r.curves),maximum(last(c.edges) for c in r.curves)
end

function release_arrival_fraction(r,l,a,b,s,e)
    r.deterministic_delay===nothing &&
        return RiverRouting.release_transfer_fraction(r.curves[l],a,b,s,e)
    d=r.deterministic_delay
    max(0.0,min(b,e+d)-max(a,s+d))/(e-s)
end

function release_transit_fraction(r,l,time,s,e)
    released=min(time,e)
    released<=s && return 0.0
    r.deterministic_delay===nothing &&
        return RiverRouting.integral_survival(r.curves[l],time-released,time-s)/(e-s)
    max(0.0,released-max(s,time-r.deterministic_delay))/(e-s)
end

# Search shifted edges directly. Inverting t=s+delay by subtraction can select
# the wrong release at an atomic jump, e.g. 1.75+0.3 == 2.05 in Float64.
function shifted_period(grid,time,delay;side=:right)
    lo,hi=0,length(grid)
    while lo<hi
        mid=(lo+hi+1)>>>1
        edge=grid[mid]+delay
        before=side==:right ? edge<=time : edge<time
        if before
            lo=mid
        else
            hi=mid-1
        end
    end
    lo
end

function RiverTransfer(r,arrival_grid,release_grid)
    RiverRouting.check_grid(arrival_grid)
    RiverRouting.check_grid(release_grid)
    lagmin,lagmax=river_delay_support(r)
    A,C,J=length(arrival_grid)-1,length(release_grid)-1,river_curve_count(r)
    release_offsets=Vector{Int}(undef,C+1)
    arrival_indices=Int[]
    values=Float64[]
    counts=zeros(Int,A)
    terminal=Matrix{Float64}(undef,J,C)
    coefficients=zeros(J)
    for j in 1:C
        release_offsets[j]=length(arrival_indices)+1
        s,e=release_grid[j],release_grid[j+1]
        lo=max(1,searchsortedlast(arrival_grid,s+lagmin))
        hi=min(A,searchsortedfirst(arrival_grid,e+lagmax)-1)
        for t in lo:hi
            for l in 1:J
                coefficients[l]=release_arrival_fraction(r,l,arrival_grid[t],arrival_grid[t+1],s,e)
            end
            any(!iszero,coefficients) || continue
            push!(arrival_indices,t)
            append!(values,coefficients)
            counts[t]+=1
        end
        for l in 1:J
            terminal[l,j]=release_transit_fraction(r,l,last(arrival_grid),s,e)
        end
    end
    release_offsets[end]=length(arrival_indices)+1
    arrival_offsets=vcat(1,1 .+ cumsum(counts))
    positions=copy(arrival_offsets[1:end-1])
    release_indices=similar(arrival_indices)
    fraction_indices=similar(arrival_indices)
    for j in 1:C,p in release_offsets[j]:(release_offsets[j+1]-1)
        t=arrival_indices[p]
        at=positions[t]
        release_indices[at]=j
        fraction_indices[at]=p
        positions[t]+=1
    end
    RiverTransfer(reshape(values,J,:),release_offsets,arrival_indices,
        arrival_offsets,release_indices,fraction_indices,terminal)
end

"""Fractions arriving in `[a,b]`, restricted to its delay support."""
function volume_transfer_row(r,grid,a,b)
    RiverRouting.check_grid(grid)
    isfinite(a) && isfinite(b) && a<b || throw(ArgumentError("ordered finite arrival window required"))
    lagmin,lagmax=river_delay_support(r)
    lo=max(1,shifted_period(grid,a,lagmax))
    hi=min(length(grid)-1,shifted_period(grid,b,lagmin;side=:left))
    release_indices=Int[];values=Float64[];J=river_curve_count(r)
    coefficients=zeros(J)
    for j in lo:hi
        for l in 1:J
            coefficients[l]=release_arrival_fraction(r,l,a,b,grid[j],grid[j+1])
        end
        any(!iszero,coefficients) || continue
        push!(release_indices,j);append!(values,coefficients)
    end
    RiverTransferRow(release_indices,reshape(values,J,:))
end

"""Instantaneous arrival coefficients, retaining either side of atomic jumps."""
function point_transfer_row(r,grid,time;side=:right)
    RiverRouting.check_grid(grid)
    isfinite(time) || throw(ArgumentError("arrival query time must be finite"))
    side in (:left,:right) || throw(ArgumentError("side must be :left or :right"))
    lagmin,lagmax=river_delay_support(r)
    lo,hi=if r.deterministic_delay===nothing
        max(1,shifted_period(grid,time,lagmax)),
            min(length(grid)-1,shifted_period(grid,time,lagmin;side=:left))
    else
        release_interval=shifted_period(grid,time,lagmin;side)
        max(1,release_interval),min(length(grid)-1,release_interval)
    end
    release_indices=Int[];values=Float64[];J=river_curve_count(r)
    coefficients=zeros(J)
    for j in lo:hi
        s,e=grid[j],grid[j+1]
        for l in 1:J
            coefficients[l]=if r.deterministic_delay===nothing
                clamp(delay_cdf(r.curves[l],time-s)-delay_cdf(r.curves[l],time-e),0.0,1.0)
            else
                d=r.deterministic_delay
                Float64(side==:right ? s+d<=time<e+d : s+d<time<=e+d)
            end
        end
        any(!iszero,coefficients) || continue
        push!(release_indices,j);append!(values,coefficients)
    end
    RiverTransferRow(release_indices,reshape(values,J,:))
end

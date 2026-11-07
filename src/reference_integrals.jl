"""Two reference-flow integrals, anchored at each supplied window start.

`C[l,1,k]` is κ∫ψ_l(q) ds through cohort k, and `C[l,2,k]` is
κ∫((s-start)/duration)ψ_l(q) ds, both in Mm³. Constant cohorts give exact
partial-cell integrals. Arrival windows then use endpoint differences instead
of sums over all overlapping cohorts. `variables` are scaled by 0.3 Mm³.
"""
struct ReferenceIntegrals
    C::Array{AffExpr,3}
    variables::Array{VariableRef,3}
    grid::Vector{Float64}
    edges::Vector{Int}
    control::Vector{Float64}
end

function reference_integrals!(m,control,r,rates,coordinates,grid,starts)
    RiverRouting.check_grid(control)
    first(control)==first(grid) && last(control)==last(grid) && all(in(grid),control) ||
        throw(ArgumentError("integral windows must partition the cohort grid"))
    J,K=river_curve_count(r),length(grid)-1
    edges=[searchsortedfirst(grid,t) for t in control]
    z=@variable(m,[1:J,1:2,1:K],base_name="reference_integral_$(r.name)")
    C=0.3.*z
    basis=Matrix{Float64}(I,J,J)
    for t in 1:length(control)-1,k in edges[t]:edges[t+1]-1
        dt=grid[k+1]-grid[k]
        centre=(grid[k]-control[t]+dt/2)/(control[t+1]-control[t])
        flux=Vector{Union{AffExpr,QuadExpr,NonlinearExpr}}(undef,J)
        for l in 1:J-1
            flux[l]=river_transfer_expression(m,r,rates[k],view(basis,:,l);
                coordinate=coordinates[k])
        end
        flux[J]=1.0*rates[k]-sum((flux[l] for l in 1:J-1);init=AffExpr(0.0))
        for l in 1:J
            prior0=k==edges[t] ? AffExpr(0.0) : C[l,1,k-1]
            prior1=k==edges[t] ? AffExpr(0.0) : C[l,2,k-1]
            @constraint(m,(C[l,1,k]-prior0-0.0036dt*flux[l])/(0.0036dt*50)==0)
            @constraint(m,(C[l,2,k]-prior1-centre*(C[l,1,k]-prior0))/(0.0036dt*50)==0)
            push!(starts,z[l,1,k]=>(prior0+0.0036dt*flux[l])/0.3)
            push!(starts,z[l,2,k]=>(prior1+centre*(C[l,1,k]-prior0))/0.3)
        end
    end
    ReferenceIntegrals(C,z,grid,edges,control)
end

function reference_primitive(data,l,n,t,time)
    left,right=data.control[t],data.control[t+1]
    time<=left && return AffExpr(0.0)
    k=clamp(searchsortedlast(data.grid,time),data.edges[t],data.edges[t+1]-1)
    prior=k==data.edges[t] ? AffExpr(0.0) : data.C[l,n,k-1]
    before=k==data.edges[t] ? AffExpr(0.0) : data.C[l,1,k-1]
    fraction=(time-data.grid[k])/(data.grid[k+1]-data.grid[k])
    weight=n==1 ? fraction : fraction*((time-left)+(data.grid[k]-left))/(2(right-left))
    prior+weight*(data.C[l,1,k]-before)
end

function reference_integral(data,l,a,b;weights=(1.0,1.0))
    b>a || return AffExpr(0.0)
    left,right=max(a,first(data.control)),min(b,last(data.control))
    right>left || return AffExpr(0.0)
    lo=clamp(searchsortedlast(data.control,left),1,length(data.control)-1)
    hi=clamp(searchsortedfirst(data.control,right)-1,lo,length(data.control)-1)
    wa,wb=weights;slope=(wb-wa)/(b-a)
    value=AffExpr(0.0)
    for t in lo:hi
        u,v=max(left,data.control[t]),min(right,data.control[t+1])
        volume=reference_primitive(data,l,1,t,v)-reference_primitive(data,l,1,t,u)
        value+=(wa+slope*(data.control[t]-a))*volume
        if slope!=0
            moment=reference_primitive(data,l,2,t,v)-reference_primitive(data,l,2,t,u)
            value+=slope*(data.control[t+1]-data.control[t])*moment
        end
    end
    value
end

function reference_arrival_volume(data,r,a,b)
    value=AffExpr(0.0)
    for (l,curve) in enumerate(r.curves),j in eachindex(curve.weights)
        w=curve.weights[j];w==0 && continue
        u,v=curve.edges[j],curve.edges[j+1]
        left,right=a-v,b-u
        x,y=minmax(a-u,b-v)
        peak=min(b-a,v-u)/(v-u)
        value+=w*peak*(reference_integral(data,l,left,x;weights=(0.0,1.0))+
            reference_integral(data,l,x,y)+
            reference_integral(data,l,y,right;weights=(1.0,0.0)))
    end
    value
end

function reference_arrival_rate(data,r,time)
    value=AffExpr(0.0)
    for (l,curve) in enumerate(r.curves),j in eachindex(curve.weights)
        w=curve.weights[j];w==0 && continue
        a,b=curve.edges[j],curve.edges[j+1]
        value+=w/(0.0036*(b-a))*reference_integral(data,l,time-b,time-a)
    end
    value
end

function reference_transit(data,r,time)
    value=AffExpr(0.0)
    for (l,curve) in enumerate(r.curves),j in eachindex(curve.weights)
        w=curve.weights[j];w==0 && continue
        a,b=curve.edges[j],curve.edges[j+1]
        value+=w*(reference_integral(data,l,time-b,time-a;weights=(0.0,1.0))+
            reference_integral(data,l,time-a,time))
    end
    value
end

"""Assign integral variables from physical cohorts for a native global start."""
function start_reference_integrals!(data,r,release,put)
    J=size(data.variables,1)
    volume=zeros(J);moment=zeros(J);flux=zeros(J)
    for t in 1:length(data.control)-1
        fill!(volume,0.0);fill!(moment,0.0)
        for k in data.edges[t]:data.edges[t+1]-1
            dt=data.grid[k+1]-data.grid[k]
            centre=(data.grid[k]-data.control[t]+dt/2)/
                (data.control[t+1]-data.control[t])
            reference_flux!(flux,r,release[k])
            for l in 1:J
                increment=0.0036dt*flux[l]
                volume[l]+=increment
                moment[l]+=centre*increment
                put(data.variables[l,1,k],volume[l]/0.3)
                put(data.variables[l,2,k],moment[l]/0.3)
            end
        end
    end
    nothing
end


# Global search needs finite domains. These are the integrals of the maximum
# admissible reference flow over each control window, including its time moment.
function bound_reference_integrals!(c,routing)
    routing.exact && return
    for (i,data) in enumerate(routing.primitives)
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

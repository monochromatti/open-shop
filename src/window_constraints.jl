# Window integrals refer to the piecewise-constant observation, not instantaneous
# river arrivals. Storage ranges refer to vertices of the piecewise-linear V(t).
window_active(c,w) = w.start_time<last(c.grid) && w.end_time>first(c.grid)

function check_window(a,b)
    isfinite(a) && isfinite(b) && b>a && isfinite(b-a) ||
        throw(ArgumentError("window endpoints must define a finite positive duration"))
end

function check_window_history(grid,values,start;vertices=false)
    isempty(grid) && isempty(values) && return nothing
    length(grid)==length(values)+(vertices ? 0 : 1) && !isempty(grid) &&
        all(isfinite,grid) && all(isfinite,values) && all(diff(grid).>0) &&
        last(grid)==start || throw(ArgumentError(
            "window history must be finite, ordered, and end at the horizon start"))
    nothing
end

function require_window_history(object,a,start)
    a>=start && return nothing
    !isempty(object.history_grid) && first(object.history_grid)<=a ||
        throw(ArgumentError("$(object.name): window needs history back to $a hours"))
    nothing
end

function check_window_ramp(r)
    r.up!==nothing || r.down!==nothing ||
        throw(ArgumentError("window ramp needs an upward or downward limit"))
    all(x->x===nothing || isfinite(x) && x>=0,(r.up,r.down)) ||
        throw(ArgumentError("window ramp differences must be finite and nonnegative"))
end

function validate_window_rules(c)
    RiverRouting.check_grid(c.grid)
    start,stop=first(c.grid),last(c.grid)
    for r in c.system.reservoirs
        check_window_history(r.history_grid,r.history_storage,start;vertices=true)
        if !isempty(r.history_storage)
            all(>=(0),r.history_storage) && last(r.history_storage)==r.v0 ||
                throw(ArgumentError("$(r.name): storage history must be nonnegative and end at v0"))
            cv=r.level_curve
            cv===nothing || first(cv.x)<=minimum(r.history_storage)<=maximum(r.history_storage)<=last(cv.x) ||
                throw(ArgumentError("$(r.name): level table does not cover storage history"))
            if cv===nothing
                a,b=extrema(vcat(r.vmin,r.vmax,r.history_storage))
                min(r.slope+2r.curvature*a,r.slope+2r.curvature*b)>0 ||
                    throw(ArgumentError("$(r.name): level law must increase over storage history"))
            end
        end
        for w in r.level_windows
            check_window(w.start_time,w.end_time)
            isfinite(w.max_change) && w.max_change>=0 ||
                throw(ArgumentError("level range must be finite and nonnegative"))
            window_active(c,w) || continue
            require_window_history(r,w.start_time,start)
            for time in (w.start_time,w.end_time)
                start<=time<=stop && time ∉ c.grid &&
                    throw(ArgumentError("grid must include reservoir level-window endpoints"))
            end
        end
    end
    for rule in c.flow_requirements
        check_window_history(rule.history_grid,rule.history_flow,start)
        for w in rule.volume_windows
            check_window(w.start_time,w.end_time)
            w.min_volume!==nothing || w.max_volume!==nothing ||
                throw(ArgumentError("volume window needs a lower or upper bound"))
            all(x->x===nothing || isfinite(x),(w.min_volume,w.max_volume)) &&
                (w.min_volume===nothing || w.max_volume===nothing || w.min_volume<=w.max_volume) ||
                throw(ArgumentError("volume window bounds must be finite and ordered"))
            window_active(c,w) || continue
            w.end_time<=stop || throw(ArgumentError("active volume window must finish within the horizon"))
            require_window_history(rule,w.start_time,start)
        end
        for r in rule.rolling_ramps
            check_window_ramp(r)
            isfinite(r.hours) && r.hours>0 && isfinite(start-r.hours) &&
                isfinite(stop+r.hours) || throw(ArgumentError("rolling duration must be finite and positive"))
            require_window_history(rule,start-r.hours,start)
        end
        for r in rule.period_ramps
            check_window_ramp(r)
            RiverRouting.check_grid(r.grid)
            length(r.grid)>=3 || throw(ArgumentError("period ramp needs at least two periods"))
            for k in 1:(length(r.grid)-1)
                a,b=r.grid[k],r.grid[k+1]
                check_window(a,b)
                a<stop && b>start || continue
                b<=stop || throw(ArgumentError("active ramp period must finish within the horizon"))
                require_window_history(rule,r.grid[max(1,k-1)],start)
            end
        end
    end
    nothing
end

# Independent numerical integration, also used for historical constants. No
# resampling, shortened windows or equal weights for unequal intervals.
function window_integral(grid,values,a,b)
    sum(max(0.0,min(b,grid[k+1])-max(a,grid[k]))*values[k]
        for k in eachindex(values);init=0.0)
end

function flow_window_volume(c,rule,q,a,b)
    0.0036*(window_integral(c.grid,q,a,b)+
        window_integral(rule.history_grid,rule.history_flow,a,b))
end

"""Cumulative observed volume in Mm³, relative to the horizon start.
Current times use two neighboring cumulative values (or one and the local
rate). Historical times are known constants; historical water is counted once.
"""
function flow_integral_at(c,rule,q,volume,time)
    time<first(c.grid) && return -0.0036*window_integral(
        rule.history_grid,rule.history_flow,time,first(c.grid))
    time==last(c.grid) && return last(volume)
    k=searchsortedlast(c.grid,time)
    previous=k==1 ? 0.0 : volume[k-1]
    previous+0.0036*(time-c.grid[k])*q[k]
end

function rolling_knots(c,rule,hours)
    start,stop=first(c.grid),last(c.grid)
    sort!(unique!(vcat(c.grid,[t+hours for t in Iterators.flatten((rule.history_grid,c.grid))
        if start<t+hours<stop])))
end

function rolling_intervals(grid,time)
    time==last(grid) && return (length(grid)-1,)
    k=searchsortedlast(grid,time)
    k>1 && grid[k]==time ? (k-1,k) : (k,)
end

function constrain_window_difference!(m,change,up,down)
    if up==0 && down==0
        @constraint(m,change==0)
    else
        up===nothing || @constraint(m,change<=up)
        down===nothing || @constraint(m,-change<=down)
    end
end

function observed_flow_bounds(c,rule,t)
    cap=sum(opinterval(c,n,:capacity,t,only(e.capacity for e in c.system.tunnels if e.name==n))
        for n in rule.tunnels;init=0.0)
    high=cap+opinterval(c,rule.name,:inflow,t,rule.inflow)
    for g in c.system.generators
        g.name in rule.generators && (high+=opinterval(c,g.name,:qmax,t,g.qmax))
    end
    for r in c.system.rivers
        if r.name in rule.rivers
            bound=opinterval(c,r.name,:capacity,t,r.capacity)
            high+=bound
            is_channel(r) && (cap+=bound)
        end
        # Arrival rates can exceed the local current release capacity after a
        # delay mixture. Bound each historical/current cohort by its released
        # volume divided by the receiving interval; routing is conservative.
        if r.name in rule.river_arrivals
            if is_channel(r)
                bound=opinterval(c,r.name,:capacity,t,r.capacity)
                high+=bound
                cap+=bound
                continue
            end
            # Upstream outfalls/confluences may carry water whose original
            # history lives on other reaches. All historical water is a safe
            # finite bound, independent of network topology.
            released=0.0036*r.capacity*(last(c.grid)-first(c.grid))+
                sum(0.0036*sum(diff(other.history_grid).*other.history_release;init=0.0)
                    for other in c.system.rivers if !is_channel(other);init=0.0)
            high+=released/(0.0036*(c.grid[t+1]-c.grid[t]))
        end
    end
    low=max(-cap+opinterval(c,rule.name,:inflow,t,rule.inflow),
        opinterval(c,rule.name,:min_flow,t,rule.min_flow))
    upper=optional_operation(c,rule.name,:max_flow,t;default=rule.max_flow)
    upper===nothing || (high=min(high,upper))
    low<=high || throw(ArgumentError("$(rule.name): flow bounds exclude every observation"))
    low,high
end

function constrain_flow_windows!(m,c,rule,q)
    isempty(rule.volume_windows) && isempty(rule.rolling_ramps) && isempty(rule.period_ramps) &&
        return nothing
    T=length(q)
    volume=VariableRef[]
    lower,upper=0.0,0.0
    for t in 1:T
        lo,hi=observed_flow_bounds(c,rule,t)
        step=0.0036*(c.grid[t+1]-c.grid[t])
        lower+=step*lo
        upper+=step*hi
        isfinite(lower) && isfinite(upper) || throw(ArgumentError("nonfinite observed-volume domain"))
        v=@variable(m,lower_bound=lower,upper_bound=upper,start=clamp(0.0,lower,upper),
            base_name="observed_volume_$(rule.name)_$t")
        previous=t==1 ? 0.0 : volume[t-1]
        @constraint(m,v==previous+step*q[t])
        push!(volume,v)
    end
    at(time)=flow_integral_at(c,rule,q,volume,time)
    for w in rule.volume_windows
        window_active(c,w) || continue
        v=at(w.end_time)-at(w.start_time)
        if w.min_volume!==nothing && w.min_volume==w.max_volume
            @constraint(m,v==w.min_volume)
        else
            w.min_volume===nothing || @constraint(m,v>=w.min_volume)
            w.max_volume===nothing || @constraint(m,v<=w.max_volume)
        end
    end
    for r in rule.rolling_ramps, time in rolling_knots(c,rule,r.hours)
        average=(at(time)-at(time-r.hours))/(0.0036*r.hours)
        for k in rolling_intervals(c.grid,time)
            constrain_window_difference!(m,q[k]-average,r.up,r.down)
        end
    end
    for r in rule.period_ramps, k in 2:(length(r.grid)-1)
        a,b,d=r.grid[k-1],r.grid[k],r.grid[k+1]
        b<last(c.grid) && d>first(c.grid) || continue
        change=(at(d)-at(b))/(0.0036*(d-b))-(at(b)-at(a))/(0.0036*(b-a))
        constrain_window_difference!(m,change,r.up,r.down)
    end
    volume
end

function level_window_history(c,r,w)
    w.start_time>=first(c.grid) && return Float64[]
    a,b=w.start_time,min(w.end_time,first(c.grid))
    times=vcat(a,[t for t in r.history_grid if a<t<b],b)
    [begin
        k=min(searchsortedlast(r.history_grid,time),length(r.history_grid)-1)
        k==0 && throw(ArgumentError("incomplete storage history"))
        fraction=(time-r.history_grid[k])/(r.history_grid[k+1]-r.history_grid[k])
        head(r,(1-fraction)*r.history_storage[k]+fraction*r.history_storage[k+1])
    end for time in times]
end

function level_window_vertices(c,w)
    a=searchsortedfirst(c.grid,max(w.start_time,first(c.grid)))
    b=searchsortedlast(c.grid,min(w.end_time,last(c.grid)))
    a:b
end

function constrain_level_windows!(m,c,level)
    bands=NamedTuple{(:reservoir,:window,:floor),Tuple{Int,Int,VariableRef}}[]
    for (i,r) in enumerate(c.system.reservoirs), (j,w) in enumerate(r.level_windows)
        window_active(c,w) || continue
        history=level_window_history(c,r,w)
        vertices=level_window_vertices(c,w)
        first(vertices)==1 && push!(history,head(r,r.v0))
        lo,hi=head(r,r.vmin)-w.max_change,head(r,r.vmax)
        if !isempty(history)
            lo=max(lo,maximum(history)-w.max_change)
            hi=min(hi,minimum(history))
        end
        lo<=hi || throw(ArgumentError("$(r.name): known history violates the level range"))
        # Use the same dimensionless level coordinate as the hydraulic heads.
        floor=@variable(m,lower_bound=lo/250,upper_bound=hi/250,
            start=clamp(head(r,r.v0),lo,hi)/250,
            base_name="level_window_$(r.name)_$j")
        for t in vertices
            t==1 && continue # Known initial surface already bounds the band.
            z=level(i,t)
            if w.max_change==0
                @constraint(m,z/250-floor==0)
            else
                @constraint(m,z/250-floor>=0)
                @constraint(m,z/250-floor<=w.max_change/250)
            end
        end
        push!(bands,(reservoir=i,window=j,floor=floor))
    end
    bands
end

function window_start_values(c,x,integrals,bands)
    pairs=Pair{VariableRef,Float64}[]
    for data in integrals
        rule=c.flow_requirements[data.index]
        q=observed_flow(c,rule,x.generator_q,x.river_release,x.arrival_volume,x.tunnel_q)
        total=0.0
        for (t,v) in enumerate(data.volume)
            total+=0.0036*(c.grid[t+1]-c.grid[t])*q[t]
            push!(pairs,v=>total)
        end
    end
    for band in bands
        r=c.system.reservoirs[band.reservoir]
        w=r.level_windows[band.window]
        levels=head.(Ref(r),x.storage[band.reservoir,level_window_vertices(c,w)])
        append!(levels,level_window_history(c,r,w))
        push!(pairs,band.floor=>minimum(levels)/250)
    end
    pairs
end

function window_residuals(c,x)
    result=Dict{String,Float64}()
    for (i,r) in enumerate(c.system.reservoirs), (j,w) in enumerate(r.level_windows)
        window_active(c,w) || continue
        z=head.(Ref(r),x.storage[i,level_window_vertices(c,w)])
        append!(z,level_window_history(c,r,w))
        result["level_window_$(r.name)_$j"]=max(0.0,maximum(z)-minimum(z)-w.max_change)
    end
    for rule in c.flow_requirements
        isempty(rule.volume_windows) && isempty(rule.rolling_ramps) && isempty(rule.period_ramps) &&
            continue
        q=observed_flow(c,rule,x.generator_q,x.river_release,x.arrival_volume,x.tunnel_q)
        for (j,w) in enumerate(rule.volume_windows)
            window_active(c,w) || continue
            v=flow_window_volume(c,rule,q,w.start_time,w.end_time)
            result["volume_window_$(rule.name)_$j"]=max(0.0,
                w.min_volume===nothing ? 0.0 : w.min_volume-v,
                w.max_volume===nothing ? 0.0 : v-w.max_volume)
        end
        for (j,r) in enumerate(rule.rolling_ramps)
            violation=0.0
            for time in rolling_knots(c,rule,r.hours)
                average=flow_window_volume(c,rule,q,time-r.hours,time)/(0.0036*r.hours)
                for k in rolling_intervals(c.grid,time)
                    change=q[k]-average
                    violation=max(violation,r.up===nothing ? 0.0 : change-r.up,
                        r.down===nothing ? 0.0 : -change-r.down)
                end
            end
            result["rolling_ramp_$(rule.name)_$j"]=violation
        end
        for (j,r) in enumerate(rule.period_ramps)
            violation=0.0
            for k in 2:(length(r.grid)-1)
                a,b,d=r.grid[k-1],r.grid[k],r.grid[k+1]
                b<last(c.grid) && d>first(c.grid) || continue
                change=flow_window_volume(c,rule,q,b,d)/(0.0036*(d-b))-
                    flow_window_volume(c,rule,q,a,b)/(0.0036*(b-a))
                violation=max(violation,r.up===nothing ? 0.0 : change-r.up,
                    r.down===nothing ? 0.0 : -change-r.down)
            end
            result["period_ramp_$(rule.name)_$j"]=violation
        end
    end
    result
end

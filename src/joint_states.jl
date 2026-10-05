function _joint_states!(
    m,
    c;
    relaxed = false,
    fixed_u = nothing,
    free_mask = nothing,
    incumbent = nothing,
)
    G=length(c.system.generators)
    T=length(c.prices)
    start=incumbent===nothing ? ones(Int, G, T) : copy(incumbent["u"])
    size(start)==(G, T) || throw(ArgumentError("incorrect incumbent shape"))
    if fixed_u!==nothing
        admissible(c, fixed_u) || throw(ArgumentError("invalid fixed commitment"))
        start=copy(fixed_u)
    end
    if free_mask!==nothing
        incumbent===nothing &&
            throw(ArgumentError("a restricted search needs an incumbent"))
        size(free_mask)==(G, T) || throw(ArgumentError("incorrect free-mask shape"))
        all(x->x==true || x==false, free_mask) ||
            throw(ArgumentError("free-mask must be boolean"))
        admissible(c, start) || throw(ArgumentError("invalid incumbent commitment"))
    end
    @variable(m, 0<=u[1:G, 1:T]<=1)
    @variable(m, 0<=su[1:G, 1:T]<=1)
    @variable(m, 0<=sd[1:G, 1:T]<=1)
    for (j, g) in enumerate(c.system.generators), t in 1:T
        relaxed || set_binary(u[j, t])
        prev=t==1 ? g.initial_on : u[j, t - 1]
        @constraint(m, u[j, t]-prev==su[j, t]-sd[j, t])
        @constraint(m, su[j, t]<=u[j, t])
        @constraint(m, su[j, t]<=1-prev)
        @constraint(m, sd[j, t]<=prev)
        @constraint(m, sd[j, t]<=1-u[j, t])
        # State alone is binary: these inequalities make startup/shutdown exact.
        for k in t:T
            c.grid[k]-c.grid[t]<g.minup-1e-8 && @constraint(m, u[j, k]>=su[j, t])
            c.grid[k]-c.grid[t]<g.mindown-1e-8 && @constraint(m, u[j, k]<=1-sd[j, t])
        end
        residual=(g.initial_on==1 ? g.minup : g.mindown)-g.initial_age
        if c.grid[t]-first(c.grid)<residual-1e-8
            fix(u[j, t], g.initial_on; force = true)
        end
        if fixed_u!==nothing || (free_mask!==nothing && !free_mask[j, t])
            fix(u[j, t], start[j, t]; force = true)
        end
        set_start_value(u[j, t], start[j, t])
        prior=t==1 ? g.initial_on : start[j, t - 1]
        set_start_value(su[j, t], max(0, start[j, t]-prior))
        set_start_value(sd[j, t], max(0, prior-start[j, t]))
    end
    (; u, su, sd, start)
end

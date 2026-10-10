# Numerical data and routing expressions for local dispatch construction.

"""Prepare capacity and gate bounds on the case's control intervals."""
function _dispatch_outlet_bounds(c)
    D, T = length(c.system.rivers), length(c.prices)
    capacity = Matrix{Float64}(undef, D, T)
    gate_lower = similar(capacity)
    gate_upper = similar(capacity)
    for (i, r) in enumerate(c.system.rivers), t in 1:T
        capacity[i, t] = opinterval(c, r.name, :capacity, t, r.capacity)
        gate_lower[i, t] = opinterval(c, r.name, :gate_min, t, r.gate_min)
        gate_upper[i, t] = opinterval(c, r.name, :gate_max, t, 1.0)
    end
    (; capacity, gate_lower, gate_upper)
end

function _dispatch_transport_expressions(c, release; transport=nothing)
    exact = all(r.deterministic_delay!==nothing for r in c.system.rivers)
    !exact && transport!==nothing && throw(ArgumentError(
        "compiled deterministic transport cannot serve distributed dispatch"))
    nd = exact ? _transport_data(c, transport) : nothing
    rd = exact ? nothing : routing_data(c)
    D, T = size(release)
    dt = diff(c.grid)
    arrivals = Matrix{Union{AffExpr,QuadExpr}}(undef, D, T)
    terminal = Union{AffExpr,QuadExpr}[]
    for (i, r) in enumerate(c.system.rivers)
        B = exact ? nothing : rd["B"][i]
        for t in 1:T
            if exact
                arrivals[i, t] = AffExpr(nd.arrival_history[i, t])+sum(
                    v*release[d, k] for (d, k, v) in nd.arrival_terms[i, t];
                    init = 0.0)
            else
                # Preserve the existing banded routing sparsity.
                ks = [k for k in 1:T if B[t, k, 1]!=0 || B[t, k, 2]!=0]
                arrivals[i, t] = AffExpr(rd["history_arrival"][i, t])+sum(
                    0.0036*dt[k]*release[i, k]*(B[t, k, 1]+
                    release[i, k]/r.capacity*(B[t, k, 2]-B[t, k, 1])) for k in ks;
                    init = 0.0)
            end
        end
        if exact
            push!(terminal, AffExpr(nd.terminal_history[i])+sum(
                v*release[d, k] for (d, k, v) in nd.terminal_terms[i]; init = 0.0))
        else
            push!(terminal, AffExpr(rd["history_terminal"][i])+sum(
                0.0036*dt[k]*release[i, k]*(1-sum(B[:, k, 1])-
                release[i, k]/r.capacity*sum(B[:, k, 2]-B[:, k, 1])) for k in 1:T;
                init = 0.0))
        end
    end
    (; exact, nd, rd, arrivals, terminal)
end

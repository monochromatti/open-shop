# Private experiment: couple water balances, passive tunnel losses and power.
# Exact hydraulic/table/commitment equations remain in the original model.
# For d = r*q*abs(q), Fenchel equality gives
#   r*abs(q)^3/3 + 2*abs(d)^(3/2)/(3*sqrt(r)) = q*d.
# Summation through incidence cancels internal node head-flow products. Storage
# uses the original midpoint head times the discrete volume change, never a
# continuous potential-energy integral. Supports below are global lower bounds.
using JuMP, OpenSHOP

const _NE_SCALE = 12500.0 # 250 m * 50 m³/s; auxiliary values stay near unity.

_ne_flow_content(q, r) = r * abs(q)^3 / 3
_ne_head_content(d, r) = 2 * abs(d)^(3 / 2) / (3 * sqrt(r))
_ne_flow_derivative(q, r) = r * q * abs(q)
_ne_head_derivative(d, r) = sign(d) * sqrt(abs(d) / r)

function _ne_variable_bounds(v)
    is_fixed(v) && return (fix_value(v), fix_value(v))
    has_lower_bound(v) && has_upper_bound(v) ||
        throw(ArgumentError("network energy needs bounded variables: $(name(v))"))
    lower_bound(v), upper_bound(v)
end

function _ne_affine_bounds(x)
    x isa Number && return (Float64(x), Float64(x))
    x isa VariableRef && return _ne_variable_bounds(x)
    x isa GenericAffExpr || return nothing
    lo = hi = Float64(x.constant)
    magnitude = abs(lo)
    free = false
    for (v, a) in x.terms
        l, h = _ne_variable_bounds(v)
        al, ah = extrema((a * l, a * h))
        lo += al
        hi += ah
        magnitude += max(abs(al), abs(ah))
        free |= l != h
    end
    # Consistent with the project's numerical bound scope, not directed rounding.
    pad = free ? 1e-10 * max(1.0, magnitude) : 0.0
    lo - pad, hi + pad
end

function _ne_incidence(b, c, i, t)
    s = c.system
    n = OpenSHOP.nodes(s)[i]
    out = AffExpr(0.0)
    for (j, e) in enumerate(s.tunnels)
        a = Int(e.source == n) - Int(e.target == n)
        a != 0 && (out += a * b.Q[j, t])
    end
    for (j, g) in enumerate(s.generators)
        p = OpenSHOP.plantof(s, g)
        a = Int(p.source == n) - Int(p.target == n)
        a != 0 && (out += a * b.GQ[j, t])
    end
    out
end

function _ne_external_injection(b, c, i, t)
    s = c.system
    n = OpenSHOP.nodes(s)[i]
    alpha = 0.0036 * (c.grid[t + 1] - c.grid[t])
    # Boundaries permit external exchange. Their head-weighted contribution is
    # exactly their hydraulic incidence; adding river water again double-counts it.
    i > length(s.reservoirs) + length(s.junctions) &&
        return _ne_incidence(b, c, i, t)
    out = AffExpr(0.0)
    if i <= length(s.reservoirs)
        r = s.reservoirs[i]
        out += OpenSHOP.opinterval(c, r.name, :inflow, t, r.inflow)
        out -= (b.V[i, t + 1] - b.V[i, t]) / alpha
    end
    for (j, r) in enumerate(s.rivers)
        # arrivals are volumes (Mm³), not release rates or same-interval releases.
        r.target == n && (out += b.arrivals[j, t] / alpha)
        r.source == n && (out -= b.RQ[j, t])
    end
    out
end

function _ne_components(b, c, t)
    s = c.system
    ix = OpenSHOP.nodeindex(s)
    adj = [Int[] for _ in OpenSHOP.nodes(s)]
    function connect(a, z)
        push!(adj[a], z)
        push!(adj[z], a)
    end
    for (j, e) in enumerate(s.tunnels)
        opening = OpenSHOP.opinterval(c, e.name, :opening, t, e.opening)
        opening == 0 && continue
        lo, hi = _ne_affine_bounds(b.Q[j, t])
        lo == hi == 0 && continue
        connect(ix[e.source], ix[e.target])
    end
    for (j, g) in enumerate(s.generators)
        lo, hi = _ne_affine_bounds(b.GQ[j, t])
        lo == hi == 0 && continue
        p = OpenSHOP.plantof(s, g)
        connect(ix[p.source], ix[p.target])
    end
    seen = falses(length(adj))
    groups = Vector{Int}[]
    for i in eachindex(adj)
        (seen[i] || isempty(adj[i])) && continue
        group = [i]
        seen[i] = true
        for node in group
            for other in adj[node]
                seen[other] && continue
                seen[other] = true
                push!(group, other)
            end
        end
        push!(groups, sort!(group))
    end
    groups
end

function _ne_curve_range(curve, lo, hi)
    vals = [OpenSHOP.table_value(curve, x; extrapolation = :linear)
            for x in OpenSHOP._global_tensor_nodes(curve.x, lo, hi)]
    extrema(vals)
end

function _ne_efficiency_upper(b, c, j, t)
    g = c.system.generators[j]
    lo, hi = _ne_affine_bounds(b.GQ[j, t])
    lo == hi == 0 && return (0.0, nothing)
    pl = OpenSHOP.plantof(c.system, g)
    qmax = sum(min(z.qmax, OpenSHOP.opinterval(c, z.name, :qmax, t, z.qmax))
               for z in c.system.generators if z.plant == pl.name; init = 0.0)
    tailmin = pl.tailwater_curve === nothing ? 0.0 :
              first(_ne_curve_range(pl.tailwater_curve, 0.0, qmax))
    dsthi = b.node_bounds[(pl.target, t)][2]
    floormin = pl.outlet_head_floor === nothing ? 0.0 :
               max(0.0, pl.outlet_head_floor - dsthi)
    tailmin + floormin >= 0 || return (0.0, "negative_extra_plant_head_loss")
    pmin = OpenSHOP.opinterval(c, g.name, :pmin, t, g.pmin)
    pmax = min(g.pmax, OpenSHOP.opinterval(c, g.name, :pmax, t, g.pmax))
    # Positive power and nonnegative efficiencies imply positive net head.
    # At zero power an allowed zero efficiency does not: q>0,h<0,P=0 is
    # possible. In that domain P/gamma <= q*network_drop is not valid.
    hd = b.shared_heads[(pl.name, t)]
    max(g.hmin, first(_ne_variable_bounds(hd))) >= 0 || pmin > 0 ||
        return (0.0, "possible_negative_running_head_at_zero_power")
    pmin <= pmax || return (0.0, "empty_running_power_domain")
    elo, ehi = g.generator_efficiency_curve === nothing ? (1.0, 1.0) :
                _ne_curve_range(g.generator_efficiency_curve, pmin, pmax)
    elo >= 0 || return (0.0, "negative_running_electrical_efficiency")
    # Running turbine efficiency is already constrained to [min_efficiency,1].
    # A tighter full-domain variable upper bound is safe even with off extrapolation.
    eta = variable_by_name(b.m, g.turbine_table === nothing ?
                           "eta_$(j)_$(t)" : "turbine_$(j)_$(t)")
    eta_hi = eta === nothing ? 1.0 : min(1.0, last(_ne_variable_bounds(eta)))
    gamma = eta_hi * ehi
    gamma > 0 && isfinite(gamma) || return (0.0, "nonpositive_efficiency_upper")
    gamma * (1 + 1e-10), nothing
end

function _ne_plan(b, c)
    s = c.system
    ix = OpenSHOP.nodeindex(s)
    records = Any[]
    skipped = Any[]
    closed = 0
    for t in eachindex(c.prices)
        closed += count(e -> OpenSHOP.opinterval(c, e.name, :opening, t, e.opening) == 0,
                        s.tunnels)
        for group in _ne_components(b, c, t)
            members = Set(group)
            gens = Any[]
            reason = nothing
            for (j, g) in enumerate(s.generators)
                p = OpenSHOP.plantof(s, g)
                ix[p.source] in members || continue
                gamma, why = _ne_efficiency_upper(b, c, j, t)
                if why !== nothing
                    reason = (unit = string(g.name), reason = why)
                    break
                end
                gamma > 0 && push!(gens, (j = j, gamma = gamma))
            end
            if reason !== nothing
                push!(skipped, Dict("interval" => t, "nodes" => string.(OpenSHOP.nodes(s)[group]),
                                   "unit" => reason.unit, "reason" => reason.reason))
                continue
            end
            tunnels = Any[]
            for (j, e) in enumerate(s.tunnels)
                ix[e.source] in members && ix[e.target] in members || continue
                opening = OpenSHOP.opinterval(c, e.name, :opening, t, e.opening)
                opening == 0 && continue
                qlo, qhi = _ne_affine_bounds(b.Q[j, t])
                qlo == qhi == 0 && continue
                d = b.H[ix[e.source], t] - b.H[ix[e.target], t]
                dlo, dhi = _ne_affine_bounds(d)
                r = e.resistance / opening
                # The exact monotone tunnel law gives a much smaller support
                # interval than the independent node head boxes in many cases.
                lawlo, lawhi = r*qlo*abs(qlo), r*qhi*abs(qhi)
                pad = 1e-10*max(1.0, abs(lawlo), abs(lawhi))
                dlo, dhi = max(dlo, lawlo-pad), min(dhi, lawhi+pad)
                dlo <= dhi || error("network energy tunnel head domain is empty")
                push!(tunnels, (j = j, r = r, d = d,
                                qlo = qlo, qhi = qhi, dlo = dlo, dhi = dhi))
            end
            injections = Any[]
            for i in group
                external = _ne_external_injection(b, c, i, t)
                incidence = _ne_incidence(b, c, i, t)
                lo, hi = _ne_affine_bounds(incidence)
                bounds = _ne_affine_bounds(external)
                if bounds !== nothing
                    lo = max(lo, bounds[1])
                    hi = min(hi, bounds[2])
                end
                lo <= hi || error("network energy injection domain is empty")
                push!(injections, (i = i, external = external, incidence = incidence,
                                   lo = lo, hi = hi, affine = bounds !== nothing))
            end
            push!(records, (t = t, nodes = group, generators = gens, tunnels = tunnels,
                            injections = injections))
        end
    end
    (; records, skipped, closed)
end

function _ne_support_epigraph!(m, x, lo, hi, r, kind, label, points, lifts)
    f = kind == :flow ? _ne_flow_content : _ne_head_content
    derivative = kind == :flow ? _ne_flow_derivative : _ne_head_derivative
    if lo == hi
        return f(lo, r), 0
    end
    maximum_value = max(f(lo, r), f(hi, r))
    z = @variable(m, lower_bound = 0.0,
                  upper_bound = maximum_value / _NE_SCALE + 1e-10,
                  base_name = label)
    push!(lifts, (z, assigned -> f(value(v -> assigned[v], x), r) / _NE_SCALE))
    count = 0
    for anchor in unique(vcat(collect(range(lo, hi; length = points)),
                              lo <= 0 <= hi ? [0.0] : Float64[]))
        anchor == 0 && continue # z >= 0 already supplies this tangent.
        fv = f(anchor, r)
        slope = derivative(anchor, r)
        pad = 1e-10 * max(1.0, abs(fv), abs(slope) * max(abs(lo), abs(hi)))
        @constraint(m, z >= (fv + slope * (x - anchor) - pad) / _NE_SCALE)
        count += 1
    end
    _NE_SCALE * z, count
end

function _ne_product!(b, c, injection, t, reference, label, lifts)
    i = injection.i
    h = b.H[i, t] - reference
    hlo, hhi = _ne_affine_bounds(h)
    slo, shi = injection.lo, injection.hi
    s = injection.external
    extra = 0
    if !injection.affine
        # Distributed routing gives a polynomial arrival expression. Keep it,
        # and share a bounded rate with the affine hydraulic incidence equation.
        z = @variable(b.m, lower_bound = slo / 50, upper_bound = shi / 50,
                      base_name = label * "_injection")
        @constraint(b.m, z == injection.external / 50)
        @constraint(b.m, z == injection.incidence / 50)
        push!(lifts, (z, assigned -> value(v -> assigned[v], injection.external) / 50))
        s = 50 * z
        extra = 2
    end
    hlo == hhi && return (hlo * s, extra, 0)
    slo == shi && return (slo * h, extra, 0)
    corners = [hlo * slo, hlo * shi, hhi * slo, hhi * shi] ./ _NE_SCALE
    pad = 1e-10 * max(1.0, maximum(abs, corners))
    w = @variable(b.m, lower_bound = minimum(corners) - pad,
                  upper_bound = maximum(corners) + pad, base_name = label)
    @constraint(b.m, w >= (hlo * s + slo * h - hlo * slo) / _NE_SCALE - pad)
    @constraint(b.m, w >= (hhi * s + shi * h - hhi * shi) / _NE_SCALE - pad)
    @constraint(b.m, w <= (hhi * s + slo * h - hhi * slo) / _NE_SCALE + pad)
    @constraint(b.m, w <= (hlo * s + shi * h - hlo * shi) / _NE_SCALE + pad)
    push!(lifts, (w, assigned -> value(v -> assigned[v], h) *
                                  value(v -> assigned[v], s) / _NE_SCALE))
    _NE_SCALE * w, extra, 4
end

"""Add valid linear-support network energy rows; exact equations remain intact."""
function add_network_energy!(b, c; support_points = 5, head_reference = nothing)
    support_points isa Integer && support_points >= 2 ||
        throw(ArgumentError("at least two support points required"))
    head_reference === nothing || isfinite(head_reference) ||
        throw(ArgumentError("finite head reference required"))
    haskey(b.m.ext, :network_energy_data) && error("network energy already installed")
    plan = _ne_plan(b, c)
    before = num_variables(b.m)
    lifts = Any[]
    rows = Any[]
    supports = products = injection_rows = 0
    for (k, record) in enumerate(plan.records)
        t = record.t
        reference = head_reference === nothing ?
                    (minimum(b.node_bounds[(OpenSHOP.nodes(c.system)[i], t)][1] for i in record.nodes) +
                     maximum(b.node_bounds[(OpenSHOP.nodes(c.system)[i], t)][2] for i in record.nodes)) / 2 :
                    Float64(head_reference)
        lhs = AffExpr(0.0)
        for tunnel in record.tunnels
            flow, n = _ne_support_epigraph!(b.m, b.Q[tunnel.j, t], tunnel.qlo, tunnel.qhi,
                tunnel.r, :flow, "energy_flow_$(k)_$(tunnel.j)", support_points, lifts)
            lhs += flow
            supports += n
            head, n = _ne_support_epigraph!(b.m, tunnel.d, tunnel.dlo, tunnel.dhi,
                tunnel.r, :head, "energy_head_$(k)_$(tunnel.j)", support_points, lifts)
            lhs += head
            supports += n
        end
        for g in record.generators
            lhs += b.P[g.j, t] / (0.00981 * g.gamma)
        end
        rhs = AffExpr(0.0)
        for injection in record.injections
            w, n, p = _ne_product!(b, c, injection, t, reference,
                                  "energy_product_$(k)_$(injection.i)", lifts)
            rhs += w
            injection_rows += n
            products += p
        end
        cut = @constraint(b.m, (lhs - rhs) / _NE_SCALE <= 1e-9)
        push!(rows, (; record, reference, lhs, rhs, cut))
    end
    b.m.ext[:network_energy_data] = (; plan, lifts, rows)
    profile = Dict("scope" => "valid redundant inequalities for original scheduling-grid equations",
        "representation" => "linear supporting epigraphs plus bounded McCormick injection products",
        "midpoint_storage" => "original H(Vmid)*(Vnext-V)/duration; no continuous potential substitution",
        "support_points" => support_points, "energy_scale" => _NE_SCALE,
        "energy_rows" => length(rows), "support_rows" => supports,
        "mccormick_rows" => products, "distributed_injection_rows" => injection_rows,
        "auxiliary_variables" => num_variables(b.m) - before,
        "closed_tunnel_intervals_excluded" => plan.closed,
        "skipped_component_intervals" => plan.skipped,
        "running_efficiency_upper" => "min(1,turbine variable upper)*electrical on-domain upper",
        "numeric_scope" => "floating-point guards; not a directed-rounding proof")
    b.m.ext[:network_energy_profile] = profile
    push!(get!(b.m.ext, :experiment_start_lifters, Any[]),
          assigned -> lift_network_energy_start!(b, assigned))
    profile
end

"""Fill every experiment auxiliary before the existing complete start audit."""
function lift_network_energy_start!(b, assigned)
    data = get(b.m.ext, :network_energy_data, nothing)
    data === nothing && return assigned
    for (v, evaluate) in data.lifts
        assigned[v] = Float64(evaluate(assigned))
        isfinite(assigned[v]) || error("nonfinite network energy lift $(name(v))")
    end
    assigned
end

"""Evaluate original energy at a root point, independently of optimistic lifts."""
function network_energy_diagnostic(b, c, at)
    data = get(b.m.ext, :network_energy_data, nothing)
    plan = data === nothing ? _ne_plan(b, c) : data.plan
    details = Any[]
    for (k, record) in enumerate(plan.records)
        t = record.t
        reference = data === nothing ?
            (minimum(b.node_bounds[(OpenSHOP.nodes(c.system)[i], t)][1] for i in record.nodes) +
             maximum(b.node_bounds[(OpenSHOP.nodes(c.system)[i], t)][2] for i in record.nodes))/2 :
            data.rows[k].reference
        flow = sum(_ne_flow_content(at(b.Q[e.j, t]), e.r) for e in record.tunnels; init = 0.0)
        head = sum(_ne_head_content(at(e.d), e.r) for e in record.tunnels; init = 0.0)
        power = sum(at(b.P[g.j, t]) / (0.00981 * g.gamma)
                    for g in record.generators; init = 0.0)
        supplied = sum((at(b.H[i.i, t]) - reference) * at(i.external)
                       for i in record.injections; init = 0.0)
        residual = flow + head + power - supplied
        row = Dict("interval" => t, "component" => k,
            "head_reference_m" => reference,
            "flow_content" => flow, "head_content" => head,
            "generation_hydraulic_lower" => power, "exact_head_injection" => supplied,
            "original_energy_residual" => residual,
            "original_energy_positive_violation" => max(0.0, residual))
        if data !== nothing
            row["lifted_support_row_residual"] = at(data.rows[k].lhs) - at(data.rows[k].rhs)
        end
        push!(details, row)
    end
    Dict("scope" => "original energy inequality at diagnostic point; not a feasible schedule",
        "energy_rows" => length(details),
        "max_original_positive_violation" => maximum(
            row["original_energy_positive_violation"] for row in details; init = 0.0),
        "sum_original_positive_violation" => sum(
            row["original_energy_positive_violation"] for row in details; init = 0.0),
        "skipped_component_intervals" => plan.skipped,
        "largest_errors" => first(sort!(details; by = row -> abs(row["original_energy_residual"]),
                                        rev = true), min(8, length(details))))
end

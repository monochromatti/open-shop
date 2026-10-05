# Exact algebraic table graphs for global MINLP solvers. No registered callbacks
# or big-M implications are used. Bounds include the numerical linear extension;
# physical operating domains/envelopes remain the caller's responsibility.
function _global_intervals(grid, lo, hi)
    all(isfinite, (lo, hi)) && lo <= hi ||
        throw(ArgumentError("global solver table coordinates need finite ordered bounds"))
    knots = unique(vcat(Float64(lo), filter(v -> lo < v < hi, grid), Float64(hi)))
    length(knots) == 1 && return [(knots[1], knots[1])]
    [(knots[k], knots[k + 1]) for k in 1:(length(knots) - 1)]
end

function _global_table!(m, c::TableCurve, x, lo, hi; name = gensym(:table))
    cells = _global_intervals(c.x, lo, hi)
    n = length(cells)
    z = n == 1 ? [1.0] : @variable(m, [1:n], binary=true, base_name="$(name)_cell")
    n > 1 && @constraint(m, sum(z) == 1)
    # Two endpoint weights per clipped/extended segment give its exact graph.
    a = @variable(m, [1:n], lower_bound=0, upper_bound=1, base_name="$(name)_left")
    b = @variable(m, [1:n], lower_bound=0, upper_bound=1, base_name="$(name)_right")
    for k in 1:n
        @constraint(m, a[k] + b[k] == z[k])
    end
    @constraint(m, x == sum(cells[k][1]*a[k] + cells[k][2]*b[k] for k in 1:n))
    yy = [
        (
            table_value(c, l; extrapolation = :linear),
            table_value(c, r; extrapolation = :linear),
        ) for (l, r) in cells
    ]
    all(isfinite, (v for pair in yy for v in pair)) ||
        throw(ArgumentError("global solver table extension bounds overflow"))
    y = @variable(
        m,
        lower_bound=minimum(min(l, r) for (l, r) in yy),
        upper_bound=maximum(max(l, r) for (l, r) in yy),
        base_name=string(name)
    )
    @constraint(m, y == sum(yy[k][1]*a[k] + yy[k][2]*b[k] for k in 1:n))
    y
end

function _global_discharge_coefficients(c, i, j, cubic)
    a = c.efficiency[i, j]
    b = c.efficiency[i + 1, j]
    cubic || return (a, b-a, 0.0, 0.0)
    d = c.discharge[i + 1]-c.discharge[i]
    u = d*c.discharge_slopes[i, j]
    v = d*c.discharge_slopes[i + 1, j]
    (a, u, 3(b-a)-2u-v, 2(a-b)+u+v)
end

function _global_turbine!(
    m,
    c::TurbineTable,
    q,
    h,
    qlo,
    qhi,
    hlo,
    hhi;
    name = gensym(:turbine),
)
    qc = _global_intervals(c.discharge, qlo, qhi)
    hc = _global_intervals(c.heads, hlo, hhi)
    cells = [(qr, hr) for qr in qc for hr in hc]
    n = length(cells)
    z = n == 1 ? [1.0] : @variable(m, [1:n], binary=true, base_name="$(name)_cell")
    n > 1 && @constraint(m, sum(z) == 1)
    ts = JuMP.VariableRef[]
    us = JuMP.VariableRef[]
    qterms = Any[]
    hterms = Any[]
    terms = Any[]
    bounds = Float64[]
    for (k, ((ql, qr), (hl, hr))) in enumerate(cells)
        qm = ql + (qr-ql)/2
        hm = hl + (hr-hl)/2
        i = _curve_segment(c.discharge, qm, :linear)
        j = _curve_segment(c.heads, hm, :linear)
        dq = c.discharge[i + 1]-c.discharge[i]
        dh = c.heads[j + 1]-c.heads[j]
        tl = (ql-c.discharge[i])/dq
        tr = (qr-c.discharge[i])/dq
        ul = (hl-c.heads[j])/dh
        ur = (hr-c.heads[j])/dh
        t = @variable(
            m,
            lower_bound=min(0, tl),
            upper_bound=max(0, tr),
            base_name="$(name)_q[$k]"
        )
        u = @variable(
            m,
            lower_bound=min(0, ul),
            upper_bound=max(0, ur),
            base_name="$(name)_h[$k]"
        )
        @constraint(m, t >= tl*z[k])
        @constraint(m, t <= tr*z[k])
        @constraint(m, u >= ul*z[k])
        @constraint(m, u <= ur*z[k])
        push!(ts, t)
        push!(us, u)
        push!(qterms, c.discharge[i]*z[k]+dq*t)
        push!(hterms, c.heads[j]*z[k]+dh*u)
        cubic =
            c.interpolation == :pchip_discharge &&
            first(c.discharge) <= qm <= last(c.discharge)
        a = _global_discharge_coefficients(c, i, j, cubic)
        b = _global_discharge_coefficients(c, i, j+1, cubic)
        # Inactive cells have z=t=u=0, including the polynomial's constant.
        # Thus no division by a selector, implication, or relaxation is involved.
        left = a[1]*z[k]+a[2]*t+a[3]*t^2+a[4]*t^3
        delta = (b[1]-a[1])+(b[2]-a[2])*t+(b[3]-a[3])*t^2+(b[4]-a[4])*t^3
        push!(terms, left+u*delta)
        tm = max(abs(tl), abs(tr))
        um = max(abs(ul), abs(ur))
        bound =
            sum(abs(a[r])*tm^(r-1) for r in 1:4) +
            um*sum(abs(b[r]-a[r])*tm^(r-1) for r in 1:4)
        isfinite(bound) ||
            throw(ArgumentError("global solver turbine polynomial bounds overflow"))
        push!(bounds, bound)
    end
    @constraint(m, q == sum(qterms))
    @constraint(m, h == sum(hterms))
    bound = maximum(bounds)
    eta = @variable(m, lower_bound=-bound, upper_bound=bound, base_name=string(name))
    @constraint(m, eta == sum(terms))
    eta
end

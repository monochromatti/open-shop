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

function _global_table!(m, c::TableCurve, x, lo, hi; name = gensym(:table), tightened = true)
    get!(m.ext,:global_tables,Dict{String,Tuple{Float64,Float64}}())[string(name)]=(lo,hi)
    cells = _global_intervals(c.x, lo, hi)
    n = length(cells)
    if tightened && n == 1
        l, r=only(cells)
        yl=table_value(c, l; extrapolation=:linear)
        yr=table_value(c, r; extrapolation=:linear)
        all(isfinite, (yl, yr)) || throw(ArgumentError("global solver table extension bounds overflow"))
        y=@variable(m, lower_bound=min(yl,yr), upper_bound=max(yl,yr), base_name=string(name))
        if l==r
            @constraint(m, x==l)
            @constraint(m, y==yl)
        else
            @constraint(m, l<=x<=r)
            i=_curve_segment(c.x,l+(r-l)/2,:linear)
            slope=(c.y[i+1]-c.y[i])/(c.x[i+1]-c.x[i])
            @constraint(m, y==c.y[i]+slope*(x-c.x[i]))
        end
        return y
    end
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

# A turbine cell is cubic in discharge and affine in head. Its range is
# attained at a head endpoint and a discharge endpoint or stationary point.
function _global_polynomial_range(a, lo, hi)
    f(x)=a[1]+x*(a[2]+x*(a[3]+x*a[4]))
    vals=Float64[f(lo),f(hi)]
    A,B,C=3a[4],2a[3],a[2]
    if A==0
        B!=0 && lo < -C/B < hi && push!(vals,f(-C/B))
    else
        d=B^2-4A*C
        if d>=0
            # The second root follows from the product, avoiding cancellation.
            z=-0.5*(B+copysign(sqrt(d),B))
            roots=z==0 ? (-B/(2A),) : (z/A,C/z)
            for x in roots
                lo<x<hi && push!(vals,f(x))
            end
        end
    end
    extrema(vals)
end

function _global_turbine_cell_range(a, b, tl, tr, ul, ur)
    ranges=[_global_polynomial_range(ntuple(k->a[k]+u*(b[k]-a[k]),4),tl,tr) for u in (ul,ur)]
    lo=minimum(first,ranges)
    hi=maximum(last,ranges)
    slack=1e-10*max(1.0,abs(lo),abs(hi))
    (lo-slack,hi+slack)
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

function _global_shift_polynomial(a, origin, scale)
    (
        a[1]+origin*(a[2]+origin*(a[3]+origin*a[4])),
        scale*(a[2]+origin*(2a[3]+3origin*a[4])),
        scale^2*(a[3]+3origin*a[4]),
        scale^3*a[4],
    )
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
    tightened = true,
    commitment = nothing,
    min_on_flow = 0.0,
    exact_bounds = false,
    range_cuts = false,
)
    bound_ranges=tightened || exact_bounds || range_cuts
    get!(m.ext,:global_turbines,Dict{String,NTuple{4,Float64}}())[string(name)]=(qlo,qhi,hlo,hhi)
    separated=tightened && commitment!==nothing && min_on_flow>0
    qc = separated ? vcat([(0.0,0.0)], qhi>=min_on_flow ? _global_intervals(c.discharge,min_on_flow,qhi) : Tuple{Float64,Float64}[]) : _global_intervals(c.discharge, qlo, qhi)
    hc = _global_intervals(c.heads, hlo, hhi)
    cells = [(qr, hr) for qr in qc for hr in hc]
    get!(m.ext,:global_turbine_cells,Dict())[string(name)]=cells
    get!(m.ext,:global_turbine_normalized,Dict{String,Bool}())[string(name)]=tightened
    n = length(cells)
    z = n == 1 ? [1.0] : @variable(m, [1:n], binary=true, base_name="$(name)_cell")
    n > 1 && @constraint(m, sum(z) == 1)
    if separated
        # Isolated q=0 cells represent off operation. On cells contain only
        # admissible discharge; their selectors are tied to commitment.
        @constraint(m, sum(z[k] for k in eachindex(cells) if cells[k][1][1]>0; init=0.0)==commitment)
    end
    ts = JuMP.VariableRef[]
    us = JuMP.VariableRef[]
    qterms = Any[]
    hterms = Any[]
    terms = Any[]
    bounds = Float64[]
    ranges = Tuple{Float64,Float64}[]
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
        cubic = c.interpolation == :pchip_discharge && first(c.discharge)<=qm<=last(c.discharge)
        a = _global_discharge_coefficients(c,i,j,cubic)
        b = _global_discharge_coefficients(c,i,j+1,cubic)
        qorigin,horigin=c.discharge[i],c.heads[j]
        qscale,hscale=dq,dh
        if tightened
            # Use local [0,1] coordinates even for a tiny clipped cell. Narrow
            # source-coordinate intervals otherwise create nearly dependent LP rows.
            aa=_global_shift_polynomial(a,tl,(qr-ql)/dq)
            bb=_global_shift_polynomial(b,tl,(qr-ql)/dq)
            a=ntuple(k->aa[k]+ul*(bb[k]-aa[k]),4)
            b=ntuple(k->aa[k]+ur*(bb[k]-aa[k]),4)
            qorigin,horigin=ql,hl
            qscale,hscale=qr-ql,hr-hl
            tl,tr=0.0,ql==qr ? 0.0 : 1.0
            ul,ur=0.0,hl==hr ? 0.0 : 1.0
        end
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
        push!(qterms, qorigin*z[k]+qscale*t)
        push!(hterms, horigin*z[k]+hscale*u)
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
        bound_ranges && push!(ranges,_global_turbine_cell_range(a,b,tl,tr,ul,ur))
    end
    @constraint(m, q == sum(qterms))
    @constraint(m, h == sum(hterms))
    bound = maximum(bounds)
    emin,emax=bound_ranges ? (minimum(first,ranges),maximum(last,ranges)) : (-bound,bound)
    eta = @variable(m, lower_bound=emin, upper_bound=emax, base_name=string(name))
    @constraint(m, eta == sum(terms))
    if range_cuts
        # At integer selection exactly one cell contributes. These linear
        # enclosure cuts strengthen its relaxation without changing the graph.
        @constraint(m, eta >= sum(ranges[k][1]*z[k] for k in 1:n))
        @constraint(m, eta <= sum(ranges[k][2]*z[k] for k in 1:n))
    end
    eta
end

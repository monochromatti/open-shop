"""A finite piecewise-linear table. Extrapolation is rejected unless requested.

Knot derivatives use the segment to the right, except at the final knot. No
smoothing, fitting, sorting, deduplication or clipping is performed.
"""
struct TableCurve
    x::Vector{Float64}
    y::Vector{Float64}
    function TableCurve(x, y)
        xx=Float64.(x)
        yy=Float64.(y)
        length(xx)==length(yy)>=2 ||
            throw(ArgumentError("a table needs at least two paired points"))
        all(isfinite, xx) && all(isfinite, yy) ||
            throw(ArgumentError("table points must be finite"))
        all(diff(xx) .> 0) ||
            throw(ArgumentError("table x coordinates must strictly increase"))
        new(xx, yy)
    end
end
TableCurve(; x, y) = TableCurve(x, y)

function _curve_segment(x, z, extrapolation)
    extrapolation in (:error, :linear) ||
        throw(ArgumentError("extrapolation must be :error or :linear"))
    isfinite(z) || throw(DomainError(z, "table coordinate must be finite"))
    if extrapolation==:error && !(first(x)<=z<=last(x))
        throw(DomainError(z, "coordinate outside table domain [$(first(x)), $(last(x))]"))
    end
    clamp(searchsortedlast(x, z), 1, length(x)-1)
end

function table_value(c::TableCurve, z; extrapolation = :error)
    i=_curve_segment(c.x, z, extrapolation)
    c.y[i]+(z-c.x[i])*(c.y[i + 1]-c.y[i])/(c.x[i + 1]-c.x[i])
end
function table_slope(c::TableCurve, z; extrapolation = :error)
    i=_curve_segment(c.x, z, extrapolation)
    (c.y[i + 1]-c.y[i])/(c.x[i + 1]-c.x[i])
end

# Shape-preserving cubic Hermite slopes (weighted harmonic interior slopes).
function _pchip_slopes(x, y)
    n=length(x)
    h=diff(x)
    delta=diff(y) ./ h
    m=zeros(n)
    n==2 && return fill(delta[1], 2)
    for k in 2:(n - 1)
        if delta[k - 1]*delta[k]>0
            w1=2h[k]+h[k - 1]
            w2=h[k]+2h[k - 1]
            m[k]=(w1+w2)/(w1/delta[k - 1]+w2/delta[k])
        end
    end
    function endpoint(h1, h2, d1, d2)
        v=((2h1+h2)*d1-h1*d2)/(h1+h2)
        sign(v)!=sign(d1) && return 0.0
        sign(d1)!=sign(d2) && abs(v)>abs(3d1) && return 3d1
        v
    end
    m[1]=endpoint(h[1], h[2], delta[1], delta[2])
    m[end]=endpoint(h[end], h[end - 1], delta[end], delta[end - 1])
    m
end
function _hermite_value(x, y, m, z, i)
    h=x[i + 1]-x[i]
    t=(z-x[i])/h
    (2t^3-3t^2+1)*y[i]+(t^3-2t^2+t)*h*m[i]+(-2t^3+3t^2)*y[i + 1]+(t^3-t^2)*h*m[i + 1]
end

"""Rectangular turbine efficiency data, with rows for discharge and columns
for head. Efficiency is bilinear; qmin/qmax are linear in head. The discharge
grid must cover the envelope, but may include non-operating points used to
define the numerical extension. Units: m, m³/s, dimensionless efficiency.
Operating head extrapolation requires `head_extrapolation=:linear`; the default
`:error` requires supplied head knots to cover the unit's operating head range.
"""
struct TurbineTable
    heads::Vector{Float64}
    discharge::Vector{Float64}
    efficiency::Matrix{Float64}
    qmin::Vector{Float64}
    qmax::Vector{Float64}
    interpolation::Symbol
    head_extrapolation::Symbol
    discharge_slopes::Matrix{Float64}
    function TurbineTable(
        heads,
        discharge,
        efficiency,
        qmin,
        qmax;
        interpolation = :bilinear,
        head_extrapolation = :error,
    )
        hh=Float64.(heads)
        qq=Float64.(discharge)
        ee=Matrix{Float64}(efficiency)
        lo=Float64.(qmin)
        hi=Float64.(qmax)
        TableCurve(hh, lo)
        TableCurve(hh, hi)
        TableCurve(qq, zeros(length(qq)))
        size(ee)==(length(qq), length(hh)) ||
            throw(ArgumentError("turbine efficiency rows must be discharge, columns head"))
        all(isfinite, ee) && all(0 .< ee .<= 1) ||
            throw(ArgumentError("turbine efficiency must be finite in (0,1]"))
        first(hh)>0 && first(qq)>=0 ||
            throw(ArgumentError("turbine heads must be positive and discharge nonnegative"))
        all(first(qq) .<= lo .<= hi .<= last(qq)) ||
            throw(ArgumentError("turbine discharge envelope exceeds table domain"))
        interpolation in (:bilinear, :pchip_discharge) ||
            throw(ArgumentError("unsupported turbine interpolation"))
        head_extrapolation in (:error, :linear) ||
            throw(ArgumentError("turbine head extrapolation must be :error or :linear"))
        slopes=hcat([_pchip_slopes(qq, ee[:, j]) for j in eachindex(hh)]...)
        new(hh, qq, ee, lo, hi, interpolation, head_extrapolation, slopes)
    end
end
TurbineTable(;
    heads,
    discharge,
    efficiency,
    qmin,
    qmax,
    interpolation = :bilinear,
    head_extrapolation = :error,
) = TurbineTable(
    heads,
    discharge,
    efficiency,
    qmin,
    qmax;
    interpolation,
    head_extrapolation,
)

function turbine_efficiency(c::TurbineTable, q, h; extrapolation = :error)
    i=_curve_segment(c.discharge, q, extrapolation)
    j=_curve_segment(c.heads, h, extrapolation)
    a=(q-c.discharge[i])/(c.discharge[i + 1]-c.discharge[i])
    b=(h-c.heads[j])/(c.heads[j + 1]-c.heads[j])
    if c.interpolation==:pchip_discharge && first(c.discharge)<=q<=last(c.discharge)
        left=_hermite_value(
            c.discharge,
            view(c.efficiency, :, j),
            view(c.discharge_slopes, :, j),
            q,
            i,
        )
        right=_hermite_value(
            c.discharge,
            view(c.efficiency, :, j+1),
            view(c.discharge_slopes, :, j+1),
            q,
            i,
        )
        return (1-b)*left+b*right
    end
    (1-a)*(1-b)*c.efficiency[i, j]+a*(1-b)*c.efficiency[i + 1, j]+(1-a)*b*c.efficiency[
        i,
        j + 1,
    ]+a*b*c.efficiency[i + 1, j + 1]
end
function _head_envelope(c::TurbineTable, values, h, extrapolation)
    j=_curve_segment(c.heads, h, extrapolation)
    values[j]+(h-c.heads[j])*(values[j + 1]-values[j])/(c.heads[j + 1]-c.heads[j])
end
turbine_qmin(c::TurbineTable, h; extrapolation = :error) =
    _head_envelope(c, c.qmin, h, extrapolation)
turbine_qmax(c::TurbineTable, h; extrapolation = :error) =
    _head_envelope(c, c.qmax, h, extrapolation)

"""Exact extrema of the bilinear, linearly extended table on a rectangle.
Used for finite off-state bounds in a mixed-integer model. Every cell is
bilinear, hence its extrema occur at corners; include all crossed knot lines.
"""
function turbine_efficiency_bounds(c::TurbineTable, qlo, qhi, hlo, hhi)
    all(isfinite, (qlo, qhi, hlo, hhi)) && qlo<=qhi && hlo<=hhi ||
        throw(ArgumentError("invalid turbine extension rectangle"))
    qq=unique(vcat(qlo, filter(x->qlo<x<qhi, c.discharge), qhi))
    hh=unique(vcat(hlo, filter(x->hlo<x<hhi, c.heads), hhi))
    extrema(turbine_efficiency(c, q, h; extrapolation = :linear) for q in qq, h in hh)
end

"""Register a JuMP operator using an explicit linear extension for trial points.
The caller must constrain and audit the physical table domain separately.
"""
table_operator(m, c::TableCurve; name = gensym(:table)) =
    JuMP.add_nonlinear_operator(m, 1, x->table_value(c, x; extrapolation = :linear); name)
turbine_operator(m, c::TurbineTable; name = gensym(:turbine)) = JuMP.add_nonlinear_operator(
    m,
    2,
    (q, h)->turbine_efficiency(c, q, h; extrapolation = :linear);
    name,
)

"""Convert shaft MW to electrical MW using a table indexed by electrical MW.

Solves P = shaft * eta(P) exactly within each linear segment. Positive segment
intercepts make P/eta(P) strictly increasing, so the physical root is unique.
"""
function electrical_power(c::TableCurve, shaft; extrapolation = :error)
    validate_electrical_curve(c)
    shaft==0 && return zero(shaft)
    required=c.x ./ c.y
    i=_curve_segment(required, shaft, extrapolation)
    slope=(c.y[i + 1]-c.y[i])/(c.x[i + 1]-c.x[i])
    intercept=c.y[i]-slope*c.x[i]
    denom=1-shaft*slope
    denom>0 || throw(
        DomainError(
            shaft,
            "generator efficiency extension has no positive electrical root",
        ),
    )
    shaft*intercept/denom
end

function validate_electrical_curve(c::TableCurve)
    first(c.x)==0 ||
        throw(ArgumentError("electrical efficiency table must start at zero MW"))
    all(0 .< c.y .<= 1) || throw(ArgumentError("generator efficiency must be in (0,1]"))
    all(c.y[i]-table_slope(c, c.x[i])*c.x[i]>0 for i in 1:(length(c.x) - 1)) ||
        throw(ArgumentError("generator P/efficiency must strictly increase"))
    true
end

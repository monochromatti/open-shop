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

"""Fixed-speed pump data: one discharge and hydraulic efficiency per head.

Heads are metres, discharge is m³/s, and efficiency is a fraction. Values are
piecewise linear in head. A single reference head requires the explicit policy
`head_extrapolation=:constant`; multiple heads support `:error` or `:linear`.
No variable-discharge efficiency surface or operating envelope is constructed.
"""
struct PumpCharacteristic
    heads::Vector{Float64}
    discharge::Vector{Float64}
    efficiency::Vector{Float64}
    head_extrapolation::Symbol
    function PumpCharacteristic(heads, discharge, efficiency; head_extrapolation=:error)
        h,q,e=Float64.(heads),Float64.(discharge),Float64.(efficiency)
        length(h)==length(q)==length(e)>0 && all(isfinite,h) &&
            all(diff(h).>0) && all(>(0),h) && all(x->isfinite(x)&&x>0,q) &&
            all(x->isfinite(x)&&0<x<=1,e) ||
            throw(ArgumentError("pump characteristic needs ordered positive heads, positive flow and efficiency in (0,1]"))
        (length(h)==1 ? head_extrapolation==:constant : head_extrapolation in (:error,:linear)) ||
            throw(ArgumentError("one pump head requires :constant; multiple heads require :error or :linear"))
        new(h,q,e,head_extrapolation)
    end
end
PumpCharacteristic(;heads,discharge,efficiency,head_extrapolation=:error) =
    PumpCharacteristic(heads,discharge,efficiency;head_extrapolation)

function _pump_value(c::PumpCharacteristic,values,h;extrapolation=:error)
    extrapolation in (:error,:linear) || throw(ArgumentError("extrapolation must be :error or :linear"))
    isfinite(h) || throw(DomainError(h,"pump head must be finite"))
    length(c.heads)==1 && return only(values)+zero(h)
    i=_curve_segment(c.heads,h,extrapolation)
    values[i]+(h-c.heads[i])*(values[i+1]-values[i])/(c.heads[i+1]-c.heads[i])
end
pump_discharge(c::PumpCharacteristic,h;kwargs...) = _pump_value(c,c.discharge,h;kwargs...)
pump_efficiency(c::PumpCharacteristic,h;kwargs...) = _pump_value(c,c.efficiency,h;kwargs...)
function _pump_slope(c::PumpCharacteristic,h)
    length(c.heads)==1 && return 0.0
    i=_curve_segment(c.heads,h,:linear)
    (c.discharge[i+1]-c.discharge[i])/(c.heads[i+1]-c.heads[i])
end

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

"""Turbine efficiency curves at ordered reference heads.

Use `curves=[TableCurve(discharge, efficiency), ...]` for one original curve per
head, including unequal discharge grids. The rectangular constructor also accepts
discharge rows and head columns. Each curve keeps its own linear or PCHIP
interpolant; the resulting efficiencies are interpolated linearly in head.

Units are m, m³/s and efficiency fractions. `qmin` and `qmax` give the operating
envelope independently of the efficiency knots. `discharge_extrapolation=:linear`
permits an envelope outside those knots, using the end secants for efficiency.
The default `:error` requires efficiency data to cover the operating discharges.

For two or more heads, `head_extrapolation=:error` requires coverage of the
unit's operating head range; `:linear` permits end-secant continuation. A single
reference head requires `:constant`, declaring efficiency and the discharge
envelope independent of head. Power still uses the actual hydraulic head.

Raw lookup APIs reject out-of-table coordinates unless `extrapolation=:linear`
is requested. A declared single-head constant response has no lookup head limit.
No calibration knots or operating limits are fitted, added or clipped.
"""
struct TurbineTable
    heads::Vector{Float64}
    curves::Vector{TableCurve}
    qmin::Vector{Float64}
    qmax::Vector{Float64}
    interpolation::Symbol
    head_extrapolation::Symbol
    discharge_extrapolation::Symbol
    _discharge_knots::Vector{Float64}
    _discharge_slopes::Vector{Vector{Float64}}
    function TurbineTable(
        heads,
        curves,
        qmin,
        qmax;
        interpolation = :bilinear,
        head_extrapolation = :error,
        discharge_extrapolation = :error,
    )
        hh=Float64.(heads)
        all(c->c isa TableCurve,curves) ||
            throw(ArgumentError("turbine curves must be TableCurve records in head order"))
        cc=[TableCurve(c.x,c.y) for c in curves]
        lo=Float64.(qmin)
        hi=Float64.(qmax)
        !isempty(hh) && length(cc)==length(lo)==length(hi)==length(hh) &&
            all(isfinite,hh) && all(i->hh[i]<hh[i+1], 1:(length(hh)-1)) &&
            all(isfinite,lo) && all(isfinite,hi) ||
            throw(ArgumentError("turbine heads and envelope need paired finite increasing data"))
        all(c->all(x->0<x<=1,c.y),cc) ||
            throw(ArgumentError("turbine efficiency must be finite in (0,1]"))
        first(hh)>0 && all(c->first(c.x)>=0,cc) ||
            throw(ArgumentError("turbine heads must be positive and discharge nonnegative"))
        all(i->0<=lo[i]<=hi[i], eachindex(lo,hi)) ||
            throw(ArgumentError("turbine discharge envelope must be nonnegative and ordered"))
        interpolation in (:bilinear, :pchip_discharge) ||
            throw(ArgumentError("unsupported turbine interpolation"))
        if length(hh)==1
            head_extrapolation==:constant ||
                throw(ArgumentError("a single turbine head requires head_extrapolation=:constant"))
        else
            head_extrapolation in (:error, :linear) ||
                throw(ArgumentError("multiple turbine heads require :error or :linear extrapolation"))
        end
        discharge_extrapolation in (:error, :linear) ||
            throw(ArgumentError("turbine discharge extrapolation must be :error or :linear"))
        discharge_extrapolation==:linear ||
            all(i->first(cc[i].x)<=lo[i]<=hi[i]<=last(cc[i].x), eachindex(lo,hi)) ||
            throw(ArgumentError("turbine discharge envelope exceeds efficiency data; declare :linear extrapolation"))
        knots=sort!(unique(vcat([c.x for c in cc]...)))
        slopes=[_pchip_slopes(c.x,c.y) for c in cc]
        all(m->all(isfinite,m),slopes) || throw(ArgumentError("turbine slopes must be finite"))
        new(hh,cc,lo,hi,interpolation,head_extrapolation,discharge_extrapolation,knots,slopes)
    end
end

function TurbineTable(heads,discharge,efficiency,qmin,qmax;kwargs...)
    values=efficiency isa AbstractVector ? reshape(Float64.(efficiency),:,1) : Matrix{Float64}(efficiency)
    size(values)==(length(discharge),length(heads)) ||
        throw(ArgumentError("turbine efficiency rows must be discharge, columns head"))
    TurbineTable(heads,[TableCurve(discharge,view(values,:,j)) for j in eachindex(heads)],qmin,qmax;kwargs...)
end

TurbineTable(;
    heads,
    discharge=nothing,
    efficiency=nothing,
    curves=nothing,
    qmin,
    qmax,
    interpolation = :bilinear,
    head_extrapolation = :error,
    discharge_extrapolation = :error,
) = if curves===nothing
    discharge!==nothing && efficiency!==nothing ||
        throw(ArgumentError("turbine data need curves or rectangular discharge and efficiency"))
    TurbineTable(heads,discharge,efficiency,qmin,qmax;interpolation,head_extrapolation,discharge_extrapolation)
else
    discharge===nothing && efficiency===nothing ||
        throw(ArgumentError("supply turbine curves or rectangular data, not both"))
    TurbineTable(heads,curves,qmin,qmax;interpolation,head_extrapolation,discharge_extrapolation)
end

function _turbine_head_weights(c::TurbineTable, h, extrapolation)
    isfinite(h) || throw(DomainError(h, "turbine head must be finite"))
    extrapolation in (:error, :linear) ||
        throw(ArgumentError("extrapolation must be :error or :linear"))
    length(c.heads)==1 && return (1, 1, zero(h))
    j=_curve_segment(c.heads,h,extrapolation)
    (j, j+1, (h-c.heads[j])/(c.heads[j+1]-c.heads[j]))
end

function turbine_efficiency(c::TurbineTable, q, h; extrapolation = :error)
    j,k,b=_turbine_head_weights(c,h,extrapolation)
    j==k && return _turbine_column(c,j,q,extrapolation)
    # A strict exact-head lookup needs only that supplied curve. Numerical
    # extensions still evaluate both columns, preserving head derivatives.
    if extrapolation==:error
        b==0 && return _turbine_column(c,j,q,extrapolation)
        b==1 && return _turbine_column(c,k,q,extrapolation)
    end
    left=_turbine_column(c,j,q,extrapolation)
    right=_turbine_column(c,k,q,extrapolation)
    left+b*(right-left)
end

function _turbine_column(c,j,q,extrapolation)
    curve=c.curves[j]
    i=_curve_segment(curve.x,q,extrapolation)
    c.interpolation==:pchip_discharge && first(curve.x)<=q<=last(curve.x) ?
        _hermite_value(curve.x,curve.y,c._discharge_slopes[j],q,i) :
        curve.y[i]+(q-curve.x[i])*(curve.y[i+1]-curve.y[i])/(curve.x[i+1]-curve.x[i])
end

function _turbine_head_columns(c,hlo,hhi)
    columns=Int[]
    for h in unique(vcat(hlo,filter(z->hlo<z<hhi,c.heads),hhi))
        j,k,b=_turbine_head_weights(c,h,:linear)
        b==1 || push!(columns,j)
        b==0 || push!(columns,k)
    end
    unique(columns)
end
function _head_envelope(c::TurbineTable, values, h, extrapolation)
    j,k,b=_turbine_head_weights(c,h,extrapolation)
    values[j]+b*(values[k]-values[j])
end
turbine_qmin(c::TurbineTable, h; extrapolation = :error) =
    _head_envelope(c, c.qmin, h, extrapolation)
turbine_qmax(c::TurbineTable, h; extrapolation = :error) =
    _head_envelope(c, c.qmax, h, extrapolation)

"""Exact extrema of the linearly extended interpolant on a rectangle.
Within each discharge cell it is cubic or linear, and it is linear in head.
Include the original polynomial's stationary points and all crossed knot lines.
"""
function turbine_efficiency_bounds(c::TurbineTable, qlo, qhi, hlo, hhi)
    all(isfinite, (qlo, qhi, hlo, hhi)) && qlo<=qhi && hlo<=hhi ||
        throw(ArgumentError("invalid turbine extension rectangle"))
    qq=unique(vcat(qlo, filter(x->qlo<x<qhi, c._discharge_knots), qhi))
    hh=unique(vcat(hlo, filter(x->hlo<x<hhi, c.heads), hhi))
    length(qq)==1 && return extrema(turbine_efficiency(c,only(qq),h;extrapolation=:linear) for h in hh)
    ranges=[_global_polynomial_range(_global_tensor_polynomial(c,a,b,h),0.0,1.0)
        for (a,b) in zip(qq[1:end-1],qq[2:end]),h in hh]
    (minimum(first,ranges),maximum(last,ranges))
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

turbine_flow_operator(m, c::TurbineTable, values; name=gensym(:turbine_flow)) =
    JuMP.add_nonlinear_operator(m,1,h->_head_envelope(c,values,h,:linear);name)

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

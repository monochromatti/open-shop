# The discrete gradient of reservoir potential is the volume-averaged level.
# Integrate piecewise-linear data segment by segment, rather than subtracting
# nearly equal potentials. The two partial derivatives are weighted mean slopes.
function _mean_level_data(c::TableCurve, a, b)
    if a>b
        value, db, da=_mean_level_data(c,b,a)
        return value,da,db
    end
    i=_curve_segment(c.x,a,:linear)
    j=_curve_segment(c.x,b,:linear)
    if i==j || a==b
        slope=(c.y[i+1]-c.y[i])/(c.x[i+1]-c.x[i])
        return c.y[i]+slope*((a-c.x[i])+(b-c.x[i]))/2,slope/2,slope/2
    end
    width=b-a
    da=zero(width)
    db=zero(width)
    left=a
    for k in i:j
        right=k==j ? b : c.x[k+1]
        theta_left=(left-a)/width
        theta_right=(right-a)/width
        weight=(right-left)/width
        centre=(theta_left+theta_right)/2
        slope=(c.y[k+1]-c.y[k])/(c.x[k+1]-c.x[k])
        da+=slope*weight*(1-centre)
        db+=slope*weight*centre
        left=right
    end
    table_value(c,a;extrapolation=:linear)+width*da,da,db
end

"""
    mean_head(reservoir, before, after)

Reservoir head in metres averaged over the storage interval from `before` to
`after` (Mm³). This is the discrete gradient of gravitational potential:
`(after-before)*mean_head = integral(level, before, after)`. At equal volumes it
equals `head(reservoir, before)`. Reversing the endpoints leaves it unchanged.

Polynomial curves are integrated algebraically; tables are integrated exactly
across their linear segments, with end-secant extension for trial points. No
subtraction of nearby potential values or small-volume cutoff is used. This
energy identity does not establish the time accuracy or uniqueness of a step.
"""
function mean_head(r::Reservoir,a,b)
    if r.level_curve===nothing
        mid=(a+b)/2
        return r.z0+r.slope*mid+r.curvature*(mid^2+(b-a)^2/12)
    end
    first(_mean_level_data(r.level_curve,a,b))
end

function _mean_head_partials(r,a,b)
    if r.level_curve===nothing
        return r.slope/2+r.curvature*(2a+b)/3,
            r.slope/2+r.curvature*(a+2b)/3
    end
    _,da,db=_mean_level_data(r.level_curve,a,b)
    da,db
end

function _mean_head_operator(m,r;name)
    c=r.level_curve
    value=(a,b)->first(_mean_level_data(c,a,b))
    gradient=(g,a,b)->begin
        _,g[1],g[2]=_mean_level_data(c,a,b)
        nothing
    end
    JuMP.add_nonlinear_operator(m,2,value,gradient;name)
end

# Integral of (level-reference_head) dV. Its small values are suitable for the
# scaled global energy equations; the arbitrary elevation datum cancels there.
function _level_potential(c::TableCurve,reference_volume,reference_head,v)
    a,b=minmax(reference_volume,v)
    i=_curve_segment(c.x,a,:linear)
    j=_curve_segment(c.x,b,:linear)
    result=zero(v)
    left=a
    for k in i:j
        right=k==j ? b : c.x[k+1]
        slope=(c.y[k+1]-c.y[k])/(c.x[k+1]-c.x[k])
        result+=(right-left)*(c.y[k]-reference_head+
            slope*((left-c.x[k])+(right-c.x[k]))/2)
        left=right
    end
    v<reference_volume ? -result : result
end

# Polynomial coefficients and exact cubic ranges for turbine interpolation.
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


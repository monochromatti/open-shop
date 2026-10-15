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
    curve=c.curves[j]
    a = curve.y[i]
    b = curve.y[i + 1]
    cubic || return (a, b-a, 0.0, 0.0)
    d = curve.x[i + 1]-curve.x[i]
    u = d*c._discharge_slopes[j][i]
    v = d*c._discharge_slopes[j][i + 1]
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

# Coefficients in the local [0,1] discharge coordinate at a supplied head.
# Each extension segment keeps its original secant, even for PCHIP tables.
function _global_column_polynomial(c,ql,qr,j)
    grid=c.curves[j].x
    # The common partition never crosses an original knot. Select by its left
    # edge: a rounded midpoint can fall into the next cell on a very narrow grid.
    i=_curve_segment(grid,ql,:linear)
    cubic=c.interpolation==:pchip_discharge && first(grid)<=ql<=qr<=last(grid)
    coefficients=_global_discharge_coefficients(c,i,j,cubic)
    dq=grid[i+1]-grid[i]
    _global_shift_polynomial(coefficients,(ql-grid[i])/dq,(qr-ql)/dq)
end

function _global_tensor_polynomial(c,ql,qr,h)
    j,k,u=_turbine_head_weights(c,h,:linear)
    a=_global_column_polynomial(c,ql,qr,j)
    (j==k || u==0) && return a
    b=_global_column_polynomial(c,ql,qr,k)
    u==1 && return b
    ntuple(i->a[i]+u*(b[i]-a[i]),4)
end


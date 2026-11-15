# Interior power extrema occur at a data knot or a root of a quadratic.
# Coefficients remain differentiable so operator gradients/Hessians include
# movement of an interior stationary head. Degenerate inflections need no root.
function cubic_extrema(c,a,b)
    value(h)=((c[4]*h+c[3])*h+c[2])*h+c[1]
    T=typeof(c[1]+c[2]+c[3]+c[4]+a+b)
    heads=T[a,b]
    A,B,C=3c[4],2c[3],c[2]
    function primal(x)
        while x isa ForwardDiff.Dual
            x=ForwardDiff.value(x)
        end
        x
    end
    if iszero(primal(A))
        iszero(primal(B)) || push!(heads,-C/B)
    else
        scale=max(abs(primal(A)),abs(primal(B)),abs(primal(C)))
        A,B,C=A/scale,B/scale,C/scale
        discriminant=B^2-4A*C
        if primal(discriminant)>0
            root=-0.5*(B+(primal(B)<0 ? -sqrt(discriminant) : sqrt(discriminant)))
            push!(heads,root/A)
            iszero(primal(root)) || push!(heads,C/root)
        end
    end
    minimum_value,maximum_value=value(a),value(a)
    at_minimum,at_maximum=a,a
    for h in heads
        primal(a)<=primal(h)<=primal(b) || continue
        p=value(h)
        if p<minimum_value
            minimum_value,at_minimum=p,h
        end
        if p>maximum_value
            maximum_value,at_maximum=p,h
        end
    end
    (;minimum=minimum_value,maximum=maximum_value,at_minimum,at_maximum)
end

"""Exact power/domain envelope across a net-head span for a held unit control.

Power inequalities are written without dividing by efficiency. Generation uses
the strictly increasing electrical relation P/e(P); pumping uses shaft work
minus P*efficiency. Within a supplied head segment these are cubic or quadratic
polynomials. Their extrema are evaluated algebraically, without time sampling.
The characteristic determines a fixed-speed pump's flow; otherwise q is held.
"""
function unit_head_envelope(g,q,head_a,head_b; pmin=g.pmin,pmax=g.pmax,qmin=g.qmin,qmax=g.qmax)
    lo,hi=minmax(head_a,head_b)
    data=g.characteristic!==nothing ? g.characteristic.heads :
        g.turbine_table!==nothing ? g.turbine_table.heads : Float64[]
    edges=vcat(lo,[h for h in data if lo<h<hi],hi)
    upper=-Inf;lower=-Inf;domain=-Inf;flow=-Inf
    at_upper=lo;at_lower=lo
    for k in 1:length(edges)-1
        a,b=edges[k],edges[k+1]
        middle=(a+b)/2
        q0,q1=q,0.0
        if g.characteristic!==nothing
            c=g.characteristic
            j=length(c.heads)==1 ? 1 : _curve_segment(c.heads,middle,:linear)
            q1=length(c.heads)==1 ? 0.0 : (c.discharge[j+1]-c.discharge[j])/(c.heads[j+1]-c.heads[j])
            q0=c.discharge[j]-q1*c.heads[j]
            e1=length(c.heads)==1 ? 0.0 : (c.efficiency[j+1]-c.efficiency[j])/(c.heads[j+1]-c.heads[j])
            e0=c.efficiency[j]-e1*c.heads[j]
            e2=0.0
        elseif g.turbine_table!==nothing
            table=g.turbine_table
            j,l,_=_turbine_head_weights(table,middle,:linear)
            left=_turbine_column(table,j,q,:linear)
            e1=j==l ? zero(left) : (_turbine_column(table,l,q,:linear)-left)/(table.heads[l]-table.heads[j])
            e0=left-e1*table.heads[j]
            e2=0.0
            for h in (a,b)
                flow=max(flow,turbine_qmin(table,h;extrapolation=:linear)-q,
                    q-turbine_qmax(table,h;extrapolation=:linear))
            end
        else
            e0=g.efficiency-g.qcurvature*((q-g.qbest)/g.qbest)^2-g.hcurvature
            e1=2g.hcurvature/g.hbest
            e2=-g.hcurvature/g.hbest^2
        end
        eta=cubic_extrema((e0,e1,e2,0.0),a,b)
        domain=max(domain,g.min_efficiency-eta.minimum,eta.maximum-1.0)
        flow=max(flow,qmin-(q0+q1*a),qmin-(q0+q1*b),q0+q1*a-qmax,q0+q1*b-qmax)
        if g.mode==:generate
            work=(0.0,0.00981*q*e0,0.00981*q*e1,0.00981*q*e2)
            extrema=cubic_extrema(work,a,b)
            up=extrema.maximum-pmax/generator_efficiency(g,pmax)
            down=pmin/generator_efficiency(g,pmin)-extrema.minimum
            up_head,down_head=extrema.at_maximum,extrema.at_minimum
        else
            up_poly=(-pmax*g.motor_efficiency*e0,
                0.00981*q0-pmax*g.motor_efficiency*e1,
                0.00981*q1-pmax*g.motor_efficiency*e2,0.0)
            down_poly=(pmin*g.motor_efficiency*e0,
                pmin*g.motor_efficiency*e1-0.00981*q0,
                pmin*g.motor_efficiency*e2-0.00981*q1,0.0)
            up_range=cubic_extrema(up_poly,a,b)
            down_range=cubic_extrema(down_poly,a,b)
            up,down=up_range.maximum,down_range.maximum
            up_head,down_head=up_range.at_maximum,down_range.at_maximum
        end
        if up>upper
            upper,at_upper=up,up_head
        end
        if down>lower
            lower,at_lower=down,down_head
        end
    end
    residual=max(upper/DISPATCH_POWER_SCALE,lower/DISPATCH_POWER_SCALE,
        domain,flow/DISPATCH_FLOW_SCALE,(g.hmin-lo)/DISPATCH_HEAD_SCALE,(hi-g.hmax)/DISPATCH_HEAD_SCALE)
    (;residual,upper,lower,domain,flow,at_upper,at_lower)
end

# Interior power extrema occur at a data knot or a root of a quadratic.
# Coefficients remain differentiable so operator gradients/Hessians include
# movement of an interior stationary head.
function cubic_extrema(c,a,b)
    value(h)=((c[4]*h+c[3])*h+c[2])*h+c[1]
    heads=(a,b,quadratic_roots(3c[4],2c[3],c[2])...)
    minimum_value,maximum_value=value(a),value(a)
    at_minimum,at_maximum=a,a
    for h in heads
        a<=h<=b || continue
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
        middle=a
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

"""Establish increasing power in head over the full declared operating box.

For generation the derivative has the sign of eta+h*eta_h; for a variable-flow
pump it has the sign of eta-h*eta_h. Each tabulated head segment is affine in
head and cubic at most in discharge, so its minimum is computed algebraically.
An uncertain sign, or a fixed-speed characteristic, retains the general span
constraints. This removes redundant inequalities without narrowing the model.
"""
function power_increases_with_head(g)
    g.hmin==g.hmax && return true
    g.characteristic!==nothing && return false
    if g.turbine_table===nothing
        lower=Inf;scale=1.0
        for q in (g.qmin,clamp(g.qbest,g.qmin,g.qmax),g.qmax)
            e0=g.efficiency-g.qcurvature*((q-g.qbest)/g.qbest)^2-g.hcurvature
            e1=2g.hcurvature/g.hbest
            e2=-g.hcurvature/g.hbest^2
            polynomial=g.mode==:generate ? (e0,2e1,3e2,0.0) : (e0,0.0,-e2,0.0)
            lower=min(lower,cubic_extrema(polynomial,g.hmin,g.hmax).minimum)
            scale=max(scale,sum(abs,polynomial))
        end
        return lower>128eps(Float64)*scale
    end
    table=g.turbine_table
    qs=_global_tensor_nodes(table._discharge_knots,g.qmin,g.qmax)
    hs=_global_tensor_nodes(table.heads,g.hmin,g.hmax)
    flows=length(qs)==1 ? [(first(qs),first(qs))] : collect(zip(qs[1:end-1],qs[2:end]))
    for (a,b) in zip(hs[1:end-1],hs[2:end]),(lo,hi) in flows
        j,k,_=_turbine_head_weights(table,a,:linear)
        left=_global_column_polynomial(table,lo,hi,j)
        slope=j==k ? ntuple(_->0.0,4) :
            (_global_column_polynomial(table,lo,hi,k).-left)./(table.heads[k]-table.heads[j])
        for h in (a,b)
            coefficient=g.mode==:generate ? 2h-table.heads[j] : -table.heads[j]
            polynomial=left.+coefficient.*slope
            lower,_=_global_polynomial_range(polynomial,0.0,1.0)
            lower>128eps(Float64)*max(1.0,sum(abs,polynomial)) || return false
        end
    end
    true
end

function unit_efficiency_bounds(g)
    if g.characteristic!==nothing
        hs=vcat(g.hmin,filter(h->g.hmin<h<g.hmax,g.characteristic.heads),g.hmax)
        return extrema(pump_efficiency(g.characteristic,h;extrapolation=:linear) for h in hs)
    elseif g.turbine_table!==nothing
        return turbine_efficiency_bounds(g.turbine_table,g.qmin,g.qmax,g.hmin,g.hmax)
    end
    polynomial_efficiency_bounds(g,g.qmin,g.qmax,g.hmin,g.hmax)
end

function discharge_envelope_monotone(g)
    monotone(values)=all(>=(0),diff(values)) || all(<=(0),diff(values))
    g.characteristic!==nothing && return monotone(g.characteristic.discharge)
    g.turbine_table===nothing && return true
    monotone(g.turbine_table.qmin) && monotone(g.turbine_table.qmax)
end

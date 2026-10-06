# Supporting inequalities strengthen the relaxation of the original power law.
# The nonlinear equality and every feasible integer schedule are retained.
# Bernstein coefficients enclose the polynomial over each discharge/head box;
# no samples or fitted curves are used to justify the inequalities.

function _power_box_upper(table, qlo, qhi, hlo, hhi, electrical_max, slope, head_slope=0.)
    a=_global_tensor_polynomial(table, qlo, qhi, hlo)
    z=_global_tensor_polynomial(table, qlo, qhi, hhi)
    d=z.-a
    dq=qhi-qlo;dh=hhi-hlo
    C=zeros(5,3)
    for k in 1:4
        for (row,flow) in ((k,qlo),(k+1,dq))
            C[row,1]+=flow*hlo*a[k]
            C[row,2]+=flow*(hlo*d[k]+dh*a[k])
            C[row,3]+=flow*dh*d[k]
        end
    end
    C .*= 0.00981*electrical_max
    C[1,1]-=slope*qlo+head_slope*hlo;C[2,1]-=slope*dq
    C[1,2]-=head_slope*dh
    upper=-Inf
    for i in 0:4,j in 0:2
        b=sum(C[k+1,l+1]*binomial(i,k)/binomial(4,k)*
            binomial(j,l)/binomial(2,l) for k in 0:i,l in 0:j)
        upper=max(upper,b)
    end
    upper+1e-9*max(1.0,sum(abs,C))
end

function _power_support(table,qlo,qhi,hlo,hhi,electrical_max,slope,head_slope=0.)
    qs=_global_tensor_nodes(table.discharge,qlo,qhi)
    hs=_global_tensor_nodes(table.heads,hlo,hhi)
    # Subdivision tightens the coefficient enclosure without changing the model.
    segments=length(qs)==1 ? [(qlo,qhi)] :
        [(a+(b-a)*k/4,a+(b-a)*(k+1)/4) for (a,b) in zip(qs[1:end-1],qs[2:end]) for k in 0:3]
    heads=length(hs)==1 ? [(hlo,hhi)] : collect(zip(hs[1:end-1],hs[2:end]))
    maximum(_power_box_upper(table,a,b,h0,h1,electrical_max,slope,head_slope)
        for (a,b) in segments,(h0,h1) in heads)
end

function _add_power_bounds!(b,c)
    count=0
    for (i,g) in enumerate(c.system.generators),t in eachindex(c.prices)
        g.turbine_table===nothing && continue
        hd=b.shared_heads[(plantof(c.system,g).name,t)]
        hlo=max(g.hmin,lower_bound(hd));hhi=min(g.hmax,upper_bound(hd))
        qlo=opinterval(c,g.name,:qmin,t,g.qmin)
        qhi=min(g.qmax,opinterval(c,g.name,:qmax,t,g.qmax))
        (hlo>hhi || qlo>qhi) && continue
        electrical=g.generator_efficiency_curve===nothing ? [1.] :
            [table_value(g.generator_efficiency_curve,p;extrapolation=:linear) for
                p in _global_tensor_nodes(g.generator_efficiency_curve.x,0.,g.pmax)]
        minimum(electrical)>=0 || continue
        electrical_max=maximum(electrical)
        table=g.turbine_table
        scale=0.00981*hhi*electrical_max
        # Slopes need no optimality property: the certified intercept supports
        # every physically feasible on-state. Multiplying it by u also covers off.
        for fraction in (0.25,0.5,0.75,1.0)
            slope=scale*fraction
            intercept=_power_support(table,qlo,qhi,hlo,hhi,electrical_max,slope)
            @constraint(b.m,(b.P[i,t]-slope*b.GQ[i,t]-intercept*b.u[i,t])/40<=0)
            count+=1
            # The head origin includes off states; a positive head slope
            # makes the off-state inequality valid throughout that interval.
            for qref in (qlo,(qlo+qhi)/2)
                hslope=.00981*electrical_max*.95*qref
                origin=lower_bound(hd)
                cut=_power_support(table,qlo,qhi,hlo,hhi,electrical_max,slope,hslope)
                @constraint(b.m,(b.P[i,t]-slope*b.GQ[i,t]-hslope*(hd-origin)-
                    (cut+hslope*origin)*b.u[i,t])/40<=0)
                count+=1
            end
        end
    end
    count
end

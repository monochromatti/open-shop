# Closed-form interval consequences of opening*head_loss = resistance*q*abs(q).
# Apply once to existing conservative head ranges; no trajectory is assumed.
function tighten_tunnel_domains!(b,c)
    byname=Dict(name(v)=>v for v in all_variables(b.m))
    tightened=0
    function narrow(v,lo,hi)
        is_fixed(v) && return
        lo=max(lower_bound(v),lo);hi=min(upper_bound(v),hi)
        lo<=hi || error("inconsistent tunnel hydraulic bounds")
        tightened+=Int(lo>lower_bound(v))+Int(hi<upper_bound(v))
        set_lower_bound(v,lo);set_upper_bound(v,hi)
    end
    signedroot(x)=copysign(sqrt(abs(x)),x)
    for (i,e) in enumerate(c.system.tunnels),t in eachindex(c.prices)
        opening=OpenSHOP.opinterval(c,e.name,:opening,t,e.opening)
        opening==0 && continue
        src=b.node_bounds[(e.source,t)];dst=b.node_bounds[(e.target,t)]
        lo=signedroot(opening*(src[1]-dst[2])/e.resistance)
        hi=signedroot(opening*(src[2]-dst[1])/e.resistance)
        slack=1e-8*max(1.,abs(lo),abs(hi))
        lo-=slack;hi+=slack
        narrow(byname["q[$i,$t]"],lo/50,hi/50)
        narrow(byname["tunnel_positive_$(i)_$(t)"],max(lo,0.),max(hi,0.))
        narrow(byname["tunnel_negative_$(i)_$(t)"],max(-hi,0.),max(-lo,0.))
        direction=byname["tunnel_direction_$(i)_$(t)"]
        if lo>0 || hi<0
            fix(direction,lo>0 ? 1. : 0.;force=true)
        end
    end
    b.m.ext[:hydraulic_bound_changes]=tightened
    tightened
end

# Hydraulic junctions have no storage. Propagating their signed flow balances
# proves one-way flow in radial intakes while preserving reversible crosslinks.
function tighten_network_domains!(b,c)
    tighten_tunnel_domains!(b,c)
    s=c.system;E=length(s.tunnels);T=length(c.prices)
    byname=Dict(name(v)=>v for v in all_variables(b.m))
    bounds(v)=is_fixed(v) ? (fix_value(v),fix_value(v)) : (lower_bound(v),upper_bound(v))
    lo=zeros(E,T);hi=similar(lo)
    for i in 1:E,t in 1:T
        a,z=bounds(byname["q[$i,$t]"]);lo[i,t]=50a;hi[i,t]=50z
    end
    plants=[OpenSHOP.plantof(s,g) for g in s.generators]
    balances=Any[]
    for node in s.junctions
        any(r->r.source==node.name || r.target==node.name,s.rivers) && continue
        edges=[(i,Int(e.target==node.name)-Int(e.source==node.name)) for
            (i,e) in enumerate(s.tunnels) if e.source==node.name || e.target==node.name]
        units=[(i,Int(p.target==node.name)-Int(p.source==node.name)) for
            (i,p) in enumerate(plants) if p.source==node.name || p.target==node.name]
        push!(balances,(edges,units))
    end
    for iteration in 1:(E+1)
        changed=false
        for (edges,units) in balances,t in 1:T
            gl=0.;gh=0.;glmag=0.;ghmag=0.;glzero=true;ghzero=true
            for (i,a) in units
                l,h=bounds(byname["gq[$i,$t]"]);l*=50;h*=50
                lowterm=min(a*l,a*h);highterm=max(a*l,a*h)
                gl+=lowterm;gh+=highterm;glmag+=abs(lowterm);ghmag+=abs(highterm)
                glzero&=lowterm==0.;ghzero&=highterm==0.
            end
            for (i,a) in edges
                l=gl;h=gh;lmag=glmag;hmag=ghmag;lzero=glzero;hzero=ghzero
                for (j,aj) in edges
                    j==i && continue
                    lowterm=min(aj*lo[j,t],aj*hi[j,t]);highterm=max(aj*lo[j,t],aj*hi[j,t])
                    l+=lowterm;h+=highterm;lmag+=abs(lowterm);hmag+=abs(highterm)
                    lzero&=lowterm==0.;hzero&=highterm==0.
                end
                low,high=minmax(-h/a,-l/a)
                slack=1e-8*max(1.,abs(low),abs(high),lmag,hmag)
                # Preserve exact zero only when every contributing endpoint is zero;
                # mixed signed cancellation must retain outward rounding.
                (a>0 ? hzero : lzero) || (low-=slack)
                (a>0 ? lzero : hzero) || (high+=slack)
                low=max(lo[i,t],low);high=min(hi[i,t],high)
                low<=high || error("inconsistent junction flow bounds")
                changed|=low>lo[i,t]+1e-7 || high<hi[i,t]-1e-7
                lo[i,t]=low;hi[i,t]=high
            end
        end
        changed || break
    end
    for i in 1:E,t in 1:T
        q=byname["q[$i,$t]"]
        is_fixed(q) && continue
        set_lower_bound(q,lo[i,t]/50);set_upper_bound(q,hi[i,t]/50)
        qp=byname["tunnel_positive_$(i)_$(t)"];qm=byname["tunnel_negative_$(i)_$(t)"]
        set_lower_bound(qp,max(0.,lo[i,t]));set_upper_bound(qp,max(0.,hi[i,t]))
        set_lower_bound(qm,max(0.,-hi[i,t]));set_upper_bound(qm,max(0.,-lo[i,t]))
        if lo[i,t]>=0 || hi[i,t]<0
            fix(byname["tunnel_direction_$(i)_$(t)"],lo[i,t]>=0 ? 1. : 0.;force=true)
        end
    end
end

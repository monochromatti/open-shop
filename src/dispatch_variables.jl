"""Optimization coordinates and their physical expressions.

Lowercase arrays are solver coordinates; uppercase arrays are physical storage
(Mm³), head (m), discharge (m³/s) and power (MW). The same coordinates serve
local dispatch and the global algebraic model. Adjustable discharge and river
gates are held on operating intervals; fixed-speed pump flows remain states.
Starts and extracted schedules are expressed in physical units; these same
scales convert them to and from solver coordinates.
"""
struct DispatchVariables
    v::Matrix{VariableRef}
    h::Matrix{VariableRef}
    q::Matrix{VariableRef}
    gq::Matrix{VariableRef}
    p::Matrix{VariableRef}
    rq::Matrix{VariableRef}
    a::Matrix{VariableRef}
    V::Matrix{AffExpr}
    H::Matrix{AffExpr}
    Q::Matrix{AffExpr}
    GQ::Matrix{AffExpr}
    P::Matrix{AffExpr}
    RQ::Matrix{AffExpr}
    storage_scale::Vector{Float64}
end

const DISPATCH_HEAD_SCALE=250.0
const DISPATCH_FLOW_SCALE=50.0
const DISPATCH_POWER_SCALE=40.0
const DISPATCH_RELEASE_SCALE=100.0

function DispatchVariables(m,s,mesh;prefix="",shared_controls=nothing,control_index=mesh.interval)
    R,G,D=length(s.reservoirs),length(s.generators),length(s.rivers)
    N=R+length(s.junctions)+length(s.boundaries)
    E,T=length(s.tunnels),length(mesh.dt)
    v=@variable(m,[1:R,1:T+1],base_name=prefix*"v")
    h=@variable(m,[1:N,1:T],base_name=prefix*"h")
    q=@variable(m,[1:E,1:T],base_name=prefix*"q")
    gq=Matrix{VariableRef}(undef,G,T)
    for t in 1:T,(i,g) in enumerate(s.generators)
        gq[i,t]=if g.characteristic!==nothing
            @variable(m,lower_bound=0.0,base_name=prefix*"gq[$i,$t]")
        elseif shared_controls!==nothing
            shared_controls.gq[i,control_index[t]]
        elseif t==first(mesh.cells[mesh.interval[t]])
            @variable(m,lower_bound=0.0,base_name=prefix*"gq[$i,$t]")
        else
            gq[i,t-1]
        end
    end
    p=@variable(m,[1:G,1:T],lower_bound=0.0,base_name=prefix*"p")
    rq=@variable(m,[1:D,1:T],lower_bound=0.0,base_name=prefix*"rq")
    a=if shared_controls===nothing
        gate=@variable(m,[1:D,1:length(mesh.cells)],lower_bound=0.0,upper_bound=1.0,
            base_name=prefix*"a")
        gate[:,mesh.interval]
    else
        shared_controls.a[:,control_index]
    end
    scale=[r.vmax for r in s.reservoirs]
    V=[scale[i]*v[i,t] for i in 1:R,t in 1:T+1]
    DispatchVariables(v,h,q,gq,p,rq,a,V,DISPATCH_HEAD_SCALE.*h,
        DISPATCH_FLOW_SCALE.*q,DISPATCH_FLOW_SCALE.*gq,
        DISPATCH_POWER_SCALE.*p,DISPATCH_RELEASE_SCALE.*rq,scale)
end

"""Assign the solver coordinates of a physical schedule, without altering it.

`put(variable,value)` also permits an audited global lift to collect assignments
before installing them in the solver. Auxiliary curve variables are separate.
"""
function set_dispatch_start!(variables::DispatchVariables,x,put=set_start_value)
    for t in axes(x.storage,2),i in axes(x.storage,1)
        put(variables.v[i,t],x.storage[i,t]/variables.storage_scale[i])
    end
    for k in eachindex(x.head)
        put(variables.h[k],x.head[k]/DISPATCH_HEAD_SCALE)
    end
    for k in eachindex(x.tunnel_q)
        put(variables.q[k],x.tunnel_q[k]/DISPATCH_FLOW_SCALE)
    end
    for k in eachindex(x.generator_q)
        put(variables.gq[k],x.generator_q[k]/DISPATCH_FLOW_SCALE)
        put(variables.p[k],x.power[k]/DISPATCH_POWER_SCALE)
    end
    for k in eachindex(x.river_release)
        put(variables.rq[k],x.river_release[k]/DISPATCH_RELEASE_SCALE)
        put(variables.a[k],x.gate[k])
    end
    nothing
end

"""Initialize dispatch from its physical reference or the standard on-state guess."""
function initialize_dispatch!(variables::DispatchVariables,c,u;warm=nothing,seed=1)
    s,T=c.system,length(c.prices)
    rng=MersenneTwister(seed)
    for (i,r) in enumerate(s.reservoirs),t in 1:T
        set_start_value(variables.h[i,t],head(r,r.v0)/DISPATCH_HEAD_SCALE)
    end
    for (i,g) in enumerate(s.generators),t in 1:T
        guess=clamp((seed==1 ? 0.6 : 0.45+0.3rand(rng))*g.qmax,g.qmin,g.qmax)*u[i,t]
        set_start_value(variables.gq[i,t],guess/DISPATCH_FLOW_SCALE)
        set_start_value(variables.p[i,t],power(g,guess,g.hbest)/DISPATCH_POWER_SCALE)
    end
    for (i,r) in enumerate(s.rivers),t in 1:T
        set_start_value(variables.rq[i,t],mean(r.history_release)/DISPATCH_RELEASE_SCALE)
        r.law==:orifice && set_start_value(variables.a[i,t],0.5)
    end
    if warm!==nothing
        set_dispatch_start!(variables,warm)
        for (i,g) in enumerate(s.generators),t in 1:T
            set_start_value(variables.gq[i,t],
                clamp(warm.generator_q[i,t],g.qmin*u[i,t],g.qmax*u[i,t])/DISPATCH_FLOW_SCALE)
            set_start_value(variables.p[i,t],
                clamp(warm.power[i,t],g.pmin*u[i,t],g.pmax*u[i,t])/DISPATCH_POWER_SCALE)
        end
    end
    nothing
end

"""Extract a physical schedule from either optimization representation."""
function _dispatch_values(b)
    c,variables=b.case,b.variables
    release=value.(variables.RQ)
    opening=_tunnel_gate_solution(b.gates)[:,b.mesh.interval]
    raw=Schedule(copy(c.grid),round.(Int,value.(b.u)),value.(variables.V),
        value.(variables.H),value.(variables.Q),Matrix{Float64}(opening),
        value.(variables.GQ),value.(variables.P),value.(variables.a),release,
        Matrix{Float64}(value.(b.routing.arrivals)),Vector{Float64}(value.(b.routing.terminal)),
        value.(b.shortfall_release),value(b.obj),value(b.gates.cost),value(b.penalty_cost),
        Dict{String,Any}("case"=>c.name,
            "fractional_tunnel_opening"=>value.(b.gates.opening)[:,b.mesh.interval]);
        reserves=_reserve_values(c,b.reserves;operating_case=b.operating_case))
    with_soft_assessment(b.operating_case,raw;included_cost=value(b.soft.cost))
end

"""Independent midpoint hydraulic step. Rates are m³/s; storage is Mm³.

`river_arrival` contains known history/prior-cohort arrivals. `current_transfer`
contains the fractions of each current release cohort arriving in this same
interval, one column per reference curve. River junction mixing and current
arrivals are evaluated inside Newton's residual, coupled to storage and heads.
Operational bounds are reported, never imposed by clipping physical states.
"""
function forward_step(
    s::HydroSystem,
    Vprev,
    dt,
    generator_q,
    gate,
    river_arrival;
    current_transfer = zeros(length(s.rivers), 2),
    current_network = nothing,
    initial_guess = nothing,
    tolerance = 1e-9,
    max_iterations = 100,
)
    nr, nj, ne, nd=length(s.reservoirs),
    length(s.junctions),
    length(s.tunnels),
    length(s.rivers)
    ng=length(s.generators)
    ix=nodeindex(s)
    nn=length(nodes(s))
    length(Vprev)==nr &&
    length(generator_q)==ng &&
    length(gate)==nd &&
    length(river_arrival)==nd || throw(DimensionMismatch("forward inputs"))
    size(current_transfer, 1)==nd &&
    size(current_transfer, 2)>=maximum((length(r.curves) for r in s.rivers); init = 1) ||
        throw(DimensionMismatch("current transfer fractions"))
    isfinite(dt) &&
    dt>0 &&
    all(isfinite, vcat(Vprev, generator_q, gate, river_arrival, vec(current_transfer))) ||
        throw(ArgumentError("finite state and positive duration required"))
    all(>=(0), generator_q) &&
    all(x->0<=x<=1, gate) &&
    all(>=(0), river_arrival) &&
    all(x->0<=x<=1, current_transfer) ||
        throw(ArgumentError("invalid controls or arrival fractions"))
    active=findall(r->r.law!=:junction, s.rivers)
    na=length(active)
    order=_forward_river_order(s)
    base=nr+nj
    nq=base+ne
    nv=nq+na
    alpha=0.0036dt
    B=zeros(base, ne)
    fixed=zeros(base)
    for (i, r) in enumerate(s.reservoirs)
        fixed[i]=r.inflow
    end
    add!(v, node, q) = haskey(ix, node) && ix[node]<=base ? (v[ix[node]]+=q) : nothing
    for (e, t) in enumerate(s.tunnels)
        add!(view(B, :, e), t.source, -1.0)
        add!(view(B, :, e), t.target, 1.0)
    end
    for (i, g) in enumerate(s.generators)
        p=plantof(s, g)
        add!(fixed, p.source, -generator_q[i])
        add!(fixed, p.target, generator_q[i])
    end
    for d in active
        haskey(ix, s.rivers[d].source) && ix[s.rivers[d].source]<=nr ||
            throw(ArgumentError("operational river law requires reservoir source"))
    end
    function heads(x)
        h=zeros(nn)
        for (i, r) in enumerate(s.reservoirs)
            h[i]=head(r, (Vprev[i]+x[i])/2)
        end
        h[(nr + 1):base].=x[(nr + 1):base]
        for (i, b) in enumerate(s.boundaries)
            h[base + i]=b.head
        end
        h
    end
    law(r, h, a) = river_law_value(r, h, a)
    function river_rates(x)
        release=zeros(nd)
        arrival=copy(river_arrival)
        release[active].=x[(nq + 1):end]
        current_network!==nothing && (arrival .+= current_network*release)
        for d in order
            r=s.rivers[d]
            if r.law==:junction
                release[d]=sum(
                    (arrival[k] for (k, e) in enumerate(s.rivers) if e.target==r.source);
                    init = 0.0,
                )
            end
            # Guards define finite trial residuals only. Accepted physical
            # releases remain untouched, and routing-domain violations fail replay.
            if current_network===nothing
                fraction=sum(
                    w*current_transfer[d, l] for (l, w) in _forward_weights(r, release[d])
                )
                arrival[d]+=max(0.0, release[d])*fraction
            end
        end
        release, arrival
    end
    function equations(x)
        h=heads(x)
        net=B*x[(base + 1):nq]+fixed
        f=zeros(nv)
        release, arrival=river_rates(x)
        for (d, r) in enumerate(s.rivers)
            add!(net, r.source, -release[d])
            add!(net, r.target, arrival[d])
        end
        for i in 1:nr
            f[i]=(x[i]-Vprev[i])/alpha-net[i]
        end
        f[(nr + 1):base].=-net[(nr + 1):base]
        for (k, t) in enumerate(s.tunnels)
            q=x[base + k]
            f[base + k]=t.opening==0 ? q :
                        t.opening*(h[ix[t.source]]-h[ix[t.target]])-t.resistance*q*abs(q)
        end
        for (k, d) in enumerate(active)
            r=s.rivers[d]
            f[nq + k]=x[nq + k]-law(r, h[ix[r.source]], gate[d])
        end
        f
    end
    x=zeros(nv)
    x[1:nr].=Vprev
    known=vcat(
        [head(r, Vprev[i]) for (i, r) in enumerate(s.reservoirs)],
        [b.head for b in s.boundaries],
    )
    x[(nr + 1):base].=isempty(known) ? 0.0 : sum(known)/length(known)
    h=heads(x)
    for (k, t) in enumerate(s.tunnels)
        dh=h[ix[t.source]]-h[ix[t.target]]
        x[base + k]=sign(dh)*sqrt(t.opening*abs(dh)/t.resistance)
    end
    for (k, d) in enumerate(active)
        x[nq + k]=law(s.rivers[d], h[ix[s.rivers[d].source]], gate[d])
    end
    if initial_guess!==nothing
        if initial_guess isa AbstractDict
            x[1:nr].=initial_guess["Vnew"]
            x[(nr + 1):base].=initial_guess["H"][(nr + 1):base]
            x[(base + 1):nq].=initial_guess["tunnel_q"]
            x[(nq + 1):end].=initial_guess["river_release"][active]
        else
            length(initial_guess)==nv || throw(DimensionMismatch("initial guess"))
            x.=initial_guess
        end
    end
    iterations=0
    for k in 1:max_iterations
        f=equations(x)
        maximum(abs, f; init = 0.0)<=tolerance && break
        J=zeros(nv, nv)
        for j in 1:nv
            epsj=1e-5*max(abs(x[j]), 1.0)
            xp=copy(x)
            xm=copy(x)
            xp[j]+=epsj
            xm[j]-=epsj
            J[:, j].=(equations(xp)-equations(xm))/(2epsj)
        end
        step=try
            -(J\f)
        catch
            -pinv(J)*f
        end
        all(isfinite, step) || (step=-pinv(J)*f)
        merit=sum(abs2, f)
        accepted=false
        scale=1.0
        for _ in 1:45
            candidate=x+scale*step
            fc=equations(candidate)
            if all(isfinite, fc) && sum(abs2, fc)<merit
                x=candidate
                accepted=true
                break
            end
            scale/=2
        end
        iterations=k
        accepted || break
    end
    h=heads(x)
    # Exact law evaluation removes Newton roundoff from a shut outlet. Its
    # storage consequence is included in the recomputed convergence residual.
    for (k, d) in enumerate(active)
        x[nq + k]=law(s.rivers[d], h[ix[s.rivers[d].source]], gate[d])
    end
    release, arrival=river_rates(x)
    gh=[
        begin
            plant=plantof(s, g)
            h[ix[plant.source]]-outlet_head(plant, h[ix[plant.target]]) - tailwater(
                plant,
                sum(
                    generator_q[j] for (j, z) in enumerate(s.generators) if z.plant==g.plant
                ),
            )
        end for g in s.generators
    ]
    pw=[power(g, generator_q[i], gh[i]) for (i, g) in enumerate(s.generators)]
    exchange=zeros(length(s.boundaries))
    function boundary!(node, q)
        j=get(ix, node, 0)-base
        j>0 && (exchange[j]+=q)
    end
    for (i, t) in enumerate(s.tunnels)
        boundary!(t.source, -x[base + i])
        boundary!(t.target, x[base + i])
    end
    for (i, g) in enumerate(s.generators)
        p=plantof(s, g)
        boundary!(p.source, -generator_q[i])
        boundary!(p.target, generator_q[i])
    end
    for (i, r) in enumerate(s.rivers)
        boundary!(r.source, -release[i])
        boundary!(r.target, arrival[i])
    end
    residual=maximum(abs, equations(x); init = 0.0)
    domain=all(isfinite, vcat(x, h, pw, release, arrival)) &&
           all(_forward_domain(r, release[d]) for (d, r) in enumerate(s.rivers))
    violations=Dict(
        "storage"=>maximum(
            (max(r.vmin-x[i], x[i]-r.vmax, 0.0) for (i, r) in enumerate(s.reservoirs));
            init = 0.0,
        ),
        "tunnel_flow"=>maximum(
            (max(abs(x[base + i])-r.capacity, 0.0) for (i, r) in enumerate(s.tunnels));
            init = 0.0,
        ),
        "river_flow"=>maximum(
            (
                max(-release[d], release[d]-r.capacity, 0.0) for
                (d, r) in enumerate(s.rivers)
            );
            init = 0.0,
        ),
    )
    Dict(
        "Vnew"=>x[1:nr],
        "H"=>h,
        "tunnel_q"=>x[(base + 1):nq],
        "river_release"=>release,
        "river_arrival"=>arrival,
        "generator_head"=>gh,
        "power"=>pw,
        "residual"=>residual,
        "converged"=>isfinite(residual) && residual<=tolerance && domain,
        "iterations"=>iterations,
        "boundary_net_inflow"=>exchange,
        "routing_domain_valid"=>domain,
        "bound_violations"=>violations,
        "law_domain_valid"=>all(
            s.rivers[d].law==:controlled ||
                s.rivers[d].allow_dry ||
                s.rivers[d].discharge_curve!==nothing ||
                h[ix[s.rivers[d].source]]>=s.rivers[d].crest for d in active
        ),
    )
end

# Topological reach evaluation also accepts boundary targets used by standalone
# integration tests. River graph/input validation belongs to the public model.
function _forward_river_order(s)
    remaining=Set(eachindex(s.rivers))
    order=Int[]
    while !isempty(remaining)
        ready=sort!([
            d for d in remaining if !any(
                k in remaining && s.rivers[k].target==s.rivers[d].source for
                k in eachindex(s.rivers)
            )
        ])
        isempty(ready) && throw(ArgumentError("river graph must be acyclic"))
        for d in ready
            r=s.rivers[d]
            r.law==:junction &&
                count(e->e.source==r.source, s.rivers)!=1 &&
                throw(ArgumentError("river junction needs exactly one outgoing reach"))
            push!(order, d)
            delete!(remaining, d)
        end
    end
    order
end
_forward_atom(r) = hasproperty(r, :deterministic_delay) && r.deterministic_delay!==nothing
function _forward_weights(r, q)
    _forward_atom(r) && return [(1, 1.0)]
    length(r.curves)==1 && return [(1, 1.0)]
    refs=[c.reference_flow for c in r.curves]
    RiverRouting.blend(r.curves, clamp(q, first(refs), last(refs)))
end
function _forward_domain(r, q)
    isfinite(q) && q>=-1e-10 && q<=r.capacity+1e-9 || return false
    _forward_atom(r) ||
        length(r.curves)==1 ||
        first(r.curves).reference_flow-1e-10<=q<=last(r.curves).reference_flow+1e-9
end
function _forward_coefficients(grid, r)
    _forward_atom(r) || return RiverRouting.coefficients(grid, grid, r.curves)
    a=RiverRouting.fixed_coefficients(grid, grid, r.deterministic_delay)
    reshape(a, size(a, 1), size(a, 2), 1)
end
function _forward_history(grid, r)
    _forward_atom(r) || return (
        RiverRouting.route(grid, r.history_grid, r.history_release, r.curves),
        RiverRouting.remaining_volume(
            first(grid),
            r.history_grid,
            r.history_release,
            r.curves,
        ),
    )
    RiverRouting.check_releases(r.history_grid, r.history_release)
    a=RiverRouting.fixed_coefficients(grid, r.history_grid, r.deterministic_delay)*(
        0.0036 .* diff(r.history_grid) .* r.history_release
    )
    time=first(grid)
    w=sum(
        0.0036*r.history_release[k]*(
            max(0.0, min(time, r.history_grid[k + 1])-r.history_grid[k])-max(
                0.0,
                min(time-r.deterministic_delay, r.history_grid[k + 1])-r.history_grid[k],
            )
        ) for k in eachindex(r.history_release)
    )
    a, w
end

"""Independent chronological replay with coupled same-step river arrivals.

Source outlet laws, midpoint storage and zero/short-delay transfers are solved
simultaneously during each forward step. Junction mixing conserves volume but
rebins the arriving water uniformly on the replay grid; refine the grid to
assess this timing approximation. Each cohort contributes exactly once.
"""
function simulate(c::ScheduleCase, generator_q, gate; grid = c.grid, transport = nothing)
    s=c.system
    ng, nd=length(s.generators), length(s.rivers)
    T=length(grid)-1
    size(generator_q)==(ng, length(c.grid)-1) && size(gate)==(nd, length(c.grid)-1) ||
        throw(DimensionMismatch("schedule dimensions"))
    RiverRouting.check_grid(grid)
    first(grid)==first(c.grid) && last(grid)==last(c.grid) ||
        throw(ArgumentError("replay grid must span original horizon"))
    all(t->any(isapprox(t, x; atol = 1e-12) for x in grid), c.grid) ||
        throw(ArgumentError("replay grid must contain original control change times"))
    indices=[
        min(searchsortedlast(c.grid, (grid[t]+grid[t + 1])/2), length(c.grid)-1) for
        t in 1:T
    ]
    gq=generator_q[:, indices]
    ga=gate[:, indices]
    nr=length(s.reservoirs)
    V=zeros(nr, T+1)
    V[:, 1].=[r.v0 for r in s.reservoirs]
    H=zeros(length(nodes(s)), T)
    Q=zeros(length(s.tunnels), T)
    P=zeros(ng, T)
    R=zeros(nd, T)
    A=zeros(nd, T)
    W=zeros(nd, T+1)
    out=zeros(T)
    steps=Any[]
    exact=all(r.deterministic_delay!==nothing for r in s.rivers)
    nc_case=_river_replace(c; grid = Float64.(grid), prices = zeros(T))
    ndat=exact ? _transport_data(nc_case, transport) : nothing
    !exact &&
        transport!==nothing &&
        throw(
            ArgumentError(
                "compiled deterministic transport cannot serve distributed replay",
            ),
        )
    B=exact ? Array{Float64,3}[] : [_forward_coefficients(grid, r) for r in s.rivers]
    # forward_step still checks the reference-curve column shape, even when
    # the exact network matrix supersedes those reference fractions.
    nc=exact ? maximum((length(r.curves) for r in s.rivers); init = 1) :
       maximum((size(b, 3) for b in B); init = 1)
    transfer=zeros(nd, nc)
    network=exact ? zeros(nd, nd) : nothing
    sizehint!(steps, T)
    if exact
        W[:, 1].=ndat.initial_transit
    else
        for (d, r) in enumerate(s.rivers)
            ah, wh=_forward_history(grid, r)
            A[d, :].=ah
            W[d, 1]=wh
        end
    end
    for t in 1:T
        dt=grid[t + 1]-grid[t]
        if !exact
            for d in 1:nd, l in axes(B[d], 3)
                transfer[d, l]=B[d][t, t, l]
            end
        end
        if exact
            fill!(network, 0.0)
            for d in 1:nd
                A[d, t]=ndat.arrival_history[d, t]
                for (i, k, v) in ndat.arrival_terms[d, t]
                    if k<t
                        A[d, t]+=v*R[i, k]
                    elseif k==t
                        network[d, i]+=v/(0.0036dt)
                    else
                        error("noncausal river transport")
                    end
                end
            end
        end
        ss=operational_system(c, grid[t], grid[t + 1])
        z=forward_step(
            ss,
            V[:, t],
            dt,
            gq[:, t],
            ga[:, t],
            A[:, t]/(0.0036dt);
            current_transfer = transfer,
            current_network = network,
            initial_guess = isempty(steps) ? nothing : last(steps),
        )
        push!(steps, z)
        V[:, t + 1].=z["Vnew"]
        H[:, t].=z["H"]
        Q[:, t].=z["tunnel_q"]
        P[:, t].=z["power"]
        R[:, t].=z["river_release"]
        out[t]=sum(z["boundary_net_inflow"])
        A[:, t].=0.0036dt .* z["river_arrival"]
        for (d, r) in enumerate(s.rivers)
            if !exact
                for (l, w) in _forward_weights(r, R[d, t])
                    A[d, (t + 1):end].+=0.0036dt*R[d, t]*w .* B[d][(t + 1):end, t, l]
                end
            end
            W[d, t + 1]=W[d, t]+0.0036dt*R[d, t]-A[d, t]
        end
    end
    Dict(
        "grid"=>collect(grid),
        "V"=>V,
        "H"=>H,
        "tunnel_q"=>Q,
        "generator_q"=>gq,
        "gate"=>ga,
        "power"=>P,
        "river_release"=>R,
        "arrival_volume"=>A,
        "transit"=>W,
        "boundary_outflow"=>out,
        "converged"=>all(z["converged"] && z["law_domain_valid"] for z in steps),
        "residual"=>maximum((z["residual"] for z in steps); init = 0.0),
        "steps"=>steps,
    )
end

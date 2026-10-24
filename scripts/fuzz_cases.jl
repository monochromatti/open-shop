module FuzzCases

using OpenSHOP, Random

const FAMILIES = ("cascade", "pressure", "rivers", "small_storage", "channels", "pumping")
const SIZES = ("small", "medium", "large")

uniform(rng, a, b) = a + (b-a)*rand(rng)
loguniform(rng, a, b) = exp(uniform(rng, log(a), log(b)))
updated(x; kwargs...) = typeof(x)(; merge(
    NamedTuple{fieldnames(typeof(x))}(Tuple(getfield(x, f) for f in fieldnames(typeof(x)))),
    (; kwargs...))...)
series(object, attribute, times, values) = OperationalSeries(;
    object, attribute, times=collect(Float64, times), values=collect(Float64, values))

function reservoir(rng, name, elevation, volume; inflow=0.0)
    lo, hi = 0.2volume, 1.8volume
    rise = uniform(rng, 1.0, 8.0)
    slope = rise/(hi-lo)
    curvature = rand(rng, (-0.1, 0.0, 0.1))*slope/hi
    z0 = elevation-slope*volume-curvature*volume^2
    curve = if rand(rng, Bool)
        knots = [lo, 0.65volume, volume, 1.35volume, hi]
        TableCurve(knots, [z0+slope*v+curvature*v^2 for v in knots])
    else
        nothing
    end
    Reservoir(; name, z0, slope, curvature, v0=volume, vmin=lo, vmax=hi,
        level_curve=curve, inflow, water_value=rand(rng, (0.0, 1000.0, 15000.0, 30000.0)))
end

function unit(rng, name, plant, flow, head, duration; mode=:generate, table=true)
    qmin = rand(rng, (0.0, 0.1, 0.25))*flow
    hmin, hmax = 0.5head, 1.5head+10.0
    efficiency = uniform(rng, 0.88, 0.94)
    turbine = if table && rand(rng, Bool)
        heads = [hmin, head, hmax]
        # Each measured curve keeps its own knots, rather than resampling PCHIP.
        curves = [TableCurve([0.0, (0.35+0.03j)*flow, 0.7flow, flow],
            [0.80, efficiency-0.02, efficiency, efficiency-0.025]) for j in 1:3]
        TurbineTable(; heads, curves, qmin=fill(qmin, 3), qmax=fill(flow, 3),
            interpolation=rand(rng, (:bilinear, :pchip_discharge)))
    else
        nothing
    end
    Generator(; name, plant, mode, qmin, qmax=flow,
        pmin=rand(rng, (0.0, 0.00981*hmin*qmin*0.75)),
        pmax=0.00981*flow*hmax/0.65, efficiency, qbest=0.7flow, hbest=head,
        hmin, hmax, qcurvature=0.05, hcurvature=0.01, turbine_table=turbine,
        motor_efficiency=mode==:pump ? 0.96 : 1.0,
        initial_on=0, initial_age=24.0, minup=0.4duration, mindown=0.4duration,
        startup=uniform(rng, 0.0, 200.0), shutdown=uniform(rng, 0.0, 30.0))
end

function controls(rng, units, T)
    u = [rand(rng, (0, 1, 1, 1)) for _ in units, _ in 1:T]
    q = [u[i,t]*rand(rng,(g.qmin==0 ? 0.0 : g.qmin+0.05g.qmax,
            uniform(rng,0.35,0.65)*g.qmax,0.95g.qmax))
        for (i,g) in enumerate(units), t in 1:T]
    u, q
end

function reach(rng, name, source, target; capacity=150.0, law=:junction,
    allocation=1.0, inflow=0.0, distributed=false, coefficient=0.0, crest=0.0)
    delay = rand(rng, (0.0, 0.125, 0.5, 1.25, 3.0))
    curves = distributed ? [RiverRouting.DelayCurve(q,
        [delay+d, delay+d+0.25, delay+d+1.0], [0.65, 0.35])
        for (q,d) in ((0.0,0.75), (capacity/2,0.25), (capacity,0.0))] :
        RiverRouting.DelayCurve[]
    River(; name, source, target, law, allocation, capacity, inflow, coefficient, crest,
        deterministic_delay=distributed ? nothing : delay, curves,
        history_grid=[-8.0,-4.0,-1.0,0.0],
        history_release=[uniform(rng, 0.0, 0.1capacity) for _ in 1:3],
        water_value=rand(rng, (0.0, 1000.0)))
end

function cascade(rng, n, grid)
    T = length(grid)-1
    volume = max(10.0, 0.4last(grid))
    lakes = [reservoir(rng, Symbol("Lake",i), 80.0*(n-i+1),
        volume*uniform(rng, 1.0, 3.0); inflow=uniform(rng, 0.0, 4.0)) for i in 1:n]
    plants = Plant[]
    units = Generator[]
    for i in 1:n
        name = Symbol("Station",i)
        push!(plants, Plant(; name, source=lakes[i].name,
            target=i==n ? :Sea : lakes[i+1].name, pmax=120.0, ramp=nothing,
            tailwater_curve=rand(rng, Bool) ? TableCurve([0.0,60.0],[0.0,1.5]) : nothing))
        for j in 1:rand(rng, 1:2)
            push!(units, unit(rng, Symbol("Unit",i,"_",j), name,
                uniform(rng, 8.0, 25.0), 80.0, minimum(diff(grid))))
        end
    end
    u, q = controls(rng, units, T)
    system = HydroSystem(; reservoirs=lakes, plants, generators=units,
        boundaries=[Boundary(name=:Sea,head=0.0)])
    (; system, u, q, gate=zeros(0,T), opening=zeros(0,T), operations=OperationalSeries[])
end

function pressure(rng, n, grid)
    T = length(grid)-1
    lakes = [reservoir(rng, Symbol("Supply",i), uniform(rng, 299.0, 301.0),
        uniform(rng, 15.0, 40.0)) for i in 1:n]
    push!(lakes, reservoir(rng, :Receiver, 150.0, 30.0))
    selectable = rand(rng, Bool)
    tunnels = [Tunnel(name=Symbol("Intake",i), source=lakes[i].name,
        target=Symbol("J",i), resistance=uniform(rng, 0.02, 0.08),
        capacity=80.0, opening=1.0, initial_opening=1.0,
        gate_positions=selectable ? [0.0,0.5,1.0] : Float64[],
        gate_adjustment_cost=uniform(rng, 0.0, 20.0),
        gate_loss_factor_curve=rand(rng, Bool) ?
            TableCurve([0.0,0.5,1.0],[0.3,0.06,0.0]) : nothing) for i in 1:n]
    # Sources and stations attach at distinct nodes of a growing pressure ring.
    for i in 1:(n==2 ? 1 : n)
        push!(tunnels,Tunnel(name=Symbol("PressureCross",i),source=Symbol("J",i),
            target=Symbol("J",mod1(i+1,n)),resistance=uniform(rng,0.02,0.08),capacity=80.0))
    end
    push!(tunnels, Tunnel(name=:LakeCross,source=lakes[1].name,target=lakes[2].name,
        resistance=0.12,capacity=20.0))
    plants = [Plant(name=Symbol("Station",i), source=Symbol("J",i),
        target=:Receiver,pmax=100.0,ramp=nothing) for i in 1:n]
    units = [unit(rng,Symbol("Unit",i),p.name,uniform(rng,5.0,10.0),140.0,
        minimum(diff(grid))) for (i,p) in enumerate(plants)]
    u,q = controls(rng,units,T)
    opening = ones(length(tunnels),T)
    for t in 1:T
        opening[1:n,t] .= rand(rng,(0.0,0.5,1.0),n)
        opening[randperm(rng,n)[1:cld(n,2)],t] .= 1.0
    end
    operations = selectable ? OperationalSeries[] :
        [series(e.name,:opening,grid[1:end-1],opening[i,:]) for (i,e) in enumerate(tunnels)]
    system = HydroSystem(; reservoirs=lakes,
        junctions=[Junction(name=Symbol("J",i),hmin=200.0,hmax=310.0) for i in 1:n],
        tunnels,plants,generators=units,
        tunnel_interlocks=[TunnelInterlock(name=:SupplyAvailability,
            tunnels=[e.name for e in tunnels[1:n]],min_open=1,max_open=n)])
    (; system,u,q,gate=zeros(0,T),opening,operations)
end

function rivers(rng,n,grid)
    T = length(grid)-1
    distributed = rand(rng,Bool)
    lakes = [reservoir(rng,Symbol("Tributary",i),uniform(rng,180.0,240.0),
        uniform(rng,15.0,40.0);inflow=uniform(rng,0.0,5.0)) for i in 1:n]
    push!(lakes,reservoir(rng,:Receiver,80.0,40.0))
    plants = [Plant(name=Symbol("Station",i),source=lakes[i].name,target=:Receiver,
        discharge_river=Symbol("Tail",i),pmax=100.0,ramp=nothing) for i in 1:n]
    units = [unit(rng,Symbol("Unit",i),p.name,uniform(rng,5.0,15.0),
        head(lakes[i],lakes[i].v0)-80.0,minimum(diff(grid))) for (i,p) in enumerate(plants)]
    reaches = [reach(rng,Symbol("Tail",i),:auto,:Fork;distributed,capacity=30.0,
        inflow=uniform(rng,0.0,1.0)) for i in 1:n]
    fraction = uniform(rng,0.1,0.9)
    append!(reaches,[
        reach(rng,:Left,:Fork,:Merge;allocation=fraction,capacity=20.0n,distributed,inflow=0.5),
        reach(rng,:Right,:Fork,:Merge;allocation=1-fraction,capacity=20.0n,distributed),
        reach(rng,:Downstream,:Merge,n>=4 ? :LowerFork : :Receiver;capacity=300.0,distributed)])
    if n>=4
        fraction = uniform(rng,0.1,0.9)
        append!(reaches,[
            reach(rng,:LowerLeft,:LowerFork,:LowerMerge;allocation=fraction,capacity=300.0,distributed),
            reach(rng,:LowerRight,:LowerFork,:LowerMerge;allocation=1-fraction,capacity=300.0,distributed),
            reach(rng,:FinalReach,:LowerMerge,:Receiver;capacity=400.0,distributed)])
    end
    push!(reaches,reach(rng,:Outlet,:Receiver,:Sea;law=:controlled,capacity=10.0,distributed))
    u,q = controls(rng,units,T)
    gate = zeros(length(reaches),T)
    gate[end,:] .= rand(rng,(0.0,0.25,0.5,0.75,1.0),T)
    # A steady prehistory obeys the same split/merge conservation laws. Its
    # station flows also agree with the supplied initial commitment ages.
    outfalls = Dict(p.discharge_river=>q[i,1] for (i,p) in enumerate(plants))
    incoming = Dict{Symbol,Float64}()
    for (i,r) in enumerate(reaches)
        flow = r.law==:junction ? r.allocation*get(incoming,r.source,0.0)+
            get(outfalls,r.name,0.0)+r.inflow : r.capacity*gate[i,1]
        reaches[i] = updated(r;history_release=fill(flow,length(r.history_release)))
        incoming[r.target] = get(incoming,r.target,0.0)+flow
    end
    system = HydroSystem(;reservoirs=lakes,plants,generators=units,rivers=reaches,
        boundaries=[Boundary(name=:Sea,head=0.0)])
    (;system,u,q,gate,opening=zeros(0,T),operations=OperationalSeries[])
end

function small_storage(rng,n,grid)
    T = length(grid)-1
    volume = loguniform(rng,1e-5,0.02) # 10 m³ to 20,000 m³ balancing pond
    lakes = [reservoir(rng,:Upper,180.0,20.0),reservoir(rng,:Pond,80.0,volume)]
    plants = [Plant(name=:UpperStation,source=:Upper,target=:Pond,pmax=50.0,ramp=nothing),
        Plant(name=:LowerStation,source=:Pond,target=:Sea,pmax=50.0,ramp=nothing)]
    flow = uniform(rng,5.0,25.0)
    units = [unit(rng,Symbol("Unit",i),p.name,flow,i==1 ? 100.0 : 80.0,
        minimum(diff(grid));table=false) for (i,p) in enumerate(plants)]
    u = ones(Int,2,T)
    q = repeat(reshape([uniform(rng,0.35,0.65)*flow for _ in 1:T],1,T),2,1)
    # Turnover can greatly exceed storage; the two stations balance each other.
    # A small independent inflow moves the level without draining the pond.
    inflow = [uniform(rng,0.0,0.05volume/(0.0036last(grid))) for _ in 1:T]
    system = HydroSystem(;reservoirs=lakes,plants,generators=units,
        boundaries=[Boundary(name=:Sea,head=0.0)])
    (;system,u,q,gate=zeros(0,T),opening=zeros(0,T),
        operations=[series(:Pond,:inflow,grid[1:end-1],inflow)])
end

function channels(rng,n,grid)
    T = length(grid)-1
    lakes = [reservoir(rng,Symbol("Lake",i),100.0+uniform(rng,-0.4,0.4),
        uniform(rng,10.0,30.0);inflow=uniform(rng,0.0,2.0)) for i in 1:n]
    reaches = River[]
    for i in 1:n-1
        gain = uniform(rng,0.3,3.0)
        curve = LevelFlowCurve([90.0,110.0],[
            TableCurve([-10.0,0.0,10.0],[-10gain,0.0,12gain]),
            TableCurve([-10.0,0.0,10.0],[-12gain,0.0,15gain])])
        push!(reaches,River(name=Symbol("Channel",i),source=lakes[i].name,
            target=lakes[i+1].name,law=:delta,delta_level_ref_up_flow_curve=curve,
            capacity=50.0,deterministic_delay=0.0,curves=RiverRouting.DelayCurve[],water_value=0.0))
    end
    law = rand(rng,(:orifice,:weir,:tabulated))
    crest = 100.0+rand(rng,(-0.5,0.0,0.5))
    outlet = reach(rng,:Outlet,lakes[end].name,:Sea;capacity=20.0,
        law=law==:tabulated ? :orifice : law,
        coefficient=uniform(rng,0.5,2.0),crest)
    outlet = updated(outlet;allow_dry=true,deterministic_delay=0.0,
        history_release=zeros(length(outlet.history_release)))
    if law==:tabulated
        outlet=updated(outlet;discharge_curve=TableCurve([95.0,100.0,105.0],[0.0,1.0,5.0]))
    end
    push!(reaches,outlet)
    plant = Plant(name=:Station,source=lakes[1].name,target=:Sea,pmax=30.0,ramp=nothing)
    units = [unit(rng,:Unit,:Station,10.0,100.0,minimum(diff(grid)))]
    u,q = controls(rng,units,T)
    gate = zeros(length(reaches),T)
    gate[end,:] .= law==:weir ? ones(T) : rand(rng,(0.0,0.25,0.5,1.0),T)
    system = HydroSystem(;reservoirs=lakes,rivers=reaches,plants=[plant],generators=units,
        boundaries=[Boundary(name=:Sea,head=0.0)])
    (;system,u,q,gate,opening=zeros(0,T),operations=OperationalSeries[])
end

function pumping(rng,n,grid)
    T = length(grid)-1
    lift = uniform(rng,80.0,300.0)
    lakes = [reservoir(rng,:Upper,50.0+lift,30.0),reservoir(rng,:Lower,50.0,30.0)]
    plants = [Plant(name=:Generating,source=:Upper,target=:Lower,pmax=300.0,ramp=nothing),
        Plant(name=:Pumping,source=:Lower,target=:Upper,pmax=300.0,ramp=nothing)]
    flow = uniform(rng,5.0,25.0)
    units = [unit(rng,:Turbine,:Generating,flow,lift,minimum(diff(grid))),
        unit(rng,:Pump,:Pumping,flow,lift,minimum(diff(grid));mode=:pump,table=false)]
    if rand(rng,Bool)
        characteristic = PumpCharacteristic([0.5lift,lift,1.5lift+10],
            [0.65flow,0.60flow,0.55flow],[0.84,0.88,0.86])
        units[2] = updated(units[2];characteristic,turbine_table=nothing)
    end
    u = zeros(Int,2,T)
    q = zeros(2,T)
    for t in 1:T
        i = rand(rng,0:2)
        i==0 && continue
        u[i,t] = 1
        q[i,t] = 0.6flow
    end
    system = HydroSystem(;reservoirs=lakes,plants,generators=units,
        commitment_exclusions=[CommitmentExclusion(name=:Reversible,groups=[[:Turbine],[:Pump]])])
    (;system,u,q,gate=zeros(0,T),opening=zeros(0,T),operations=OperationalSeries[])
end

"""Construct a synthetic basin and controls before adding witness-compatible restrictions.

There is no rejection sampling: a failed hydraulic reference remains a finding.
The saved JSON, rather than the seed alone, is the cross-version reproduction.
"""
function generate_case(seed::Integer, family::AbstractString, size::AbstractString)
    family in FAMILIES || throw(ArgumentError("unknown family: $family"))
    size in SIZES || throw(ArgumentError("unknown size: $size"))
    rng = Xoshiro(seed)
    nrange, trange = size=="small" ? (2:3,3:5) : size=="medium" ? (4:6,6:10) : (8:12,12:24)
    n,T = rand(rng,nrange),rand(rng,trange)
    grid = vcat(0.0,cumsum(rand(rng,(0.25,0.5,1.0,1.5),T)))
    build = getfield(@__MODULE__,Symbol(family))
    data = build(rng,n,grid)
    units = [updated(g;initial_on=data.u[i,1]) for (i,g) in enumerate(data.system.generators)]
    system = updated(data.system;generators=units)
    prices = [uniform(rng,-60.0,160.0) for _ in 1:T]
    prices[rand(rng,1:T)] = rand(rng,(-150.0,0.0,350.0))
    operations = copy(data.operations)
    for r in system.reservoirs
        r.name==:Pond && family=="small_storage" && continue
        base = r.inflow
        push!(operations,series(r.name,:inflow,grid[1:end-1],
            [base*uniform(rng,0.0,2.0) for _ in 1:T]))
    end
    case = ScheduleCase(;name="fuzz $family seed=$seed size=$size",system,grid,prices,operations)
    (;case,u=data.u,q=data.q,gate=data.gate,opening=data.opening)
end

"""Add operating restrictions around an independently reconstructed feasible schedule."""
function restrict_case(seed, case, witness)
    rng = Xoshiro(seed ⊻ 0x1c4b7a21)
    s,grid = case.system,case.grid
    dt = diff(grid)
    operations = copy(case.operations)
    for (i,r) in enumerate(s.reservoirs)
        width = r.vmax-r.vmin
        margin = rand(rng,(0.002,0.03,0.15))*width
        v = witness.storage[i,:]
        # A held restriction covers both ends of its interval, not just its start.
        lower = vcat(min.(v[1:end-1],v[2:end]),last(v)).-margin
        upper = vcat(max.(v[1:end-1],v[2:end]),last(v)).+margin
        push!(operations,series(r.name,:vmin,grid,max.(r.vmin,lower)))
        push!(operations,series(r.name,:vmax,grid,min.(r.vmax,upper)))
    end
    for (i,g) in enumerate(s.generators)
        q = witness.generator_q[i,:]
        p = witness.power[i,:]
        push!(operations,series(g.name,:qmax,grid[1:end-1],
            max.(g.qmin+0.02g.qmax,min.(g.qmax,q.+0.2g.qmax))))
        push!(operations,series(g.name,:pmax,grid[1:end-1],min.(g.pmax,p.+0.2g.pmax)))
        forced = [rand(rng)<0.25 ? Float64(witness.commitment[i,t]) : -1.0 for t in eachindex(dt)]
        push!(operations,series(g.name,:forced_on,grid[1:end-1],forced))
        push!(operations,series(g.name,:maintenance,grid[1:end-1],
            [witness.commitment[i,t]==0 && rand(rng)<0.3 ? 1.0 : 0.0 for t in eachindex(dt)]))
        # A hard zero-discharge schedule requests an off state in the input API.
        # Retain idle committed units by leaving their flow unscheduled.
        if rand(rng)<0.3 && g.characteristic===nothing &&
           all(t->q[t]>0 || witness.commitment[i,t]==0,eachindex(dt))
            push!(operations,series(g.name,:discharge,grid[1:end-1],q))
        end
        # Midpoint-to-midpoint rates, as required by the public ramp convention.
        elapsed = [(dt[t-1]+dt[t])/2 for t in 2:length(dt)]
        qrate = maximum(abs.(diff(q))./elapsed;init=0.0)+0.02g.qmax
        prate = maximum(abs.(diff(p))./elapsed;init=0.0)+0.02g.pmax
        for (attribute,value) in ((:discharge_ramp_up,qrate),(:discharge_ramp_down,qrate),
            (:ramp_up,prate),(:ramp_down,prate))
            push!(operations,series(g.name,attribute,[first(grid)],[value]))
        end
    end
    for (i,e) in enumerate(s.tunnels)
        if !isempty(e.gate_positions) && rand(rng)<0.5
            # A late schedule leaves the earlier gate choices free.
            k = rand(rng,2:length(dt))
            push!(operations,series(e.name,:opening,grid[k:end-1],witness.tunnel_opening[i,k:end]))
        end
        cap = min.(e.capacity,abs.(witness.tunnel_q[i,:]).+0.2e.capacity)
        push!(operations,series(e.name,:capacity,grid[1:end-1],cap))
    end
    for (i,r) in enumerate(s.rivers)
        r.law==:controlled || continue
        lower = max.(r.gate_min,witness.gate[i,:].-0.15)
        upper = min.(1.0,witness.gate[i,:].+0.15)
        push!(operations,series(r.name,:gate_min,grid[1:end-1],lower))
        push!(operations,series(r.name,:gate_max,grid[1:end-1],upper))
        push!(operations,series(r.name,:min_release,grid[1:end-1],0.8witness.river_release[i,:]))
    end
    requirements = isempty(s.generators) ? FlowRequirement[] : [FlowRequirement(
        name=:ObservedStationFlow,generators=[first(s.generators).name],
        min_flow=0.7minimum(witness.generator_q[1,:]),
        max_flow=maximum(witness.generator_q[1,:])+0.1first(s.generators).qmax,
        max_volume=1.2*0.0036sum(dt.*witness.generator_q[1,:]))]
    if !isempty(s.rivers)
        arrivals=witness.arrival_volume[1,:]
        rates=arrivals./(0.0036.*dt)
        push!(requirements,FlowRequirement(name=:ObservedArrival,
            river_arrivals=[first(s.rivers).name],min_flow=minimum(rates)-0.1maximum(abs,rates)-1e-5,
            max_flow=maximum(rates)+0.1maximum(abs,rates)+1e-5,
            min_volume=sum(arrivals)-0.1abs(sum(arrivals))-1e-6,
            max_volume=sum(arrivals)+0.1abs(sum(arrivals))+1e-6))
    end
    rivers=[begin
        initial=something(r.initial_release,isempty(r.history_release) ? 0.0 : last(r.history_release))
        positive=filter(>(0.0),vcat(witness.river_release[i,:],initial))
        minimum_flow=isempty(positive) ? min(0.1,0.01r.capacity) : 0.5minimum(positive)
        change_cost=uniform(rng,0.0,5.0)
        if r.law==:controlled && minimum_flow>=1e-5
            updated(r;flow_change_cost=change_cost,min_allowed_flow=minimum_flow,
                flow_start_cost=uniform(rng,0.0,100.0),flow_stop_cost=uniform(rng,0.0,30.0))
        else
            updated(r;flow_change_cost=change_cost)
        end
    end for (i,r) in enumerate(s.rivers)]
    references=Set(vcat([r.name for r in s.reservoirs],[b.name for b in s.boundaries]))
    surface(name)=begin
        i=findfirst(r->r.name==name,s.reservoirs)
        i===nothing ? fill(only(b.head for b in s.boundaries if b.name==name),length(grid)) :
            head.(Ref(s.reservoirs[i]),witness.storage[i,:])
    end
    ix=nodeindex(s)
    generators=[begin
        p=plantof(s,g)
        if p.source in references && p.target in references
            sign=g.mode==:pump ? -1.0 : 1.0
            heads=sign.*(witness.head[ix[p.source],:]-witness.head[ix[p.target],:])
            vertices=sign.*(surface(p.source)-surface(p.target))
            minimum(heads)>=0 && minimum(vertices)>=0 ? updated(g;
                gross_head_min=0.9minimum(heads),gross_head_max=1.1maximum(heads)+1e-6,
                startup_gross_head_min=0.9minimum(vertices)) : g
        else
            g
        end
    end for g in s.generators]
    r = first(s.reservoirs)
    allowance = 0.15(r.vmax-r.vmin)
    soft = [SoftConstraint(name=:TerminalStorage,object=r.name,quantity=:volume,
        time=last(grid),lower=min(r.vmax,witness.storage[1,end]+0.3allowance),
        max_lower_slack=allowance,lower_cost=uniform(rng,100.0,30000.0))]
    updated(case;system=updated(s;rivers,generators),operations,flow_requirements=requirements,soft_constraints=soft)
end

end

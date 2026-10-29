using OpenSHOP

"""Synthetic integrated basin, with steady reference controls for reproduction.

Volumes are Mm³, heads metres and flows m³/s. The default 2,000 m³ forebay
turns over about 43 times an hour at its reference flow. The response curves
are illustrative input data, not a calibration of a particular watercourse.
"""
# Upper reservoirs share a pressure ring. Station outfalls merge, split and
# rejoin before reaching a small forebay feeding the downstream station.
function basin_case(;hours=48,branches=4,pond_volume=0.002,distributed=true)
    grid=collect(0.0:1.0:hours)
    T=length(grid)-1
    reservoirs=[Reservoir(name=Symbol("Lake",i),v0=40.0,vmin=15.0,vmax=65.0,
        z0=380.0+0.3i,slope=0.25,inflow=6.0,water_value=50000.0,
        level_curve=TableCurve([15.0,30.0,40.0,50.0,65.0],
            [384.0,387.0,390.0,394.0,401.0].+0.3i)) for i in 1:branches]
    push!(reservoirs,Reservoir(name=:Forebay,v0=pond_volume,
        vmin=0.2pond_volume,vmax=1.8pond_volume,
        z0=99.0,slope=1.0/pond_volume,water_value=12000.0))
    tunnels=[Tunnel(name=Symbol("Intake",i),source=Symbol("Lake",i),
        target=Symbol("J",i),resistance=0.035,capacity=40.0,
        gate_positions=i<=2 ? [0.0,0.5,1.0] : Float64[],initial_opening=1.0,
        gate_adjustment_cost=10.0,
        gate_loss_factor_curve=TableCurve([0.0,0.5,1.0],[0.2,0.04,0.0])) for i in 1:branches]
    append!(tunnels,[Tunnel(name=Symbol("Cross",i),source=Symbol("J",i),
        target=Symbol("J",mod1(i+1,branches)),resistance=0.08,capacity=30.0)
        for i in 1:branches])
    plants=[Plant(name=Symbol("Upper",i),source=Symbol("J",i),target=:Forebay,
        discharge_river=Symbol("Tail",i),pmax=60.0,ramp=25.0,
        discharge_ramp_up=8.0,discharge_ramp_down=8.0) for i in 1:branches]
    push!(plants,Plant(name=:Lower,source=:Forebay,target=:Sea,pmax=100.0,
        ramp=50.0,discharge_ramp_up=20.0branches,discharge_ramp_down=20.0branches))
    units=[Generator(name=Symbol("Unit",i),plant=p.name,qmin=2.0,qmax=16.0,
        pmin=1.0,pmax=60.0,hmin=240.0,hmax=330.0,qbest=11.0,hbest=290.0,
        minup=3.0,mindown=2.0,initial_on=1,initial_age=12.0,startup=180.0,
        turbine_table=TurbineTable(heads=[240.0,290.0,330.0],qmin=fill(2.0,3),qmax=fill(16.0,3),
            curves=[TableCurve([0.0,4.0,10.0,16.0],[0.78,0.87,0.93-0.002i,0.9])
                for _ in 1:3])) for (i,p) in enumerate(plants[1:branches])]
    push!(units,Generator(name=:LowerUnit,plant=:Lower,qmin=2.0,qmax=16.0branches,
        pmin=0.5,pmax=100.0,hmin=95.0,hmax=105.0,qbest=10.0branches,hbest=100.0,
        minup=2.0,mindown=1.0,initial_on=1,initial_age=12.0,startup=120.0,
        efficiency=0.92,qcurvature=0.05,hcurvature=0.0))
    function reach(name,source,target,capacity,flow;allocation=1.0,offset=0.0)
        curves=distributed ? [RiverRouting.DelayCurve(q,
            [offset+lag,offset+lag+0.5,offset+lag+1.5],[0.7,0.3])
            for (q,lag) in ((0.0,1.5),(capacity/2,0.7),(capacity,0.2))] :
            RiverRouting.DelayCurve[]
        River(;name,source,target,capacity,allocation,law=:junction,
            curves,deterministic_delay=distributed ? nothing : offset+0.5,
            history_grid=[-24.0,0.0],history_release=[flow],initial_release=flow,
            water_value=12000.0)
    end
    rivers=[reach(Symbol("Tail",i),:auto,:Fork,20.0,6.0;offset=0.2i)
        for i in 1:branches]
    append!(rivers,[reach(:West,:Fork,:Merge,12.0branches,2.4branches;
            allocation=0.4,offset=0.4),
        reach(:East,:Fork,:Merge,18.0branches,3.6branches;allocation=0.6,offset=0.8),
        reach(:Main,:Merge,:Forebay,30.0branches,6.0branches;offset=0.3)])
    system=HydroSystem(;reservoirs,tunnels,plants,generators=units,rivers,
        junctions=[Junction(name=Symbol("J",i),hmin=300.0,hmax=410.0) for i in 1:branches],
        boundaries=[Boundary(name=:Sea,head=0.0)])
    operations=[OperationalSeries(object=Symbol("Lake",i),attribute=:inflow,
        times=grid[1:end-1],values=[6.0*(1+0.25sin((t+2i)/12)) for t in 0:(T-1)])
        for i in 1:branches]
    # Environmental release and finite lookback constraints are independent of
    # the reference trajectory. Neither fixes the economic dispatch schedule.
    requirements=[FlowRequirement(name=:DownstreamFlow,generators=[:LowerUnit],
        min_flow=2.0,max_flow=16.0branches,
        history_grid=[-3.0,0.0],history_flow=[6.0branches],
        rolling_ramps=[RollingRamp(hours=3.0,up=12.0branches,down=12.0branches)])]
    prices=[55.0+35sin(2pi*(t-7)/24)+15sin(2pi*t/8) for t in 0:(T-1)]
    prices[7:24:end].=-20.0
    prices[19:24:end].=170.0
    c=ScheduleCase(;name="synthetic integrated basin",system,grid,prices,operations,
        flow_requirements=requirements)
    u=ones(Int,length(units),T)
    q=fill(6.0,length(units),T)
    q[end,:].=6.0branches
    (;case=c,u,q,gate=zeros(length(rivers),T),opening=ones(length(tunnels),T))
end

"""Shift the previous numerical iterate; hold its final controls in the new tail.

This is a starting point, not a physical continuation or an accepted schedule.
The optimizer must recompute states and flows under the new forecasts.
"""
function shifted_start(x,grid)
    profile=schedule_dict(x)
    intervals=[clamp(searchsortedlast(x.grid,t),1,length(x.grid)-1) for t in grid[1:end-1]]
    vertices=[clamp(searchsortedlast(x.grid,t),1,length(x.grid)) for t in grid]
    for key in ("u","H","tunnel_q","tunnel_opening","generator_q","power",
        "gate","river_release","arrival_volume","shortfall_release")
        haskey(profile,key) && (profile[key]=profile[key][:,intervals])
    end
    profile["V"]=profile["V"][:,vertices]
    profile["grid"]=copy(grid)
    profile
end

"""Replan a fixed-length forecast window after each executed prefix.

Restart carries the executed state and observations. A new case supplies future
prices and inflows; absolute time is retained in every forecast and history.
Only numerical guesses are shifted. No past decision is optimized again.
"""
function rolling_horizon(;horizon=48,advance=6,replans=2,branches=4,time_limit=180.0)
    advance<horizon || throw(ArgumentError("advance must be shorter than the horizon"))
    template=basin_case(;hours=horizon+advance*replans,branches).case
    function window(state,requirements,start,revision)
        grid=collect(start:1.0:(start+horizon))
        indices=[searchsortedfirst(template.grid,t) for t in grid[1:end-1]]
        # Each revision changes forecasts for the unexecuted future only.
        operations=[OperationalSeries(object=z.object,attribute=z.attribute,times=z.times,
            values=z.values.*(1+0.03sin(revision))) for z in template.operations]
        ScheduleCase(;name="rolling synthetic basin",system=state,grid,
            prices=template.prices[indices].+5sin(revision),operations,
            flow_requirements=requirements)
    end
    c=window(template.system,template.flow_requirements,0.0,0)
    warm=nothing
    results=NamedTuple[]
    for revision in 0:replans
        result=schedule_case(c;warm,time_limit,proposal_time_limit=5.0,
            nlp_time_limit=time_limit,max_refinements=1)
        result["accepted"] || error("no accepted schedule for window starting $(first(c.grid))")
        x=result["solution"]
        push!(results,(;case=c,schedule=x,seconds=result["seconds"]))
        println("window $(first(c.grid))–$(last(c.grid)) h: objective=$(x.objective), seconds=$(result["seconds"])")
        revision==replans && break
        start=first(c.grid)+advance
        continued=restart_case(c,x,start)
        c=window(continued.system,continued.flow_requirements,start,revision+1)
        warm=shifted_start(x,c.grid)
    end
    results
end

if abspath(PROGRAM_FILE)==@__FILE__
    rolling_horizon()
end

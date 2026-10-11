# Connected twelve-reservoir benchmark with six units, two cross-linked intakes,
# selectable intake/bypass gates, and delayed many-to-one river connections.
include(joinpath(@__DIR__,"..","examples","mixed_source_network.jl"))
using .MixedSourceNetworkExample: mixed_source_fixture
using OpenSHOP, JuMP, JSON3

function gate_benchmark_copy(x;kwargs...)
    typeof(x)(;merge((;(k=>getfield(x,k) for k in fieldnames(typeof(x)))...),
        (;kwargs...))...)
end

function gate_network_case(hours)
    hours isa Integer && hours>=2 || throw(ArgumentError("at least two hourly intervals required"))
    base=mixed_source_fixture().case.system
    scale=20.0
    reservoirs=Reservoir[];junctions=Junction[];boundaries=Boundary[]
    tunnels=Tunnel[];plants=Plant[];generators=Generator[];rivers=River[]
    river_junctions=RiverJunction[];interlocks=TunnelInterlock[]
    flow_requirements=FlowRequirement[]
    for basin in 1:2
        rename(n)=n==:auto ? n : Symbol(n,"_",basin)
        for r in base.reservoirs
            push!(reservoirs,gate_benchmark_copy(r;name=rename(r.name),
                z0=r.z0+0.4*(basin-1)+(r.name==:UpperA ? 1.0 : r.name==:UpperB ? -16.0 : 0.0),
                slope=r.slope/scale,
                v0=scale*r.v0,vmin=scale*r.vmin,vmax=scale*r.vmax))
        end
        append!(junctions,[gate_benchmark_copy(j;name=rename(j.name),
            hmin=240.0+0.4*(basin-1),hmax=j.hmax+0.4*(basin-1)) for j in base.junctions])
        append!(boundaries,[gate_benchmark_copy(b;name=rename(b.name)) for b in base.boundaries])
        append!(river_junctions,[RiverJunction(name=rename(j.name)) for j in base.river_junctions])
        for e in base.tunnels
            positions=e.name in (:SupplyA,:SupplyB) ? [0.0,1.0] :
                e.name==:RoutedCanal ? [0.0,0.5,1.0] : Float64[]
            push!(tunnels,gate_benchmark_copy(e;name=rename(e.name),source=rename(e.source),
                target=rename(e.target),opening=e.name==:SupplyB ? 0.0 : e.opening,
                gate_positions=positions,gate_adjustment_cost=isempty(positions) ? 0.0 : 10.0,
                discharge_river=e.discharge_river===nothing ? nothing : rename(e.discharge_river)))
        end
        push!(interlocks,TunnelInterlock(name=rename(:SupplyInterlock),
            tunnels=rename.([:SupplyA,:SupplyB]),min_open=1,max_open=1))
        push!(flow_requirements,FlowRequirement(name=rename(:RequiredHighRelease),
            generators=rename.([:UnitA,:UnitB]),min_flow=8.0))
        append!(plants,[gate_benchmark_copy(p;name=rename(p.name),source=rename(p.source),
            target=rename(p.target),discharge_river=rename(p.discharge_river)) for p in base.plants])
        append!(generators,[gate_benchmark_copy(g;name=rename(g.name),plant=rename(g.plant))
            for g in base.generators])
        append!(rivers,[gate_benchmark_copy(r;name=rename(r.name),source=rename(r.source),
            target=rename(r.target)) for r in base.rivers])
    end
    response=LevelFlowCurve([70.0],[TableCurve([-10.0,0.0,10.0],[-40.0,0.0,40.0])])
    push!(rivers,River(name=:LowerConnection,source=:Lower_1,target=:Lower_2,law=:delta,
        capacity=20.0,water_value=0.0,deterministic_delay=0.0,
        curves=RiverRouting.DelayCurve[],delta_level_ref_up_flow_curve=response))
    system=HydroSystem(;reservoirs,junctions,boundaries,tunnels,plants,generators,rivers,
        river_junctions,tunnel_interlocks=interlocks)
    daily=[30.0,30.0,25.0,20.0,35.0,60.0,80.0,100.0,90.0,70.0,50.0,40.0,
        35.0,30.0,25.0,40.0,70.0,100.0,120.0,100.0,60.0,-20.0,-40.0,10.0]
    case=ScheduleCase(;name="twelve-reservoir gates $(hours)h",system=system,
        grid=Float64.(0:hours),prices=[daily[mod1(t,24)] for t in 1:hours],flow_requirements)
    OpenSHOP.validate_inputs(case)
    case
end

function run_gate_network(hours;proposal_seconds=5.0,nlp_seconds=20.0)
    case=gate_network_case(hours)
    prepared=schedule_case(case;proposal_time_limit=proposal_seconds,
        nlp_time_limit=nlp_seconds,max_refinements=0,operational_margin=0.02)
    prepared["accepted"] || error("$(hours)h gate network did not produce an accepted schedule")
    proposal_count=count(p->haskey(p,"u"),prepared["proposals"])
    proposal_count>0 || error("gate network has no feasible joint proposal")
    actual=case_from_dict(prepared["case"]);x=prepared["solution"]
    audit=replay_audit(actual,x;grid=OpenSHOP.refined_grid(actual.grid;factor=4))
    audit["valid"] || error("fourfold replay failed: $(audit["errors"])")
    global_build=@elapsed model=OpenSHOP._build_global_dispatch(actual;warm=x,joint=true)
    start=OpenSHOP._lift_start!(model,actual,x)
    start["valid"] || error("complete global-model start failed: $(start)")
    changes=sum(count(!iszero,diff(vcat(OpenSHOP.tunnel_initial_opening(case.system.tunnels[i]),
        x["tunnel_opening"][i,:]))) for i in axes(x["tunnel_opening"],1))
    hours>=24 && changes==0 && error("daily negative prices did not induce any gate switches")
    return Dict("hours"=>hours,"reservoirs"=>length(case.system.reservoirs),
        "generators"=>length(case.system.generators),"tunnels"=>length(case.system.tunnels),
        "rivers"=>length(case.system.rivers),
        "gate_state_binaries"=>count(v->is_binary(v) && startswith(name(v),"gate_select_"),all_variables(model.m)),
        "global_model_variables"=>num_variables(model.m),
        "global_construction_seconds"=>global_build,"global_start_valid"=>start["valid"],
        "global_start_assigned"=>start["assigned"],
        "accepted"=>prepared["accepted"],"objective"=>x["objective"],
        "feasible_joint_proposals"=>proposal_count,
        "schedule_seconds"=>prepared["seconds"],"gate_adjustment_cost"=>x["gate_adjustment_cost"],
        "position_changes"=>changes,
        "fourfold_replay"=>audit["valid"],
        "storage_replay_difference_Mm3"=>audit["max_storage_difference_Mm3"],
        "max_water_residual_Mm3"=>audit["max_water_residual_Mm3"],
        "global_optimality_proven"=>false,
        "proposal_seconds"=>proposal_seconds,"nlp_seconds"=>nlp_seconds,
        "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),
        "blas_threads"=>OpenSHOP.BLAS.get_num_threads(),
        "cpu_name"=>Sys.CPU_NAME,"kernel"=>string(Sys.KERNEL))
end

if abspath(PROGRAM_FILE)==@__FILE__
    # Exclude compilation and solver loading from measured scheduling calls.
    run_gate_network(2;proposal_seconds=10.0,nlp_seconds=45.0)
    hours=isempty(ARGS) ? [48,168] : parse.(Int,ARGS)
    for horizon in hours
        println(JSON3.write(run_gate_network(horizon)))
    end
end

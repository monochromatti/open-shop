module FuzzWorker

using OpenSHOP, JSON3, LinearAlgebra
import Ipopt, SCIP
include("fuzz_cases.jl")

function checkpoint(path, value)
    writejson(path*".tmp",value)
    mv(path*".tmp",path;force=true)
end

function finding!(report, kind, stage, detail; severity="error")
    push!(report["findings"],Dict("kind"=>kind,"stage"=>stage,
        "severity"=>severity,"detail"=>detail))
end

function matrix(rows, n, T; integer=false)
    length(rows)==n && all(length(r)==T for r in rows) ||
        throw(DimensionMismatch("saved control dimensions"))
    A = zeros(integer ? Int : Float64,n,T)
    for i in 1:n, t in 1:T
        A[i,t] = rows[i][t]
    end
    A
end

function load_input(path)
    input = JSON3.read(read(path,String),Dict{String,Any})
    c = case_from_dict(input["case"])
    controls = input["controls"]
    s,T = c.system,length(c.prices)
    (;case=c,u=matrix(controls["commitment"],length(s.generators),T;integer=true),
        q=matrix(controls["generator_q"],length(s.generators),T),
        gate=matrix(controls["gate"],length(s.rivers),T),
        opening=matrix(controls["tunnel_opening"],length(s.tunnels),T))
end

function input_dict(data)
    Dict("case"=>case_dict(data.case),"controls"=>Dict(
        "commitment"=>data.u,"generator_q"=>data.q,"gate"=>data.gate,
        "tunnel_opening"=>data.opening))
end

reconstruct(data) = dispatch_from_controls(data.case,data.u,data.q,data.gate;
    tunnel_opening=data.opening)

# Independently integrate a held input, including changes inside a control window.
function input_integral(c, object, attribute, a, b, default)
    j = findfirst(z->z.object==object && z.attribute==attribute,c.operations)
    j===nothing && return (b-a)*default
    z = c.operations[j]
    edges = vcat(a,[t for t in z.times if a<t<b],b)
    sum(1:(length(edges)-1);init=0.0) do i
        k = searchsortedlast(z.times,edges[i])
        (edges[i+1]-edges[i])*(k==0 ? default : z.values[k])
    end
end

# This ledger uses the returned physical boundary exchange and public input
# series, rather than the optimization constraints or validation residuals.
function water_error(c, trajectory)
    inventory = vec(sum(trajectory.storage;dims=1)+sum(trajectory.transit;dims=1))
    expected = first(inventory)
    error = 0.0
    for t in eachindex(c.prices)
        natural = 0.0
        for r in Iterators.flatten((c.system.reservoirs,c.system.rivers))
            natural += input_integral(c,r.name,:inflow,c.grid[t],c.grid[t+1],r.inflow)
        end
        expected += 0.0036*(natural-(c.grid[t+1]-c.grid[t])*trajectory.boundary_outflow[t])
        error = max(error,abs(inventory[t+1]-expected))
    end
    error
end

function check_schedule!(report,c,x,stage,replay_factor)
    audit = validate(c,x)
    audit["valid"] || finding!(report,"invalid_schedule",stage,audit)
    components = objective_components(c,x)
    tolerance = 1e-6+1e-7max(1.0,abs(x.objective))
    abs(components.total-x.objective)<=tolerance ||
        finding!(report,"objective_accounting",stage,Dict(
            "reported"=>x.objective,"recomputed"=>components.total))
    trajectory = simulate(c,x.generator_q,x.gate;commitment=x.commitment,
        tunnel_opening=x.tunnel_opening,head_reference=x.head)
    trajectory.converged || finding!(report,"hydraulic_nonconvergence",stage,
        Dict("raw_residual"=>trajectory.residual))
    ledger = water_error(c,trajectory)
    ledger<=1e-6 || finding!(report,"water_conservation",stage,ledger)
    fine = replay_audit(c,x;grid=OpenSHOP.refined_grid(c.grid;factor=replay_factor))
    if !fine["valid"]
        finding!(report,"replay_sensitivity",stage,fine;severity="warning")
    end
    Dict("audit"=>audit,"water_error_Mm3"=>ledger,"objective_components"=>Dict(
        string(k)=>v for (k,v) in pairs(components)),"replay"=>fine)
end

function continuation_check!(report,c,x)
    # The first complete interval is executed; the remainder retains its cohorts,
    # gate positions, output history and outstanding commitment obligations.
    k = max(2,cld(length(c.grid),2))
    try
        tail = restart_case(c,x,c.grid[k])
        y = dispatch_from_controls(tail,x.commitment[:,k:end],x.generator_q[:,k:end],
            x.gate[:,k:end];tunnel_opening=x.tunnel_opening[:,k:end],
            head_reference=x.head[:,k:end])
        audit = validate(tail,y)
        difference = maximum(abs,y.storage-x.storage[:,k:end];init=0.0)
        audit["valid"] && difference<=1e-6 || finding!(report,"restart_difference",
            "reference",Dict("audit"=>audit,"storage_difference_Mm3"=>difference))
    catch error
        finding!(report,"restart_exception","reference",sprint(showerror,error))
    end
end

function topology_features(c)
    s = c.system
    Dict("reservoirs"=>length(s.reservoirs),"pressure_junctions"=>length(s.junctions),
        "tunnels"=>length(s.tunnels),"plants"=>length(s.plants),
        "units"=>length(s.generators),"rivers"=>length(s.rivers),
        "intervals"=>length(c.prices),"horizon_hours"=>last(c.grid)-first(c.grid),
        "minimum_reservoir_Mm3"=>isempty(s.reservoirs) ? nothing : minimum(r.v0 for r in s.reservoirs),
        "tabulated_reservoirs"=>count(r->r.level_curve!==nothing,s.reservoirs),
        "tabulated_turbines"=>count(g->g.turbine_table!==nothing,s.generators),
        "fixed_speed_pumps"=>count(g->g.characteristic!==nothing,s.generators),
        "selectable_gates"=>count(e->!isempty(e.gate_positions),s.tunnels),
        "distributed_rivers"=>count(r->r.deterministic_delay===nothing,s.rivers),
        "river_laws"=>sort!(unique(string(r.law) for r in s.rivers)),
        "operation_attributes"=>sort!(unique(string(o.attribute) for o in c.operations)),
        "soft_constraints"=>length(c.soft_constraints))
end

function check_result!(report,c,result,stage,reference,options)
    audit = nothing
    if result.solution!==nothing
        # A rejected iterate may be useful evidence, but is not a promised schedule.
        original = validate(c,result.solution)
        if result.accepted || original["valid"]
            audit = check_schedule!(report,c,result.solution,stage,options["replay_factor"])
        end
        result.accepted && !original["valid"] &&
            finding!(report,"accepted_invalid_schedule",stage,original)
        if stage=="local" && result.accepted
            result.solution.commitment==reference.commitment &&
                result.solution.tunnel_opening==reference.tunnel_opening ||
                finding!(report,"fixed_controls_changed",stage,"local solve changed fixed controls")
        end
    elseif result.accepted
        finding!(report,"accepted_without_schedule",stage,result.status)
    end
    known = reference.objective
    for previous in values(get(report,"stages",Dict()))
        checks = previous["checks"]
        if checks!==nothing && checks["audit"]["valid"]
            known=max(known,previous["objective"])
        end
    end
    if result.solution!==nothing && validate(c,result.solution)["valid"]
        known = max(known,result.solution.objective)
    end
    tolerance = 1e-5+1e-7max(1.0,abs(known))
    function check_bound(name,bound,source)
        if bound!==nothing && bound<known-tolerance &&
           !any(f->f["kind"]=="invalid_upper_bound" && f["stage"]==name,report["findings"])
            finding!(report,"invalid_upper_bound",name,Dict(
                "source"=>source,"bound"=>bound,"known_feasible_objective"=>known))
        end
    end
    for (name,bound) in (("global_bound",result.global_bound),
        ("rejected_global_bound",get(result.diagnostics,"rejected_global_bound",nothing)))
        check_bound(stage,bound,name)
    end
    # A later solve may disprove a bound that enclosed all earlier incumbents.
    for (name,previous) in get(report,"stages",Dict())
        check_bound(name,previous["global_bound"],"global_bound")
    end
    if stage!="local" && result.status in ("INFEASIBLE","INFEASIBLE_OR_UNBOUNDED")
        finding!(report,"false_infeasibility",stage,result.status)
    end
    start = get(result.diagnostics,"start_audit",nothing)
    stage=="warm" && start!==nothing && !start["valid"] &&
        finding!(report,"reference_rejected_by_model",stage,start)
    if result.global_certificate
        if result.global_bound===nothing || result.solution===nothing
            finding!(report,"certificate_without_bounds",stage,result.status)
        else
            gap = max(0.0,result.global_bound-result.solution.objective)
            gap<=options["relative_gap"]*max(1.0,abs(result.solution.objective))+tolerance ||
                finding!(report,"certificate_gap",stage,gap)
        end
    end
    if !result.accepted
        kind = occursin("LIMIT",result.status) || occursin("BUDGET_EXHAUSTED",result.status) ?
            "search_inconclusive" : "solver_failure"
        finding!(report,kind,stage,result.status;severity="warning")
    elseif stage=="local" && result.objective<reference.objective-tolerance
        finding!(report,"local_objective_regression",stage,
            Dict("reference"=>reference.objective,"returned"=>result.objective);severity="warning")
    end
    Dict("status"=>result.status,"accepted"=>result.accepted,"objective"=>result.objective,
        "global_bound"=>result.global_bound,"relative_gap"=>result.relative_gap,
        "global_certificate"=>result.global_certificate,"total_seconds"=>result.total_seconds,
        "native_incumbent_available"=>get(get(result.diagnostics,"scip_diagnostics",Dict()),
            "native_incumbent_available",nothing),
        "incumbent_source"=>get(result.diagnostics,"incumbent_source",nothing),"checks"=>audit)
end

function run(request, directory)
    BLAS.set_num_threads(1)
    # JSON numbers need not preserve all 64 integer bits through generic readers.
    seed = request["seed"] isa AbstractString ? parse(Int,request["seed"]) : Int(request["seed"])
    report = Dict{String,Any}("seed"=>string(seed),"family"=>request["family"],
        "size"=>request["size"],"findings"=>Any[],"stages"=>Dict{String,Any}(),
        "versions"=>Dict("Julia"=>string(VERSION),"OpenSHOP"=>string(pkgversion(OpenSHOP)),
            "Ipopt.jl"=>string(pkgversion(Ipopt)),"SCIP.jl"=>string(pkgversion(SCIP)),
            "JuMP"=>string(pkgversion(OpenSHOP.JuMP)),
            "ForwardDiff"=>string(pkgversion(OpenSHOP.ForwardDiff))))
    started = time()
    phase = "generation"
    progress() = checkpoint(joinpath(directory,"progress.json"),Dict(
        "phase"=>phase,"elapsed_seconds"=>time()-started))
    try
        progress()
        data = if request["replay"]===nothing
            FuzzCases.generate_case(seed,request["family"],request["size"])
        else
            load_input(request["replay"])
        end
        # Persist before even constructing a reference, including failure cases.
        checkpoint(joinpath(directory,"input.json"),input_dict(data))
        report["features"] = topology_features(data.case)
        phase = "input_validation"
        progress()
        data = merge(data,(;case=case_from_dict(case_dict(data.case))))
        phase = "reference"
        progress()
        reference = reconstruct(data)
        original = validate(data.case,reference)
        if !original["valid"]
            finding!(report,"reference_failure",phase,original)
        else
            if request["replay"]===nothing
                c = FuzzCases.restrict_case(seed,data.case,reference)
                data = merge(data,(;case=c))
                checkpoint(joinpath(directory,"input.json"),input_dict(data))
                phase = "restricted_input_validation"
                progress()
                data = merge(data,(;case=case_from_dict(case_dict(c))))
                phase = "reference"
                progress()
                reference = reconstruct(data)
            end
            c = data.case
            report["features"] = topology_features(c)
            writejson(joinpath(directory,"reference.json"),reference)
            report["reference"] = check_schedule!(report,c,reference,phase,request["replay_factor"])
            checkpoint(joinpath(directory,"partial.json"),report)
            if report["reference"]["audit"]["valid"]
                continuation_check!(report,c,reference)
                for stage in request["stages"]
                    phase = stage
                    progress()
                    try
                        result = if stage=="local"
                            solve_case(c;u=data.u,tunnel_opening=data.opening,warm=reference,
                                time_limit=request["time_limit"])
                        else
                            solve(c;initial=stage=="warm" ? reference : nothing,
                                time_limit=request["time_limit"],relative_gap=request["relative_gap"],
                                diagnostics_path=joinpath(directory,stage*".scip.log"))
                        end
                        writejson(joinpath(directory,stage*".json"),result)
                        report["stages"][stage] = check_result!(report,c,result,stage,reference,request)
                    catch error
                        finding!(report,"exception",stage,sprint(showerror,error,catch_backtrace()))
                    end
                    checkpoint(joinpath(directory,"partial.json"),report)
                    GC.gc()
                end
            end
        end
    catch error
        finding!(report,"exception",phase,sprint(showerror,error,catch_backtrace()))
    end
    report["elapsed_seconds"] = time()-started
    report["outcome"] = any(f->f["severity"]=="error",report["findings"]) ? "fail" :
        isempty(report["findings"]) ? "pass" : "review"
    checkpoint(joinpath(directory,"result.json"),report)
    report
end

function main(args)
    length(args)==1 || error("usage: fuzz_worker.jl REQUEST.json")
    request = JSON3.read(read(only(args),String),Dict{String,Any})
    run(request,dirname(abspath(only(args))))
end

end

if abspath(PROGRAM_FILE)==@__FILE__
    FuzzWorker.main(ARGS)
end

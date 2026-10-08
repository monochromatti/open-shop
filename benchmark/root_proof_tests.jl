using Test
include("root_proof_profile.jl")

@testset "Root proof measurement guards" begin
    points=[(seconds=1.0,upper=120.0,native_seconds=0.5,run=1),
            (seconds=3.0,upper=107.0,native_seconds=2.5,run=1)]
    thresholds=proof_thresholds(points,100.0,5.0)
    @test only(filter(t->t["relative_gap"]==0.07,thresholds))["first_observed_seconds"]==3.0
    @test only(filter(t->t["gap_percent"]==5.0,thresholds))["right_censored"]
    @test only(filter(t->t["gap_percent"]==5.0,thresholds))["first_observed_seconds"]===nothing
    @test PROOF_PROFILES["current"].parameters==Pair{String,Any}[]
    @test PROOF_PROFILES["obbt_1"].parameters[1].second==1.0
    @test PROOF_PROFILES["obbt_quarter"].parameters[1].second==0.25
    @test PROOF_PROFILES["obbt_cap30"].cap==30.0
    @test_throws ErrorException install_root_progress!(SCIP.Optimizer(),time_ns();obbt_cap_seconds=-1.0)
    @test_throws ErrorException install_root_progress!(SCIP.Optimizer(),time_ns();obbt_cap_seconds=NaN)
end

@testset "Completed OBBT boundary stop" begin
    c=readcase(joinpath(@__DIR__,"cases","turbine-tables.json"))
    prepared=schedule_case(c;proposal_time_limit=5.0,nlp_time_limit=20.0,
        max_refinements=0,operational_margin=0.1)
    @test prepared["accepted"]
    seed=prepared["solution"]
    q=copy(seed["generator_q"]);q[seed["u"].==0].=0.0
    initial=dispatch_from_controls(c,seed["u"],q,seed["gate"])
    observed=Any[]
    setup=b->append!(observed,install_root_progress!(unsafe_backend(b.m),time_ns();obbt_cap_seconds=0.0))
    # First call warms the physical builder and native callback compilation.
    OpenSHOP._solve(c;initial,time_limit=10.0,relative_gap=1e-8,optimizer_setup=setup)
    empty!(observed)
    result=OpenSHOP._solve(c;initial,time_limit=10.0,relative_gap=1e-8,optimizer_setup=setup)
    stats=root_progress_statistics(observed)
    @test isempty(stats["errors"])
    @test result["accepted"]
    @test result["global_bound"]>=result["feasible_lower_bound"]-1e-6
    @test !isempty(stats["points"])
    @test stats["obbt_disabled_at_boundary"]
    @test stats["obbt_calls_at_disable"]>=1
    changed=[p for p in stats["points"] if p["obbt_frequency"]==-1]
    @test !isempty(changed)
    @test all(p->p["obbt_calls"]==stats["obbt_calls_at_disable"],changed)
    prop=SCIP.SCIPfindProp(first(observed).optimizer,"obbt")
    @test SCIP.SCIPpropGetFreq(prop)==-1
    @test SCIP.SCIPpropGetNCalls(prop)==stats["obbt_calls_at_disable"]
    @test all(p->p["depth"]==0,stats["points"])
    @test all(p->p["local_root_lp_objective"]!==nothing && isfinite(p["local_root_lp_objective"]),stats["points"])
    # The observation is separate from the global-bound certificate path.
    @test occursin("diagnostic only",stats["scope"])
    @test all(e->e.skipped>=0,observed)
end

@testset "Root LP objective converts maximization sign and offset" begin
    m=Model(SCIP.Optimizer);set_silent(m)
    set_optimizer_attribute(m,"presolving/maxrounds",0)
    @variable(m,0<=x<=2)
    @variable(m,0<=y<=2)
    @constraint(m,x+2y<=3)
    @constraint(m,2x+y<=3)
    @objective(m,Max,(7+3x+3y)/10000)
    JuMP.MOI.Utilities.attach_optimizer(backend(m))
    observers=install_root_progress!(unsafe_backend(m),time_ns())
    optimize!(m)
    stats=root_progress_statistics(observers)
    @test isempty(stats["errors"])
    @test objective_value(m)≈13/10000
    @test !isempty(stats["points"])
    @test all(p->p["local_root_lp_objective"]≈13.0,stats["points"])
end

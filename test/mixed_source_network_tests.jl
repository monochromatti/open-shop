include(joinpath(@__DIR__,"..","examples","mixed_source_network.jl"))
using .MixedSourceNetworkExample: mixed_source_fixture, mixed_source_small_fixture

# Independent piecewise rate integration and historical inventory calculations;
# these deliberately do not use the routing/operating helpers under test.
function mixed_source_integral(c,object,initial,a,b)
    operation=findfirst(z->z.object==object && z.attribute==:inflow,c.operations)
    operation===nothing && return initial*(b-a)
    series=c.operations[operation]
    cuts=sort(unique(vcat(a,[t for t in series.times if a<t<b],b)))
    sum(begin
        k=searchsortedlast(series.times,cuts[i])
        rate=k==0 ? initial : series.values[k]
        rate*(cuts[i+1]-cuts[i])
    end for i in 1:(length(cuts)-1))
end

function mixed_source_history_inventory(r)
    rate=first(r.history_release)
    @assert all(==(rate),r.history_release)
    if r.deterministic_delay!==nothing
        return 0.0036*rate*r.deterministic_delay
    end
    refs=[z.reference_flow for z in r.curves]
    index=clamp(searchsortedlast(refs,rate),1,length(refs)-1)
    fraction=(rate-refs[index])/(refs[index+1]-refs[index])
    means=[sum(curve.weights .* (curve.edges[1:end-1]+curve.edges[2:end])/2)
        for curve in r.curves]
    0.0036*rate*((1-fraction)*means[index]+fraction*means[index+1])
end

function mixed_source_conservation(c,z)
    initial_inventory=sum(mixed_source_history_inventory(r) for r in c.system.rivers)
    @test sum(z["transit"][:,1])≈initial_inventory atol=1e-10
    initial=sum(r.v0 for r in c.system.reservoirs)+initial_inventory
    for edge in 2:length(c.grid)
        horizon=c.grid[edge]
        incoming=sum(mixed_source_integral(c,r.name,r.inflow,first(c.grid),horizon)
            for r in c.system.reservoirs)+
            sum(mixed_source_integral(c,r.name,r.inflow,first(c.grid),horizon)
                for r in c.system.rivers)
        outgoing=sum(diff(c.grid)[1:edge-1] .* z["boundary_outflow"][1:edge-1])
        @test sum(z["V"][:,edge])+sum(z["transit"][:,edge])≈
            initial+0.0036*(incoming-outgoing) atol=2e-8
    end
end

@testset "Mixed plant outfalls, natural river inflows and exact objective" begin
    fixture=mixed_source_small_fixture();c=fixture.case
    @test OpenSHOP.validate_inputs(c)
    @test length(c.system.river_junctions)==3
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    public=case_dict(c)
    @test isempty(public["river_junctions"])
    @test public["plants"][1]["discharge_river"]==:Tail1
    x=dispatch_from_controls(c,fixture.u,fixture.q,fixture.gates)
    @test x["validation"]["valid"]
    @test replay_audit(c,x)["valid"]
    # The Unit1 off pulse [.5,1.5] reaches the receiver at[1.05,2.05].
    expected_arrival=[5.75,4.85,4.65]
    @test x["arrival_volume"][3,:]≈0.0036 .* diff(c.grid) .* expected_arrival atol=1e-10
    receiver=[5.0,5.00495,5.01161,5.01755]
    @test vec(x["V"][3,:])≈receiver atol=1e-10
    @test x["V"][1:2,:]≈fill(5.0,2,4) atol=1e-10
    receiver_heads=49.0 .+ 0.2 .* (receiver[1:end-1]+receiver[2:end])/2
    expected_power=0.00981 .* fixture.q .* [0.9*(100-h) for h in receiver_heads]'
    expected_power[2,:].=0.00981*0.9 .* fixture.q[2,:] .* (120 .- receiver_heads)
    @test x["power"]≈expected_power atol=1e-9
    expected_objective=sum(c.prices[t]*(c.grid[t+1]-c.grid[t])*expected_power[i,t]
        for i in 1:2,t in 1:3)
    @test x["objective"]≈expected_objective atol=1e-8
    z=simulate(c,fixture.q,fixture.gates)
    mixed_source_conservation(c,z)

    # Restart retains delayed cohorts from before the commitment change.
    continued=restart_case(c,x,0.5)
    remainder=dispatch_from_controls(continued,fixture.u[:,2:end],
        fixture.q[:,2:end],fixture.gates[:,2:end])
    @test remainder["validation"]["valid"]
    @test replay_audit(continued,remainder)["valid"]
    @test remainder["V"]≈x["V"][:,2:end] atol=1e-9
    @test remainder["terminal_transit"]≈x["terminal_transit"] atol=1e-10
    @test case_dict(case_from_dict(case_dict(continued)))==case_dict(continued)

    local_result=solve_case(c;u=fixture.u,warm=x,time_limit=45.0)
    @test get(get(local_result,"validation",Dict()),"valid",false)
    if haskey(local_result,"objective")
        @test local_result["objective"]≈expected_objective atol=1e-3
        @test replay_audit(c,local_result)["valid"]
    end
    proposal=OpenSHOP.propose_commitment(c;reference=x,time_limit=10.0)
    @test haskey(proposal,"u")
    if haskey(proposal,"u")
        @test proposal["u"]==fixture.u
        @test admissible(c,proposal["u"])
    end
    result=solve(c;initial=x,time_limit=30.0,relative_gap=1e-4)
    @test result["accepted"]
    @test result["global_certificate"]
    @test !result["commitment_fixed"]
    if result["solution"]!==nothing
        @test result["solution"]["objective"]≈expected_objective atol=1e-3
    end
end

@testset "Connected mixed-source network with pressurized loop and nested merges" begin
    for distributed in (false,true)
        fixture=mixed_source_fixture(;distributed);c=fixture.case
        @test OpenSHOP.validate_inputs(c)
        @test length(c.system.reservoirs)==6
        @test length(c.system.plants)==3
        @test length(c.system.tunnels)==5
        @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
        if distributed
            @test all(length(r.curves)==3 for r in c.system.rivers)
        else
            @test all(r.deterministic_delay!==nothing for r in c.system.rivers)
        end
        x=dispatch_from_controls(c,fixture.u,fixture.q,fixture.gates)
        @test x["validation"]["valid"]
        coarse_replay=replay_audit(c,x)
        @test coarse_replay["valid"]==!distributed
        z=simulate(c,fixture.q,fixture.gates)
        @test z["converged"]
        mixed_source_conservation(c,z)
        names=Dict(r.name=>i for (i,r) in enumerate(c.system.rivers))
        @test x["tunnel_q"][3,1]<-0.5
        @test x["tunnel_q"][3,4]>0.5 # A shutdown reverses the intake cross-flow.
        @test all(x["tunnel_q"][5,:].>0)
        @test x["river_release"][names[:CanalRiver],:]≈x["tunnel_q"][5,:] atol=1e-9
        @test x["river_release"][names[:TailB],:]≈fixture.q[2,:] atol=1e-9
        @test x["river_release"][names[:TailA],4]≈0.8 atol=1e-9 # Off unit, natural inflow remains.
        # Pond receives routed arrivals, with no instantaneous credit from the
        # canal or HighA/HighB despite being their hydraulic target.
        pond_net=x["arrival_volume"][names[:UpperMerge],1]/(0.0036*diff(c.grid)[1])-
            fixture.q[3,1]-x["river_release"][names[:PondOutlet],1]-
            x["river_release"][names[:PondSpill],1]
        @test x["V"][3,2]-x["V"][3,1]≈0.0036*diff(c.grid)[1]*pond_net atol=1e-9
        @test fixture.q[1,1]+fixture.q[2,1]+x["tunnel_q"][5,1]>20.0
        for t in eachindex(c.prices)
            source=x["H"][1,t];reference=x["H"][3,t]
            opening=c.grid[t]<2.0 || c.grid[t]>=4.2 ? 1.0 : 0.8
            @test x["tunnel_q"][5,t]≈sqrt(opening*(source-reference)/1.5) atol=1e-8
        end
        # A restart after several distinct delay edges preserves all inventories.
        continued=restart_case(c,x,2.0)
        remainder=dispatch_from_controls(continued,fixture.u[:,4:end],
            fixture.q[:,4:end],fixture.gates[:,4:end])
        @test remainder["validation"]["valid"]
        @test remainder["V"]≈x["V"][:,4:end] atol=1e-8
        @test remainder["terminal_transit"]≈x["terminal_transit"] atol=1e-8
        @test replay_audit(continued,remainder)["valid"]==!distributed

        local_model=OpenSHOP._build_dispatch(c;u=fixture.u,warm=x,feasibility_only=true)
        @test local_model.m isa JuMP.Model
        global_model=OpenSHOP._build_global_dispatch(c;joint=true,warm=x)
        lift=OpenSHOP._lift_start!(global_model,c,x)
        @test lift["valid"]
        local_result=solve_case(c;u=fixture.u,warm=x,feasibility_only=true,time_limit=20.0)
        @test get(get(local_result,"validation",Dict()),"valid",false)
        if !distributed && haskey(local_result,"validation") && local_result["validation"]["valid"]
            @test replay_audit(c,local_result)["valid"]
        end
        if distributed
            # Coarse confluence averages are a discretization, not exact waves.
            # Refine physical integration while retaining the original operation
            # windows, exactly as replay_audit does; no physical data or tolerance
            # is changed and no optimization is repeated on this refined grid.
            grid=OpenSHOP.refined_grid(c.grid;factor=8)
            indices=[searchsortedlast(c.grid,(grid[t]+grid[t+1])/2)
                for t in 1:length(grid)-1]
            physical_case=OpenSHOP.without_interval_controls(with_grid(c,grid))
            physical=dispatch_from_controls(physical_case,fixture.u[:,indices],
                fixture.q[:,indices],fixture.gates[:,indices])
            @test physical["validation"]["valid"]
            windows=OpenSHOP.operating_window_values(c,physical,grid,indices,fixture.u)
            @test maximum(values(OpenSHOP.operating_residuals(c,windows)))<=1e-4
            fine_replay=replay_audit(physical_case,physical)
            @test fine_replay["valid"]
            refined_restart=restart_case(physical_case,physical,2.0)
            edge=findfirst(==(2.0),grid)
            refined_remainder=dispatch_from_controls(refined_restart,
                physical["u"][:,edge:end],physical["generator_q"][:,edge:end],
                physical["gate"][:,edge:end])
            @test refined_remainder["validation"]["valid"]
            @test refined_remainder["V"]≈physical["V"][:,edge:end] atol=1e-8
            @test refined_remainder["terminal_transit"]≈physical["terminal_transit"] atol=1e-8
            @test replay_audit(refined_restart,refined_remainder)["valid"]
        end
    end
end

@testset "Neighboring distributed references coexist with deterministic reaches" begin
    fixture=mixed_source_small_fixture(;mixed=true);c=fixture.case
    @test length(c.system.rivers[1].curves)==3
    @test c.system.rivers[1].deterministic_delay===nothing
    @test all(r.deterministic_delay!==nothing for r in c.system.rivers[2:end])
    x=dispatch_from_controls(c,fixture.u,fixture.q,fixture.gates)
    @test x["validation"]["valid"]
    @test x["river_release"][1,:]≈[2.25,0.25,2.25] atol=1e-10
    @test replay_audit(c,x)["valid"]
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    z=simulate(c,fixture.q,fixture.gates)
    mixed_source_conservation(c,z)
    continued=restart_case(c,x,0.5)
    remainder=dispatch_from_controls(continued,fixture.u[:,2:end],
        fixture.q[:,2:end],fixture.gates[:,2:end])
    @test remainder["validation"]["valid"]
    @test remainder["V"]≈x["V"][:,2:end] atol=1e-9
    @test remainder["terminal_transit"]≈x["terminal_transit"] atol=1e-9
    @test replay_audit(continued,remainder)["valid"]
    local_result=solve_case(c;u=fixture.u,warm=x,time_limit=45.0)
    @test get(get(local_result,"validation",Dict()),"valid",false)
    if haskey(local_result,"validation") && local_result["validation"]["valid"]
        @test local_result["objective"]≈x["objective"] atol=1e-3
        @test replay_audit(c,local_result)["valid"]
    end
    proposal=OpenSHOP.propose_commitment(c;reference=x,time_limit=10.0)
    @test haskey(proposal,"u")
    if haskey(proposal,"u")
        @test proposal["u"]==fixture.u
        @test admissible(c,proposal["u"])
    end
    result=solve(c;initial=x,time_limit=30.0,relative_gap=1e-4)
    @test result["accepted"]
    @test result["global_certificate"]
    @test !result["transport_exact"]
    @test !result["commitment_fixed"]
end

@testset "Mixed-source dispatch actually optimizes a free unit" begin
    fixture=mixed_source_small_fixture()
    c=OpenSHOP._river_replace(fixture.case;name="optimized mixed source unit",
        operations=[z for z in fixture.case.operations
            if !(z.object==:Unit2 && z.attribute==:discharge)])
    initial=dispatch_from_controls(c,fixture.u,fixture.q,fixture.gates)
    @test initial["validation"]["valid"]
    maximum_q=copy(fixture.q);maximum_q[2,:].=5.0
    best=dispatch_from_controls(c,fixture.u,maximum_q,fixture.gates)
    @test best["validation"]["valid"]
    @test replay_audit(c,best)["valid"]
    # Extra Unit2 discharge reaches Receiver after .75+.3=1.05h.
    receiver=[5.0,5.00495,5.01485,5.02799]
    upper_b=[5.0,4.9964,4.9892,4.982]
    @test best["arrival_volume"][3,:]≈0.0036 .* diff(c.grid) .* [5.75,5.75,6.65] atol=1e-10
    @test vec(best["V"][3,:])≈receiver atol=1e-10
    @test vec(best["V"][2,:])≈upper_b atol=1e-10
    receiver_heads=49.0 .+0.2 .* (receiver[1:end-1]+receiver[2:end])/2
    source_b=119.0 .+0.2 .* (upper_b[1:end-1]+upper_b[2:end])/2
    expected_power=zeros(2,3)
    expected_power[1,:].=0.00981*0.9 .* maximum_q[1,:] .* (100 .-receiver_heads)
    expected_power[2,:].=0.00981*0.9 .* maximum_q[2,:] .* (source_b .-receiver_heads)
    expected_objective=sum(c.prices[t]*diff(c.grid)[t]*expected_power[i,t]
        for i in 1:2,t in 1:3)
    @test best["objective"]≈expected_objective atol=1e-8
    @test best["objective"]>initial["objective"]+100.0
    local_result=solve_case(c;u=fixture.u,warm=initial,time_limit=20.0)
    @test get(get(local_result,"validation",Dict()),"valid",false)
    local_replay=nothing
    if haskey(local_result,"generator_q")
        @test vec(local_result["generator_q"][2,:])≈fill(5.0,3) atol=2e-4
        @test local_result["objective"]≈expected_objective atol=1e-3
        local_replay=replay_audit(c,local_result)
        @test local_replay["valid"]
    end
    result=solve(c;initial,time_limit=30.0,relative_gap=1e-4)
    @test result["accepted"]
    if result["global_certificate"]
        @test result["global_bound"]!==nothing
        @test result["global_bound"]>=expected_objective-1e-6
    else
        # An unresolved first root LP bound must remain diagnostic, not a proof.
        @test result["global_bound"]===nothing
        @test result["bound_rejection"]==
            "SCIP terminated without a finite first root LP bound despite LP iterations; native bound is retained only in diagnostics"
        diagnostics=result["scip_diagnostics"]
        @test result["status"]=="OPTIMAL"
        @test diagnostics["lp_iterations"]>0
        @test diagnostics["first_root_lp_upper_bound"]===nothing
        @test isfinite(diagnostics["final_upper_bound"])
        @test diagnostics["final_upper_bound"]>=expected_objective-1e-6
    end
    @test !result["commitment_fixed"]
    if result["solution"]!==nothing
        @test vec(result["solution"]["generator_q"][2,:])≈fill(5.0,3) atol=2e-4
        @test all(result["solution"]["u"][2,:].==1)
        @test result["solution"]["objective"]≈expected_objective atol=1e-3
    end
end

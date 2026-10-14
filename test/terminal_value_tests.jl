include("../examples/terminal_value.jl")
using .TerminalValueExample: terminal_value_case
tv_copy(x;kwargs...)=typeof(x)(;merge((; (k=>getfield(x,k) for k in fieldnames(typeof(x)))...),
    (;kwargs...))...)

@testset "Energy-to-volume valuation uses explicit reference conditions" begin
    @test energy_conversion_factor(120.0,0.90)≈294.3 atol=1e-12
    @test energy_conversion_factor(80.0,0.92)≈200.56 atol=1e-12
    @test energy_conversion_factor(0.0,1.0)==0.0
    @test volume_water_value(30.0,294.3;downstream_value=6016.8)≈14845.8 atol=1e-10
    @test volume_water_value(-10.0,200.56)≈-2005.6 atol=1e-10
    for (h,eta) in ((-1.0,0.9),(100.0,0.0),(100.0,1.1),(Inf,0.9),(100.0,NaN))
        @test_throws ArgumentError energy_conversion_factor(h,eta)
    end
    @test_throws ArgumentError volume_water_value(10.0,-1.0)
    @test_throws ArgumentError volume_water_value(Inf,200.56)
    @test_throws ArgumentError volume_water_value(30.0,200.56;downstream_value=NaN)
end

@testset "Cascade and historical transit objective agree with hand accounting" begin
    c=terminal_value_case();u=ones(Int,2,2);q=[4.0 4.0;2.0 2.0]
    x=dispatch_from_controls(c,u,q,zeros(1,2))
    storage=[2.0 1.9856 1.9712;2.0 2.0036 2.0108]
    heads=[199.9928 199.9784;80.0018 80.0072]
    expected_power=[0.00981*4.0*0.90*(heads[1,t]-heads[2,t]) for t in 1:2]'
    expected_power=vcat(expected_power,[0.00981*2.0*0.92*heads[2,t] for t in 1:2]')
    @test x.storage≈storage atol=1e-10
    @test x.head[1:2,:]≈heads atol=1e-10
    @test x.power≈expected_power atol=1e-9
    @test only(x.terminal_transit)≈0.0072 atol=1e-12
    @test validate(c,x)["valid"]
    @test replay_audit(c,x)["valid"]
    parts=objective_components(c,x)
    @test parts.energy_revenue≈100.0*sum(expected_power) atol=1e-8
    @test parts.reservoir_value≈14845.8*(-0.0288)+6016.8*0.0108 atol=1e-8
    @test parts.transit_value≈6016.8*(0.0072-0.0036) atol=1e-8
    @test parts.transition_cost==8.0
    @test parts.gate_cost==parts.release_penalty==0.0
    @test parts.total≈x.objective atol=1e-8
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)

    # Initial values are subtracted constants, including historical transit.
    absolute=parts.energy_revenue-8.0+
        14845.8*storage[1,end]+6016.8*storage[2,end]+6016.8*0.0072
    offset=14845.8*2.0+6016.8*2.0+6016.8*0.0036
    @test absolute-offset≈parts.total atol=1e-8
    data=schedule_dict(x);data["objective"]+=1.0
    @test !validate(c,data)["valid"]
    @test objective_components(c,data).total≈parts.total atol=1e-8

    for result in (solve_case(c;u,warm=x,time_limit=30.0),
        solve(c;initial=x,time_limit=60.0,relative_gap=1e-5))
        @test result.accepted
        if result.solution!==nothing
            @test result.solution.objective≈parts.total atol=2e-3
            @test objective_components(c,result.solution).total≈result.solution.objective atol=1e-5
            @test replay_audit(c,result.solution)["valid"]
        end
    end
end

@testset "Valuation is additive across restart and currency normalization" begin
    c=terminal_value_case();u=ones(Int,2,2);q=[4.0 4.0;2.0 2.0]
    x=dispatch_from_controls(c,u,q,zeros(1,2))
    continued=restart_case(c,x,1.0)
    rest=dispatch_from_controls(continued,u[:,2:end],q[:,2:end],zeros(1,1))
    prefix=tv_copy(c;grid=[0.0,1.0],prices=[100.0])
    past=dispatch_from_controls(prefix,u[:,1:1],q[:,1:1],zeros(1,1))
    @test x.objective≈past.objective+rest.objective atol=1e-8
    @test objective_components(continued,rest).total≈rest.objective atol=1e-8
    @test objective_components(prefix,past).total≈past.objective atol=1e-8
    @test rest.storage≈x.storage[:,2:end] atol=1e-10
    @test replay_audit(continued,rest)["valid"]

    # Converting every monetary input by the same fixed rate preserves controls
    # and scales the complete objective; no exchange rate is inferred internally.
    rate=0.1
    system=tv_copy(c.system;
        reservoirs=[tv_copy(r;water_value=rate*r.water_value) for r in c.system.reservoirs],
        rivers=[tv_copy(r;water_value=rate*r.water_value) for r in c.system.rivers],
        generators=[tv_copy(g;startup=rate*g.startup,shutdown=rate*g.shutdown) for g in c.system.generators])
    converted=tv_copy(c;system,prices=rate.*c.prices)
    y=dispatch_from_controls(converted,u,q,zeros(1,2))
    @test y.storage==x.storage
    @test y.power==x.power
    @test y.objective≈rate*x.objective atol=1e-8
    original=objective_components(c,x);scaled=objective_components(converted,y)
    @test all(isapprox(scaled[k],rate*original[k];atol=1e-8) for k in keys(original))
    @test volume_water_value(rate*30.0,energy_conversion_factor(120.0,0.90);
        downstream_value=rate*6016.8)≈rate*14845.8 atol=1e-10

    # Startup and shutdown charges are one event each, even with zero output later.
    cycle=[1 0;1 0]
    stopped=dispatch_from_controls(c,cycle,q.*cycle,zeros(1,2))
    @test objective_components(c,stopped).transition_cost==19.0
    @test objective_components(c,stopped).total≈stopped.objective atol=1e-8
end

@testset "Distributed travel-time inventory has the same independent accounting" begin
    base=terminal_value_case()
    river=tv_copy(only(base.system.rivers);deterministic_delay=nothing,
        curves=[RiverRouting.DelayCurve(4.0,[0.0,1.0],[1.0])])
    c=tv_copy(base;system=tv_copy(base.system;rivers=[river]))
    x=dispatch_from_controls(c,ones(Int,2,2),[4.0 4.0;2.0 2.0],zeros(1,2))
    # Uniform travel times in [0,1] have mean .5 hour. Historical inflow is
    # 2 m³/s and final inflow is 4 m³/s, so the two inventories follow directly.
    @test only(x.terminal_transit)≈0.0036*4.0*0.5 atol=1e-12
    @test x.storage≈[2.0 1.9856 1.9712;2.0 2.0036 2.0108] atol=1e-10
    parts=objective_components(c,x)
    @test parts.transit_value≈6016.8*0.0036*(4.0-2.0)*0.5 atol=1e-9
    @test parts.total≈x.objective atol=1e-8
    @test validate(c,x)["valid"]
end

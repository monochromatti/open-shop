include(joinpath(@__DIR__,"..","examples","shared_tailrace.jl"))
using .SharedTailraceExample: shared_tailrace_case

@testset "Distinct penstocks with one shared quadratic tailrace" begin
    c=shared_tailrace_case()
    u=ones(Int,2,2);q=[5.0 5.0;3.0 3.0];gates=zeros(0,2)
    x=dispatch_from_controls(c,u,q,gates)
    @test x["validation"]["valid"]
    # This reference follows storage conservation and the three pressure laws
    # directly, without calling the model's constitutive functions.
    volumes=[100.0,99.9856,99.9568]
    heads=149.0 .+ 0.01 .* (volumes[1:end-1]+volumes[2:end])/2
    tail=0.03*8.0^2
    power=[0.00981*0.90*5*(h-0.08*5^2-tail) for h in heads]'
    power=vcat(power,[0.00981*0.92*3*(h-0.20*3^2-tail) for h in heads]')
    @test vec(x["V"])≈volumes atol=1e-9
    @test x["tunnel_q"]≈[q;8.0 8.0] atol=1e-8
    @test x["power"]≈power atol=1e-8
    expected=sum(diff(c.grid).*c.prices.*vec(sum(power;dims=1)))-5.0
    @test x["objective"]≈expected atol=1e-7
    @test replay_audit(c,x)["valid"]

    # Increasing B's flow reduces A's head even with A's flow held fixed.
    other=dispatch_from_controls(c,u,[5.0 5.0;4.0 4.0],gates)
    @test other["validation"]["valid"]
    drop=0.03*(9^2-8^2)
    # Include the additional midpoint drawdown caused by B's extra release.
    extra_drawdown=0.01 .* [0.0009,0.0036]
    @test vec(x["power"][1,:]-other["power"][1,:])≈
        0.00981*0.90*5 .* (drop .+extra_drawdown) atol=1e-8
    @test replay_audit(c,other)["valid"]

    restored=case_from_dict(case_dict(c))
    copy_dispatch=dispatch_from_controls(restored,u,q,gates)
    @test copy_dispatch["power"]≈power atol=1e-8

    # At optimum the 10 m³/s common tailrace capacity binds. With fixed total
    # flow, storage and common loss are fixed; differentiating the strictly
    # concave allocation objective gives this independent quadratic equation.
    optimal_volumes=[100.0,99.982,99.946]
    optimal_heads=149.0 .+0.01 .*
        (optimal_volumes[1:end-1]+optimal_volumes[2:end])/2
    a=3*0.92*0.20-3*0.90*0.08
    b=-6*0.92*0.20*10
    optimal_a=[(-b-sqrt(b^2-4*a*(3*0.92*0.20*100+
        (0.90-0.92)*(h-3))))/(2*a) for h in optimal_heads]
    optimal_b=10 .-optimal_a
    optimal_power=0.00981 .* (0.90 .*optimal_a .*
        (optimal_heads .-3 .-0.08 .*optimal_a.^2)+0.92 .*optimal_b .*
        (optimal_heads .-3 .-0.20 .*optimal_b.^2))
    optimum=sum(diff(c.grid).*c.prices.*optimal_power)-5.0
    @test optimum>expected+100.0

    local_result=solve_case(c;u,warm=x,time_limit=30.0)
    @test get(get(local_result,"validation",Dict()),"valid",false)
    if haskey(local_result,"objective")
        @test local_result["objective"]≈optimum atol=5e-3
        @test vec(local_result["generator_q"][1,:])≈optimal_a atol=0.02
        @test vec(local_result["generator_q"][2,:])≈optimal_b atol=0.02
        @test replay_audit(c,local_result)["valid"]
    end
    result=solve(c;initial=x,time_limit=60.0,relative_gap=1e-6)
    @test result.accepted
    @test result["global_bound"]===nothing || result["global_bound"]>=optimum-1e-5
    if result.solution!==nothing
        @test result.solution["objective"]≈optimum atol=5e-3
        @test result.solution["tunnel_q"][3,:]≈fill(10.0,2) atol=1e-4
        @test vec(result.solution["generator_q"][1,:])≈optimal_a atol=0.02
        @test replay_audit(c,result.solution)["valid"]
    end
end

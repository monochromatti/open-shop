include(joinpath(@__DIR__,"..","examples","unequal_turbine_grids.jl"))
using .UnequalTurbineGridsExample: unequal_grid_case

@testset "Unequal headwise grids preserve supplied knots and hand interpolation" begin
    for interpolation in (:bilinear,:pchip_discharge)
        table=only(unequal_grid_case(;interpolation).system.generators).turbine_table
        for (j,h) in enumerate(table.heads),(q,eta) in zip(table.curves[j].x,table.curves[j].y)
            @test OpenSHOP.turbine_efficiency(table,q,h)≈eta atol=1e-14
        end
        expected=if interpolation==:bilinear
            (0.875+(0.87+0.5*(0.94-0.87)/3))/2
        else
            # Original first-column slopes are [.08,0,-.03]. At q=2.5 its
            # first Hermite cell gives .891875. Second-column slopes are
            # [17/300,21/620,1/300]; q=2.5 is t=1/6 in its second cell.
            second=(87//100)*(25//27)+3*(21//620)*(25//216)+
                (94//100)*(2//27)+3*(1//300)*(-5//216)
            (0.891875+Float64(second))/2
        end
        @test OpenSHOP.turbine_efficiency(table,2.5,150.0)≈expected atol=1e-14
        @test OpenSHOP.turbine_efficiency(table,0.0,150.0;extrapolation=:linear)≈0.76 atol=1e-14
        @test OpenSHOP.turbine_efficiency(table,6.0,150.0;extrapolation=:linear)≈
            (0.87+(0.94+(0.94-0.87)/3))/2 atol=1e-14
        @test OpenSHOP.turbine_qmin(table,150.0)==1.0
        @test OpenSHOP.turbine_qmax(table,150.0)==5.0
        @test_throws DomainError OpenSHOP.turbine_efficiency(table,0.0,150.0)
    end
    original=only(unequal_grid_case().system.generators).turbine_table
    refitted=TurbineTable([100.0],[1.0,2.0,3.0,5.0],[0.80,0.87,0.90,0.88],
        [1.0],[5.0];interpolation=:pchip_discharge,head_extrapolation=:constant)
    @test OpenSHOP.turbine_efficiency(original,2.5,100.0)≈0.891875 atol=1e-14
    @test OpenSHOP.turbine_efficiency(refitted,2.5,100.0)≈0.89025 atol=1e-14
end

@testset "Unequal efficiency coverage is distinct from operating limits" begin
    curves=[TableCurve([1.0,3.0,5.0],[0.8,0.9,0.88]),
        TableCurve([2.0,4.0,6.0],[0.82,0.9,0.94])]
    for interpolation in (:bilinear,:pchip_discharge)
        table=TurbineTable([100.0,200.0],curves,[1.0,2.0],[5.0,6.0];interpolation)
        @test OpenSHOP.turbine_efficiency(table,1.0,100.0)==0.8
        @test OpenSHOP.turbine_efficiency(table,6.0,200.0)==0.94
        @test_throws DomainError OpenSHOP.turbine_efficiency(table,1.0,150.0)
        @test_throws DomainError OpenSHOP.turbine_efficiency(table,6.0,150.0)
        @test OpenSHOP.turbine_efficiency(table,1.0,150.0;extrapolation=:linear)≈0.79
        @test OpenSHOP.turbine_efficiency(table,6.0,150.0;extrapolation=:linear)≈0.905
        @test OpenSHOP.turbine_qmin(table,150.0)==1.5
        @test OpenSHOP.turbine_qmax(table,150.0)==5.5
        c=unequal_grid_case(;interpolation)
        g=OpenSHOP._river_replace(only(c.system.generators);turbine_table=table)
        c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;generators=[g]))
        @test_throws ArgumentError OpenSHOP.validate_inputs(c)
        admitted=OpenSHOP._river_replace(g;qmin=2.0)
        common=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;generators=[admitted]))
        @test OpenSHOP.validate_inputs(common)
        # A singleton operating head participates only in its own supplied curve.
        fixed=OpenSHOP._river_replace(g;hmin=100.0,hmax=100.0)
        @test OpenSHOP.validate_curve_data(OpenSHOP._river_replace(c.system;generators=[fixed]))
        # Below the first head, linear head extrapolation involves both columns.
        extended=TurbineTable([100.0,200.0],curves,[1.0,2.0],[5.0,6.0];
            interpolation,head_extrapolation=:linear)
        wide=OpenSHOP._river_replace(g;turbine_table=extended,hmin=50.0,hmax=90.0)
        @test_throws ArgumentError OpenSHOP.validate_curve_data(OpenSHOP._river_replace(c.system;generators=[wide]))
    end
    table=TurbineTable([100.0,200.0],curves,[1.0,2.0],[5.0,6.0])
    curves[1].x[2]=2.0;curves[1].y[1]=0.1
    @test table.curves[1].x==[1.0,3.0,5.0]
    @test OpenSHOP.turbine_efficiency(table,1.0,100.0)==0.8
end

@testset "Headwise turbine data round trip without synthesized calibration" begin
    c=unequal_grid_case()
    dictionary=case_dict(c)
    raw=dictionary["generators"][1]["turbine_table"]
    @test raw["curves"][1]["x"]==[1.0,3.0,5.0]
    @test raw["curves"][2]["x"]==[1.0,2.0,5.0]
    @test !haskey(raw,"discharge")
    @test case_dict(case_from_dict(dictionary))==dictionary
    mktempdir() do dir
        path=joinpath(dir,"unequal.json")
        writejson(path,dictionary)
        @test case_dict(readcase(path))==dictionary
    end
    ambiguous=deepcopy(dictionary)
    ambiguous["generators"][1]["turbine_table"]["discharge"]=[1.0,3.0,5.0]
    @test_throws ArgumentError case_from_dict(ambiguous)
    @test_throws ArgumentError TurbineTable(heads=[100.0],curves=[TableCurve([1.0,3.0],[0.8,0.9])],
        discharge=[1.0,3.0],efficiency=[0.8,0.9],qmin=[1.0],qmax=[3.0],head_extrapolation=:constant)
end

@testset "Global turbine equations retain original headwise polynomials" begin
    for interpolation in (:bilinear,:pchip_discharge)
        table=only(unequal_grid_case(;interpolation).system.generators).turbine_table
        second=Float64((87//100)*(25//27)+3*(21//620)*(25//216)+
            (94//100)*(2//27)+3*(1//300)*(-5//216))
        left=interpolation==:bilinear ? 0.875 : 0.891875
        right=interpolation==:bilinear ? 0.87+0.5*(0.94-0.87)/3 : second
        fourth=interpolation==:bilinear ? 0.87+2*(0.94-0.87)/3 :
            Float64((87//100)*(7//27)+3*(21//620)*(2//27)+
                (94//100)*(20//27)+3*(1//300)*(-4//27))
        m=Model();@variable(m,0.0<=q<=6.0);@variable(m,75.0<=h<=225.0)
        graph=OpenSHOP._global_tensor_turbine_graph!(m,table,q,h,0.0,6.0,75.0,225.0)
        for (qv,hv,expected) in ((2.5,100.0,left),(2.5,150.0,(left+right)/2),
            (2.5,225.0,left+1.25*(right-left)),(4.0,200.0,fourth),
            (0.0,150.0,0.76),(6.0,150.0,(0.87+0.94+(0.94-0.87)/3)/2))
            point=OpenSHOP._global_tensor_turbine_values(graph,qv,hv)
            point[q]=qv;point[h]=hv
            @test point[graph.eta]≈expected atol=1e-12
            @test isempty(primal_feasibility_report(m,point;atol=1e-9))
        end
        bounds=OpenSHOP.turbine_efficiency_bounds(table,0.0,6.0,75.0,225.0)
        @test all(bounds[1]-1e-12<=OpenSHOP.turbine_efficiency(table,q,h;extrapolation=:linear)<=bounds[2]+1e-12
            for q in 0.0:0.05:6.0,h in 75.0:2.5:225.0)
    end
end

@testset "Unequal grids retain forward, nonlinear and global schedule economics" begin
    for interpolation in (:bilinear,:pchip_discharge)
        c=unequal_grid_case(;interpolation)
        x=dispatch_from_controls(c,ones(Int,1,2),fill(5.0,1,2),zeros(0,2))
        heads=[149.9955,149.9865]
        eta=0.88 .+ (heads .- 100.0).*(0.94-0.88)/100
        expected=0.00981 .* 5.0 .* heads .* eta
        @test x.storage≈reshape([2.0,1.991,1.982],1,:) atol=1e-12
        @test vec(x.power)≈expected atol=1e-10
        @test x.objective≈sum(c.prices.*expected) atol=1e-8
        @test validate(c,x)["valid"]
        @test replay_audit(c,x)["valid"]
        next=restart_case(c,x,1.0)
        continued=dispatch_from_controls(next,ones(Int,1,1),fill(5.0,1,1),zeros(0,1))
        @test continued.storage[:,end]≈x.storage[:,end] atol=1e-12
        @test continued.power[:,1]≈x.power[:,2] atol=1e-10
        b=OpenSHOP._build_global_dispatch(c;joint=true)
        @test OpenSHOP._lift_start!(b,c,x)["valid"]
        for result in (solve_case(c;u=ones(Int,1,2),warm=x,time_limit=60.0),
            solve(c;initial=x,time_limit=60.0,relative_gap=1e-4))
            @test result.accepted
            if result.solution!==nothing
                @test result.solution.generator_q≈fill(5.0,1,2) atol=1e-4
                @test result.solution.objective≈x.objective atol=1e-3
                @test replay_audit(c,result.solution)["valid"]
            end
        end
    end
end

@testset "Zero unit minima use the same unequal-grid turbine equations" begin
    c=unequal_grid_case()
    table=TurbineTable(heads=[100.0,200.0],curves=[
        TableCurve([0.0,3.0,5.0],[0.8,0.9,0.88]),
        TableCurve([0.0,2.0,5.0],[0.82,0.87,0.94])],
        qmin=[0.0,0.0],qmax=[5.0,5.0],interpolation=:pchip_discharge)
    g=OpenSHOP._river_replace(only(c.system.generators);qmin=0.0,pmin=0.0,turbine_table=table)
    c=OpenSHOP._river_replace(c;system=OpenSHOP._river_replace(c.system;generators=[g]))
    @test OpenSHOP.validate_inputs(c)
    idle=dispatch_from_controls(c,ones(Int,1,2),zeros(1,2),zeros(0,2))
    @test validate(c,idle)["valid"]
    @test idle.power==zeros(1,2)
    for result in (solve_case(c;u=ones(Int,1,2),warm=idle,time_limit=60.0),
        solve(c;initial=idle,time_limit=60.0,relative_gap=1e-4))
        @test result.accepted
        if result.solution!==nothing
            @test result.solution.generator_q≈fill(5.0,1,2) atol=1e-4
            @test replay_audit(c,result.solution)["valid"]
        end
    end
end

@testset "Polynomial restriction selects the original cell on adjacent Float64 knots" begin
    a=nextfloat(1.0);b=nextfloat(a)
    table=TurbineTable([100.0,200.0],[
        TableCurve([1.0,a,b,2.0],[0.8,0.85,0.9,0.92]),
        TableCurve([1.0,1.5,2.0],[0.82,0.87,0.94])],
        [1.0,1.0],[2.0,2.0];interpolation=:pchip_discharge)
    polynomial=OpenSHOP._global_tensor_polynomial(table,a,b,100.0)
    @test polynomial[1]≈0.85 atol=1e-14
    @test sum(polynomial)≈0.9 atol=1e-14
    @test OpenSHOP.turbine_efficiency(table,a,100.0)==0.85
end

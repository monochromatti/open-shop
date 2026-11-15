using JuMP, Ipopt, ForwardDiff
const HC=OpenSHOP

@testset "Explicit table interpolation and domains" begin
    c=HC.TableCurve([0.0, 1.0, 3.0], [2.0, 4.0, 5.0])
    @test HC.table_value(c, 0.0)==2.0
    @test HC.table_value(c, 1.0)==4.0
    @test HC.table_value(c, 2.0)==4.5
    @test HC.table_value(c, 3.0)==5.0
    @test HC.table_slope(c, 1.0)==0.5
    @test_throws DomainError HC.table_value(c, -0.1)
    @test_throws DomainError HC.table_value(c, 3.1)
    @test_throws DomainError HC.table_value(c, NaN)
    @test HC.table_value(c, -1.0; extrapolation = :linear)==0.0
    @test HC.table_value(c, 4.0; extrapolation = :linear)==5.5
    @test_throws ArgumentError HC.table_value(c, 1.0; extrapolation = :clamp)
    @test_throws ArgumentError HC.TableCurve([0.0, 0.0], [1.0, 2.0])
    @test_throws ArgumentError HC.TableCurve([0.0, 1.0], [1.0, NaN])
    @test_throws ArgumentError HC.TableCurve([0.0], [1.0])
    @test_throws ArgumentError HC.TableCurve([0.0, 1.0], [1.0])
    @test_throws ArgumentError HC.TableCurve([1.0, 0.0], [1.0, 2.0])
end

@testset "Turbine family and implicit electrical output" begin
    t=HC.TurbineTable(
        [50.0, 100.0],
        [0.0, 10.0, 20.0],
        [0.7 0.72; 0.9 0.94; 0.8 0.84],
        [4.0, 5.0],
        [18.0, 20.0],
    )
    @test HC.turbine_efficiency(t, 10.0, 50.0)≈0.9
    @test HC.turbine_efficiency(t, 5.0, 75.0)≈0.815
    @test HC.turbine_qmin(t, 75.0)≈4.5
    @test HC.turbine_qmax(t, 75.0)≈19.0
    @test all(
        isapprox.(HC.turbine_efficiency_bounds(t, 0.0, 20.0, 50.0, 100.0), (0.7, 0.94)),
    )
    ext=HC.turbine_efficiency_bounds(t, 0.0, 20.0, 25.0, 120.0)
    @test all(
        ext[1]-1e-12<=HC.turbine_efficiency(t, q, h; extrapolation = :linear)<=ext[2]+1e-12
        for q in 0.0:0.5:20.0, h in 25.0:1.0:120.0
    )
    @test_throws DomainError HC.turbine_efficiency(t, 21.0, 75.0)
    @test_throws DomainError HC.turbine_qmax(t, 101.0)
    @test_throws ArgumentError HC.TurbineTable(
        [50.0, 100.0],
        [0.0, 20.0],
        fill(1.1, 2, 2),
        [4.0, 5.0],
        [18.0, 20.0],
    )
    @test_throws ArgumentError HC.TurbineTable(
        [50.0, 100.0],
        [0.0, 20.0],
        fill(0.9, 2, 2),
        [4.0, 5.0],
        [18.0, 21.0],
    )
    eg=HC.TableCurve([0.0, 10.0, 20.0], [0.9, 0.95, 0.96])
    @test HC.electrical_power(eg, 0.0)==0.0
    @test HC.electrical_power(eg, 10/0.95)≈10.0
    for p in [1.0, 7.0, 10.0, 14.0, 20.0]
        shaft=p/HC.table_value(eg, p)
        @test HC.electrical_power(eg, shaft)≈p
    end
    @test_throws DomainError HC.electrical_power(eg, 22.0)
    @test_throws DomainError HC.electrical_power(eg, -1.0)
    @test_throws ArgumentError HC.electrical_power(
        HC.TableCurve([1.0, 2.0], [0.9, 0.95]),
        1.0,
    )
    @test_throws ArgumentError HC.electrical_power(
        HC.TableCurve([0.0, 1.0, 2.0], [0.4, 0.4, 1.0]),
        1.0,
    )
end

@testset "JuMP table operators differentiate and solve" begin
    c=HC.TableCurve([0.0, 1.0, 3.0], [2.0, 4.0, 5.0])
    m=Model(Ipopt.Optimizer)
    set_silent(m)
    @variable(m, 0<=x<=3, start=1.5)
    f=HC.table_operator(m, c; name = :level_table_test)
    @constraint(m, f(x)==4.5)
    @objective(m, Min, x^2)
    optimize!(m)
    @test termination_status(m)==HC.MOI.LOCALLY_SOLVED
    @test value(x)≈2.0 atol=1e-6
    t=HC.TurbineTable(
        [50.0, 100.0],
        [0.0, 10.0, 20.0],
        [0.7 0.72; 0.9 0.94; 0.8 0.84],
        [4.0, 5.0],
        [18.0, 20.0],
    )
    m2=Model(Ipopt.Optimizer)
    set_silent(m2)
    @variable(m2, 1<=q<=9, start=4.0)
    @variable(m2, 60<=h<=90, start=75.0)
    eta=HC.turbine_operator(m2, t; name = :turbine_table_test)
    @constraint(m2, h==75.0)
    @constraint(m2, eta(q, h)==0.815)
    @objective(m2, Min, q^2)
    optimize!(m2)
    @test termination_status(m2)==HC.MOI.LOCALLY_SOLVED
    @test value(q)≈5.0 atol=1e-6
end

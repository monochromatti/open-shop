using Test
include("root_profile.jl")

@testset "Certified power envelopes" begin
    for interpolation in (:bilinear,:pchip_discharge)
        table=TurbineTable([50.,90.,140.],[2.,5.,10.,16.],
            [.70 .75 .79;.91 .95 .93;.85 .92 .96;.73 .81 .86],
            [2.,2.,2.],[16.,16.,16.];interpolation)
        for (qlo,qhi,hlo,hhi) in ((0.,19.,30.,170.),(3.,14.,60.,130.),
                (5.,5.+1e-8,89.,90.),(7.,7.,100.,100.))
            for slope in (-.5,0.,.5,1.5)
                intercept=power_support(table,qlo,qhi,hlo,hhi,.98,slope)
                residual=maximum(0.00981*q*h*.98*
                    turbine_efficiency(table,q,h;extrapolation=:linear)-slope*q-intercept
                    for q in range(qlo,qhi;length=101),h in range(hlo,hhi;length=37))
                @test residual<=1e-8
                @test isfinite(intercept)
            end
        end
    end
end

@testset "Native statistics and root capture" begin
    mktempdir() do output
        rows=root_profile(joinpath(@__DIR__,"cases","turbine-tables.json"),output;
            seconds=30.,profiles=["baseline"],capture=true)
        row=only(rows)
        @test row["accepted"]
        @test row["diagnostics_close_error"]===nothing
        @test row["solver_error"]===nothing
        @test haskey(row["scip_statistics"],"lp")
        @test row["first_root_lp"]!==nothing
        @test row["first_root_lp"]["fractional_commitment_max"]>=0
    end
end

# Independent high-precision primitive of the sum of uniform release time and
# each uniform delay bin. It is deliberately different from the Float64 kernel.
function reference_cohort_fraction(c,a,b,s,e)
    setprecision(256) do
        a,b,s,e=BigFloat.((a,b,s,e))
        primitive(x)=sum(BigFloat(c.weights[k])*(
            max(zero(x),x-BigFloat(c.edges[k]))^2-
            max(zero(x),x-BigFloat(c.edges[k+1]))^2)/
            (2*(BigFloat(c.edges[k+1])-BigFloat(c.edges[k]))) for k in eachindex(c.weights))
        Float64((primitive(b-s)-primitive(b-e)-primitive(a-s)+primitive(a-e))/(e-s))
    end
end

@testset "Cohort fractions remain physical on narrow and unequal intervals" begin
    curve=RiverRouting.DelayCurve(0.0,[0.1,0.5,3.0],[0.2,0.8])
    for (s,e) in ((0.8,nextfloat(0.8)),(0.8,0.9),(0.0,1.0),(-3.0,-1.1))
        times=sort!(unique!([s-1,s+0.1,s+0.5,s+1.0,e+3.0,e+4.0]))
        fractions=[RiverRouting.cohort_transfer_fraction(curve,a,b,s,e)
            for (a,b) in zip(times[1:end-1],times[2:end])]
        expected=[reference_cohort_fraction(curve,a,b,s,e)
            for (a,b) in zip(times[1:end-1],times[2:end])]
        @test all(0 .<= fractions .<= 1)
        @test fractions≈expected atol=1e-14
        @test sum(fractions)≈1.0 atol=1e-14
    end
    # Broad arrivals, tiny releases and the reverse occur at different places
    # in a network. Exercise support crossings and fully arrived cohorts too.
    for width in (eps(0.8),1e-12,0.01,1.0,10.0), delay in (0.01,0.4,8.0)
        c=RiverRouting.DelayCurve(0.0,[0.1,0.1+delay],[1.0])
        for (a,b) in ((0.0,0.5),(0.9,1.0),(0.8,2.3),(1.0,nextfloat(1.0)),(10.0,20.0))
            actual=RiverRouting.cohort_transfer_fraction(c,a,b,0.8,0.8+width)
            expected=reference_cohort_fraction(c,a,b,0.8,0.8+width)
            @test actual≈expected atol=1e-13
        end
    end
end

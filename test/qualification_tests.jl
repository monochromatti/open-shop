include("../scripts/qualify.jl")

@testset "Qualification measures loss only against qualified same-input solutions" begin
    rows=[Dict("label"=>"good","accepted"=>true,"qualified"=>true,"seconds"=>2.,"objective"=>-120.),
        Dict("label"=>"unfinished","accepted"=>false,"qualified"=>false,"seconds"=>3.,"objective"=>nothing)]
    reference=(;objective=-100.)
    comparison=Qualification.compare(rows,reference)
    @test comparison[1]["value_shortfall"]==20.
    @test comparison[1]["relative_shortfall"]≈0.2
    @test comparison[2]["value_shortfall"]===nothing
    @test all(x["relative_shortfall"]===nothing for x in Qualification.compare(rows,nothing))
end

@testset "Operating revisions retain the executed hydraulic state and histories" begin
    data=Qualification.input_data(Dict("seed"=>"17","family"=>"integrated",
        "size"=>"medium","hours"=>12,"branches"=>3,"pond_volume"=>0.0002))
    c=data.case
    x=Qualification.FuzzWorker.reconstruct(data)
    @test validate(c,x)["valid"]
    tail=restart_case(c,x,2.)
    for event in ("forecast","outage","restriction")
        next=Qualification.revise(tail,event)
        @test OpenSHOP.validate_inputs(next)
        @test first(next.grid)==2.
        @test [r.v0 for r in next.system.reservoirs]≈x.storage[:,3]
        @test [r.history_release for r in next.system.rivers]==
            [r.history_release for r in tail.system.rivers]
        @test [g.initial_age for g in next.system.generators]==
            [g.initial_age for g in tail.system.generators]
        if event=="forecast"
            for r in next.system.reservoirs
                @test OpenSHOP.opvalue(next,r.name,:inflow,2.,r.inflow)≈
                    0.7OpenSHOP.opvalue(tail,r.name,:inflow,2.,r.inflow)
            end
        else
            g=first(next.system.generators)
            attribute=event=="outage" ? :maintenance : :qmax
            default=event=="outage" ? 0. : g.qmax
            @test OpenSHOP.opvalue(next,g.name,attribute,2.,default)==default
            @test OpenSHOP.opvalue(next,g.name,attribute,4.,default)≈
                (event=="outage" ? 1. : max(g.qmin,0.6g.qmax))
            @test OpenSHOP.opvalue(next,g.name,attribute,8.,default)==default
        end
    end
end

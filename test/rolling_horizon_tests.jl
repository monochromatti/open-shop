include(joinpath(@__DIR__,"..","examples","rolling_horizon.jl"))

@testset "Integrated basin restart retains storage, transit and operating history" begin
    data=basin_case(;hours=8,branches=3,pond_volume=0.0002)
    c=data.case
    x=dispatch_from_controls(c,data.u,data.q,data.gate;tunnel_opening=data.opening)
    @test validate(c,x)["valid"]
    @test replay_audit(c,x)["valid"]
    time=2.0
    edge=searchsortedfirst(c.grid,time)
    continued=restart_case(c,x,time)
    rest=dispatch_from_controls(continued,data.u[:,edge:end],data.q[:,edge:end],
        data.gate[:,edge:end];tunnel_opening=data.opening[:,edge:end])
    @test validate(continued,rest)["valid"]
    @test rest.storage≈x.storage[:,edge:end] atol=1e-9
    @test rest.arrival_volume≈x.arrival_volume[:,edge:end] atol=1e-9
    @test rest.terminal_transit≈x.terminal_transit atol=1e-9
    @test all(g.initial_age==14.0 for g in continued.system.generators)

    grid=collect(time:1.0:(time+8.0))
    updated=ScheduleCase(;name="new forecast",system=continued.system,grid,
        prices=fill(75.0,8),operations=continued.operations,
        flow_requirements=continued.flow_requirements)
    guess=shifted_start(x,grid)
    y=dispatch_from_controls(updated,guess["u"],guess["generator_q"],guess["gate"];
        tunnel_opening=guess["tunnel_opening"])
    @test validate(updated,y)["valid"]
    @test replay_audit(updated,y)["valid"]
    @test y.storage[:,1]==x.storage[:,edge]
    old=simulate(c,data.q,data.gate;commitment=data.u,tunnel_opening=data.opening)
    new=simulate(updated,y.generator_q,y.gate;commitment=y.commitment,
        tunnel_opening=y.tunnel_opening)
    @test new.transit[:,1]≈old.transit[:,edge] atol=1e-10
    # Whole-basin conservation uses only input inflows and the sea-bound unit.
    inventory=vec(sum(new.storage;dims=1)+sum(new.transit;dims=1))
    expected=first(inventory)
    for t in eachindex(updated.prices)
        inflow=sum(z.values[searchsortedlast(z.times,grid[t])] for z in updated.operations)
        expected+=0.0036*(grid[t+1]-grid[t])*(inflow-y.generator_q[end,t])
        @test inventory[t+1]≈expected atol=1e-9
    end
end

@testset "River integrals reproduce independent cohort transport" begin
    grid=[0.0,0.2,0.8,1.5,2.0]
    base=graph_case(;grid,history=true)
    controls=[0.0 17.0 12.0 26.0;5.0 1.0 4.0 3.0;2.0 8.0 1.0 5.0]
    responses=[RiverRouting.DelayCurve(q,[0.1,0.5,3.0],w) for (q,w) in
        ((0.0,[0.2,0.8]),(10.0,[0.5,0.5]),(30.0,[0.8,0.2]))]
    for curves in (responses[1:1],responses[[1,3]],responses),
        mode in (:distributed,:mixed,:atomic),
        integral_grid in (grid,[0.0,2.0],OpenSHOP.refined_grid(grid))
        rivers=[begin
            atomic=mode==:atomic || mode==:mixed && i==1
            OpenSHOP._river_replace(r;
                curves=atomic ? RiverRouting.DelayCurve[] : curves,
                deterministic_delay=atomic ? 0.35 : nothing)
        end for (i,r) in enumerate(base.system.rivers)]
        c=OpenSHOP._river_replace(replace_graph(base;rivers);
            operations=[OperationalSeries(object=:MC,attribute=:capacity,
                times=[0.0,0.2],values=[10.0,30.0])])
        model=Model()
        @variable(model,q[1:3,1:4])
        set_start_value.(q,controls)
        routing=OpenSHOP.RiverTransport(model,c,q;integral_grid)
        expected=OpenSHOP.routed_schedule(c,controls)
        set_start_value.(q,expected.release)
        OpenSHOP.start_transport!(routing)
        OpenSHOP.bound_transport!(c,routing)
        @test value.(start_value,routing.arrivals)≈expected.arrival_volume atol=1e-12
        @test value.(start_value,routing.terminal)≈expected.transit[:,end] atol=1e-12
        for (i,r) in enumerate(rivers)
            cohort=expected.cohorts[i]
            for (a,b) in ((0.0,2.0),(0.13,1.63),(0.8,1.1))
                expression=OpenSHOP.routing_volume_expression(model,c,routing,i,a,b)
                volume=only(OpenSHOP.route_volumes(r,[a,b],cohort.grid,cohort.release))+
                    only(OpenSHOP.route_volumes(r,[a,b],r.history_grid,r.history_release))
                @test value(start_value,expression)≈volume atol=1e-12
            end
            for time in (0.0,0.13,grid[2]+0.35,0.8,1.31,2.0),side in (:left,:right)
                expression=OpenSHOP.routing_point_expression(model,c,routing,i,time;side)
                @test value(start_value,expression)≈
                    point_arrival(r,cohort.grid,cohort.release,time;side) atol=1e-10
            end
        end
        point=Dict(v=>start_value(v) for v in all_variables(model))
        @test isempty(primal_feasibility_report(model,point;atol=1e-9))
    end
end

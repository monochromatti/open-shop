function low_flow_delay_case()
    curves=[RiverRouting.DelayCurve(0.0,[2.0,3.0],[1.0]),
        RiverRouting.DelayCurve(10.0,[0.5,1.5],[1.0])]
    lake=Reservoir(;name=:Lake,z0=100.0,v0=2.0,vmin=1.0,vmax=3.0,water_value=1.0)
    receiver=Reservoir(;name=:Receiver,z0=20.0,v0=2.0,vmin=1.0,vmax=3.0,water_value=1.0)
    rivers=[
        River(;name=:Outlet,source=:Lake,target=:Merged,curves,capacity=10.0,
            law=:controlled,water_value=1.0,history_grid=[-6.0,0.0],history_release=[2.5]),
        River(;name=:Natural,target=:Merged,
            curves=[RiverRouting.DelayCurve(5.0,[0.5,1.5],[1.0])],capacity=10.0,
            inflow=0.25,water_value=1.0,history_grid=[-6.0,0.0],history_release=[0.25]),
        River(;name=:Merged,target=:Receiver,curves,capacity=10.0,inflow=0.5,
            water_value=1.0,history_grid=[-6.0,0.0],history_release=[3.25]),
    ]
    flows=[0.0,2.5,1.0,0.0]
    operations=[OperationalSeries(;object=:Outlet,attribute=attr,times=[0.0,1.0,2.0,3.0],
        values=flows ./ 10) for attr in (:gate_min,:gate_max)]
    ScheduleCase(;name="low_flow_confluence",system=HydroSystem(;reservoirs=[lake,receiver],rivers),
        grid=collect(0.0:4.0),prices=zeros(4),operations)
end

@testset "Strict delay coverage, including zero and below-reference flow" begin
    curves=[RiverRouting.DelayCurve(0.0,[2.0,3.0],[1.0]),
        RiverRouting.DelayCurve(10.0,[0.5,1.5],[1.0])]
    @test collect(RiverRouting.blend(curves,0.0))==[(1,1.0),(2,0.0)]
    @test collect(RiverRouting.blend(curves,2.5))==[(1,0.75),(2,0.25)]
    @test RiverRouting.route([0.0,2.0],[0.0,1.0],[0.0],curves)==[0.0]
    @test RiverRouting.remaining_volume(2.0,[0.0,1.0],[0.0],curves)==0.0
    # At two hours the slow bin has not arrived. The fast bin has delivered
    # 7/8 of its one-hour release. Its mixture weight at q=2.5 is 1/4.
    arrived=only(RiverRouting.route([0.0,2.0],[0.0,1.0],[2.5],curves))
    @test arrived≈0.0036*2.5*0.25*0.875 atol=1e-12
    @test RiverRouting.remaining_volume(2.0,[0.0,1.0],[2.5],curves)≈
        0.0036*2.5-arrived atol=1e-12
    @test RiverRouting.remaining_volume(0.0,[-2.0,-1.0],[2.5],curves)≈
        0.0036*2.5*(0.75+0.25*0.125) atol=1e-12
    @test sum(RiverRouting.route([0.0,1.0,2.0,5.0],[0.0,1.0,2.0],[0.0,2.5],curves))≈
        0.0036*2.5 atol=1e-12

    positive=[RiverRouting.DelayCurve(5.0,[1.0,2.0],[1.0]),
        RiverRouting.DelayCurve(10.0,[0.0,1.0],[1.0])]
    @test isnothing(RiverRouting.check_curves(positive))
    @test_throws ArgumentError RiverRouting.check_curves(positive;capacity=10.0)
    @test_throws ArgumentError RiverRouting.route([0.0,3.0],[0.0,1.0],[0.0],positive)
    @test_throws ArgumentError RiverRouting.route([0.0,3.0],[0.0,1.0],[2.5],positive)
    @test sum(RiverRouting.route([0.0,3.0],[0.0,1.0],[5.0],positive))≈0.0036*5
    fixed=[first(positive)]
    @test isnothing(RiverRouting.check_curves(fixed;capacity=20.0))
    for q in (0.0,0.5,5.0,20.0)
        @test collect(RiverRouting.blend(fixed,q))==[(1,1.0)]
        @test sum(RiverRouting.route([0.0,3.0],[-1.0,0.0],[q],fixed))+
            only(RiverRouting.route([-3.0,0.0],[-1.0,0.0],[q],fixed))≈0.0036*q atol=1e-12
    end
    # Bin weights are probabilities, not densities. Bin widths differ, so a
    # density interpretation would not give this mean delay of 2.375 hours.
    unequal=[RiverRouting.DelayCurve(5.0,[0.0,1.0,5.0],[0.25,0.75])]
    @test RiverRouting.remaining_volume(0.0,[-6.0,0.0],[0.5],unequal)≈
        0.0036*0.5*2.375 atol=1e-12
    @test_throws ArgumentError RiverRouting.DelayCurve(5.0,[0.0,1.0,4.0],[0.25,0.25])

    c=low_flow_delay_case()
    gates=zeros(3,4);gates[1,:]=[0.0,0.25,0.1,0.0]
    x=dispatch_from_controls(c,zeros(Int,0,4),zeros(0,4),gates)
    @test OpenSHOP.validate_inputs(c)
    @test x["validation"]["valid"]
    @test x.river_release[1,:]≈[0.0,2.5,1.0,0.0] atol=1e-10
    @test x.river_release[2,:]≈fill(0.25,4) atol=1e-10
    @test x.river_release[3,:]≈0.5 .+ vec(sum(x.arrival_volume[1:2,:];dims=1))./0.0036 atol=1e-10
    initial_transit=routing_data(c)["history_initial"]
    @test sum(x.storage[:,end])+sum(x.terminal_transit)≈
        sum(r.v0 for r in c.system.reservoirs)+sum(initial_transit)+0.0036*4*0.75 atol=1e-9
    @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
    restored=restart_case(c,x,2.0)
    rest=dispatch_from_controls(restored,zeros(Int,0,2),zeros(0,2),gates[:,3:4])
    @test rest["validation"]["valid"]
    @test rest.storage≈x.storage[:,3:end] atol=1e-9
    @test rest.terminal_transit≈x.terminal_transit atol=1e-10
    # Distributed confluences rebin arrivals. Report this grid-dependent error
    # separately instead of treating volume conservation as timing accuracy.
    transport=transport_audit(c,x;schedule_tolerance=0.01,absolute_tolerance=1e-5)
    @test transport["converged"]
    @test haskey(transport,"cumulative_schedule_discrepancy_Mm3")
    @test OpenSHOP._lift_start!(OpenSHOP._build_global_dispatch(c;joint=true),c,x)["valid"]
    local_result=solve_case(c;u=zeros(Int,0,4),warm=x,time_limit=10.0)
    @test local_result["accepted"]
    @test local_result["terminal_transit"]≈x.terminal_transit atol=1e-6

    # Full model admission is strict even when every selected release lies
    # inside a supplied positive subdomain. Raw routing above remains available.
    good_river=River(;name=:Reach,source=:Lake,target=:Receiver,curves=fixed,
        capacity=10.0,law=:controlled,water_value=1.0,
        history_grid=[-6.0,0.0],history_release=[5.0])
    good=ScheduleCase(;name="positive_reference",system=HydroSystem(
        reservoirs=c.system.reservoirs,rivers=[good_river]),grid=[0.0,1.0],prices=[0.0])
    witness=dispatch_from_controls(good,zeros(Int,0,1),zeros(0,1),fill(0.5,1,1))
    @test witness["validation"]["valid"]
    bad_river=OpenSHOP._river_replace(good_river;curves=positive)
    bad=OpenSHOP._river_replace(good;system=OpenSHOP._river_replace(good.system;rivers=[bad_river]))
    @test_throws ArgumentError OpenSHOP.validate_inputs(bad)
    @test_throws ArgumentError OpenSHOP._build_dispatch(bad;u=zeros(Int,0,1))
    @test_throws ArgumentError OpenSHOP._build_global_dispatch(bad;joint=true)
    @test_throws ArgumentError simulate(bad,zeros(0,1),fill(0.5,1,1))
    @test_throws ArgumentError dispatch_from_controls(bad,zeros(Int,0,1),zeros(0,1),fill(0.5,1,1))
    @test_throws ArgumentError forward_step(bad.system,[2.0,2.0],1.0,Float64[],[0.5],[0.0])
    @test_throws ArgumentError case_from_dict(case_dict(bad))
    @test !validate(bad,witness)["valid"]
    @test validate(bad,witness)["residuals"]["curve_input"]==Inf
end

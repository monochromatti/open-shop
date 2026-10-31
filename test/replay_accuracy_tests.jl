function passive_replay_case(scale; tabulated=false)
    lake=Reservoir(name=:Forebay,z0=1.8,slope=1000.0/scale,
        level_curve=tabulated ? TableCurve(scale.*[0.0001,0.0003],[1.9,2.1]) : nothing,
        v0=0.00022scale,vmin=0.0001scale,vmax=0.0003scale,
        inflow=40.0scale,water_value=0.0)
    system=HydroSystem(reservoirs=[lake],boundaries=[Boundary(name=:Tail,head=0.0)],
        tunnels=[Tunnel(name=:Outlet,source=:Forebay,target=:Tail,
            resistance=0.00125/scale^2,capacity=60.0scale)])
    ScheduleCase(name="physically similar passive reservoir",system=system,
        grid=[0.0,1.0],prices=[0.0])
end

@testset "Replay accuracy respects physical similarity and level changes" begin
    for scale in (0.01,1.0,100.0), tabulated in (false,true)
        c=passive_replay_case(scale;tabulated)
        x=dispatch_from_controls(c,zeros(Int,0,1),zeros(0,1),zeros(0,1))
        @test validate(c,x)["valid"]
        audit=replay_audit(c,x)
        @test audit["physically_valid"]
        @test !audit["numerically_converged"]
        @test !audit["valid"]
        # Each rescaling has the same head trajectory and capacity-relative
        # mesh discrepancy. Conservation cannot establish time accuracy.
        @test audit["max_level_difference_m"]≈0.030694129 rtol=1e-5
        @test audit["max_storage_difference_Mm3"]/(0.0002scale)≈0.153470645 rtol=1e-5
        head_only=replay_audit(c,x;storage_convergence=Inf,storage_relative_convergence=Inf)
        @test !head_only["valid"]
        @test replay_audit(c,x;storage_convergence=Inf,
            storage_relative_convergence=Inf,level_convergence=Inf)["valid"]
    end
    c=with_grid(passive_replay_case(1.0),OpenSHOP.refined_grid([0.0,1.0];factor=64))
    x=dispatch_from_controls(c,zeros(Int,0,64),zeros(0,64),zeros(0,64))
    @test replay_audit(c,x)["valid"]
    @test_throws ArgumentError replay_audit(c,x;storage_relative_convergence=-1.0)
    @test_throws ArgumentError replay_audit(c,x;level_convergence=NaN)
end

@testset "Flat levels do not conceal storage error in a small receiving pond" begin
    reservoirs=[Reservoir(name=:Source,z0=100.0,slope=0.0,v0=1.0,
        vmin=0.0,vmax=2.0,water_value=0.0),
        Reservoir(name=:Pond,z0=50.0,slope=0.0,v0=0.00005,
            vmin=0.0,vmax=0.0001,water_value=0.0)]
    curve=[RiverRouting.DelayCurve(0.0,[0.0,1.0],[1.0])]
    rivers=[River(name=:Tributary,source=:Source,target=:Main,law=:controlled,
        capacity=0.02,water_value=0.0,curves=curve,
        history_grid=[-2.0,0.0],history_release=[0.0]),
        River(name=:Main,source=:auto,target=:Pond,law=:junction,
            capacity=0.02,water_value=0.0,curves=curve,
            history_grid=[-2.0,0.0],history_release=[0.0])]
    c=ScheduleCase(name="small receiving pond",system=HydroSystem(;reservoirs,rivers),
        grid=[0.0,1.0],prices=[0.0])
    for mesh in (Float64[],[0.0,1.0])
        case=OpenSHOP._river_replace(c;routing_grid=mesh)
        x=dispatch_from_controls(case,zeros(Int,0,1),zeros(0,1),reshape([1.0,0.0],2,1))
        audit=replay_audit(case,x)
        @test validate(case,x)["valid"]
        # Two independent uniform [0,1] hour responses deliver one sixth
        # of the first hour's 72 m³ release by the one-hour horizon.
        @test abs(x.arrival_volume[2,1]-12e-6)>0.1e-6
        @test audit["physically_valid"]
        @test audit["max_level_difference_m"]==0.0
        @test !audit["numerically_converged"]
        @test !audit["valid"]
        @test replay_audit(case,x;storage_relative_convergence=Inf)["valid"]
    end
end

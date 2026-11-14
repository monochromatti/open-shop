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
        validation=validate_replay(c,x)
        @test validation["physically_valid"]
        @test !validation["numerically_converged"]
        @test !validation["valid"]
        # Each rescaling has the same head trajectory and capacity-relative
        # mesh discrepancy. Conservation cannot establish time accuracy.
        @test validation["max_level_difference_m"]≈0.030694129 rtol=1e-5
        @test validation["max_storage_difference_Mm3"]/(0.0002scale)≈0.153470645 rtol=1e-5
        head_only=validate_replay(c,x;storage_convergence=Inf,storage_relative_convergence=Inf)
        @test !head_only["valid"]
        @test validate_replay(c,x;storage_convergence=Inf,
            storage_relative_convergence=Inf,level_convergence=Inf)["valid"]
    end
    c=with_grid(passive_replay_case(1.0),OpenSHOP.refined_grid([0.0,1.0];factor=64))
    x=dispatch_from_controls(c,zeros(Int,0,64),zeros(0,64),zeros(0,64))
    @test validate_replay(c,x)["valid"]
    @test_throws ArgumentError validate_replay(c,x;storage_relative_convergence=-1.0)
    @test_throws ArgumentError validate_replay(c,x;level_convergence=NaN)
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
        validation=validate_replay(case,x)
        @test validate(case,x)["valid"]
        # Two independent uniform [0,1] hour responses deliver one sixth
        # of the first hour's 72 m³ release by the one-hour horizon.
        @test abs(x.arrival_volume[2,1]-12e-6)>0.1e-6
        @test validation["physically_valid"]
        @test validation["max_level_difference_m"]==0.0
        @test !validation["numerically_converged"]
        @test !validation["valid"]
        @test validate_replay(case,x;storage_relative_convergence=Inf)["valid"]
    end
end

@testset "Transport refinement solves a prescribed inflow on the original control grid" begin
    reservoirs=[Reservoir(name=:Source,z0=100.0,slope=0.01,v0=1.0,
        vmin=0.0,vmax=2.0,water_value=0.0),
        Reservoir(name=:Pond,z0=50.0,slope=0.1,v0=0.00005,
            vmin=0.0,vmax=0.0001,water_value=0.0)]
    curve=[RiverRouting.DelayCurve(0.0,[0.0,1.0],[1.0])]
    rivers=[River(name=:Tributary,source=:Source,target=:Main,law=:controlled,
        capacity=0.02,gate_min=1.0,water_value=0.0,curves=curve,
        history_grid=[-2.0,0.0],history_release=[0.0]),
        River(name=:Main,source=:auto,target=:Pond,law=:junction,
            capacity=0.02,water_value=0.0,curves=curve,
            history_grid=[-2.0,0.0],history_release=[0.0])]
    c=ScheduleCase(name="prescribed delayed inflow",system=HydroSystem(;reservoirs,rivers),
        grid=[0.0,1.0],prices=[0.0])
    result=solve_verified(c;u=zeros(Int,0,1),time_limit=15.0,max_refinements=1,
        feasibility_only=true)
    @test result["accepted"]
    if result["accepted"]
        actual,x=result["case"],result["solution"]
        @test x.grid==c.grid
        @test isempty(c.routing_grid)
        @test validate(actual,x)["valid"]
        @test validate_replay(actual,x)["valid"]
        @test x.arrival_volume[2,1]≈12e-6 atol=0.15e-6
        # Preparation's finite problem must survive a budget-exhausted public
        # solve and serialization, rather than reverting to the original mesh.
        elapsed=Ref(0.0)
        prepare=(case;kwargs...)->begin
            elapsed[]=31.0
            Dict("accepted"=>true,"solution"=>x,"problem"=>actual)
        end
        solved=OpenSHOP._solve(c,()->elapsed[],prepare;time_limit=30.0,
            preparation_time_limit=10.0)
        @test solved.accepted
        @test validate(solved.problem,solved.solution)["valid"]
        restored=SolveResult(OpenSHOP.jsonready(solved))
        @test validate(restored.problem,restored.solution)["valid"]
        native=OpenSHOP._solve(c,time,prepare;
            time_limit=30.0,preparation_time_limit=10.0)
        @test native.accepted
        @test validate(native.problem,native.solution)["valid"]
        @test native.global_certificate
    end
end

@testset "Transport refinement meets receiving-storage accuracy" begin
    for scale in (0.01,1.0,100.0), tributaries in (1,2)
        reservoirs=[Reservoir(name=:Source,z0=100.0,v0=scale,vmin=0.0,vmax=2scale,
            water_value=0.0),
            Reservoir(name=:Pond,z0=50.0,v0=0.00005scale,vmin=0.0,vmax=0.0001scale,
                water_value=0.0)]
        curve=[RiverRouting.DelayCurve(0.0,[0.0,1.0],[1.0])]
        flow=0.02scale
        rivers=[River(name=Symbol(:Tributary,i),source=:Source,target=:Main,
            law=:controlled,capacity=flow,water_value=0.0,curves=curve,
            history_grid=[-2.0,0.0],history_release=[0.0]) for i in 1:tributaries]
        push!(rivers,River(name=:Main,source=:auto,target=:Pond,law=:junction,
            capacity=tributaries*flow,water_value=0.0,curves=curve,
            history_grid=[-2.0,0.0],history_release=[0.0]))
        c=ScheduleCase(name="receiving-volume accuracy",system=HydroSystem(;reservoirs,rivers),
            grid=[0.0,1.0],prices=[0.0],routing_grid=[0.0,1.0])
        u=zeros(Int,0,1)
        q=zeros(0,1)
        gate=ones(length(rivers),1)
        gate[end,1]=0.0
        x=dispatch_from_controls(c,u,q,gate)
        released=0.0036tributaries*flow
        # Two uniform delays deliver one sixth of the release by hour one.
        # On the initial mesh the error is 1/12 of the release, but the raw
        # coarse/fine difference is only 1/16. A 7% tolerance separates them.
        tolerance=0.07released
        grid=OpenSHOP.refine_transport_grid(c,x,[Inf,tolerance],1.0)
        resolved=OpenSHOP._river_replace(c;routing_grid=grid)
        y=dispatch_from_controls(resolved,u,q,gate)
        @test abs(y.arrival_volume[end,1]-released/6)<=tolerance
        @test validate(resolved,y)["valid"]
    end
end

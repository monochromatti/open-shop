using JSON3

function operating_replay_case(; distributed=false)
    reservoir=Reservoir(name=:Lake,z0=100.0,slope=1.0,v0=5.0,
        vmin=1.0,vmax=10.0,inflow=2.0,water_value=0.0)
    curves=distributed ? [
        RiverRouting.DelayCurve(0.0,[0.0,0.5],[1.0]),
        RiverRouting.DelayCurve(10.0,[0.5,1.0],[1.0]),
    ] : RiverRouting.DelayCurve[]
    river=River(name=:Outlet,source=:Lake,target=:Sea,law=:controlled,
        inflow=1.0,capacity=10.0,water_value=0.0,curves=curves,
        deterministic_delay=distributed ? nothing : 0.5,
        history_grid=[-1.0,0.0],history_release=[0.0])
    system=HydroSystem(reservoirs=[reservoir],
        boundaries=[Boundary(name=:Sea,head=0.0)],rivers=[river])
    operations=[
        OperationalSeries(object=:Lake,attribute=:inflow,
            times=[0.0,1.0,2.0],values=[2.0,4.0,1.0]),
        OperationalSeries(object=:Outlet,attribute=:inflow,
            times=[0.0,1.0,2.0],values=[1.0,2.0,0.0]),
        OperationalSeries(object=:Outlet,attribute=:capacity,
            times=[0.0,1.0,2.0],values=[10.0,6.0,4.0]),
    ]
    ScheduleCase(name="operating hydraulic replay",system=system,grid=[0.0,1.0,2.0,3.0],
        prices=zeros(3),operations=operations)
end

@testset "Physical replay JSON retains public array names and shapes" begin
    c=operating_replay_case()
    trajectory=simulate(c;generator_q=zeros(0,3),gate=fill(0.5,1,3))
    step=forward_step(c.system,[5.0],1.0,Float64[],[0.5],[0.0])
    mktempdir() do directory
        path=joinpath(directory,"trajectory.json")
        writejson(path,trajectory)
        saved=JSON3.read(read(path,String))
        @test Set(String.(keys(saved))) == Set(("grid","V","H","tunnel_q",
            "tunnel_opening","generator_q","gate","power","river_release",
            "arrival_volume","transit","boundary_outflow","converged","residual","steps"))
        @test length(saved.V) == size(trajectory.storage,1)
        @test length(first(saved.V)) == length(c.grid)
        @test collect(first(saved.V)) == vec(trajectory.storage)
        @test length(saved.H) == size(trajectory.head,1)
        @test length(first(saved.H)) == length(c.prices)
        @test length(saved.steps) == length(c.prices)
        @test collect(saved.steps[1].Vnew) == trajectory.steps[1].storage
        @test collect(saved.steps[1].H) == trajectory.steps[1].head
        @test saved.steps[1].bound_violations.river_flow == 0.0
        step_path=joinpath(directory,"step.json")
        writejson(step_path,step)
        saved_step=JSON3.read(read(step_path,String))
        @test Set(String.(keys(saved_step))) == Set(("Vnew","H","tunnel_q",
            "tunnel_opening","river_release","river_injection","river_arrival",
            "generator_head","power","residual","converged","iterations",
            "boundary_net_inflow","routing_domain_valid","bound_violations","law_domain_valid"))
        @test collect(saved_step.Vnew) == step.storage
        @test collect(saved_step.H) == step.head
        @test saved_step.iterations == step.iterations
        bounds_path=joinpath(directory,"bounds.json")
        writejson(bounds_path,step.bound_violations)
        saved_bounds=JSON3.read(read(bounds_path,String))
        @test Set(String.(keys(saved_bounds))) == Set(("storage","tunnel_flow","river_flow"))
        @test saved_bounds.storage == step.bound_violations.storage
    end
end

@testset "Operating tunnel limits do not clip hydraulic flow" begin
    reservoir=Reservoir(name=:Lake,z0=100.0,slope=0.0,v0=5.0,
        vmin=1.0,vmax=10.0,water_value=0.0)
    tunnel=Tunnel(name=:Outlet,source=:Lake,target=:Sea,resistance=1.0,capacity=10.0)
    system=HydroSystem(reservoirs=[reservoir],tunnels=[tunnel],
        boundaries=[Boundary(name=:Sea,head=0.0)])
    operations=[
        OperationalSeries(object=:Lake,attribute=:inflow,
            times=[0.0,1.0,2.0],values=[2.0,1.0,0.0]),
        OperationalSeries(object=:Outlet,attribute=:capacity,
            times=[0.0,1.0,2.0],values=[10.0,5.0,0.0]),
        OperationalSeries(object=:Outlet,attribute=:opening,
            times=[0.0,1.0,2.0],values=[1.0,1.0,0.0]),
    ]
    c=ScheduleCase(name="operating tunnel capacity",system=system,operations=operations,
        grid=[0.0,1.0,2.0,3.0],prices=zeros(3))
    trajectory=simulate(c;generator_q=zeros(0,3),gate=zeros(0,3))
    @test trajectory.converged
    @test vec(trajectory.tunnel_q) ≈ [10.0,10.0,0.0] atol=1e-12
    @test trajectory.tunnel_q[1,3] == 0.0
    @test vec(diff(trajectory.storage;dims=2)) ≈ 0.0036 .* [-8.0,-9.0,0.0] atol=1e-10
    @test trajectory.boundary_outflow ≈ [10.0,10.0,0.0] atol=1e-12
    @test [step.bound_violations.tunnel_flow for step in trajectory.steps] ≈ [0.0,5.0,0.0]
    @test tunnel.capacity == 10.0 && tunnel.opening == 1.0
end

@testset "Typed replay retains operating inflows, capacities and conservation" begin
    for distributed in (false,true)
        c=operating_replay_case(;distributed)
        q=zeros(0,3)
        gate=reshape([0.5,0.25,1.0],1,:)
        trajectory=simulate(c;generator_q=q,gate)
        @test trajectory isa OpenSHOP.HydraulicTrajectory
        @test trajectory.steps isa Vector{OpenSHOP.HydraulicStep}
        @test trajectory.converged
        @test trajectory.storage === trajectory["V"]
        @test trajectory.head === trajectory["H"]
        @test @inferred(sum(trajectory.power)) == 0.0
        @test_throws KeyError trajectory["missing"]
        @test_throws MethodError setindex!(trajectory, zeros(1,4), "V")
        # The controlled withdrawal excludes each reach's natural inflow.
        withdrawals=[(10.0-1.0)*0.5,(6.0-2.0)*0.25,4.0]
        releases=withdrawals+[1.0,2.0,0.0]
        @test vec(trajectory.river_release) ≈ releases atol=1e-10
        @test vec(diff(trajectory.storage;dims=2)) ≈
            0.0036 .* ([2.0,4.0,1.0]-withdrawals) atol=1e-10
        for t in 1:3
            inventory_change=sum(trajectory.storage[:,t+1]-trajectory.storage[:,t])+
                sum(trajectory.transit[:,t+1]-trajectory.transit[:,t])
            @test inventory_change ≈ 0.0036*(
                [2.0,4.0,1.0][t]+[1.0,2.0,0.0][t]-trajectory.boundary_outflow[t]
            ) atol=2e-10
            @test trajectory.steps[t].bound_violations.river_flow == 0.0
        end
        if distributed
            # Operating limits do not replace the original calibration domain:
            # the q=4 cohort still interpolates against reference flows 0 and 10.
            river=only(c.system.rivers)
            expected=RiverRouting.route(c.grid,c.grid,releases,river.curves)
            @test vec(trajectory.arrival_volume) ≈ expected atol=1e-12
            @test last(river.curves).reference_flow == 10.0
            @test river.capacity == 10.0
        end
        shut=simulate(c;generator_q=q,gate=zeros(1,3))
        @test shut.converged
        @test vec(shut.river_release) == [1.0,2.0,0.0]
        @test shut.storage[:,end] ≈ [5.0+0.0036*7.0] atol=1e-10
        first_step=trajectory.steps[1]
        last_step=trajectory.steps[end]
        @test first_step.storage !== last_step.storage
        @test first_step.head !== last_step.head
        @test first_step.river_release !== last_step.river_release
        @test first_step.river_arrival !== last_step.river_arrival
        # Outputs saved in the trajectory and in each step have independent
        # ownership; changing one does not overwrite a later step or a column.
        saved=copy(trajectory.head)
        saved_last=copy(last_step.head)
        first_step.head[1]+=1.0
        @test trajectory.head == saved
        @test last_step.head == saved_last
    end
end

@testset "Operating capacity reports signed delta limits without clipping" begin
    curve=LevelFlowCurve([100.0],[TableCurve([-10.0,0.0,10.0],[-40.0,0.0,40.0])])
    for direction in (-1.0,1.0)
        reservoirs=[
            Reservoir(name=:Upper,z0=100.0,slope=1.0,v0=3.0+direction,
                vmin=1.0,vmax=6.0,water_value=0.0),
            Reservoir(name=:Lower,z0=100.0,slope=1.0,v0=3.0,
                vmin=1.0,vmax=6.0,water_value=0.0),
        ]
        river=River(name=:Channel,source=:Upper,target=:Lower,law=:delta,
            capacity=50.0,water_value=0.0,deterministic_delay=0.0,
            curves=RiverRouting.DelayCurve[],delta_level_ref_up_flow_curve=curve)
        system=HydroSystem(;reservoirs,rivers=[river])
        operations=[OperationalSeries(object=:Channel,attribute=:capacity,
            times=[0.0,1.0],values=[10.0,1.0])]
        c=ScheduleCase(name="operating delta capacity",system=system,operations=operations,
            grid=[0.0,1.0,2.0],prices=zeros(2))
        trajectory=simulate(c;generator_q=zeros(0,2),gate=zeros(1,2))
        q1=4direction/(1+0.0036*4)
        q2=4*(direction-2*0.0036*q1)/(1+0.0036*4)
        @test vec(trajectory.river_release) ≈ [q1,q2] atol=1e-9
        @test trajectory.steps[1].converged
        @test !trajectory.steps[2].routing_domain_valid
        @test !trajectory.converged
        @test trajectory.steps[2].bound_violations.river_flow ≈ abs(q2)-1 atol=1e-9
        @test trajectory.steps[2].residual <= 1e-9
        @test maximum(abs,trajectory.transit) < 1e-12
        @test vec(sum(trajectory.storage;dims=1)) ≈ fill(6.0+direction,3) atol=1e-10
        @test river.capacity == 50.0
    end
end

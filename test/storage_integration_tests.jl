function energy_cycle_case(;tabulated=false,efficiency=1.0,prescribed=false,concave=false)
    slope,curvature=concave ? (3.0,-0.2) : (1.0,0.2)
    lake=Reservoir(name=:Lake,z0=100.0,slope=slope,curvature=curvature,
        level_curve=tabulated ? TableCurve([0.5,1.0,2.0,3.0,3.5],
            [100+slope*v+curvature*v^2 for v in (0.5,1.0,2.0,3.0,3.5)]) : nothing,
        v0=1.0,vmin=0.5,vmax=3.5,water_value=0.0)
    units=[Generator(name=name,plant=plant,mode=mode,qmin=0.0,qmax=600.0,
        pmin=0.0,pmax=1000.0,hmin=90.0,hmax=120.0,hbest=100.0,qbest=300.0,
        efficiency=efficiency,qcurvature=0.0,hcurvature=0.0,initial_on=0,
        minup=0.0,mindown=0.0,startup=0.0) for (name,plant,mode) in
        ((:Turbine,:Generation,:generate),(:Pump,:Pumping,:pump))]
    system=HydroSystem(reservoirs=[lake],boundaries=[Boundary(name=:Sea,head=0.0)],
        plants=[Plant(name=:Generation,source=:Lake,target=:Sea,pmax=1000.0,ramp=nothing),
            Plant(name=:Pumping,source=:Sea,target=:Lake,pmax=1000.0,ramp=nothing)],generators=units)
    q=[0.0 1.0 1.0;2.0 0.0 0.0]/0.0036
    operations=prescribed ? [OperationalSeries(object=units[i].name,attribute=:discharge,
        times=[0.0,1.0,2.0],values=collect(q[i,:])) for i in 1:2] : OperationalSeries[]
    c=ScheduleCase(name="energy-consistent cycle",system=system,grid=[0.0,1.0,2.0,3.0],
        prices=ones(3),operations=operations)
    c,[0 1 1;1 0 0],q
end

@testset "Volume-averaged heads include crossings, reversals and zero changes" begin
    polynomial=Reservoir(name=:Polynomial,z0=100.0,slope=1.0,curvature=0.2,
        v0=1.0,vmin=0.0,vmax=5.0,water_value=0.0)
    table=TableCurve([0.0,1.0,2.0,4.0],[100.0,101.0,104.0,108.0])
    tabulated=Reservoir(name=:Table,level_curve=table,v0=1.0,vmin=0.0,vmax=4.0,water_value=0.0)
    @test mean_head(polynomial,1.0,3.0)≈100+2+0.2*13/3
    @test mean_head(tabulated,0.0,4.0)≈103.75
    @test mean_head(tabulated,0.5,3.0)≈103.15
    @test mean_head(tabulated,-1.0,0.0)≈99.5
    @test mean_head(tabulated,4.0,5.0)≈109.0
    for r in (polynomial,tabulated), (a,b) in ((0.5,3.0),(1.0,1.0),
        (1.0,nextfloat(1.0)),(prevfloat(1.0),nextfloat(1.0)))
        h=mean_head(r,a,b)
        @test isfinite(h)
        @test h≈mean_head(r,b,a) atol=1e-13
        @test min(head(r,a),head(r,b))<=h+1e-13
        @test h<=max(head(r,a),head(r,b))+1e-13
        a==b && @test h==head(r,a)
    end
    # Independent exact primitive for this particular piecewise-linear table.
    primitive(v)=v<=1 ? 100v+v^2/2 : v<=2 ?
        100.5+101*(v-1)+1.5*(v-1)^2 : 203+104*(v-2)+(v-2)^2
    for (a,b) in ((0.0,4.0),(0.5,3.0),(3.0,0.5),(1.0,2.0),(-1.0,5.0))
        @test (b-a)*mean_head(tabulated,a,b)≈primitive(b)-primitive(a) atol=1e-12
    end
    for (a,b) in ((0.2,3.2),(3.2,0.2),(0.3,0.6),(1.3,2.6))
        da,db=OpenSHOP._mean_head_partials(tabulated,a,b)
        delta=1e-5
        @test da≈(mean_head(tabulated,a+delta,b)-mean_head(tabulated,a-delta,b))/(2delta) atol=1e-8
        @test db≈(mean_head(tabulated,a,b+delta)-mean_head(tabulated,a,b-delta))/(2delta) atol=1e-8
    end
end

@testset "Connected passive hydraulics account for storage, electrical work and losses" begin
    reservoirs=[Reservoir(name=name,z0=z,slope=s,curvature=c,v0=v,
        vmin=0.1,vmax=5.0,water_value=0.0) for (name,z,s,c,v) in
        ((:A,100.0,1.0,0.2,3.0),(:B,98.0,2.0,0.1,1.0),(:C,95.0,1.0,0.05,1.0))]
    system=HydroSystem(reservoirs=reservoirs,junctions=[Junction(name=:J,hmin=90.0,hmax=110.0)],
        tunnels=[Tunnel(name=name,source=source,target=target,resistance=0.1,capacity=20.0)
            for (name,source,target) in ((:AJ,:A,:J),(:JB,:J,:B),(:BC,:B,:C))],
        rivers=[River(name=:Channel,source=:A,target=:C,law=:delta,
            curves=RiverRouting.DelayCurve[],deterministic_delay=0.0,capacity=20.0,water_value=0.0,
            delta_level_ref_up_flow_curve=LevelFlowCurve([0.0],
                [TableCurve([-20.0,0.0,20.0],[-10.0,0.0,10.0])]))],
        plants=[Plant(name=:Station,source=:B,target=:C,pmax=100.0,ramp=nothing)],
        generators=[Generator(name=:Unit,plant=:Station,qmin=0.0,qmax=10.0,
            pmin=0.0,pmax=100.0,hmin=0.0,hmax=20.0,qbest=3.0,hbest=5.0,
            efficiency=0.9,qcurvature=0.0,hcurvature=0.0,minup=0.0,mindown=0.0)])
    dt=0.1
    before=[r.v0 for r in reservoirs]
    step=forward_step(system,before,dt,[3.0],[0.0],[0.0])
    @test step.converged
    @test sum(step.storage)≈sum(before) atol=1e-10
    potential(r,v)=r.z0*v+r.slope*v^2/2+r.curvature*v^3/3
    change=(0.00981/0.0036)*sum(potential(r,step.storage[i])-potential(r,before[i])
        for (i,r) in enumerate(reservoirs))
    tunnel_loss=0.00981dt*sum(t.resistance*abs(step.tunnel_q[i])^3
        for (i,t) in enumerate(system.tunnels))
    channel_loss=0.00981dt*(step.head[1]-step.head[3])*only(step.river_release)
    hydraulic_work=only(step.power)*dt/0.9
    @test tunnel_loss>0 && channel_loss>0 && hydraulic_work>0
    @test change+tunnel_loss+channel_loss+hydraulic_work≈0.0 atol=1e-9
    @test change+only(step.power)*dt<0
end

@testset "Free discharge optimization uses the integrated reservoir work" begin
    for tabulated in (false,true), concave in (false,true)
        original,_,_=energy_cycle_case(;tabulated,concave)
        source=only(original.system.reservoirs)
        lake=Reservoir(name=:Lake,z0=source.z0,slope=source.slope,curvature=source.curvature,
            level_curve=source.level_curve,v0=3.0,vmin=1.0,vmax=3.5,water_value=0.0)
        unit=Generator(name=:Unit,plant=:Station,qmin=0.0,qmax=1/0.0036,
            pmin=0.0,pmax=1000.0,hmin=90.0,hmax=120.0,qbest=300.0,hbest=100.0,
            efficiency=1.0,qcurvature=0.0,hcurvature=0.0,initial_on=1,minup=0.0,
            mindown=0.0,startup=0.0)
        c=ScheduleCase(name="free discharge energy maximum",grid=[0.0,1.0],prices=[1.0],
            system=HydroSystem(reservoirs=[lake],boundaries=[Boundary(name=:Sea,head=0.0)],
                plants=[Plant(name=:Station,source=:Lake,target=:Sea,pmax=1000.0,ramp=nothing)],
                generators=[unit]))
        u=ones(Int,1,1)
        warm=dispatch_from_controls(c,u,fill(unit.qmax,1,1),zeros(0,1))
        # Positive elevation makes work strictly increasing with withdrawn volume;
        # the discharge cap therefore fixes the optimum at V=2 Mm³.
        work=tabulated ? (head(lake,2.0)+head(lake,3.0))/2 :
            100+lake.slope*2.5+lake.curvature*19/3
        expected=(0.00981/0.0036)*work
        for result in (solve_case(c;u,warm,time_limit=30.0),
            solve(c;initial=warm,fixed_u=u,time_limit=45.0))
            @test result.accepted
            @test result.objective≈expected atol=1e-4
            @test result.solution.storage[1,end]≈2.0 atol=1e-6
            @test result.global_bound===nothing || result.global_bound>=expected-1e-5
        end
    end
end

@testset "Tabulated and polynomial cycles account for energy and conversion losses" begin
    for tabulated in (false,true), concave in (false,true), eta in (1.0,0.9)
        c,u,q=energy_cycle_case(;tabulated,concave,efficiency=eta)
        coarse=dispatch_from_controls(c,u,q,zeros(0,3))
        @test validate(c,coarse)["valid"]
        @test coarse.storage≈[1.0 3.0 2.0 1.0] atol=1e-10
        # Independent potential integrals between 1 and 3 Mm³, in m*Mm³.
        lake=only(c.system.reservoirs)
        work=(0.00981/0.0036)*(200+4lake.slope+lake.curvature*(tabulated ? 9 : 26/3))
        expected=work*(eta-1/eta)
        @test sum(net_power(c,coarse))≈expected atol=1e-9
        for factor in (2,4)
            fine=simulate(c,q,zeros(0,3);commitment=u,
                grid=OpenSHOP.refined_grid(c.grid;factor))
            @test fine.converged
            @test fine.storage[1,end]≈1.0 atol=1e-10
            @test sum(net_power(c,fine).*diff(fine.grid))≈expected atol=1e-9
        end
    end
end

@testset "Both optimizers retain energy-consistent prescribed cycles" begin
    for tabulated in (false,true)
        c,u,q=energy_cycle_case(;tabulated,prescribed=true)
        warm=dispatch_from_controls(c,u,q,zeros(0,3))
        for result in (solve_case(c;u,warm,time_limit=30.0),
            solve(c;initial=warm,fixed_u=u,time_limit=45.0))
            @test result.accepted
            @test result.solution!==nothing
            if result.solution!==nothing
                @test validate(c,result.solution)["valid"]
                @test sum(net_power(c,result.solution))≈0.0 atol=1e-7
                @test result.solution.storage≈warm.storage atol=1e-8
            end
            if result["solver"]=="SCIP"
                @test result["start_audit"]["valid"]
                @test result.global_bound===nothing || result.global_bound>=-1e-6
            end
        end
    end
end

@testset "Zero storage change retains the physical head and a tight global bound" begin
    r=Reservoir(name=:Lake,level_curve=TableCurve([0.5,1.0,2.0,3.5],
        [100.55,101.2,102.8,105.95]),v0=1.0,vmin=0.5,vmax=3.5,inflow=10.0,water_value=0.0)
    unit=Generator(name=:Unit,plant=:Station,qmin=0.0,qmax=10.0,pmin=0.0,pmax=100.0,
        hmin=95.0,hmax=110.0,qbest=10.0,hbest=100.0,efficiency=1.0,
        qcurvature=0.0,hcurvature=0.0,initial_on=1,minup=0.0,mindown=0.0,startup=0.0)
    c=ScheduleCase(name="constant storage throughflow",grid=[0.0,1.0],prices=[1.0],
        system=HydroSystem(reservoirs=[r],boundaries=[Boundary(name=:Sea,head=0.0)],
            plants=[Plant(name=:Station,source=:Lake,target=:Sea,pmax=100.0,ramp=nothing)],
            generators=[unit]),
        operations=[OperationalSeries(object=:Unit,attribute=:discharge,times=[0.0],values=[10.0])])
    u=ones(Int,1,1)
    warm=dispatch_from_controls(c,u,fill(10.0,1,1),zeros(0,1))
    expected=0.00981*10*101.2
    for result in (solve_case(c;u,warm,time_limit=30.0),
        solve(c;initial=warm,fixed_u=u,time_limit=45.0))
        @test result.accepted
        @test result.solution.storage≈[1.0 1.0] atol=1e-10
        @test result.solution.head[1,1]≈101.2 atol=1e-8
        @test result.objective≈expected atol=1e-7
        if result["solver"]=="SCIP"
            @test result.global_certificate
            @test result.global_bound≈expected atol=1e-5
        end
    end
end

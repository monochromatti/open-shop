@testset "Revenue iterates retain independently reconstructed controls" begin
    lake=Reservoir(name=:Lake,v0=2.0,vmin=1.0,vmax=3.0,z0=100.0,
        slope=1.0,water_value=0.0)
    unit=Generator(name=:Unit,plant=:Station,qmin=1.0,qmax=5.0,
        pmin=0.1,pmax=5.0,qbest=3.0,hbest=100.0,hmin=90.0,hmax=110.0,
        efficiency=0.9,qcurvature=0.0,hcurvature=0.0,minup=0.0,mindown=0.0)
    system=HydroSystem(;reservoirs=[lake],generators=[unit],
        plants=[Plant(name=:Station,source=:Lake,target=:Sea,pmax=5.0)],
        boundaries=[Boundary(name=:Sea,head=0.0)])
    c=ScheduleCase(;name="retained revenue iterate",system,
        grid=[0.0,1.0,2.0],prices=[20.0,80.0])
    u=ones(Int,1,2)
    warm=dispatch_from_controls(c,u,fill(2.5,1,2),zeros(0,2))
    b=OpenSHOP._build_dispatch(c;u,warm)
    OpenSHOP._local_optimizer!(b.m,"Ipopt",30.0)
    # A short solve must still preserve the physical starting plan. Neither
    # the final native state nor a particular iteration count is an oracle.
    set_optimizer_attribute(b.m,"max_iter",1)
    retained=OpenSHOP.retain_feasible_dispatch!(b)
    optimize!(b.m)
    x=retained.solution[]
    @test x!==nothing
    if x!==nothing
        @test validate(c,x)["valid"]
        @test x.objective>=warm.objective-1e-8
        @test x.commitment==u
        @test vec(diff(x.storage;dims=2))≈-0.0036vec(x.generator_q) atol=1e-10
        @test validate_replay(c,x;grid=OpenSHOP.refined_grid(c.grid;factor=4))["valid"]
    end
end

@testset "Retained controls include spinning capacity" begin
    c=spinning_case()
    u=ones(Int,1,2)
    allocation=ReserveAllocation(fill(3.0,1,1,2),fill(1.0,1,1,2),
        fill(5.0+3/0.8829,1,2),fill(5.0-1/0.8829,1,2),zeros(1,2),zeros(1,2))
    warm=dispatch_from_controls(c,u,fill(5.0,1,2),zeros(0,2);reserves=allocation)
    @test validate(c,warm)["valid"]
    b=OpenSHOP._build_dispatch(c;u,warm)
    OpenSHOP._local_optimizer!(b.m,"Ipopt",30.0)
    set_optimizer_attribute(b.m,"max_iter",0)
    retained=OpenSHOP.retain_feasible_dispatch!(b)
    optimize!(b.m)
    x=retained.solution[]
    @test x!==nothing
    if x!==nothing
        @test validate(c,x)["valid"]
        @test OpenSHOP.validate_reserves(c,x)["valid"]
        @test x.reserves.up≈allocation.up atol=1e-8
        @test x.reserves.down≈allocation.down atol=1e-8
        @test x.objective≈warm.objective atol=1e-6
    end
end

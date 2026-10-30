@testset "A zero-power equality does not prevent dispatch of another unit" begin
    lake=Reservoir(name=:Lake,z0=100.0,slope=1.0,v0=2.0,vmin=1.9,vmax=3.0,water_value=0.0)
    common=(;qmin=0.0,qmax=5.0,pmin=0.0,pmax=10.0,qbest=3.0,
        efficiency=0.9,qcurvature=0.0,hcurvature=0.0,minup=0.0,mindown=0.0,
        initial_on=1,initial_age=24.0,startup=0.0)
    system=HydroSystem(reservoirs=[lake],
        boundaries=[Boundary(name=:FlatUp,head=100.0),Boundary(name=:FlatDown,head=100.0),
            Boundary(name=:Sea,head=0.0)],
        plants=[Plant(name=:Flat,source=:FlatUp,target=:FlatDown,pmax=10.0,ramp=nothing),
            Plant(name=:Productive,source=:Lake,target=:Sea,pmax=10.0,ramp=nothing)],
        generators=[Generator(;name=:FlatUnit,plant=:Flat,hmin=0.0,hmax=10.0,hbest=1.0,common...),
            Generator(;name=:Unit,plant=:Productive,hmin=90.0,hmax=110.0,hbest=100.0,common...)])
    c=ScheduleCase(name="bound-active conversion",system=system,grid=[0.0,1.0,2.0],prices=[50.0,70.0])
    u=ones(Int,2,2)
    # The flat-head unit must produce zero. In the other plant, the positive
    # head makes full discharge profitable and storage never binds.
    expected=dispatch_from_controls(c,u,[2.0 2.0;5.0 5.0],zeros(0,2))
    for result in (solve_case(c;u,time_limit=20.0),solve_verified(c;u,time_limit=20.0,max_refinements=0))
        @test result["accepted"]
        if result["accepted"]
            x=result["solution"]
            @test all(iszero,x.power[1,:])
            @test x.generator_q[2,:]≈[5.0,5.0] atol=1e-5
            @test all(q->0.0<=q<=5.0,x.generator_q)
            @test x.objective≈expected.objective atol=1e-3
            @test validate(c,x)["valid"]
            @test replay_audit(c,x)["valid"]
        end
    end
end

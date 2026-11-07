@testset "One dispatch respects capacities on both time meshes" begin
    lake=Reservoir(name=:Lake,z0=100.0,slope=0.1,v0=5.0,vmin=1.0,vmax=10.0,
        water_value=0.0)
    plant=Plant(name=:Station,source=:Lake,target=:Sea,pmax=10.0)
    unit=Generator(name=:Unit,plant=:Station,qmin=0.1,qmax=2.0,pmin=0.01,pmax=10.0,
        hmin=90.0,hmax=110.0,minup=0.0,mindown=0.0,
        qbest=1.0,hbest=100.0,efficiency=1.0,qcurvature=0.0,hcurvature=0.0)
    system=HydroSystem(reservoirs=[lake],boundaries=[Boundary(name=:Sea,head=0.0)],
        plants=[plant],generators=[unit])
    capacity=OperationalSeries(object=:Unit,attribute=:qmax,times=[0.0,1.0],values=[2.0,1.0])
    c=ScheduleCase(name="shared discharge capacity",system=system,grid=[0.0,1.0,2.0],
        prices=[100.0,100.0],operations=[capacity])
    result=solve_case(c;u=ones(Int,1,2),replay_grid=[0.0,0.5,1.0,1.5,2.0],time_limit=60.0)
    @test result.accepted
    if result.solution!==nothing
        @test vec(result.solution.generator_q) ≈ [2.0,1.0] atol=1e-4
        @test validate(c,result.solution)["valid"]
        @test replay_audit(c,result.solution)["valid"]
    end
end

@testset "A prescribed fully open outlet stays open on both meshes" begin
    lake=Reservoir(name=:Lake,z0=100.0,slope=0.1,v0=5.0,vmin=1.0,vmax=10.0,
        water_value=0.0)
    outlet=River(name=:Outlet,source=:Lake,target=:Sea,law=:controlled,capacity=10.0,
        gate_min=1.0,water_value=0.0,deterministic_delay=0.5,curves=RiverRouting.DelayCurve[],
        history_grid=[-1.0,0.0],history_release=[10.0])
    system=HydroSystem(reservoirs=[lake],boundaries=[Boundary(name=:Sea,head=0.0)],
        rivers=[outlet])
    c=ScheduleCase(name="shared fully open outlet",system=system,grid=[0.0,2.0],prices=[0.0])
    result=solve_case(c;u=zeros(Int,0,1),replay_grid=[0.0,1.0,2.0],time_limit=30.0)
    @test result.accepted
    if result.solution!==nothing
        @test result.solution.gate[1,1] == 1.0
        @test result.solution.river_release[1,1] ≈ 10.0 atol=1e-8
        @test result.solution.storage[1,end] ≈ 5.0-0.0036*2*10 atol=1e-9
        @test replay_audit(c,result.solution)["valid"]
    end
end

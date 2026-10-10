function commitment_control_case(; grid=[0.0,0.5,1.5,2.0,3.5],
    operations=OperationalSeries[], plant_minup=0.0, plant_mindown=0.0,
    plant_initial_on=nothing, plant_initial_age=nothing,
    unit_initial_on=0, unit_initial_age=10.0, unit_minup=0.0, unit_mindown=0.0)
    lake=Reservoir(name=:Lake,z0=100.0,slope=1.0,v0=2.0,vmin=1.0,vmax=3.0,water_value=0.0)
    plant=Plant(name=:Plant,source=:Lake,target=:Sea,pmax=10.0,
        minup=plant_minup,mindown=plant_mindown,
        initial_on=plant_initial_on,initial_age=plant_initial_age)
    units=[Generator(name=name,plant=:Plant,qmin=1.0,qmax=5.0,pmin=0.1,pmax=5.0,
        hmin=90.0,hmax=110.0,qbest=3.0,hbest=100.0,
        qcurvature=0.0,hcurvature=0.0,initial_on=unit_initial_on,
        initial_age=unit_initial_age,minup=unit_minup,mindown=unit_mindown)
        for name in (:UnitA,:UnitB)]
    system=HydroSystem(reservoirs=[lake],junctions=Junction[],
        boundaries=[Boundary(name=:Sea,head=0.0)],tunnels=Tunnel[],plants=[plant],
        generators=units,rivers=River[])
    ScheduleCase(name="commitment controls",system=system,grid=Float64.(grid),
        prices=zeros(length(grid)-1),operations=operations)
end

@testset "Hard commitment controls intersect rather than overwrite" begin
    maintenance=OperationalSeries(object=:UnitA,attribute=:maintenance,
        times=[0.0,1.5],values=[1.0,0.0])
    c=commitment_control_case(;operations=[maintenance])
    bounds=OpenSHOP._commitment_bounds(c)
    @test bounds.unit_upper==[0 0 1 1;1 1 1 1]
    @test admissible(c,[0 0 1 1;1 1 1 1])
    @test !admissible(c,ones(Int,2,4))
    plant_maintenance=OperationalSeries(object=:Plant,attribute=:maintenance,
        times=[0.5,2.0],values=[1.0,0.0])
    c=commitment_control_case(;operations=[plant_maintenance])
    @test OpenSHOP._commitment_bounds(c).unit_upper==[1 0 0 1;1 0 0 1]
    @test admissible(c,[1 0 0 1;0 0 0 1])
    forced=OperationalSeries(object=:UnitA,attribute=:forced_on,
        times=[0.0],values=[1.0])
    c=commitment_control_case(;operations=[maintenance,forced])
    @test_throws ArgumentError OpenSHOP._commitment_bounds(c)
    @test_throws ArgumentError OpenSHOP._joint_states!(JuMP.Model(),c)
    @test !admissible(c,ones(Int,2,4))
    outage=OperationalSeries(object=:UnitA,attribute=:forced_on,
        times=[0.0],values=[0.0])
    residual=commitment_control_case(;operations=[outage],unit_initial_on=1,
        unit_initial_age=0.0,unit_minup=2.0)
    @test_throws ArgumentError OpenSHOP._commitment_bounds(residual)
    @test_throws ArgumentError OpenSHOP._joint_states!(JuMP.Model(),residual)
    @test !admissible(residual,zeros(Int,2,4))
end

@testset "Scheduled power and discharge imply the same physical state" begin
    power=OperationalSeries(object=:UnitA,attribute=:power,times=[0.5,2.0],values=[1.0,0.0])
    discharge=OperationalSeries(object=:Plant,attribute=:discharge,
        times=[0.0,1.5],values=[0.0,3.0])
    c=commitment_control_case(;operations=[power])
    b=OpenSHOP._commitment_bounds(c)
    @test b.unit_lower[1,:]==[0,1,1,0]
    @test b.unit_upper[1,:]==[1,1,1,0]
    c=commitment_control_case(;operations=[discharge])
    b=OpenSHOP._commitment_bounds(c)
    @test b.unit_upper==[0 0 1 1;0 0 1 1]
    @test b.plant_lower==[0 0 1 1]
    @test !admissible(c,zeros(Int,2,4))
    @test admissible(c,[0 0 1 0;0 0 0 1])
    m=JuMP.Model(OpenSHOP.HiGHS.Optimizer)
    JuMP.set_silent(m)
    states=OpenSHOP._joint_states!(m,c)
    @test all(isnothing,states.plant_u) # No unnecessary dwell-state auxiliaries.
    JuMP.optimize!(m)
    @test JuMP.termination_status(m)==OpenSHOP.MOI.OPTIMAL
    @test admissible(c,round.(Int,JuMP.value.(states.u)))
    positive=OperationalSeries(object=:Plant,attribute=:power,times=[0.0],values=[2.0])
    closure=OperationalSeries(object=:Plant,attribute=:maintenance,times=[0.0],values=[1.0])
    c=commitment_control_case(;operations=[positive,closure])
    @test_throws ArgumentError OpenSHOP._commitment_bounds(c)
    # Positive power and zero discharge cannot be reconciled by commitment.
    zero=OperationalSeries(object=:UnitA,attribute=:discharge,times=[0.5],values=[0.0])
    c=commitment_control_case(;operations=[power,zero])
    @test_throws ArgumentError OpenSHOP._commitment_bounds(c)
end

@testset "Plant dwell acts on aggregate state and permits unit handover" begin
    c=commitment_control_case(;plant_minup=2.0,plant_mindown=1.0,
        plant_initial_on=0,plant_initial_age=10.0)
    handover=[1 0 1 0;0 1 0 0]
    @test admissible(c,handover)
    too_short=[1 0 0 0;0 0 0 0]
    @test !admissible(c,too_short)
    too_soon=[1 0 0 1;0 1 0 0]
    @test !admissible(c,too_soon) # Off at1.5, back on at2.0: only0.5h off.
    # Exhaust all schedules against independent switch-time duration checks,
    # then compare the actual linear MIP feasibility for each assignment.
    function expected(u)
        on=vec(any(u.==1;dims=1))
        switches=findall(t->on[t]!=(t==1 ? false : on[t-1]),eachindex(on))
        for k in 1:(length(switches)-1)
            t,next=switches[k],switches[k+1]
            c.grid[next]-c.grid[t] >= (on[t] ? 2.0 : 1.0) || return false
        end
        true
    end
    m=JuMP.Model(OpenSHOP.HiGHS.Optimizer)
    JuMP.set_silent(m)
    states=OpenSHOP._joint_states!(m,c)
    @test all(x->x isa JuMP.VariableRef,states.plant_u)
    for code in 0:255
        u=reshape([Int((code>>k)&1) for k in 0:7],2,4)
        accepted=expected(u)
        @test admissible(c,u)==accepted
        for j in 1:2,t in 1:4
            JuMP.fix(states.u[j,t],u[j,t];force=true)
        end
        JuMP.optimize!(m)
        status=JuMP.termination_status(m)
        @test status==(accepted ? OpenSHOP.MOI.OPTIMAL : OpenSHOP.MOI.INFEASIBLE)
    end
end

@testset "Plant initial and terminal histories remain explicit" begin
    missing=commitment_control_case(;plant_minup=2.0)
    @test_throws ArgumentError OpenSHOP._commitment_bounds(missing)
    inconsistent=commitment_control_case(;plant_initial_on=1,plant_initial_age=10.0)
    @test_throws ArgumentError OpenSHOP._commitment_bounds(inconsistent)
    c=commitment_control_case(;plant_minup=3.0,plant_mindown=2.0,
        plant_initial_on=0,plant_initial_age=10.0)
    u=[0 0 1 0;0 0 0 1]
    @test admissible(c,u) # The plant stays on across the handover.
    history=only(OpenSHOP._plant_commitment_history(c,u))
    @test history.state==1
    @test history.age==2.0
    @test history.residual==1.0
    at_edge=only(OpenSHOP._plant_commitment_history(c,u,2.0))
    @test at_edge.state==1 && at_edge.age==0.5 && at_edge.residual==2.5
    @test_throws ArgumentError OpenSHOP._plant_commitment_history(c,u,2.1)
    residual=commitment_control_case(;plant_minup=3.0,plant_initial_on=1,
        plant_initial_age=1.0,unit_initial_on=1,unit_initial_age=1.0)
    @test !admissible(residual,[0 0 0 0;0 0 0 0])
    @test admissible(residual,[1 0 0 0;0 1 1 0])
end

@testset "Legacy size and plant minimum-down constraints" begin
    legacy=commitment_control_case()
    m=JuMP.Model()
    states=OpenSHOP._joint_states!(m,legacy)
    @test JuMP.num_variables(m)==3*2*4
    @test all(isnothing,states.plant_u)
    c=commitment_control_case(;plant_mindown=1.0,plant_initial_on=1,
        plant_initial_age=10.0,unit_initial_on=1)
    @test !admissible(c,[0 1 1 1;0 0 0 0])
    @test admissible(c,[0 0 1 1;0 0 0 0])
    residual=commitment_control_case(;plant_mindown=1.0,plant_initial_on=0,
        plant_initial_age=0.0,unit_initial_age=0.0)
    @test OpenSHOP._commitment_bounds(residual).unit_upper==[0 0 1 1;0 0 1 1]
    @test !admissible(residual,[1 1 1 1;0 0 0 0])
    @test admissible(residual,[0 0 1 1;0 0 0 0])
    forced=OperationalSeries(object=:UnitA,attribute=:forced_on,times=[0.0],values=[1.0])
    conflict=commitment_control_case(;plant_mindown=1.0,plant_initial_on=0,
        plant_initial_age=0.0,unit_initial_age=0.0,operations=[forced])
    @test_throws ArgumentError OpenSHOP._joint_states!(JuMP.Model(),conflict)
end

@testset "Aggregate initial age respects known unit histories" begin
    on=commitment_control_case(;unit_initial_on=1,unit_initial_age=10.0)
    @test OpenSHOP._commitment_bounds(on).plant_initial_age==[10.0]
    too_young=commitment_control_case(;unit_initial_on=1,unit_initial_age=10.0,
        plant_initial_age=5.0)
    @test_throws ArgumentError OpenSHOP._commitment_bounds(too_young)
    handover=commitment_control_case(;unit_initial_on=1,unit_initial_age=10.0,
        plant_initial_age=12.0)
    @test OpenSHOP._commitment_bounds(handover).plant_initial_age==[12.0]
    off=commitment_control_case(;unit_initial_age=10.0)
    units=copy(off.system.generators)
    units[2]=OpenSHOP._river_replace(units[2];initial_age=4.0)
    off=OpenSHOP._river_replace(off;
        system=OpenSHOP._river_replace(off.system;generators=units))
    @test OpenSHOP._commitment_bounds(off).plant_initial_age==[4.0]
    inconsistent=OpenSHOP._river_replace(off;
        system=OpenSHOP._river_replace(off.system;plants=[OpenSHOP._river_replace(
            only(off.system.plants);initial_age=10.0)]))
    @test_throws ArgumentError OpenSHOP._commitment_bounds(inconsistent)
    continuation=OpenSHOP._river_replace(off;
        system=OpenSHOP._river_replace(off.system;plants=[OpenSHOP._river_replace(
            only(off.system.plants);initial_age=4.0,mindown=5.0)]))
    @test OpenSHOP._commitment_bounds(continuation).unit_upper==[0 0 1 1;0 0 1 1]
end

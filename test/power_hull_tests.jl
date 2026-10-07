using Test, JuMP, OpenSHOP, HiGHS

function power_hull_test_residual(m,assigned)
    maximum(all_constraints(m;include_variable_in_set_constraints=true);init=0.0) do ref
        obj=constraint_object(ref)
        x=obj.func isa Number ? obj.func : JuMP.value(v->assigned[v],obj.func)
        set=obj.set
        set isa JuMP.MOI.EqualTo ? abs(x-set.value) :
            set isa JuMP.MOI.LessThan ? max(0.0,x-set.upper) :
            set isa JuMP.MOI.GreaterThan ? max(0.0,set.lower-x) :
            set isa JuMP.MOI.Interval ? max(0.0,set.lower-x,x-set.upper) :
            set isa JuMP.MOI.ZeroOne ? max(abs(x-round(x)),max(0.0,-x,x-1)) :
            error("unsupported product-hull test set")
    end
end

@testset "Conditional on-state efficiency ranges" begin
    analytic=Generator(name=:Unit,plant=:Plant,qmin=1.0,qmax=20.0,pmin=0.01,pmax=100.0,hmin=1.0,hmax=200.0,efficiency=0.95,qbest=10.0,hbest=100.0,
        qcurvature=0.2,hcurvature=0.1)
    lo,hi=OpenSHOP._power_hull_eta_bounds(analytic,(8.0,12.0),(95.0,105.0))
    @test lo<0.94175 && hi>0.95
    @test lo≈0.94175 atol=2e-10
    @test hi≈0.95 atol=2e-10
    # Negative curvature remains supported by the analytic extrema helper.
    convex=Generator(name=:Unit,plant=:Plant,qmin=1.0,qmax=20.0,pmin=0.01,pmax=100.0,hmin=1.0,hmax=200.0,efficiency=0.8,qbest=10.0,hbest=100.0,
        qcurvature=-0.2,hcurvature=0.0)
    lo,hi=OpenSHOP._power_hull_eta_bounds(convex,(8.0,12.0),(100.0,100.0))
    @test lo<0.8 && hi>0.808
    for interpolation in (:bilinear,:pchip_discharge)
        table=TurbineTable([50.0,90.0,140.0],[2.0,5.0,10.0,16.0],
            [0.70 0.75 0.79;0.91 0.95 0.93;0.85 0.92 0.96;0.73 0.81 0.86],
            fill(2.0,3),fill(16.0,3);interpolation)
        g=Generator(name=:Unit,plant=:Plant,qmin=1.0,qmax=20.0,pmin=0.01,pmax=100.0,hmin=1.0,hmax=200.0,turbine_table=table)
        for (qbox,hbox) in (((0.0,20.0),(20.0,180.0)),((3.1,14.3),(60.0,130.0)),
            ((5.0,5.0+1e-8),(89.0,90.0)),((7.0,7.0),(100.0,100.0)),
            ((10.0,10.0),(20.0,180.0)))
            lo,hi=OpenSHOP._power_hull_eta_bounds(g,qbox,hbox)
            sampled=[OpenSHOP.turbine_efficiency(table,q,h;extrapolation=:linear)
                for q in range(qbox...;length=41),h in range(hbox...;length=19)]
            @test lo<=minimum(sampled)<=maximum(sampled)<=hi
            @test hi-lo<maximum(sampled)-minimum(sampled)+0.001
        end
        qbox=(3.1,14.3);hbox=(60.0,130.0)
        hbox,etabox=OpenSHOP._power_hull_on_boxes(g,qbox,hbox,(0.0,1.0),0.01,1.0)
        m=Model();@variable(m,0<=q<=20);@variable(m,-20<=h<=200)
        @variable(m,-2<=eta<=2);@variable(m,0<=p<=100);@variable(m,u,Bin)
        record=OpenSHOP._power_hull_unit!(m,q,h,eta,p,u,qbox,hbox,etabox,(-20.0,200.0),
            (-2.0,2.0),1.0,1.0;name=:table_tight)
        m.ext[:global_power_hulls]=Dict("table_tight"=>record)
        for (uv,qv,hv) in ((1.0,3.1,60.0),(1.0,7.25,101.0),(1.0,14.3,130.0),
            (0.0,0.0,-20.0))
            ev=OpenSHOP.turbine_efficiency(table,qv,hv;extrapolation=:linear)
            assigned=Dict(q=>qv,h=>hv,eta=>ev,p=>0.00981*qv*hv*ev,u=>uv)
            OpenSHOP._lift_power_hulls!(m,assigned)
            @test length(assigned)==num_variables(m)
            @test power_hull_test_residual(m,assigned)<1e-11
        end
    end
    # Head extrapolation can defeat the original PCHIP monotonicity: the
    # derivative extrema, rather than only efficiency values at q knots, matter.
    table=TurbineTable([10.0,20.0],[1.0,2.0,3.0,4.0],
        [0.6 0.9;0.7 0.8;0.8 0.9;0.9 0.8],fill(1.0,2),fill(4.0,2);
        interpolation=:pchip_discharge)
    g=Generator(name=:Unit,plant=:Plant,qmin=1.0,qmax=20.0,pmin=0.01,pmax=100.0,hmin=1.0,hmax=200.0,turbine_table=table)
    lo,hi=OpenSHOP._power_hull_eta_bounds(g,(2.0,3.0),(30.0,30.0))
    values=[OpenSHOP.turbine_efficiency(table,q,30.0;extrapolation=:linear)
        for q in range(2.0,3.0;length=301)]
    @test lo<=minimum(values) && maximum(values)<=hi
    @test hi>maximum(OpenSHOP.turbine_efficiency(table,q,30.0;extrapolation=:linear)
        for q in (2.0,3.0))+0.001
end

@testset "Tight boxes preserve exact on/off lifts" begin
    g=Generator(name=:Unit,plant=:Plant,qmin=1.0,qmax=2.0,pmin=0.05,pmax=1.0,
        hmin=0.1,hmax=10.0,min_efficiency=0.0,efficiency=0.9,qcurvature=0.0,hcurvature=0.0)
    qbox=(1.0,2.0);hbox=(0.1,10.0);etabox=(0.0,1.0)
    tighthead,tighteta=OpenSHOP._power_hull_on_boxes(g,qbox,hbox,etabox,0.05,1.0)
    limitinghead=0.05/(0.00981*2*0.9)
    @test hbox[1]<tighthead[1]<limitinghead
    @test tighthead[1]≈limitinghead rtol=1e-9
    @test tighteta[1]<0.9<tighteta[2]
    @test tighteta[2]-tighteta[1]<3e-10
    @test OpenSHOP._power_hull_on_boxes(g,qbox,hbox,etabox,0.0,1.0)[1]==hbox
    @test OpenSHOP._power_hull_on_boxes(g,qbox,hbox,etabox,0.05,0.0)[1]==hbox
    m=Model();@variable(m,0<=q<=2);@variable(m,-20<=h<=10)
    eta=@variable(m,lower_bound=-2,upper_bound=1.2,base_name="eta_1_1")
    @variable(m,0<=p<=1);@variable(m,u,Bin)
    system=HydroSystem(reservoirs=Reservoir[],junctions=Junction[],boundaries=Boundary[],
        tunnels=Tunnel[],rivers=River[],plants=[Plant(name=:Plant,source=:Source,target=:Target,pmax=1.0)],
        generators=[g])
    c=ScheduleCase(name="tight_hull",system=system,grid=[0.0,1.0],prices=[1.0])
    b=(m=m,u=reshape([u],1,1),GQ=reshape([q],1,1),P=reshape([p],1,1),shared_heads=Dict((:Plant,1)=>h))
    profile=OpenSHOP._add_power_hull!(b,c)
    @test profile["variables_added"]==10
    @test only(values(m.ext[:global_power_hulls])).axes[2][1]==tighthead[1]
    for (uv,qv,hv,ev) in ((1.0,2.0,limitinghead,0.9),(1.0,1.5,5.0,0.9),
        (0.0,0.0,-20.0,-2.0),(0.0,0.0,-1.0,1.2))
        assigned=Dict(q=>qv,h=>hv,eta=>ev,p=>0.00981*qv*hv*ev,u=>uv)
        OpenSHOP._lift_power_hulls!(m,assigned)
        @test length(assigned)==num_variables(m)
        @test power_hull_test_residual(m,assigned)<1e-12
    end
end

@testset "Joint power hull contains physical on/off points" begin
    for boxes in (
        ((1.0,2.0),(2.0,5.0),(0.3,0.9)),
        ((1.0,1.0),(2.0,5.0),(0.3,0.9)),
        ((1.0,2.0),(3.0,3.0),(0.7,0.7)),
        ((1.0,1.0),(3.0,3.0),(0.7,0.7)),
        ((0.0,-0.0),(0.0,-0.0),(0.0,-0.0)))
        qbox,hbox,etabox=boxes
        m=Model();@variable(m,0<=q<=qbox[2]);@variable(m,-20<=h<=10)
        @variable(m,-2<=eta<=1.2);@variable(m,0<=p<=1);@variable(m,u,Bin)
        before=count(is_binary,all_variables(m))
        record=OpenSHOP._power_hull_unit!(m,q,h,eta,p,u,qbox,hbox,etabox,(-20.0,10.0),(-2.0,1.2),
            0.8,0.95;name=:test)
        m.ext[:global_power_hulls]=Dict("test"=>record)
        @test count(is_binary,all_variables(m))==before
        @test length(record.weights)==prod(box[1]==box[2] ? 1 : 2 for box in boxes)
        @test num_variables(m)-5<=10
        for qv in unique([qbox[1],sum(qbox)/2,qbox[2]]),
            hv in unique([hbox[1],sum(hbox)/2,hbox[2]]),
            ev in unique([etabox[1],sum(etabox)/2,etabox[2]]),electrical in (0.8,0.9,0.95)
            assigned=Dict(q=>qv,h=>hv,eta=>ev,p=>0.00981*qv*hv*ev*electrical,u=>1.0)
            OpenSHOP._lift_power_hulls!(m,assigned)
            @test length(assigned)==num_variables(m)
            @test power_hull_test_residual(m,assigned)<1e-12
            @test sum(assigned[w]*prod(corner) for (w,corner) in zip(record.weights,record.corners))≈qv*hv*ev
        end
        # Full off-inclusive product bounds preserve heads/eta outside the on box.
        for hv in (-20.0,0.0,10.0),ev in (-2.0,0.0,1.2)
            assigned=Dict(q=>0.0,h=>hv,eta=>ev,p=>0.0,u=>0.0)
            OpenSHOP._lift_power_hulls!(m,assigned)
            @test power_hull_test_residual(m,assigned)==0.0
            @test all(w->assigned[w]==0.0,record.weights)
        end
    end
end

function power_hull_test_mccormick!(m,z,x,y,lx,ux,ly,uy)
    @constraint(m,z>=lx*y+ly*x-lx*ly)
    @constraint(m,z>=ux*y+uy*x-ux*uy)
    @constraint(m,z<=ux*y+ly*x-ux*ly)
    @constraint(m,z<=lx*y+uy*x-lx*uy)
end

@testset "Complete hull rejects sequential-product optimism" begin
    # Asymmetric boxes and differently positioned factors create a strict gap.
    # Actual q*h*eta=3.0; recursive McCormick permits4.11; joint hull permits3.705.
    m=Model(HiGHS.Optimizer);set_silent(m)
    @variable(m,1<=q<=2);@variable(m,2<=h<=5);@variable(m,0.3<=eta<=0.9)
    @constraint(m,q==1.25);@constraint(m,h==3.2);@constraint(m,eta==0.75)
    @variable(m,0.3<=flow<=1.8);@variable(m,0<=p<=1)
    power_hull_test_mccormick!(m,flow,q,eta,1.0,2.0,0.3,0.9)
    power_hull_test_mccormick!(m,p/0.00981,flow,h,0.3,1.8,2.0,5.0)
    @objective(m,Max,p)
    optimize!(m)
    @test termination_status(m)==JuMP.MOI.OPTIMAL
    @test objective_value(m)/0.00981≈4.11 atol=1e-10
    record=OpenSHOP._power_hull_unit!(m,q,h,eta,p,1.0,(1.0,2.0),(2.0,5.0),(0.3,0.9),
        (2.0,5.0),(0.3,0.9),1.0,1.0;name=:joint)
    optimize!(m)
    @test termination_status(m)==JuMP.MOI.OPTIMAL
    @test objective_value(m)/0.00981≈3.705 atol=1e-10
    @test objective_value(m)>=0.00981*1.25*3.2*0.75
    @test length(record.weights)==8
end

function power_hull_test_fixture(;off=false,negative_electrical=false,empty_box=false)
    m=Model();@variable(m,0<=q<=2);@variable(m,-20<=h<=5)
    eta=@variable(m,lower_bound=-2,upper_bound=0.9,base_name="eta_1_1")
    @variable(m,0<=p<=1);@variable(m,u,Bin)
    off && fix(u,0;force=true)
    curve=negative_electrical ? TableCurve([0.0,0.1],[0.9,0.8]) : nothing
    generator=Generator(name=:Unit,plant=:Plant,qmin=1.0,qmax=2.0,pmin=0.01,pmax=1.0,
        hmin=empty_box ? 10.0 : 2.0,hmax=20.0,min_efficiency=0.3,efficiency=0.75,qcurvature=0.0,hcurvature=0.0,generator_efficiency_curve=curve)
    system=HydroSystem(reservoirs=Reservoir[],junctions=Junction[],boundaries=Boundary[],
        tunnels=Tunnel[],rivers=River[],plants=[Plant(name=:Plant,source=:Source,target=:Target,pmax=1.0)],
        generators=[generator])
    c=ScheduleCase(name="hull_fixture",system=system,grid=[0.0,1.0],prices=[1.0])
    b=(m=m,u=reshape([u],1,1),GQ=reshape([q],1,1),P=reshape([p],1,1),shared_heads=Dict((:Plant,1)=>h))
    b,c,(q=q,h=h,eta=eta,p=p,u=u)
end

@testset "Power hull domain coverage and complete start lifting" begin
    b,c,v=power_hull_test_fixture()
    profile=OpenSHOP._add_power_hull!(b,c)
    @test profile["units_added"]==1
    @test profile["variables_added"]==10
    assigned=Dict(v.q=>1.25,v.h=>3.2,v.eta=>0.75,v.p=>0.00981*3.0,v.u=>1.0)
    OpenSHOP._lift_power_hulls!(b.m,assigned)
    @test length(assigned)==num_variables(b.m)
    @test power_hull_test_residual(b.m,assigned)<1e-12
    @test_throws ArgumentError OpenSHOP._add_power_hull!(b,c)
    for (kwargs,reason) in (((off=true,),"fixed_off"),((negative_electrical=true,),"electrical"),((empty_box=true,),"empty_on_box"))
        b,c,_=power_hull_test_fixture(;kwargs...)
        profile=OpenSHOP._add_power_hull!(b,c)
        @test profile["units_added"]==0
        @test profile["skipped"][reason]==1
    end
end

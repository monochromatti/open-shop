using Test, JuMP, OpenSHOP, HiGHS
isdefined(@__MODULE__,:ProductHull) || include(joinpath(@__DIR__,"product_hull.jl"))
using .ProductHull

function hull_residual(m,assigned)
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

@testset "Joint power hull contains physical on/off points" begin
    for lower_power in (false,true),boxes in (
        ((1.0,2.0),(2.0,5.0),(0.3,0.9)),
        ((1.0,1.0),(2.0,5.0),(0.3,0.9)),
        ((1.0,2.0),(3.0,3.0),(0.7,0.7)),
        ((1.0,1.0),(3.0,3.0),(0.7,0.7)),
        ((0.0,-0.0),(0.0,-0.0),(0.0,-0.0)))
        qbox,hbox,etabox=boxes
        m=Model();@variable(m,0<=q<=qbox[2]);@variable(m,-20<=h<=10)
        @variable(m,-2<=eta<=1.2);@variable(m,0<=p<=1);@variable(m,u,Bin)
        before=count(is_binary,all_variables(m))
        record=ProductHull._add_unit!(m,q,h,eta,p,u,qbox,hbox,etabox,(-20.0,10.0),(-2.0,1.2),
            0.8,0.95;name=:test,lower_power)
        m.ext[:product_hull]=Dict("test"=>record)
        @test count(is_binary,all_variables(m))==before
        @test length(record.weights)==prod(box[1]==box[2] ? 1 : 2 for box in boxes)
        @test num_variables(m)-5<=10
        for qv in unique([qbox[1],sum(qbox)/2,qbox[2]]),
            hv in unique([hbox[1],sum(hbox)/2,hbox[2]]),
            ev in unique([etabox[1],sum(etabox)/2,etabox[2]]),electrical in (0.8,0.9,0.95)
            assigned=Dict(q=>qv,h=>hv,eta=>ev,p=>0.00981*qv*hv*ev*electrical,u=>1.0)
            lift_product_hull!(m,assigned)
            @test length(assigned)==num_variables(m)
            @test hull_residual(m,assigned)<1e-12
            @test sum(assigned[w]*prod(corner) for (w,corner) in zip(record.weights,record.corners))≈qv*hv*ev
        end
        # Full off-inclusive product bounds preserve heads/eta outside the on box.
        for hv in (-20.0,0.0,10.0),ev in (-2.0,0.0,1.2)
            assigned=Dict(q=>0.0,h=>hv,eta=>ev,p=>0.0,u=>0.0)
            lift_product_hull!(m,assigned)
            @test hull_residual(m,assigned)==0.0
            @test all(w->assigned[w]==0.0,record.weights)
        end
    end
end

function mccormick!(m,z,x,y,lx,ux,ly,uy)
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
    mccormick!(m,flow,q,eta,1.0,2.0,0.3,0.9)
    mccormick!(m,p/0.00981,flow,h,0.3,1.8,2.0,5.0)
    @objective(m,Max,p)
    optimize!(m)
    @test termination_status(m)==JuMP.MOI.OPTIMAL
    @test objective_value(m)/0.00981≈4.11 atol=1e-10
    record=ProductHull._add_unit!(m,q,h,eta,p,1.0,(1.0,2.0),(2.0,5.0),(0.3,0.9),
        (2.0,5.0),(0.3,0.9),1.0,1.0;name=:joint)
    optimize!(m)
    @test termination_status(m)==JuMP.MOI.OPTIMAL
    @test objective_value(m)/0.00981≈3.705 atol=1e-10
    @test objective_value(m)>=0.00981*1.25*3.2*0.75
    @test length(record.weights)==8
end

function wrapper_fixture(;off=false,negative_electrical=false,empty_box=false)
    m=Model();@variable(m,0<=q<=2);@variable(m,-20<=h<=5)
    eta=@variable(m,lower_bound=-2,upper_bound=0.9,base_name="eta_1_1")
    @variable(m,0<=p<=1);@variable(m,u,Bin)
    off && fix(u,0;force=true)
    curve=negative_electrical ? TableCurve([0.0,0.1],[0.9,0.8]) : nothing
    generator=Generator(name=:Unit,plant=:Plant,qmin=1.0,qmax=2.0,pmin=0.01,pmax=1.0,
        hmin=empty_box ? 10.0 : 2.0,hmax=20.0,min_efficiency=0.3,generator_efficiency_curve=curve)
    system=HydroSystem(reservoirs=Reservoir[],junctions=Junction[],boundaries=Boundary[],
        tunnels=Tunnel[],rivers=River[],plants=[Plant(name=:Plant,source=:Source,target=:Target,pmax=1.0)],
        generators=[generator])
    c=ScheduleCase(name="hull_fixture",system=system,grid=[0.0,1.0],prices=[1.0])
    b=(m=m,u=reshape([u],1,1),GQ=reshape([q],1,1),P=reshape([p],1,1),shared_heads=Dict((:Plant,1)=>h))
    b,c,(q=q,h=h,eta=eta,p=p,u=u)
end

@testset "Power hull coverage and automatic start extension" begin
    b,c,v=wrapper_fixture()
    profile=add_product_hull!(b,c)
    @test profile["units_added"]==1
    @test profile["variables_added"]==10
    assigned=Dict(v.q=>1.25,v.h=>3.2,v.eta=>0.75,v.p=>0.00981*3.0,v.u=>1.0)
    foreach(lift!->lift!(assigned),b.m.ext[:experiment_start_lifters])
    @test length(assigned)==num_variables(b.m)
    @test hull_residual(b.m,assigned)<1e-12
    @test_throws ArgumentError add_product_hull!(b,c)
    for (kwargs,reason) in (((off=true,),"fixed_off"),((negative_electrical=true,),"electrical"),((empty_box=true,),"empty_on_box"))
        b,c,_=wrapper_fixture(;kwargs...)
        profile=add_product_hull!(b,c)
        @test profile["units_added"]==0
        @test profile["skipped"][reason]==1
        @test !haskey(b.m.ext,:experiment_start_lifters)
    end
    b,c,_=wrapper_fixture()
    @test add_product_hull!(b,c;max_units=0)["skipped"]["size"]==1
    b,c,_=wrapper_fixture()
    @test add_product_hull!(b,c;selected=["absent"])["units_added"]==0
end

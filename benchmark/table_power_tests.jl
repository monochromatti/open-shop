using Test, JuMP, OpenSHOP, HiGHS
isdefined(@__MODULE__,:TablePower) || include(joinpath(@__DIR__,"table_power.jl"))
using .TablePower

@testset "Table-coordinate certificates cover interpolation and extensions" begin
    support,evaluations,hits,seconds=TablePower._support_cache()
    for interpolation in (:bilinear,:pchip_discharge)
        table=TurbineTable([1.0,2.0,3.0],[1.0,1.5,2.0],
            [.7 .5 .4;.95 .8 .6;.8 .7 .4],fill(1.0,3),fill(2.0,3);interpolation)
        qnodes=OpenSHOP._global_tensor_nodes(table.discharge,0.0,3.0)
        hnodes=OpenSHOP._global_tensor_nodes(table.heads,-1.0,5.0)
        for (qbox,hbox) in (((1.2,1.8),(1.4,2.6)),((0.25,2.5),(0.7,3.4)),
            ((1.3,1.3),(2.2,2.2)))
            a=.00981*.95*hbox[2]*.5;b=.00981*.95*qbox[2]*.75
            hc=TablePower._head_coefficients(table,qbox,hbox,hnodes,.95,a,support)
            qc=TablePower._discharge_coefficients(table,qbox,hbox,qnodes,.95,b,support)
            head_error=-Inf;discharge_error=-Inf
            for q in range(qbox...;length=17),h in range(hbox...;length=19)
                power=.00981*.95*q*h*OpenSHOP.turbine_efficiency(table,q,h;extrapolation=:linear)
                mu=OpenSHOP._global_tensor_coordinate_weights(hnodes,h)
                lambda=OpenSHOP._global_tensor_coordinate_weights(qnodes,q)
                head_error=max(head_error,power-a*q-sum(hc.*mu))
                discharge_error=max(discharge_error,power-b*h-sum(qc.*lambda))
            end
            @test head_error<=1e-10
            @test discharge_error<=1e-10
            @test all(isfinite,hc)&&all(isfinite,qc)
        end
        old_evaluations=evaluations[];old_hits=hits[]
        first=support(table,1.2,1.8,1.4,2.6,.95,0.1,0.2)
        @test support(table,1.2,1.8,1.4,2.6,.95,0.1,0.2)==first
        @test evaluations[]==old_evaluations+1 && hits[]==old_hits+1
    end
    @test seconds[]>=0
    # Decreasing head efficiency makes h*eta concave: raw nodal supports fail.
    table=TurbineTable([1.0,3.0],[1.0,2.0],[.9 .3;.9 .3],fill(1.0,2),fill(2.0,2))
    raw=[support(table,1.0,1.0,h,h,1.0,0.0) for h in (1.0,3.0)]
    exact=.00981*1.0*2.0*.6
    @test exact>sum(raw)/2+1e-3
    corrected=TablePower._head_coefficients(table,(1.0,1.0),(1.0,3.0),[1.0,3.0],1.0,0.0,support)
    @test sum(corrected)/2>=exact
end

function table_power_gate_residual(m,assigned)
    maximum(all_constraints(m;include_variable_in_set_constraints=true);init=0.0) do ref
        obj=constraint_object(ref)
        x=obj.func isa Number ? obj.func : JuMP.value(v->assigned[v],obj.func)
        s=obj.set
        s isa JuMP.MOI.EqualTo ? abs(x-s.value) : s isa JuMP.MOI.LessThan ? max(0.0,x-s.upper) :
            s isa JuMP.MOI.GreaterThan ? max(0.0,s.lower-x) :
            s isa JuMP.MOI.ZeroOne ? abs(x-round(x)) : error("unexpected gate test set")
    end
end

@testset "Signed intercept gates and constant coordinates" begin
    for coefficients in ([-2.0,3.0],[-2.0],[-2.0,-2.0])
        m=Model();@variable(m,u,Bin)
        weights=length(coefficients)==1 ? [1.0] : @variable(m,[1:2],lower_bound=0,upper_bound=1)
        length(coefficients)>1 && @constraint(m,sum(weights)==1)
        records=Any[]
        z=TablePower._intercept!(m,coefficients,weights,u,true,"test",records)
        @test length(records)==(coefficients==[-2.0,3.0] ? 1 : 0)
        for uv in (0.0,1.0),t in (0.0,.25,1.0)
            assigned=Dict(u=>uv)
            if length(weights)>1
                assigned[weights[1]]=1-t;assigned[weights[2]]=t
            end
            TablePower._lift_gates!(records,assigned)
            expression=sum(coefficients[j]*(weights[j] isa Number ? weights[j] : assigned[weights[j]]) for j in eachindex(weights))
            @test TablePower._value(z,assigned)≈uv*expression
            @test table_power_gate_residual(m,assigned)<1e-12
            @test length(assigned)==num_variables(m)
        end
    end
end

function table_power_fixture(;interpolation=:pchip_discharge,analytic=false,negative_electrical=false,positive_electrical=false,operating=false)
    table=analytic ? nothing : TurbineTable([1.0,50.0,100.0],[5.0,10.0,20.0],
        [.15 .2 .25;.85 .9 .95;.75 .8 .85],fill(5.0,3),fill(20.0,3);
        interpolation,head_extrapolation=:linear)
    electrical=negative_electrical ? TableCurve([0.0,1.0],[.9,.8]) :
        positive_electrical ? TableCurve([0.0,2.0,3.0,5.0,25.0],[.8,.85,.97,.86,.8]) : nothing
    generators=[Generator(name=name,plant=plant,qmin=5.0,qmax=20.0,
        pmin=i==1 ? 1.0 : 1e-5,pmax=25.0,efficiency=.9,qbest=10.0,qcurvature=i==1 ? 0.0 : 2.0,
        hbest=50.0,hcurvature=0.0,hmin=i==1 ? 1.0 : .001,hmax=105.0,
        initial_on=0,minup=0.0,mindown=0.0,startup=0.0,turbine_table=table,generator_efficiency_curve=electrical)
        for (i,(name,plant)) in enumerate(((:OnUnit,:OnPlant),(:OffUnit,:OffPlant)))]
    system=HydroSystem(reservoirs=[Reservoir(name=:Lake,z0=40.0,slope=20.0,v0=2.001,
        vmin=1.0,vmax=4.0,inflow=6.0,water_value=10.0)],junctions=Junction[],
        boundaries=[Boundary(name=:LowerTail,head=20.0),Boundary(name=:UpperTail,head=80.0)],
        tunnels=Tunnel[],rivers=River[],plants=[Plant(name=:OnPlant,source=:Lake,target=:LowerTail,pmax=25.0),
        Plant(name=:OffPlant,source=:Lake,target=:UpperTail,pmax=25.0)],generators=generators)
    operations=operating ? [OperationalSeries(object=:OnUnit,attribute=attribute,times=[0.0],values=[value])
        for (attribute,value) in ((:qmin,7.0),(:qmax,9.0),(:pmin,2.0),(:pmax,4.0))] : OperationalSeries[]
    ScheduleCase(name="table_power",system=system,grid=[0.0,1.0],prices=[100.0],operations=operations)
end

@testset "Coordinate power supports retain complete physical starts" begin
    for interpolation in (:bilinear,:pchip_discharge)
        c=table_power_fixture(;interpolation)
        seed=dispatch_from_controls(c,reshape([1,0],2,1),reshape([8.0,0.0],2,1),zeros(0,1))
        @test seed["validation"]["valid"]
        for (axes,gate) in ((:head,false),(:head,true),(:discharge,false),(:discharge,true),(:both,false),(:both,true)),fixed in (false,true)
            b=OpenSHOP._build_global_dispatch(c;joint=true,fixed_u=fixed ? seed["u"] : nothing)
            before=count(is_binary,all_variables(b.m))
            profile=add_table_power!(b,c;axes,gate)
            @test profile["units_added"]==(fixed ? 1 : 2)
            @test count(is_binary,all_variables(b.m))==before
            @test profile["variables_added"]<=profile["head_rows"]+profile["discharge_rows"]
            !gate && @test profile["variables_added"]==0
            @test 0<=profile["certificate_seconds"]<=profile["setup_seconds"]
            audit=OpenSHOP._lift_start!(b,c,seed)
            @test audit["valid"]
            @test audit["assigned"]==audit["variables"]==num_variables(b.m)
            @test audit["objective"]≈seed["objective"] atol=1e-6
            if !fixed
                off=b.m.ext[:global_power_hulls]["power_hull_2_1"]
                @test start_value(off.head)<0 && start_value(off.eta)<0
            end
        end
    end
    c=table_power_fixture(;analytic=true)
    b=OpenSHOP._build_global_dispatch(c;joint=true)
    @test add_table_power!(b,c;axes=:both,gate=true)["skipped"]["analytic"]==2
    @test_throws ArgumentError add_table_power!(b,c)
    @test_throws ArgumentError add_table_power!(b,c;axes=:wrong)
    c=table_power_fixture(;negative_electrical=true)
    @test_throws ArgumentError OpenSHOP._build_global_dispatch(c;joint=true)
    # Defensive helper coverage for malformed curves rejected by the public builder.
    b=OpenSHOP._build_global_dispatch(table_power_fixture();joint=true)
    @test add_table_power!(b,c)["skipped"]["electrical"]==2
end

@testset "Head table support removes joint-product relaxation optimism" begin
    table=TurbineTable([1.0,3.0],[1.0,2.0],[.9 .3;.9 .3],fill(1.0,2),fill(2.0,2))
    m=Model(HiGHS.Optimizer);set_silent(m)
    @variable(m,1<=q<=2);@variable(m,1<=h<=3)
    eta=@variable(m,lower_bound=.3,upper_bound=.9,base_name="eta_1_1")
    @variable(m,0<=p<=.1)
    @constraint(m,q==1.5);@constraint(m,h==2.0)
    qc=OpenSHOP._global_tensor_coordinate!(m,q,[1.0,2.0];name=:q)
    hc=OpenSHOP._global_tensor_coordinate!(m,h,[1.0,3.0];name=:h)
    @constraint(m,eta==.9*hc.weights[1]+.3*hc.weights[2])
    record=OpenSHOP._power_hull_unit!(m,q,h,eta,p,1.0,(1.0,2.0),(1.0,3.0),(.3,.9),
        (1.0,3.0),(.3,.9),1.0,1.0;name=:joint)
    m.ext[:global_power_hulls]=Dict("power_hull_1_1"=>record)
    m.ext[:global_tensor_turbines]=Dict("turbine_1_1"=>(qcoordinate=qc,hcoordinate=hc))
    generator=Generator(name=:Unit,plant=:Plant,qmin=1.0,qmax=2.0,pmin=.001,pmax=.1,
        hmin=1.0,hmax=3.0,turbine_table=table)
    system=HydroSystem(reservoirs=Reservoir[],junctions=Junction[],boundaries=Boundary[],tunnels=Tunnel[],rivers=River[],
        plants=[Plant(name=:Plant,source=:Source,target=:Target,pmax=.1)],generators=[generator])
    c=ScheduleCase(name="lp_witness",system=system,grid=[0.0,1.0],prices=[1.0])
    b=(m=m,u=reshape([1.0],1,1),GQ=reshape([q],1,1),P=reshape([p],1,1))
    @objective(m,Max,p)
    optimize!(m)
    @test termination_status(m)==JuMP.MOI.OPTIMAL
    @test objective_value(m)/.00981≈2.85 atol=1e-9
    profile=add_table_power!(b,c;axes=:head,gate=true)
    optimize!(m)
    @test termination_status(m)==JuMP.MOI.OPTIMAL
    @test objective_value(m)/.00981≈2.25 atol=2e-6
    @test objective_value(m)>=.00981*1.5*2*.6
    @test profile["variables_added"]==0
end

@testset "Interior electrical maximum and tightened operations" begin
    for interpolation in (:bilinear,:pchip_discharge)
        c=table_power_fixture(;interpolation,positive_electrical=true,operating=true)
        seed=dispatch_from_controls(c,reshape([1,0],2,1),reshape([8.0,0.0],2,1),zeros(0,1))
        @test seed["validation"]["valid"]
        g=first(c.system.generators)
        endpoints=max(OpenSHOP.generator_efficiency(g,2.0),OpenSHOP.generator_efficiency(g,4.0))
        @test OpenSHOP.generator_efficiency(g,seed["power"][1,1])>endpoints
        for gate in (false,true),fixed in (false,true)
            b=OpenSHOP._build_global_dispatch(c;joint=true,fixed_u=fixed ? seed["u"] : nothing)
            prior=Set(all_constraints(b.m;include_variable_in_set_constraints=false))
            profile=add_table_power!(b,c;axes=:head,gate)
            @test profile["units_added"]==2-fixed
            # The new physical power row must use .97 at the interior knot,
            # rather than the smaller values at the operating endpoints.
            power=variable_by_name(b.m,"p[1,1]")
            flow=variable_by_name(b.m,"gq[1,1]")
            onheadhi=last(b.m.ext[:global_power_hulls]["power_hull_1_1"].axes[2])
            expected=-.00981*.97*onheadhi*.25*50/40
            cuts=[constraint_object(ref).func for ref in all_constraints(b.m;include_variable_in_set_constraints=false)
                if !(ref in prior) && constraint_object(ref).func isa AffExpr &&
                   coefficient(constraint_object(ref).func,power)!=0.0]
            @test length(cuts)==4
            @test minimum(abs(coefficient(expr,flow)-expected) for expr in cuts)<1e-12
            audit=OpenSHOP._lift_start!(b,c,seed)
            @test audit["valid"]
            @test audit["assigned"]==audit["variables"]==num_variables(b.m)
            @test audit["objective"]≈seed["objective"] atol=1e-6
        end
    end
end

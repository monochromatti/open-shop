using Test, JuMP, OpenSHOP, HiGHS
isdefined(@__MODULE__,:TableConsistency) || include(joinpath(@__DIR__,"table_consistency.jl"))
using .TableConsistency

function consistency_residual(m,assigned)
    residual=0.0
    for ref in all_constraints(m;include_variable_in_set_constraints=true)
        obj=constraint_object(ref)
        if obj.set isa JuMP.MOI.SOS2
            values=[JuMP.value(v->assigned[v],x) for x in obj.func]
            e=OpenSHOP._sos2_residual(values,obj.set.weights)
        else
            x=obj.func isa Number ? obj.func : JuMP.value(v->assigned[v],obj.func)
            set=obj.set
            e=set isa JuMP.MOI.EqualTo ? abs(x-set.value) :
                set isa JuMP.MOI.LessThan ? max(0.0,x-set.upper) :
                set isa JuMP.MOI.GreaterThan ? max(0.0,set.lower-x) :
                set isa JuMP.MOI.Interval ? max(0.0,set.lower-x,x-set.upper) :
                error("unsupported consistency test set")
        end
        residual=max(residual,e)
    end
    residual
end

@testset "Transportation consistency preserves exact table points" begin
    for interpolation in (:bilinear,:pchip_discharge)
        table=TurbineTable([50.0,100.0,150.0],[2.0,5.0,9.0,14.0],
            [0.75 0.80 0.85;0.95 0.90 0.85;0.80 0.90 0.98;0.70 0.75 0.80],
            [2.0,2.0,2.0],[14.0,14.0,14.0];interpolation)
        # Extended negative discharge/head domains are numerical graph tests,
        # not operating domains. Coefficient signs must not affect validity.
        for domain in ((-5.0,17.0,-200.0,180.0),(0.0,14.0,60.0,130.0),
                       (5.0,5.0+1e-8,99.0,100.0),(7.0,7.0,80.0,120.0),(0.0,14.0,100.0,100.0))
            qlo,qhi,hlo,hhi=domain
            m=Model();@variable(m,qlo<=q<=qhi);@variable(m,hlo<=h<=hhi)
            eta=OpenSHOP._global_tensor_turbine!(m,table,q,h,domain...;name=:test_table)
            data=m.ext[:global_tensor_turbines]["test_table"]
            before=count(is_binary,all_variables(m))
            profile=add_table_consistency!((m=m,),nothing;mode=:all)
            @test count(is_binary,all_variables(m))==before
            @test profile["tables_added"]==Int(qlo!=qhi && hlo!=hhi)
            points=unique(vcat(qlo,qhi,(qlo+qhi)/2,filter(x->qlo<=x<=qhi,table.discharge),
                vec([(1-t)*data.qcoordinate.nodes[i]+t*data.qcoordinate.nodes[i+1]
                 for i in 1:(length(data.qcoordinate.nodes)-1),t in (1/3,2/3)])))
            heads=unique(vcat(hlo,hhi,(hlo+hhi)/2,filter(x->hlo<=x<=hhi,table.heads)))
            maximum_residual=0.0;maximum_error=0.0
            for qv in points,hv in heads
                assigned=OpenSHOP._global_tensor_turbine_values(data,qv,hv)
                assigned[q]=qv;assigned[h]=hv
                for lifter in get(m.ext,:experiment_start_lifters,[])
                    lifter(assigned)
                end
                @test length(assigned)==num_variables(m)
                maximum_residual=max(maximum_residual,consistency_residual(m,assigned))
                maximum_error=max(maximum_error,abs(assigned[eta]-OpenSHOP.turbine_efficiency(table,qv,hv;extrapolation=:linear)))
            end
            @test maximum_residual<1e-9
            @test maximum_error<1e-11
            @test all(v->has_lower_bound(v)&&has_upper_bound(v),all_variables(m))
        end
    end
end

@testset "A coupled table hull excludes independent-product optimism" begin
    # The additive bilinear surface has eta=.6+.1*lambda2+.2*mu2.
    # At both weight vectors(.5,.5), its exact value and joint hull value are.75.
    # Independent column/head McCormick products instead permit eta=.80.
    m=Model(HiGHS.Optimizer);set_silent(m)
    @variable(m,lambda[1:2]>=0);@variable(m,mu[1:2]>=0)
    foreach(v->fix(v,0.5;force=true),lambda)
    foreach(v->fix(v,0.5;force=true),mu)
    @variable(m,0.6<=eta<=0.9)
    @variable(m,products[1:2]>=0)
    column_values=(0.65,0.85)
    for (j,(lo,hi)) in enumerate(((0.6,0.7),(0.8,0.9)))
        value=column_values[j]
        @constraint(m,products[j]>=lo*mu[j])
        @constraint(m,products[j]>=value-hi*(1-mu[j]))
        @constraint(m,products[j]<=hi*mu[j])
        @constraint(m,products[j]<=value-lo*(1-mu[j]))
    end
    @constraint(m,eta==sum(products))
    @objective(m,Max,eta)
    optimize!(m)
    @test termination_status(m)==JuMP.MOI.OPTIMAL
    @test objective_value(m)≈0.8 atol=1e-10
    data=(nodal=[0.6 0.8;0.7 0.9],qcoordinate=(weights=lambda,),hcoordinate=(weights=mu,),
          eta=eta,r=Dict(),s=Dict(),A=zeros(1,2),B=zeros(1,2))
    m.ext[:global_tensor_turbines]=Dict("additive"=>data)
    add_table_consistency!((m=m,),nothing;mode=:bilinear)
    optimize!(m)
    @test termination_status(m)==JuMP.MOI.OPTIMAL
    @test objective_value(m)≈0.75 atol=1e-10
end

@testset "Cubic basis products preserve their shared head mass" begin
    table=TurbineTable([50.0,100.0],[0.0,1.0,2.0],
        [0.70 0.75;0.90 0.95;0.80 0.85],[0.0,0.0],[2.0,2.0];interpolation=:pchip_discharge)
    m=Model();@variable(m,0<=q<=2);@variable(m,50<=h<=100)
    OpenSHOP._global_tensor_turbine!(m,table,q,h,0.0,2.0,50.0,100.0;name=:cubic)
    data=m.ext[:global_tensor_turbines]["cubic"]
    add_table_consistency!((m=m,),nothing)
    assigned=OpenSHOP._global_tensor_turbine_values(data,0.5,75.0)
    assigned[q]=0.5;assigned[h]=75.0
    lift_table_consistency!(m,assigned)
    @test consistency_residual(m,assigned)<1e-10
    record=m.ext[:table_consistency]["cubic"]
    i=first(sort!(collect(keys(data.r))))
    rv=assigned[data.r[i]];upper=upper_bound(data.r[i])
    # Each independent McCormick envelope permits upper/2 at mu=.5.
    optimistic=upper/2
    @test optimistic<=rv+1e-12
    @test optimistic>=rv-upper/2-1e-12
    @test 2optimistic>rv+1e-3
    for j in 1:2
        assigned[record.R[(i,j)]]=optimistic
    end
    @test consistency_residual(m,assigned)>1e-3
end

@testset "Consistency coverage and size controls" begin
    function size_fixture()
        m=Model();@variable(m,0<=q<=2);@variable(m,50<=h<=100)
        for (name,interpolation) in ((:linear,:bilinear),(:cubic,:pchip_discharge))
            table=TurbineTable([50.0,100.0],[0.0,1.0,2.0],
                [0.70 0.75;0.90 0.95;0.80 0.85],[0.0,0.0],[2.0,2.0];interpolation)
            OpenSHOP._global_tensor_turbine!(m,table,q,h,0.0,2.0,50.0,100.0;name)
        end
        m
    end
    m=size_fixture()
    profile=add_table_consistency!((m=m,),nothing;mode=:bilinear)
    @test profile["table_names"]==["linear"]
    @test profile["skipped_mode"]==1
    @test profile["variables_added"]==6
    @test_throws ArgumentError add_table_consistency!((m=m,),nothing)
    selected=add_table_consistency!((m=size_fixture(),),nothing;selected=[:cubic])
    @test selected["table_names"]==["cubic"]
    @test selected["pchip_tables"]==1
    limited=add_table_consistency!((m=size_fixture(),),nothing;max_added_variables=6)
    @test limited["table_names"]==["linear"]
    @test limited["skipped_size"]==1
    @test add_table_consistency!((m=size_fixture(),),nothing;max_tables=0)["tables_added"]==0
    empty_model=Model()
    @test_throws ArgumentError add_table_consistency!((m=empty_model,),nothing;mode=:unknown)
    @test_throws ArgumentError add_table_consistency!((m=empty_model,),nothing;max_tables=-1)
    @test add_table_consistency!((m=empty_model,),nothing;max_added_variables=0)["tables_added"]==0
    @test !haskey(empty_model.ext,:experiment_start_lifters)
end

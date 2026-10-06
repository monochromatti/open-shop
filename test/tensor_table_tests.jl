@testset "Independent tensor table graphs" begin
    function tensor_residual(m,values)
        residual=0.0
        for ref in all_constraints(m;include_variable_in_set_constraints=true)
            obj=constraint_object(ref)
            if obj.set isa OpenSHOP.MOI.SOS2
                x=[values[v] for v in obj.func]
                order=sortperm(obj.set.weights)
                x=x[order]
                violation=minimum(maximum(abs(x[j]) for j in eachindex(x) if j!=i && j!=i+1;init=0.0) for i in 1:(length(x)-1))
            else
                x=obj.func isa Number ? obj.func : JuMP.value(v->values[v],obj.func)
                set=obj.set
                violation=set isa OpenSHOP.MOI.EqualTo ? abs(x-set.value) :
                    set isa OpenSHOP.MOI.LessThan ? max(0.0,x-set.upper) :
                    set isa OpenSHOP.MOI.GreaterThan ? max(0.0,set.lower-x) :
                    set isa OpenSHOP.MOI.Interval ? max(0.0,set.lower-x,x-set.upper) :
                    error("unsupported tensor test constraint")
            end
            residual=max(residual,violation)
        end
        residual
    end
    for interpolation in (:bilinear,:pchip_discharge)
        table=TurbineTable([50.0,90.0,140.0],[2.0,5.0,10.0,16.0],
            [0.70 0.75 0.79;0.91 0.95 0.93;0.85 0.92 0.96;0.73 0.81 0.86],
            [2.0,2.0,2.0],[16.0,16.0,16.0];interpolation)
        # Extension, clipped original cells, narrow cells, and singleton axes.
        for (qlo,qhi,hlo,hhi) in ((-1.0,19.0,30.0,170.0),(3.0,14.0,60.0,130.0),
                                  (5.0,5.0+1e-8,89.0,90.0),(7.0,7.0,100.0,100.0),
                                  (0.0,16.0,90.0,90.0),(10.0,10.0,30.0,170.0))
            m=Model()
            @variable(m,qlo<=q<=qhi)
            @variable(m,hlo<=h<=hhi)
            eta=OpenSHOP._global_tensor_turbine!(m,table,q,h,qlo,qhi,hlo,hhi;name=:eta)
            data=m.ext[:global_tensor_turbines]["eta"]
            qs=unique(vcat(qlo,qhi,(qlo+qhi)/2,filter(x->qlo<=x<=qhi,table.discharge)))
            hs=unique(vcat(hlo,hhi,(hlo+hhi)/2,filter(x->hlo<=x<=hhi,table.heads)))
            append!(qs,[(1-t)*data.qcoordinate.nodes[i]+t*data.qcoordinate.nodes[i+1]
                for i in 1:(length(data.qcoordinate.nodes)-1),t in (1/3,2/3)])
            max_error=0.0
            max_residual=0.0
            for qv in qs,hv in hs
                values=OpenSHOP._global_tensor_turbine_values(data,qv,hv)
                values[q]=qv;values[h]=hv
                @test length(values)==num_variables(m)
                max_error=max(max_error,abs(values[eta]-OpenSHOP.turbine_efficiency(table,qv,hv;extrapolation=:linear)))
                max_residual=max(max_residual,tensor_residual(m,values))
            end
            @test max_error<1e-11
            @test max_residual<1e-10
            @test all(v->has_lower_bound(v)&&has_upper_bound(v),all_variables(m))
            @test count(is_binary,all_variables(m))==0
            interpolation==:bilinear && @test isempty(data.r)
            @test length(data.w)==length(data.r)
        end
        # A common head variable and identical clipped nodes reuse one SOS2 axis.
        m=Model();@variable(m,0<=q<=16);@variable(m,60<=h<=130)
        OpenSHOP._global_tensor_turbine!(m,table,q,h,0.0,16.0,60.0,130.0;name=:one)
        first_head=m.ext[:global_tensor_turbines]["one"].hcoordinate
        envelope=TableCurve(table.heads,table.qmin)
        OpenSHOP._global_tensor_table!(m,envelope,h,60.0,130.0;name=:envelope)
        @test m.ext[:global_tensor_tables]["envelope"].coordinate===first_head
        @test length(m.ext[:global_tensor_coordinates])==2
    end
    # Irregular, nonmonotone efficiency tables exercise independent polynomial
    # ranges and exact interpolation away from the hand-chosen fixture knots.
    rng=OpenSHOP.MersenneTwister(1407)
    for interpolation in (:bilinear,:pchip_discharge), trial in 1:12
        table=TurbineTable([35.0,80.0,130.0],[1.0,4.0,9.0,17.0],
            0.65 .+ 0.33 .* rand(rng,4,3),fill(1.0,3),fill(17.0,3);interpolation)
        for (qlo,qhi,hlo,hhi) in ((-1.0,21.0,10.0,180.0),(3.1,14.3,55.0,115.0),
                                  (8.0,8.0,90.0,90.0),(2.0,15.0,90.0,90.0))
            m=Model();@variable(m,qlo<=q<=qhi);@variable(m,hlo<=h<=hhi)
            eta=OpenSHOP._global_tensor_turbine!(m,table,q,h,qlo,qhi,hlo,hhi;name=:random_eta)
            data=m.ext[:global_tensor_turbines]["random_eta"]
            points=vcat([(qlo,hlo),(qhi,hhi)],
                [(qlo+rand(rng)*(qhi-qlo),hlo+rand(rng)*(hhi-hlo)) for _ in 1:6])
            for (qv,hv) in points
                values=OpenSHOP._global_tensor_turbine_values(data,qv,hv)
                values[q]=qv;values[h]=hv
                original=OpenSHOP.turbine_efficiency(table,qv,hv;extrapolation=:linear)
                @test values[eta]≈original atol=1e-11
                @test tensor_residual(m,values)<1e-9
                @test length(values)==num_variables(m)
                @test lower_bound(eta)<=original<=upper_bound(eta)
            end
        end
    end
    @test_throws DomainError OpenSHOP._global_tensor_coordinate_weights([0.0,1.0],1.1)
    @test OpenSHOP._global_tensor_coordinate_weights([3.0],3.0)==[1.0]
end

@testset "Tensor dispatch, shared heads and complete starts" begin
    copywith(x;kwargs...)=typeof(x)(;merge(NamedTuple{fieldnames(typeof(x))}(Tuple(getfield(x,n) for n in fieldnames(typeof(x)))),(;kwargs...))...)
    for interpolation in (:bilinear,:pchip_discharge)
        fixture=analytic_global_fixture(:on)
        c=fixture.case
        table=TurbineTable([35.0,80.0,110.0],[5.0,10.0,20.0],
            [0.82 0.86 0.88;0.89 0.93 0.94;0.84 0.90 0.92],
            [5.0,5.0,5.0],[20.0,20.0,20.0];interpolation)
        g=copywith(only(c.system.generators);turbine_table=table)
        sys=copywith(c.system;generators=[g,copywith(g;name=:SecondUnit)])
        c=copywith(c;system=sys)
        seed=dispatch_from_controls(c,reshape([1,0],2,1),reshape([8.0,0.0],2,1),zeros(0,1))
        @test seed["validation"]["valid"]
        b=OpenSHOP._build_global_dispatch(c;joint=true)
        audit=OpenSHOP._lift_start!(b,c,seed)
        @test audit["valid"]
        @test audit["assigned"]==audit["variables"]
        @test audit["objective"]≈seed["objective"] atol=1e-6
        tables=b.m.ext[:global_tensor_turbines]
        @test tables["turbine_1_1"].hcoordinate===tables["turbine_2_1"].hcoordinate
        fixed=OpenSHOP._build_global_dispatch(c;joint=true,fixed_u=seed["u"])
        lifted=OpenSHOP._lift_start!(fixed,c,seed)
        @test lifted["valid"]
        @test lifted["assigned"]==lifted["variables"]
        @test lifted["objective"]≈seed["objective"] atol=1e-6
        @test fixed.m.ext[:global_tensor_turbines]["turbine_2_1"].qcoordinate.nodes==[0.0]
        result=solve(c;initial=seed,time_limit=20.0,relative_gap=1e-3)
        @test result["accepted"]
        @test result["start_audit"]["valid"]
        @test result["model_profile"]["sos2_constraints"]>0
        @test result["scip_diagnostics"]["available"]
        @test result["raw_sos2_residual"]<=1e-6
        @test result["feasible_lower_bound"]>=seed["objective"]-1e-6
    end
    @test OpenSHOP._sos2_residual([0.5,0.5,0.0],[1.,2.,3.])==0
    @test OpenSHOP._sos2_residual([0.5,0.0,0.5],[1.,2.,3.])==0.5
    @test OpenSHOP._sos2_residual([0.5,0.0,0.5],[1.,3.,2.])==0
    @test OpenSHOP._sos2_residual([1.0,0.0,1e-9],[1.,2.,3.])==1e-9
end

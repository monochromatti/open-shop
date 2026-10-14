@testset "Delta-level reference discharge curves" begin
    a=TableCurve([-3.0,0.0,2.0,4.0],[-9.0,0.0,8.0,10.0])
    b=TableCurve([-4.0,-1.0,0.0,1.0,3.0],[-12.0,-2.0,0.0,6.0,14.0])
    c=OpenSHOP.LevelFlowCurve(;reference_levels=[100.0,110.0],curves=[a,b])
    f(d,r)=OpenSHOP.level_flow_value(c,d,r;extrapolation=:linear)
    @test OpenSHOP.level_flow_value(c,0.5,105.0)==2.5
    @test f(0.5,95.0)==2.0
    @test f(0.5,115.0)==3.0
    @test all(f(0.0,r)==0.0 for r in (90.0,100.0,105.0,110.0,120.0))
    @test f(-2.0,105.0)≈-17/3
    @test f(-2.0,105.0)!=-f(2.0,105.0)
    @test f(5.0,110.0)==22.0
    @test_throws DomainError OpenSHOP.level_flow_value(c,5.0,105.0)
    # An exact reference only needs that curve's delta domain.
    @test OpenSHOP.level_flow_value(c,4.0,100.0)==10.0
    @test OpenSHOP.level_flow_value(c,-4.0,110.0)==-12.0
    @test_throws DomainError f(NaN,105.0)
    @test_throws DomainError f(1.0,Inf)
    @test_throws ArgumentError OpenSHOP.level_flow_value(c,1.0,105.0;extrapolation=:clamp)
    single=OpenSHOP.LevelFlowCurve([100.0],[a])
    @test OpenSHOP.level_flow_value(single,1.0,-1000.0)==4.0
    @test OpenSHOP.level_flow_derivatives(single,1.0,1000.0)==
        (delta=4.0,reference=0.0,delta_reference=0.0)
    @test_throws ArgumentError OpenSHOP.LevelFlowCurve(Float64[],TableCurve[])
    @test_throws ArgumentError OpenSHOP.LevelFlowCurve([100.0,100.0],[a,b])
    @test_throws ArgumentError OpenSHOP.LevelFlowCurve([100.0,Inf],[a,b])
    @test_throws ArgumentError OpenSHOP.LevelFlowCurve([100.0],[a,b])
    @test_throws ArgumentError OpenSHOP.LevelFlowCurve([100.0],
        [TableCurve([-1.0,1.0],[-1.0,1.0])])
    @test_throws ArgumentError OpenSHOP.LevelFlowCurve([100.0],
        [TableCurve([-1.0,0.0,1.0],[-1.0,0.1,1.0])])
    @test_throws ArgumentError OpenSHOP.LevelFlowCurve([100.0],
        [TableCurve([-1.0,0.0,1.0,2.0],[-1.0,0.0,2.0,1.0])])
    # Compare analytic gradients and the cross derivative with independently
    # evaluated secants away from knot lines, including delta extrapolation.
    eps=1e-4
    for (delta,reference) in ((0.5,105.0),(-2.0,103.0),(5.0,107.0),(0.5,95.0),(0.5,115.0))
        d=OpenSHOP.level_flow_derivatives(c,delta,reference;extrapolation=:linear)
        @test d.delta≈(f(delta+eps,reference)-f(delta-eps,reference))/(2eps) atol=1e-8
        @test d.reference≈(f(delta,reference+eps)-f(delta,reference-eps))/(2eps) atol=1e-8
        cross=(f(delta+eps,reference+eps)-f(delta+eps,reference-eps)-
               f(delta-eps,reference+eps)+f(delta-eps,reference-eps))/(4eps^2)
        @test d.delta_reference≈cross atol=2e-6
    end
    @test OpenSHOP.level_flow_derivatives(c,0.5,105.0)==
        (delta=5.0,reference=0.1,delta_reference=0.2)
    @test OpenSHOP.level_flow_derivatives(c,2.0,100.0).delta==1.0
    m=Model()
    op=OpenSHOP.level_flow_operator(m,c;name=:delta_reference_test)
    @variable(m,delta)
    @variable(m,reference)
    expr=op(delta,reference)
    @test JuMP.value(v->v==delta ? 0.5 : 105.0,expr)==2.5
    registered=OpenSHOP.MOI.get(m,OpenSHOP.MOI.UserDefinedFunction(:delta_reference_test,2))
    gradient=zeros(2)
    registered[2](gradient,0.5,105.0)
    @test gradient==[5.0,0.1]
    hessian=fill(NaN,2,2)
    registered[3](hessian,0.5,105.0)
    @test hessian[1,1]==hessian[2,2]==0.0
    @test hessian[2,1]==0.2
end

@testset "Level-flow equations reproduce uneven calibration grids" begin
    c=OpenSHOP.LevelFlowCurve([100.0,110.0,120.0],[
        TableCurve([-3.0,0.0,2.0,4.0],[-9.0,0.0,8.0,10.0]),
        TableCurve([-4.0,-1.0,0.0,1.0,3.0],[-12.0,-2.0,0.0,6.0,14.0]),
        TableCurve([-2.0,0.0,0.5,5.0],[-5.0,0.0,4.0,22.0])])
    function level_flow_graph_residual(m,values)
        residual=0.0
        for ref in all_constraints(m;include_variable_in_set_constraints=true)
            obj=constraint_object(ref)
            if obj.set isa OpenSHOP.MOI.SOS2
                x=[values[v] for v in obj.func][sortperm(obj.set.weights)]
                violation=minimum(maximum(abs(x[j]) for j in eachindex(x)
                    if j!=i && j!=i+1;init=0.0) for i in 1:(length(x)-1))
            else
                x=obj.func isa Number ? obj.func : JuMP.value(v->values[v],obj.func)
                set=obj.set
                violation=set isa OpenSHOP.MOI.EqualTo ? abs(x-set.value) :
                    set isa OpenSHOP.MOI.LessThan ? max(0.0,x-set.upper) :
                    set isa OpenSHOP.MOI.GreaterThan ? max(0.0,set.lower-x) :
                    set isa OpenSHOP.MOI.Interval ? max(0.0,set.lower-x,x-set.upper) :
                    error("unsupported level-flow graph constraint")
            end
            residual=max(residual,violation)
        end
        residual
    end
    for (dlo,dhi,rlo,rhi) in ((-6.0,7.0,90.0,130.0),(-1.3,2.8,103.0,117.0),
                              (0.0,0.0,90.0,130.0),(-2.0,3.0,105.0,105.0),
                              (0.7,0.7,115.0,115.0))
        m=Model()
        @variable(m,dlo<=delta<=dhi)
        @variable(m,rlo<=reference<=rhi)
        flow=OpenSHOP._global_level_flow!(m,c,delta,reference,dlo,dhi,rlo,rhi;name=:flow)
        data=m.ext[:global_level_flow_tables]["flow"]
        deltas=unique(vcat(data.dcoordinate.nodes,
            [(1-t)*data.dcoordinate.nodes[i]+t*data.dcoordinate.nodes[i+1]
             for i in 1:(length(data.dcoordinate.nodes)-1) for t in (1/3,2/3)]))
        refs=unique(vcat(data.rcoordinate.nodes,
            [(1-t)*data.rcoordinate.nodes[i]+t*data.rcoordinate.nodes[i+1]
             for i in 1:(length(data.rcoordinate.nodes)-1) for t in (1/3,2/3)]))
        max_error=0.0;max_residual=0.0;bound_violation=0.0
        for d in deltas,r in refs
            values=OpenSHOP._global_level_flow_values(data,d,r)
            values[delta]=d
            values[reference]=r
            actual=OpenSHOP.level_flow_value(c,d,r;extrapolation=:linear)
            max_error=max(max_error,abs(values[flow]-actual))
            max_residual=max(max_residual,level_flow_graph_residual(m,values))
            bound_violation=max(bound_violation,lower_bound(flow)-actual,actual-upper_bound(flow))
        end
        @test max_error<1e-11
        @test max_residual<1e-10
        @test bound_violation==0.0
    end
    # A single reference gives the same flow at every reference level.
    single=OpenSHOP.LevelFlowCurve([105.0],[first(c.curves)])
    for (dlo,dhi) in ((-5.0,6.0),(0.0,0.0))
        m=Model()
        @variable(m,dlo<=delta<=dhi)
        @variable(m,90<=reference<=130)
        flow=OpenSHOP._global_level_flow!(m,single,delta,reference,
            dlo,dhi,90.0,130.0;name=:single)
        data=m.ext[:global_level_flow_tables]["single"]
        for d in unique(vcat(data.dcoordinate.nodes,(dlo+dhi)/2)),r in (90.0,105.0,130.0)
            values=OpenSHOP._global_level_flow_values(data,d,r)
            values[delta]=d
            values[reference]=r
            @test level_flow_graph_residual(m,values)<1e-10
            @test values[flow]≈OpenSHOP.level_flow_value(single,d,r;extrapolation=:linear) atol=1e-11
        end
    end
end

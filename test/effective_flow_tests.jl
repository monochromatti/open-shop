@testset "Effective flow product lift" begin
    function effective_residual(m,values)
        residual=0.0
        for ref in all_constraints(m;include_variable_in_set_constraints=true)
            obj=constraint_object(ref)
            x=obj.func isa Number ? obj.func : JuMP.value(v->values[v],obj.func)
            set=obj.set
            e=set isa OpenSHOP.MOI.EqualTo ? abs(x-set.value) :
                set isa OpenSHOP.MOI.LessThan ? max(0.0,x-set.upper) :
                set isa OpenSHOP.MOI.GreaterThan ? max(0.0,set.lower-x) :
                set isa OpenSHOP.MOI.Interval ? max(0.0,set.lower-x,x-set.upper) :
                error("unsupported effective-flow test constraint")
            residual=max(residual,e)
        end
        residual
    end
    for (eta_lower,eta_upper,min_efficiency) in ((-0.5,1.2,0.6),(0.4,0.9,0.4),(-1.0,-0.2,0.0))
        m=Model();@variable(m,0<=q<=10);@variable(m,eta_lower<=eta<=eta_upper)
        flow=OpenSHOP._global_effective_flow!(m,q,eta,10.0,min_efficiency;name=:effective)
        data=m.ext[:global_effective_flows]["effective"]
        @test data.q===q
        @test data.eta===eta
        @test data.flow===flow
        @test lower_bound(flow)==0.0
        @test upper_bound(flow)==10*max(0.0,min(1.0,eta_upper))
        @test (lower_bound(eta),upper_bound(eta))==(eta_lower,eta_upper)
        # All off-state continuations, including negative and >1, remain valid.
        for ev in (eta_lower,(eta_lower+eta_upper)/2,eta_upper)
            @test effective_residual(m,Dict(q=>0.0,eta=>ev,flow=>0.0))==0.0
        end
        if eta_upper>=min_efficiency
            for qv in (0.0,3.0,10.0),ev in (min_efficiency,(min_efficiency+min(1.0,eta_upper))/2,min(1.0,eta_upper))
                values=Dict(q=>qv,eta=>ev,flow=>qv*ev)
                @test effective_residual(m,values)<1e-12
            end
        end
        # A wrong lifted value violates the product equality even when Q is off.
        @test effective_residual(m,Dict(q=>0.0,eta=>eta_lower,flow=>0.1))>=0.1
    end
    m=Model();@variable(m,0<=q<=10);@variable(m,-1<=eta<=1.2)
    flow=OpenSHOP._global_effective_flow!(m,50*(q/50),eta,10.0,0.6)
    @test effective_residual(m,Dict(q=>5.0,eta=>0.9,flow=>4.5))<1e-12
    # These operating values violate declared physical efficiency bounds; cuts
    # exclude them without imposing those bounds on the off-state continuation.
    @test effective_residual(m,Dict(q=>5.0,eta=>1.1,flow=>5.5))>=0.5
    @test effective_residual(m,Dict(q=>5.0,eta=>0.5,flow=>2.5))>=0.5
    @test_throws ArgumentError OpenSHOP._global_effective_flow!(m,q,eta,-1.0,0.6)
    @test_throws ArgumentError OpenSHOP._global_effective_flow!(m,q,eta,10.0,1.1)
    @test_throws ArgumentError OpenSHOP._global_effective_flow!(m,q,eta,Inf,0.6)
    unbounded=@variable(m)
    @test_throws ArgumentError OpenSHOP._global_effective_flow!(m,q,unbounded,10.0,0.6)
end

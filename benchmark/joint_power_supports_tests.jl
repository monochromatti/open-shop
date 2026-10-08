using Test, JuMP, OpenSHOP, HiGHS, SCIP
if !isdefined(Main,:JointPowerOracle)
    include("joint_power_supports.jl")
end
if !isdefined(Main,:table_power_bounds_case)
    include(joinpath(@__DIR__,"..","test","table_power_bounds_tests.jl"))
end

joint_test_value(oracle,qa,hb,q,h)=sum(qa.*OpenSHOP._global_tensor_coordinate_weights(oracle.qnodes,q))+
    sum(hb.*OpenSHOP._global_tensor_coordinate_weights(oracle.hnodes,h))

@testset "Joint polynomial transform retains the exact physical graph" begin
    for interpolation in (:bilinear,:pchip_discharge)
        table=TurbineTable([1.,1.5,2.],[1.,2.,3.],
            [.7 .95 .8;.5 .8 .7;.4 .6 .4],fill(1.,3),fill(3.,3);interpolation)
        for (qlo,qhi,hlo,hhi) in ((.25,.75,.6,.9),(1.,2.,1.,1.5),(2.,3.,1.5,2.),(3.,4.,2.,2.5))
            B=_joint_power_bernstein(table,qlo,qhi,hlo,hhi,.95)
            @test maximum(B)≈OpenSHOP._power_box_upper(table,qlo,qhi,hlo,hhi,.95,0.,0.) atol=1e-13
            # Bernstein polynomial identities check the implementation; the
            # complete coefficient enclosure (not these points) certifies cuts.
            for x in (.0,.17,.53,1.),y in (.0,.31,.79,1.)
                q=qlo+(qhi-qlo)*x;h=hlo+(hhi-hlo)*y
                physical=.00981*.95*q*h*OpenSHOP.turbine_efficiency(table,q,h;extrapolation=:linear)
                reconstructed=sum(B[i+1,j+1]*binomial(4,i)*x^i*(1-x)^(4-i)*
                    binomial(2,j)*y^j*(1-y)^(2-j) for i in 0:4,j in 0:2)
                @test reconstructed>=physical-1e-12
                @test reconstructed-physical<1e-7
            end
        end
    end
end

@testset "Joint supports certify original concavity, PCHIP and extrapolation" begin
    for interpolation in (:bilinear,:pchip_discharge)
        table=TurbineTable([1.,1.5,2.],[1.,2.,3.],
            [.7 .95 .8;.5 .8 .7;.4 .6 .4],fill(1.,3),fill(3.,3);interpolation)
        for (qbox,hbox) in (((1.2,1.8),(1.1,1.7)),((.25,4.),(.6,2.5)),((1.3,1.3),(1.6,1.6)))
            qnodes=OpenSHOP._global_tensor_nodes(table.discharge,0.,4.)
            hnodes=OpenSHOP._global_tensor_nodes(table.heads,-2.,3.)
            oracle=JointPowerOracle(table,qbox,hbox,qnodes,hnodes,.95)
            q=(qbox[1]+qbox[2])/2;h=(hbox[1]+hbox[2])/2
            qw=OpenSHOP._global_tensor_coordinate_weights(qnodes,q)
            hw=OpenSHOP._global_tensor_coordinate_weights(hnodes,h)
            support=joint_power_support(oracle,qw,hw)
            @test support!==nothing&&support.certified
            support===nothing&&continue
            @test support.residual_upper<=0
            @test joint_power_certificate(oracle,support.qcoefficients,support.hcoefficients)<=1e-10
            @test oracle.bernstein_rows==60*length(oracle.cells)
            # Physical schedule points supplement the global certificate.
            for q in range(qbox...;length=13),h in range(hbox...;length=11)
                exact=.00981*.95*q*h*OpenSHOP.turbine_efficiency(table,q,h;extrapolation=:linear)
                @test joint_test_value(oracle,support.qcoefficients,support.hcoefficients,q,h)>=exact-1e-10
                # Zero-to-maximum electrical efficiency remains safe even
                # where a linear extension has negative hydraulic power.
                @test joint_test_value(oracle,support.qcoefficients,support.hcoefficients,q,h)>=max(0.,exact)-1e-10
            end
            second=joint_power_support(oracle,0.4 .* qw,0.4 .* hw)
            @test second!==nothing&&second.certified
            @test second.value≈.4*support.value atol=1e-7
            @test oracle.solves==2
            @test joint_power_support(oracle,qw,hw;time_limit=0.)===nothing
        end
    end
    # Head efficiency decreases linearly, so h*eta is strictly concave.
    table=TurbineTable([1.,3.],[1.,2.],[.9 .3;.9 .3],fill(1.,2),fill(2.,2))
    oracle=JointPowerOracle(table,(1.,1.),(1.,3.),[1.],[-2.,1.,3.],1.)
    rawq=[0.];rawh=[0.,.00981*.9,.00981*.9]
    @test joint_power_certificate(oracle,rawq,rawh)>1e-3
    @test .00981*2*.6>joint_test_value(oracle,rawq,rawh,1.,2.)+1e-3
    corrected=joint_power_support(oracle,[1.],[0.,.5,.5])
    @test joint_power_certificate(oracle,corrected.qcoefficients,corrected.hcoefficients)<=1e-10
    @test corrected.value>=.00981*2*.6
    @test joint_power_certificate(oracle,[-10.],fill(-10.,3))>10.
    # Arbitrary signed coefficient decompositions retain exactly the same
    # physical inequality, including a negative off-domain nodal intercept.
    qa=corrected.qcoefficients.-3.;hb=corrected.hcoefficients.+3.
    @test minimum(qa)<0
    @test joint_power_certificate(oracle,qa,hb)<=1e-9
    @test joint_test_value(oracle,qa,hb,1.,2.)≈corrected.value atol=1e-10
    @test_throws ArgumentError JointPowerOracle(table,(1.,2.),(1.,3.),[1.],[1.,3.],1.)
    @test_throws ArgumentError JointPowerOracle(table,(1.,2.),(1.,3.),[1.,2.],[1.,3.],1.;coefficient_bound=.001)
    @test_throws ArgumentError joint_power_support(oracle,[NaN],[0.,.5,.5])
    @test_throws ArgumentError joint_power_certificate(oracle,[0.],fill(NaN,3))
end

@testset "Joint native affine rows preserve full starts and signed off heads" begin
    for interpolation in (:bilinear,:pchip_discharge),fixed in (false,true)
        c=table_power_bounds_case(;interpolation)
        seed=dispatch_from_controls(c,reshape([1,0],2,1),reshape([8.,0.],2,1),zeros(0,1))
        @test seed["validation"]["valid"]
        b=OpenSHOP._build_global_dispatch(c;joint=true,fixed_u=fixed ? seed["u"] : nothing)
        @test OpenSHOP._lift_start!(b,c,seed)["valid"]
        before=(num_variables(b.m),num_constraints(b.m;count_variable_in_set_constraints=false))
        set_optimizer(b.m,SCIP.Optimizer);JuMP.MOI.Utilities.attach_optimizer(backend(b.m))
        separator=install_joint_power_supports!(b,c;certificate_budget=0.)
        @test before==(num_variables(b.m),num_constraints(b.m;count_variable_in_set_constraints=false))
        assigned=Dict(v=>Float64(start_value(v)) for v in all_variables(b.m))
        @test length(assigned)==num_variables(b.m)
        byreference=Dict(SCIP.VarRef(optimizer_index(v).value)=>assigned[v] for v in all_variables(b.m))
        native=[byreference[r] for r in separator.references]
        at(j)=native[j]
        for coordinate in separator.coordinates
            q=coordinate.discharge.source;h=coordinate.head.source
            oracle=JointPowerOracle(q.table,q.qbox,q.hbox,q.nodes,h.nodes,q.electrical_max)
            qw=[OpenSHOP._power_cut_value(x,at) for x in coordinate.discharge.onweights]
            hw=[OpenSHOP._power_cut_value(x,at) for x in coordinate.head.onweights]
            candidate=joint_power_support(oracle,qw,hw)
            @test candidate!==nothing&&candidate.certified
            # A gauge shift deliberately exercises negative affine discharge
            # constants and signed coefficients in native row assembly.
            qa=candidate.qcoefficients.-3.;hb=candidate.hcoefficients.+3.
            row=_joint_power_row(coordinate,qa,hb)
            original=h.power-sum(qa[j]*q.onweights[j] for j in eachindex(qa))-
                sum(hb[j]*h.onweights[j] for j in eachindex(hb))
            residual=40*OpenSHOP._power_cut_value(row,at)
            @test residual≈OpenSHOP._table_power_value(original,assigned) atol=1e-9
            @test residual<=1e-7
            @test sum(qw)≈OpenSHOP._table_power_value(q.commitment,assigned) atol=1e-12
            @test sum(hw)≈sum(qw) atol=1e-12
            if OpenSHOP._table_power_value(q.commitment,assigned)==0
                @test abs(residual)<=1e-12
                @test all(iszero,qw)&&all(iszero,hw)
            end
        end
        @test isempty(joint_power_supports_statistics(separator)["errors"])
        if !fixed
            off=b.m.ext[:global_power_hulls]["power_hull_2_1"]
            @test start_value(off.head)<0
        end
    end
    c=table_power_bounds_case(;analytic=true)
    b=OpenSHOP._build_global_dispatch(c;joint=true)
    @test install_joint_power_supports!(b,c)===nothing
end

function joint_native_witness(;certificate_budget=10.,max_cuts=64,corrupt=false)
    table=TurbineTable([1.,3.],[1.,3.],fill(.9,2,2),[1.,1.],[3.,3.])
    qbox=(1.,3.);hbox=(1.25,1.25)
    support,_,_,_=OpenSHOP._table_power_support_cache()
    slopes=[.00981*3*f for f in OpenSHOP._TABLE_POWER_FRACTIONS]
    coefficients=[OpenSHOP._table_power_head_coefficients(table,qbox,(1.,3.),[1.25],1.,a,support) for a in slopes]
    m=Model(SCIP.Optimizer);set_silent(m)
    @variable(m,u,Bin);@constraint(m,u==1)
    @variable(m,1<=q<=3);@variable(m,h_on);@constraint(m,h_on==1.25)
    @variable(m,0<=p<=1)
    @variable(m,0<=lambda[1:3]<=1)
    @constraint(m,lambda[1]==0);@constraint(m,sum(lambda)==1)
    @constraint(m,q==lambda[2]+3*lambda[3])
    for j in eachindex(slopes)
        @constraint(m,p<=slopes[j]*q+coefficients[j][1]*u)
    end
    @variable(m,dummy,Bin);@variable(m,0<=auxiliary<=2)
    @constraint(m,auxiliary<=2*dummy);@constraint(m,auxiliary<=3-2*dummy)
    flow_cost=.00981*.9*1.25
    @objective(m,Max,p-flow_cost*q+.1*auxiliary)
    set_optimizer_attribute(m,"presolving/maxrounds",0)
    set_optimizer_attribute(m,"misc/usesymmetry",0)
    set_optimizer_attribute(m,"limits/time",10.)
    JuMP.MOI.Utilities.attach_optimizer(backend(m))
    optimizer=unsafe_backend(m)
    @test SCIP.SCIPsetHeuristics(optimizer,SCIP.SCIP_PARAMSETTING_OFF,true)==SCIP.SCIP_OKAY
    @test SCIP.SCIPsetSeparating(optimizer,SCIP.SCIP_PARAMSETTING_OFF,true)==SCIP.SCIP_OKAY
    set_optimizer_attribute(m,"separating/minefficacyroot",0.)
    set_optimizer_attribute(m,"separating/poolfreq",1)
    set_optimizer_attribute(m,"cutselection/hybrid/minorthoroot",0.)
    m.ext[:global_table_power_supports]=[
        OpenSHOP._TablePowerSupportCoordinate(1,1,:head,table,qbox,hbox,1.,[1.25],
            OpenSHOP._TablePowerTerm[u],q,p,u,q,h_on),
        OpenSHOP._TablePowerSupportCoordinate(1,1,:discharge,table,qbox,hbox,1.,[0.,1.,3.],
            OpenSHOP._TablePowerTerm[lambda[1]+u-1;lambda[2:3]],h_on,p,u,q,h_on)]
    b=(m=m,);c=(prices=[100.],grid=[0.,1.])
    separator=install_joint_power_supports!(b,c;certificate_budget,max_cuts)
    corrupt&&(table.efficiency[1,1]=NaN)
    optimize!(m)
    @test termination_status(m)==JuMP.MOI.OPTIMAL
    (;m,separator,flow_cost,q,p)
end

@testset "Joint root separator actually applies certified cuts and respects caps" begin
    witness=joint_native_witness()
    statistics=joint_power_supports_statistics(witness.separator)
    @test isempty(statistics["errors"])
    @test statistics["calls"]>0
    @test statistics["cuts_added_to_global_pool"]>0
    @test statistics["native_cuts_applied"]>0
    @test statistics["rounds"]<=statistics["max_rounds"]
    @test statistics["coordinate_checks"]<=statistics["max_coordinate_checks"]
    @test statistics["infeasible_flags"]==0
    @test objective_value(witness.m)≈.1 atol=1e-7
    @test value(witness.p)≈witness.flow_cost*value(witness.q) atol=1e-7
    @test all(cut->cut.candidate.certified&&cut.candidate.residual_upper<=0,witness.separator.certified)
    for keyword in ((certificate_budget=0.,),(max_cuts=0,))
        capped=joint_native_witness(;keyword...)
        stats=joint_power_supports_statistics(capped.separator)
        @test isempty(stats["errors"])
        @test stats["cuts_added_to_global_pool"]==0
        @test stats["coordinate_checks"]==0
        @test objective_value(capped.m)>.101
    end
    failed=joint_native_witness(;corrupt=true)
    stats=joint_power_supports_statistics(failed.separator)
    @test !isempty(stats["errors"])
    @test stats["cuts_added_to_global_pool"]==0
    @test stats["native_cuts_applied"]==0
end

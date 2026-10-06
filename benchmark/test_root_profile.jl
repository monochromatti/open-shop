using Test
include("root_profile.jl")

@testset "Certified power envelopes" begin
    for interpolation in (:bilinear,:pchip_discharge)
        table=TurbineTable([50.,90.,140.],[2.,5.,10.,16.],
            [.70 .75 .79;.91 .95 .93;.85 .92 .96;.73 .81 .86],
            [2.,2.,2.],[16.,16.,16.];interpolation)
        for (qlo,qhi,hlo,hhi) in ((0.,19.,30.,170.),(3.,14.,60.,130.),
                (5.,5.0+1e-8,89.,90.),(7.,7.,100.,100.))
            for slope in (-.5,0.,.5,1.5)
                intercept=power_support(table,qlo,qhi,hlo,hhi,.98,slope)
                residual=maximum(0.00981*q*h*.98*
                    OpenSHOP.turbine_efficiency(table,q,h;extrapolation=:linear)-slope*q-intercept
                    for q in range(qlo,qhi;length=101),h in range(hlo,hhi;length=37))
                @test residual<=1e-8
                @test isfinite(intercept)
            end
        end
    end
end

@testset "Tunnel domains enclose signed hydraulic flows" begin
    c=readcase(joinpath(@__DIR__,"cases","distributed-rivers.json"))
    b=OpenSHOP._build_global_dispatch(c;joint=true)
    before=(num_variables(b.m),num_constraints(b.m;count_variable_in_set_constraints=false))
    @test tighten_tunnel_domains!(b,c)>=0
    @test before==(num_variables(b.m),num_constraints(b.m;count_variable_in_set_constraints=false))
    byname=Dict(name(v)=>v for v in all_variables(b.m))
    function bounds(v)
        is_fixed(v) ? (fix_value(v),fix_value(v)) : (lower_bound(v),upper_bound(v))
    end
    for (i,e) in enumerate(c.system.tunnels),t in eachindex(c.prices)
        opening=OpenSHOP.opinterval(c,e.name,:opening,t,e.opening)
        cap=OpenSHOP.opinterval(c,e.name,:capacity,t,e.capacity)
        src=b.node_bounds[(e.source,t)];dst=b.node_bounds[(e.target,t)]
        for h0 in range(src...;length=13),h1 in range(dst...;length=13)
            flow=copysign(sqrt(abs(opening*(h0-h1)/e.resistance)),h0-h1)
            abs(flow)>cap && continue
            for (n,x) in (("q[$i,$t]",flow/50),
                    ("tunnel_positive_$(i)_$(t)",max(flow,0.)),
                    ("tunnel_negative_$(i)_$(t)",max(-flow,0.)))
                lo,hi=bounds(byname[n])
                @test lo-1e-9<=x<=hi+1e-9
            end
        end
    end
end


@testset "Junction propagation preserves forward trajectories" begin
    for reversed in (false,true)
        s=HydroSystem(reservoirs=[Reservoir(name=:lake,z0=100.,slope=.001,
                v0=2.,vmin=0.,vmax=4.,water_value=1.)],
            junctions=[Junction(name=:intake,hmin=50.,hmax=150.)],
            boundaries=[Boundary(name=:tail,head=0.)],
            tunnels=[Tunnel(name=:feed,source=reversed ? :intake : :lake,
                target=reversed ? :lake : :intake,resistance=.001,capacity=60.)],
            plants=[Plant(name=:plant,source=:intake,target=:tail,pmax=50.)],
            generators=[Generator(name=Symbol("g",i),plant=:plant,qmin=.001,qmax=25.,
                pmin=.001,pmax=30.,hmin=1.,hmax=160.,qcurvature=0.,hcurvature=0.,
                initial_on=0,minup=0.,mindown=0.) for i in 1:2],rivers=River[])
        c=ScheduleCase(name="radial",system=s,grid=[0.,1.],prices=[100.])
        b=OpenSHOP._build_global_dispatch(c;joint=true)
        tighten_network_domains!(b,c)
        add_plant_energy_bounds!(b,c)
        q=variable_by_name(b.m,"q[1,1]")
        @test reversed ? upper_bound(q)==0. : lower_bound(q)==0.
        for rate in (0.,.1,20.,40.)
            initial=dispatch_from_controls(c,fill(Int(rate>0),2,1),fill(rate/2,2,1),zeros(0,1))
            @test initial["validation"]["valid"]
            @test OpenSHOP._lift_start!(b,c,initial)["valid"]
        end
    end
end

@testset "Native statistics and root capture" begin
    mktempdir() do output
        rows=root_profile(joinpath(@__DIR__,"cases","turbine-tables.json"),output;
            seconds=30.,profiles=["baseline","hydraulic_domains","plant_energy","static_symmetry"],capture=true)
        for row in rows
            @test row["accepted"]
            @test row["diagnostics_close_error"]===nothing
            @test row["solver_error"]===nothing
            @test haskey(row["scip_statistics"],"lp")
            @test row["first_root_lp"]!==nothing
            @test row["first_root_lp"]["fractional_commitment_max"]>=0
        end
    end
end

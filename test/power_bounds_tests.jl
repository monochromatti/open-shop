@testset "Certified power envelopes" begin
    for interpolation in (:bilinear,:pchip_discharge)
        table=TurbineTable([50.,90.,140.],[2.,5.,10.,16.],
            [.70 .75 .79;.91 .95 .93;.85 .92 .96;.73 .81 .86],
            [2.,2.,2.],[16.,16.,16.];interpolation)
        for (qlo,qhi,hlo,hhi) in ((0.,19.,30.,170.),(3.,14.,60.,130.),
                (5.,5.0+1e-8,89.,90.),(7.,7.,100.,100.))
            for slope in (-.5,0.,.5,1.5), hslope in (-.02,0.,.05)
                intercept=OpenSHOP._power_support(table,qlo,qhi,hlo,hhi,.98,slope,hslope)
                residual=maximum(0.00981*q*h*.98*
                    OpenSHOP.turbine_efficiency(table,q,h;extrapolation=:linear)-slope*q-hslope*h-intercept
                    for q in range(qlo,qhi;length=101),h in range(hlo,hhi;length=37))
                @test residual<=1e-8
                @test isfinite(intercept)
            end
        end
    end
end


@testset "Power planes cover off head and tightened operating domains" begin
    for interpolation in (:bilinear,:pchip_discharge), electrical in
            (nothing,TableCurve([0.,20.],[.9,.98]))
        table=TurbineTable([50.,100.],[2.,5.,10.,16.],
            [.2 .8;.6 .9;.7 .95;.3 .85],[2.,2.],[16.,16.];interpolation)
        g=Generator(name=:g,plant=:plant,qmin=2.,qmax=16.,pmin=.001,pmax=20.,
            hmin=10.,hmax=110.,turbine_table=table,generator_efficiency_curve=electrical)
        s=HydroSystem(reservoirs=Reservoir[],junctions=Junction[],
            boundaries=[Boundary(name=:source,head=100.),Boundary(name=:tail,head=0.)],
            tunnels=Tunnel[],plants=[Plant(name=:plant,source=:source,target=:tail,pmax=20.)],
            generators=[g],rivers=River[])
        c=ScheduleCase(name="plane_domain",system=s,grid=[0.,1.],prices=[100.],
            operations=[OperationalSeries(object=:g,attribute=:qmin,times=[0.],values=[5.]),
                OperationalSeries(object=:g,attribute=:qmax,times=[0.],values=[10.])])
        m=Model(OpenSHOP.SCIP.Optimizer)
        set_silent(m)
        @variable(m,0<=q<=16)
        @variable(m,0<=p<=20)
        @variable(m,u,Bin)
        @variable(m,-60<=h<=150)
        b=(m=m,P=reshape([p],1,1),GQ=reshape([q],1,1),u=reshape([u],1,1),
            shared_heads=Dict((:plant,1)=>h))
        OpenSHOP._add_power_bounds!(b,c)
        cuts=all_constraints(m;include_variable_in_set_constraints=false)
        function violation(power,flow,head,on)
            point=Dict(p=>power,q=>flow,h=>head,u=>on)
            maximum(cuts;init=-Inf) do cut
                object=constraint_object(cut)
                value(v->point[v],object.func)-object.set.upper
            end
        end
        @test OpenSHOP.turbine_efficiency(table,0.,-50.;extrapolation=:linear)<0
        off_violation=maximum(violation(0.,0.,head,0.) for head in range(-60.,150.;length=15))
        on_violation=-Inf
        for flow in range(5.,10.;length=11),head in range(10.,110.;length=15)
            eta=OpenSHOP.turbine_efficiency(table,flow,head;extrapolation=:linear)
            0<=eta<=1 || continue
            shaft=.00981*flow*head*eta
            power=electrical===nothing ? shaft : .9*shaft/(1-.004*shaft)
            on_violation=max(on_violation,violation(power,flow,head,1.))
        end
        @test off_violation<=1e-8
        @test on_violation<=1e-8
        # Bounds must also exclude a loose 20 MW relaxation point without
        # excluding the independently calculated physical output.
        fix(q,7.5;force=true);fix(h,75.0;force=true);fix(u,1.0;force=true)
        @objective(m,Max,p)
        optimize!(m)
        shaft=.00981*7.5*75.0*OpenSHOP.turbine_efficiency(table,7.5,75.0)
        physical=electrical===nothing ? shaft : .9*shaft/(1-.004*shaft)
        @test termination_status(m)==OpenSHOP.MOI.OPTIMAL
        @test physical-1e-8<=objective_value(m)<20.0-1e-3
    end
end

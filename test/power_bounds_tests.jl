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
            (nothing,TableCurve([0.,20.],[.9,.98]),TableCurve([0.,20.],[-.1,.9]))
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
        m=Model()
        @variable(m,0<=q<=16)
        @variable(m,0<=p<=20)
        @variable(m,u,Bin)
        @variable(m,-60<=h<=150)
        b=(m=m,P=reshape([p],1,1),GQ=reshape([q],1,1),u=reshape([u],1,1),
            shared_heads=Dict((:plant,1)=>h))
        unsupported=electrical!==nothing && minimum(electrical.y)<0
        @test OpenSHOP._add_power_bounds!(b,c)==(unsupported ? 0 : 12)
        cuts=all_constraints(m;include_variable_in_set_constraints=false)
        if unsupported
            @test isempty(cuts)
            continue
        end
        function check(power,flow,head,on)
            point=Dict(p=>power,q=>flow,h=>head,u=>on)
            for cut in cuts
                object=constraint_object(cut)
                @test value(v->point[v],object.func)<=object.set.upper+1e-8
            end
        end
        @test OpenSHOP.turbine_efficiency(table,0.,-50.;extrapolation=:linear)<0
        for head in range(-60.,150.;length=15)
            check(0.,0.,head,0.)
        end
        for flow in range(5.,10.;length=11),head in range(10.,110.;length=15)
            eta=OpenSHOP.turbine_efficiency(table,flow,head;extrapolation=:linear)
            0<=eta<=1 || continue
            shaft=.00981*flow*head*eta
            power=electrical===nothing ? shaft : .9*shaft/(1-.004*shaft)
            check(power,flow,head,1.)
        end
    end
end

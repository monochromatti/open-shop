using Test, JuMP, OpenSHOP

function power_hull_integration_case(interpolation)
    table=interpolation===nothing ? nothing : TurbineTable(
        [1.0,50.0,100.0],[5.0,10.0,20.0],
        [0.15 0.20 0.25;0.85 0.90 0.95;0.75 0.80 0.85],
        fill(5.0,3),fill(20.0,3);interpolation,head_extrapolation=:linear)
    generators=[Generator(name=name,plant=plant,qmin=5.0,qmax=20.0,pmin=i==1 ? 1.0 : 1e-5,pmax=25.0,
        efficiency=0.9,qbest=10.0,qcurvature=i==1 ? 0.0 : 2.0,hbest=50.0,hcurvature=0.0,
        hmin=i==1 ? 1.0 : 0.001,hmax=105.0,initial_on=0,initial_age=8.0,minup=0.0,mindown=0.0,startup=0.0,
        turbine_table=table) for (i,(name,plant)) in enumerate(((:OnUnit,:OnPlant),(:OffUnit,:OffPlant)))]
    system=HydroSystem(
        reservoirs=[Reservoir(name=:Lake,z0=40.0,slope=20.0,v0=2.001,vmin=1.0,vmax=4.0,inflow=6.0,water_value=10.0)],
        junctions=Junction[],boundaries=[Boundary(name=:LowerTail,head=20.0),Boundary(name=:UpperTail,head=80.0)],
        tunnels=Tunnel[],rivers=River[],
        plants=[Plant(name=:OnPlant,source=:Lake,target=:LowerTail,pmax=25.0),
                Plant(name=:OffPlant,source=:Lake,target=:UpperTail,pmax=25.0)],
        generators=generators)
    ScheduleCase(name="power_hull_$(interpolation)",system=system,grid=[0.0,1.0],prices=[100.0])
end

@testset "Power hull integrates with complete physical starts" begin
    for interpolation in (nothing,:bilinear,:pchip_discharge)
        c=power_hull_integration_case(interpolation)
        seed=dispatch_from_controls(c,reshape([1,0],2,1),reshape([8.0,0.0],2,1),zeros(0,1))
        @test seed["validation"]["valid"]
        for fixed in (false,true)
            b=OpenSHOP._build_global_dispatch(c;joint=true,fixed_u=fixed ? seed["u"] : nothing)
            @test haskey(b.m.ext,:global_power_hulls)
            records=b.m.ext[:global_power_hulls]
            @test length(records)==(fixed ? 1 : 2)
            audit=OpenSHOP._lift_start!(b,c,seed)
            @test audit["valid"]
            @test audit["assigned"]==audit["variables"]==num_variables(b.m)
            @test audit["objective"]≈seed["objective"] atol=1e-6
            @test all(v->start_value(v)!==nothing,all_variables(b.m))
            for record in values(records)
                mass=record.u isa Number ? record.u : start_value(record.u)
                @test sum(start_value,record.weights)≈mass atol=1e-12
                @test start_value(record.onhead)≈mass*start_value(record.head) atol=1e-12
                @test start_value(record.oneta)≈mass*start_value(record.eta) atol=1e-12
            end
            if !fixed
                off=records["power_hull_2_1"]
                # The off plant has negative actual head and negative efficiency
                # continuation; neither value belongs to the on-state box.
                @test start_value(off.head)<0 && lower_bound(off.head)<0
                @test start_value(off.eta)<0 && lower_bound(off.eta)<0
                @test start_value(off.onhead)==start_value(off.oneta)==0.0
                @test all(w->start_value(w)==0.0,off.weights)
            end
        end
    end
end

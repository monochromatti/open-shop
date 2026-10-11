function tunnel_loss_case(; curve=TableCurve([0.0, 0.5, 1.0], [8.0, 3.0, 0.0]),
    openings=[0.25, 0.5, 0.0], boundary_head=0.0)
    lake=Reservoir(name=:GateLake, z0=90.0, slope=1.0, v0=10.0,
        vmin=8.0, vmax=12.0, water_value=10.0)
    tunnel=Tunnel(name=:GatedTunnel, source=:GateLake, target=:GateTail,
        resistance=1.0, capacity=20.0, opening=0.5, gate_loss_factor_curve=curve)
    system=HydroSystem(reservoirs=[lake], junctions=Junction[],
        boundaries=[Boundary(name=:GateTail, head=boundary_head)], tunnels=[tunnel],
        plants=Plant[], generators=Generator[], rivers=River[])
    ScheduleCase(name="prescribed gate loss factors", system=system,
        grid=collect(0.0:length(openings)), prices=zeros(length(openings)),
        operations=[OperationalSeries(object=:GatedTunnel, attribute=:opening,
            times=collect(0.0:(length(openings)-1)), values=Float64.(openings))])
end

# Solve R*q*|q| = H0 - alpha*q/2 analytically; closure is a separate state.
function tunnel_loss_expected(; openings=[0.25, 0.5, 0.0], boundary_head=0.0)
    volume=10.0
    flows=Float64[]
    volumes=[volume]
    for opening in openings
        if opening==0
            flow=0.0
        else
            resistance=opening<=0.5 ? 9.0-10.0opening : 7.0-6.0opening
            delta=90.0+volume-boundary_head
            flow=sign(delta)*(-0.0018+sqrt(0.0018^2+4resistance*abs(delta)))/(2resistance)
        end
        volume-=0.0036flow
        push!(flows,flow)
        push!(volumes,volume)
    end
    flows,volumes
end

@testset "Tunnel gate losses, reversal, closure and input round trip" begin
    for tail in (0.0,200.0)
        c=tunnel_loss_case(; boundary_head=tail)
        @test OpenSHOP.validate_inputs(c)
        @test case_dict(case_from_dict(case_dict(c)))==case_dict(c)
        x=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3))
        expected,storage=tunnel_loss_expected(; boundary_head=tail)
        @test x["validation"]["valid"]
        @test vec(x["tunnel_q"])≈expected atol=1e-8
        @test vec(x["V"])≈storage atol=1e-10
        @test x["tunnel_q"][1,3]≈0.0 atol=1e-9
        @test replay_audit(c,x)["valid"]
        continued=restart_case(c,x,1.0)
        remainder=dispatch_from_controls(continued,zeros(Int,0,2),zeros(0,2),zeros(0,2))
        @test remainder["validation"]["valid"]
        @test remainder["tunnel_q"]≈x["tunnel_q"][:,2:end] atol=1e-8
        @test remainder["V"]≈x["V"][:,2:end] atol=1e-10
        bad=deepcopy(x)
        bad["tunnel_q"][1,2]*=0.5
        @test !validate(c,bad)["valid"]
    end
    for curve in (TableCurve([0.1,1.0],[1.0,0.0]),
        TableCurve([0.0,0.9],[1.0,0.0]), TableCurve([0.0,1.0],[1.0,-0.1]))
        @test_throws ArgumentError OpenSHOP.validate_inputs(tunnel_loss_case(;curve))
    end
    # Legacy opening scaling is preserved when no gate-loss data is supplied.
    tunnel=only(tunnel_loss_case(;curve=nothing).system.tunnels)
    @test OpenSHOP.tunnel_flow(tunnel,100.0,0.25)==5.0
    @test OpenSHOP.tunnel_flow(tunnel,-100.0,0.25)==-5.0
    @test OpenSHOP.tunnel_flow(tunnel,100.0,0.0)==0.0
    supplied=only(tunnel_loss_case().system.tunnels)
    @test OpenSHOP.tunnel_flow(supplied,100.0,0.5)==5.0
    @test OpenSHOP.tunnel_flow(supplied,100.0,0.0)==0.0
end

@testset "Tunnel gate losses in local and SCIP optimization" begin
    for tail in (0.0,200.0)
        c=tunnel_loss_case(;boundary_head=tail)
        expected,storage=tunnel_loss_expected(;boundary_head=tail)
        initial=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3))
        local_result=solve_case(c;warm=initial,time_limit=10.0)
        @test local_result["validation"]["valid"]
        @test vec(local_result["tunnel_q"])≈expected atol=1e-6
        @test vec(local_result["V"])≈storage atol=1e-8
        b=OpenSHOP._build_global_dispatch(c;warm=initial)
        JuMP.set_optimizer(b.m,OpenSHOP.SCIP.Optimizer)
        JuMP.set_silent(b.m)
        JuMP.set_optimizer_attribute(b.m,"limits/time",15.0)
        JuMP.set_optimizer_attribute(b.m,"numerics/feastol",1e-8)
        JuMP.optimize!(b.m)
        @test JuMP.has_values(b.m)
        if JuMP.has_values(b.m)
            raw=OpenSHOP._dispatch_values(b)
            raw["u"]=zeros(Int,0,3)
            @test validate(c,raw)["valid"]
            @test vec(raw["tunnel_q"])≈expected atol=1e-5
            @test vec(raw["V"])≈storage atol=1e-7
        end
    end
end

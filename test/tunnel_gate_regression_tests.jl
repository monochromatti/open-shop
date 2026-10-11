@testset "Discrete gate extraction preserves mode identity" begin
    c=tg_flow_case(positions=[0.0,1e-12,1.0],opening=0.0,
        curve=TableCurve([0.0,1.0],[8.0,0.0]))
    data=OpenSHOP._tunnel_gate_states!(Model(),c)
    assigned=Dict(w[k]=>4e-7 for w in values(data.weights) for k in 2:length(w))
    evaluate(x)=JuMP.value(v->assigned[v],x)
    # Both selectors represent zero within solver tolerance. The opening
    # expression is nearer the positive 1e-12 mode, but the state is closed.
    @test evaluate(data.opening[1,1])>1e-12
    @test OpenSHOP._tunnel_gate_solution(data,evaluate)==zeros(1,3)
    for weights in values(data.weights)
        assigned[weights[2]]=1-4e-7
        assigned[weights[3]]=4e-7
    end
    @test OpenSHOP._tunnel_gate_solution(data,evaluate)==fill(1e-12,1,3)
    x=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3);
        tunnel_opening=fill(1e-12,1,3))
    @test x["validation"]["valid"]
    @test all(>(3.0),x["tunnel_q"])
    for weights in values(data.weights)
        assigned[weights[2]]=4e-7
        assigned[weights[3]]=1-4e-7
    end
    @test OpenSHOP._tunnel_gate_solution(data,evaluate)==ones(1,3)
    assigned[first(data.weights)[2][2]]=0.25
    @test_throws ArgumentError OpenSHOP._tunnel_gate_solution(data,evaluate)
    for weights in values(data.weights),k in 2:length(weights)
        assigned[weights[k]]=1.0
    end
    @test_throws ArgumentError OpenSHOP._tunnel_gate_solution(data,evaluate)
    @test_throws ArgumentError OpenSHOP._tunnel_gate_solution(data,evaluate;tolerance=0.5)

    ordinary=tg_flow_case(positions=[0.0,1.0],opening=0.0)
    @test_throws ArgumentError validate_tunnel_openings(ordinary,fill(1e-9,1,3))
    @test_throws ArgumentError OpenSHOP._tunnel_gate_states!(Model(),ordinary;fixed=fill(1e-9,1,3))
    @test_throws ArgumentError solve_case(ordinary;tunnel_opening=fill(1e-9,1,3))
    @test_throws ArgumentError simulate(ordinary,zeros(0,3),zeros(0,3);tunnel_opening=fill(1e-9,1,3))
    fixed=OpenSHOP._tunnel_gate_states!(Model(),c;fixed=fill(1e-12,1,3))
    @test fixed.open==ones(1,3)
    @test isempty(fixed.weights)
    valuable=tg_copy(c;system=tg_copy(c.system;
        reservoirs=[tg_copy(only(c.system.reservoirs);water_value=1000.0)]))
    closed=dispatch_from_controls(valuable,zeros(Int,0,3),zeros(0,3),zeros(0,3);
        tunnel_opening=zeros(1,3))
    joint=OpenSHOP.solve(valuable;initial=closed,time_limit=20.0,relative_gap=1e-4)
    @test joint["accepted"]
    @test joint["global_certificate"]
    @test joint["solution"]["tunnel_opening"]==zeros(1,3)
    @test joint["solution"]["objective"]==0.0
    prescribed=tg_flow_case(positions=Float64[],opening=0.1)
    prescribed=tg_copy(prescribed;grid=[0.0,0.3,0.7,1.0])
    literal=dispatch_from_controls(prescribed,zeros(Int,0,3),zeros(0,3),zeros(0,3);
        tunnel_opening=fill(0.1,1,3))
    @test literal["validation"]["valid"]
    @test replay_audit(prescribed,literal)["valid"]
    shut=tg_copy(prescribed;system=tg_copy(prescribed.system;
        tunnels=[tg_copy(only(prescribed.system.tunnels);opening=0.0)]))
    @test_throws ArgumentError validate_tunnel_openings(shut,fill(1e-9,1,3))
    @test_throws ArgumentError OpenSHOP._tunnel_gate_states!(Model(),shut;fixed=fill(1e-9,1,3))
end

function floating_gate_case(;junctions=[Junction(name=:Island,hmin=200.0,hmax=250.0)],
    tunnels=nothing,grid=[0.0,1.0])
    if tunnels===nothing
        tunnels=[
            Tunnel(name=:ClosedIntake,source=:Lake,target=first(junctions).name,
                resistance=1.0,capacity=10.0,opening=0.0,gate_positions=[0.0,1.0]),
            Tunnel(name=:ClosedOutlet,source=last(junctions).name,target=:Tail,
                resistance=1.0,capacity=10.0,opening=0.0,gate_positions=[0.0,1.0]),
        ]
    end
    s=HydroSystem(;reservoirs=[Reservoir(name=:Lake,z0=100.0,slope=0.1,
        v0=10.0,vmin=0.0,vmax=20.0,water_value=0.0)],junctions,
        boundaries=[Boundary(name=:Tail,head=0.0)],tunnels,
        plants=Plant[],generators=Generator[],rivers=River[])
    ScheduleCase(;name="floating pressure component",system=s,grid,prices=zeros(length(grid)-1))
end

@testset "Closed gates preserve floating pressure references and bounds" begin
    c=floating_gate_case()
    @test OpenSHOP.validate_inputs(c)
    opening=zeros(2,1)
    x=dispatch_from_controls(c,zeros(Int,0,1),zeros(0,1),zeros(0,1);tunnel_opening=opening)
    @test x["validation"]["valid"]
    @test x["H"][2,1]==225.0
    @test x["V"]==fill(10.0,1,2)
    @test x["tunnel_q"]==zeros(2,1)
    @test replay_audit(c,x)["valid"]
    for reference in (200.0,250.0)
        shifted=deepcopy(x);shifted["H"][2,1]=reference
        @test validate(c,shifted)["valid"]
        @test replay_audit(c,shifted)["valid"]
    end
    bad=deepcopy(x);bad["H"][2,1]=199.0
    @test !validate(c,bad)["valid"]
    result=OpenSHOP.solve(c;fixed_tunnel_opening=opening,time_limit=20.0)
    @test result["accepted"]
    @test result["global_certificate"]
    @test result["raw_solver_solution"]["validation"]["valid"]
    @test result["solution"]["H"][2,1]==225.0
    local_result=solve_case(c;tunnel_opening=opening,time_limit=10.0)
    @test local_result["validation"]["valid"]

    pair=floating_gate_case(junctions=[Junction(name=:Left,hmin=200.0,hmax=240.0),
        Junction(name=:Right,hmin=220.0,hmax=260.0)])
    bridge=Tunnel(name=:Bridge,source=:Left,target=:Right,resistance=1.0,capacity=10.0)
    pair=tg_copy(pair;system=tg_copy(pair.system;tunnels=[pair.system.tunnels;bridge]))
    opening=reshape([0.0,0.0,1.0],3,1)
    y=dispatch_from_controls(pair,zeros(Int,0,1),zeros(0,1),zeros(0,1);tunnel_opening=opening)
    @test y["validation"]["valid"]
    @test y["H"][2:3,1]≈[230.0,230.0] atol=1e-9
    @test y["tunnel_q"]==zeros(3,1)
    shifted=deepcopy(y);shifted["H"][2:3,1].+=5.0
    @test validate(pair,shifted)["valid"]
    relative=deepcopy(y);relative["H"][2,1]+=1.0
    @test !validate(pair,relative)["valid"]
    @test replay_audit(pair,y)["valid"]
    for ref in (-100.0,1000.0)
        guess=Dict("Vnew"=>[10.0],"H"=>[101.0,ref,ref+5.0,0.0],
            "tunnel_q"=>[1.0,-1.0,2.0],"river_release"=>Float64[])
        z=forward_step(pair.system,[10.0],1.0,Float64[],Float64[],Float64[];
            tunnel_opening=vec(opening),initial_guess=guess)
        @test z["converged"]
        @test z["H"][2:3]≈[230.0,230.0] atol=1e-8
        @test maximum(abs,z["tunnel_q"])<1e-8
    end
    incompatible=tg_copy(pair;system=tg_copy(pair.system;
        junctions=[Junction(name=:Left,hmin=200.0,hmax=210.0),
            Junction(name=:Right,hmin=220.0,hmax=230.0)]))
    invalid=dispatch_from_controls(incompatible,zeros(Int,0,1),zeros(0,1),zeros(0,1);
        tunnel_opening=opening)
    @test !invalid["validation"]["valid"]

    islands=floating_gate_case(junctions=[Junction(name=:Left,hmin=200.0,hmax=240.0),
        Junction(name=:Right,hmin=300.0,hmax=400.0)])
    @test OpenSHOP._floating_pressure_components(islands.system,[0.0,0.0],Float64[])==[[2],[3]]
    separate=dispatch_from_controls(islands,zeros(Int,0,1),zeros(0,1),zeros(0,1);
        tunnel_opening=zeros(2,1))
    @test separate["validation"]["valid"]
    @test separate["H"][2:3,1]==[220.0,350.0]
    independent=deepcopy(separate);independent["H"][2:3,1]=[235.0,305.0]
    @test validate(islands,independent)["valid"]
end

@testset "Pressure references follow closure and reconnection" begin
    c=floating_gate_case(junctions=[Junction(name=:Island,hmin=0.0,hmax=150.0)],
        grid=[0.0,1.0,2.0,3.0])
    opening=[0.0 1.0 0.0;0.0 0.0 0.0]
    x=dispatch_from_controls(c,zeros(Int,0,3),zeros(0,3),zeros(0,3);tunnel_opening=opening)
    @test x["validation"]["valid"]
    @test vec(x["H"][2,:])≈[75.0,101.0,75.0] atol=1e-9
    @test replay_audit(c,x)["valid"]
    shifted=deepcopy(x);shifted["H"][2,1]=140.0;shifted["H"][2,3]=10.0
    @test validate(c,shifted)["valid"]
    shifted["H"][2,2]=110.0
    @test !validate(c,shifted)["valid"]
    restart=restart_case(c,x,1.0)
    continuation=dispatch_from_controls(restart,zeros(Int,0,2),zeros(0,2),zeros(0,2);
        tunnel_opening=opening[:,2:3])
    @test continuation["validation"]["valid"]
    @test continuation["H"]≈x["H"][:,2:3] atol=1e-9
end

@testset "Pressure reference removal retains all water balances" begin
    c=floating_gate_case(junctions=[Junction(name=:Left,hmin=200.0,hmax=240.0),
        Junction(name=:Right,hmin=200.0,hmax=240.0)])
    river=River(name=:Outfall,source=:Mix,target=:Tail,capacity=10.0,
        curves=RiverRouting.DelayCurve[],deterministic_delay=0.0,water_value=0.0)
    bridge=Tunnel(name=:Bridge,source=:Left,target=:Right,resistance=1.0,capacity=10.0)
    unit=Generator(name=:Unit,plant=:Station,qmin=1.0,qmax=6.0,pmin=0.1,pmax=10.0,
        hmin=0.0,hmax=300.0,qbest=3.0,hbest=50.0)
    system=tg_copy(c.system;tunnels=[c.system.tunnels;bridge],
        plants=[Plant(name=:Station,source=:Left,target=:Right,pmax=10.0,
            discharge_river=:Outfall)],generators=[unit],rivers=[river],
        river_junctions=[RiverJunction(name=:Mix)])
    opening=[0.0,0.0,1.0]
    @test OpenSHOP._floating_pressure_components(system,opening,[5.0])==[[2,3]]
    # Five m³/s leaves an isolated component with no supply. Least squares may
    # produce a trial iterate, but the full original residual must reject it.
    z=forward_step(system,[10.0],1.0,[5.0],[0.0],[0.0];tunnel_opening=opening,
        current_transfer=ones(1,2))
    @test !z["converged"]
    @test z["residual"]>=2.5-1e-6
    # An active plant joining separate pressure components has observable
    # head differences and must not acquire independent pressure references.
    separate=tg_copy(system;tunnels=c.system.tunnels)
    @test isempty(OpenSHOP._floating_pressure_components(separate,[0.0,0.0],[5.0]))
    floor=tg_copy(system;plants=[tg_copy(only(system.plants);outlet_head_floor=210.0)])
    @test isempty(OpenSHOP._floating_pressure_components(floor,opening,[5.0]))
end

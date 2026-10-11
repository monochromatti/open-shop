include(joinpath(@__DIR__,"..","examples","tunnel_gates.jl"))
using .TunnelGatesExample: tunnel_gate_case

tg_copy(x;kwargs...)=typeof(x)(;merge((; (k=>getfield(x,k) for k in fieldnames(typeof(x)))...),
    (;kwargs...))...)

function tg_flow_case(;positions=[0.0,0.5,1.0],opening=0.5,initial=nothing,cost=0.0,
    curve=TableCurve([0.0,0.5,1.0],[8.0,3.0,0.0]),tail=0.0,
    interlocks=TunnelInterlock[],operations=OperationalSeries[])
    lake=Reservoir(name=:GateLake,z0=90.0,slope=1.0,v0=10.0,
        vmin=8.0,vmax=12.0,water_value=0.0)
    tunnel=Tunnel(;name=:GateLink,source=:GateLake,target=:GateTail,resistance=1.0,
        capacity=20.0,opening,gate_positions=Float64.(positions),initial_opening=initial,
        gate_adjustment_cost=cost,gate_loss_factor_curve=curve)
    system=HydroSystem(reservoirs=[lake],junctions=Junction[],
        boundaries=[Boundary(name=:GateTail,head=tail)],tunnels=[tunnel],plants=Plant[],
        generators=Generator[],rivers=River[],tunnel_interlocks=interlocks)
    ScheduleCase(;name="discrete pressure gate",system,grid=[0.0,1.0,2.0,3.0],
        prices=zeros(3),operations)
end

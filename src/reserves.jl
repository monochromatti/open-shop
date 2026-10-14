# Spinning generation capacity at the scheduled hydraulic node heads.
# Capacity is a commitment, not an activated flow in the nominal water balance.

_reserve_input(x::ReserveAllocation,T) = x
function _reserve_array(x)
    converted=Array{Float64,3}(x)
    all(isfinite,converted) &&
    all(i->iszero(x[i])==iszero(converted[i]),eachindex(converted)) ||
        throw(ArgumentError("reserve controls must remain finite in Float64 without changing zero states"))
    converted
end
function _reserve_input(x::AbstractDict,T)
    cubes(key)=if x[key] isa AbstractArray{<:Real,3}
        _reserve_array(x[key])
    else
        groups=[_schedule_matrix(group,T;control=true) for group in x[key]]
        permutedims(cat(groups...;dims=3),(1,3,2))
    end
    ReserveAllocation(cubes("up"),cubes("down"),_schedule_matrix(x["up_q"],T),
        _schedule_matrix(x["down_q"],T),_schedule_matrix(x["up_shortfall"],T),
        _schedule_matrix(x["down_shortfall"],T))
end
function jsonready(x::ReserveAllocation)
    jsonready(Dict("up"=>[x.up[:,k,:] for k in axes(x.up,2)],
        "down"=>[x.down[:,k,:] for k in axes(x.down,2)],
        "up_q"=>x.up_q,"down_q"=>x.down_q,
        "up_shortfall"=>x.up_shortfall,"down_shortfall"=>x.down_shortfall))
end
function _reserve_slice(x::ReserveAllocation,indices)
    ReserveAllocation(x.up[:,:,indices],x.down[:,:,indices],x.up_q[:,indices],
        x.down_q[:,indices],x.up_shortfall[:,indices],x.down_shortfall[:,indices])
end
_reserve_slice(::Nothing,indices) = nothing

function reserve_members(c)
    index=Dict(g.name=>i for (i,g) in enumerate(c.system.generators))
    [[index[name] for name in group.generators] for group in c.reserve_groups]
end
function validate_reserve_inputs(c)
    physical=Set(x.name for collection in (c.system.reservoirs,c.system.junctions,
        c.system.boundaries,c.system.tunnels,c.system.tunnel_interlocks,c.system.plants,
        c.system.generators,c.system.commitment_exclusions,c.system.rivers,
        c.system.river_junctions,c.flow_requirements,c.soft_constraints)
        for x in collection)
    for group in c.reserve_groups
        group.name in physical && throw(ArgumentError("duplicate reserve group name $(group.name)"))
        push!(physical,group.name)
        !isempty(group.generators) && length(unique(group.generators))==length(group.generators) ||
            throw(ArgumentError("reserve group requires distinct generating units"))
        all(n->any(g.name==n && g.mode==:generate for g in c.system.generators),group.generators) ||
            throw(ArgumentError("spinning reserves support generating units only"))
        for field in (:up_requirement,:down_requirement,:activation_hours,:up_max,:down_max)
            value=getfield(group,field)
            value===nothing || (isfinite(value) && value>=0) ||
                throw(ArgumentError("nonnegative finite reserve $field required"))
        end
        isfinite(group.up_price) && isfinite(group.down_price) &&
        isfinite(group.response_hours) && group.response_hours>0 ||
            throw(ArgumentError("finite prices and positive reserve response hours required"))
        group.shortfall_cost===nothing ||
            (isfinite(group.shortfall_cost) && group.shortfall_cost>0) ||
            throw(ArgumentError("reserve shortfall cost must be positive or nothing"))
        if group.activation_hours>0
            for name in group.generators
                g=only(g for g in c.system.generators if g.name==name)
                g.reserve_up_energy!==nothing && g.reserve_down_energy!==nothing ||
                    throw(ArgumentError("sustained reserve capacity needs qualified unit MWh budgets"))
            end
        end
    end
    for g in c.system.generators, field in (:reserve_up_energy,:reserve_down_energy)
        value=getfield(g,field)
        value===nothing || (isfinite(value) && value>=0) ||
            throw(ArgumentError("qualified reserve energy must be nonnegative and finite"))
    end
    nothing
end

function reserve_economics(c,r::ReserveAllocation)
    revenue=0.0;cost=0.0
    for (k,group) in enumerate(c.reserve_groups),t in eachindex(c.prices)
        dt=c.grid[t+1]-c.grid[t]
        up=sum(@view r.up[:,k,t]);down=sum(@view r.down[:,k,t])
        revenue+=dt*(opinterval(c,group.name,:up_price,t,group.up_price)*up+
            opinterval(c,group.name,:down_price,t,group.down_price)*down)
        if group.shortfall_cost!==nothing
            up_shortfall=max(0.0,opinterval(c,group.name,:up_requirement,t,group.up_requirement)-up)
            down_shortfall=max(0.0,opinterval(c,group.name,:down_requirement,t,group.down_requirement)-down)
            cost+=dt*group.shortfall_cost*(up_shortfall+down_shortfall)
        end
    end
    (;revenue,cost,value=revenue-cost)
end
reserve_economics(c,::Nothing) = (;revenue=0.0,cost=0.0,value=0.0)

function _zero_reserves(c,q)
    isempty(c.reserve_groups) && return nothing
    G,K,T=length(c.system.generators),length(c.reserve_groups),length(c.prices)
    up=zeros(G,K,T);down=copy(up)
    su=[opinterval(c,g.name,:up_requirement,t,g.up_requirement) for g in c.reserve_groups,t in 1:T]
    sd=[opinterval(c,g.name,:down_requirement,t,g.down_requirement) for g in c.reserve_groups,t in 1:T]
    ReserveAllocation(up,down,copy(q),copy(q),su,sd)
end

# This operation contains only shared linear allocation, response and budget rows.
function _reserve_variables!(m,c,u,P,Q)
    isempty(c.reserve_groups) && return nothing
    G,K,T=length(c.system.generators),length(c.reserve_groups),length(c.prices)
    members=reserve_members(c)
    @variable(m,reserve_up[i=1:G,k=1:K,t=1:T],lower_bound=0,upper_bound=c.system.generators[i].pmax)
    @variable(m,reserve_down[i=1:G,k=1:K,t=1:T],lower_bound=0,upper_bound=c.system.generators[i].pmax)
    @variable(m,reserve_shortfall_up[1:K,1:T]>=0)
    @variable(m,reserve_shortfall_down[1:K,1:T]>=0)
    @variable(m,reserve_q[i=1:G,t=1:T,d=1:2],lower_bound=0,upper_bound=c.system.generators[i].qmax/50)
    @variable(m,reserve_p[i=1:G,t=1:T,d=1:2],lower_bound=0,upper_bound=c.system.generators[i].pmax/40)
    Qpoint=50 .* reserve_q;Ppoint=40 .* reserve_p
    for (k,group) in enumerate(c.reserve_groups),t in 1:T
        ceiling=sum(c.system.generators[i].pmax for i in members[k])
        for (direction,R,shortfall) in ((:up,reserve_up,reserve_shortfall_up),
            (:down,reserve_down,reserve_shortfall_down))
            requirement=opinterval(c,group.name,Symbol(direction,:_requirement),t,getfield(group,Symbol(direction,:_requirement)))
            cap=optional_operation(c,group.name,Symbol(direction,:_max),t;
                default=getfield(group,Symbol(direction,:_max)))
            cap===nothing && (cap=ceiling)
            set_upper_bound(shortfall[k,t],requirement)
            group.shortfall_cost===nothing && fix(shortfall[k,t],0;force=true)
            @constraint(m,sum(R[i,k,t] for i in 1:G)+shortfall[k,t]>=requirement)
            @constraint(m,sum(R[i,k,t] for i in 1:G)<=cap)
            for (i,g) in enumerate(c.system.generators)
                if !(i in members[k])
                    fix(R[i,k,t],0;force=true)
                else
                    available=opinterval(c,g.name,:reserve_available,t,g.reserve_available)
                    state=u[i,t] isa Number ? u[i,t] : is_fixed(u[i,t]) ? fix_value(u[i,t]) : nothing
                    if available==0 || state==0
                        fix(R[i,k,t],0;force=true)
                    else
                        @constraint(m,R[i,k,t]<=g.pmax*available*u[i,t])
                    end
                end
            end
        end
    end
    for (i,g) in enumerate(c.system.generators),t in 1:T
        up=sum(reserve_up[i,k,t] for k in 1:K)
        down=sum(reserve_down[i,k,t] for k in 1:K)
        @constraint(m,Ppoint[i,t,1]==P[i,t]+up)
        @constraint(m,Ppoint[i,t,2]==P[i,t]-down)
        @constraint(m,Qpoint[i,t,1]>=Q[i,t])
        @constraint(m,Qpoint[i,t,2]<=Q[i,t])
        if g.mode==:pump
            @constraint(m,Qpoint[i,t,1]==Q[i,t])
            @constraint(m,Qpoint[i,t,2]==Q[i,t])
        end
        plant_units=findall(z->z.plant==g.plant,c.system.generators)
        groups=findall(k->!isempty(intersect(plant_units,members[k])),1:K)
        if !isempty(groups)
            response=minimum(c.reserve_groups[k].response_hours for k in groups)
            a=(c.grid[t]+c.grid[t+1])/2
            for (direction,delta_p,delta_q) in ((:up,up,Qpoint[i,t,1]-Q[i,t]),
                (:down,down,Q[i,t]-Qpoint[i,t,2]))
                rate=ramp_allowance(c,g,Symbol(:ramp_,direction),a,a+response)
                flow=ramp_allowance(c,g,Symbol(:discharge_ramp_,direction),a,a+response)
                rate===nothing || @constraint(m,delta_p<=rate)
                flow===nothing || @constraint(m,delta_q<=flow)
                budget=getfield(g,Symbol(:reserve_,direction,:_energy))
                budget===nothing || @constraint(m,sum(c.reserve_groups[k].activation_hours*
                    (direction==:up ? reserve_up[i,k,t] : reserve_down[i,k,t]) for k in groups)<=budget)
            end
        end
        for d in 1:2
            @constraint(m,Qpoint[i,t,d]>=opinterval(c,g.name,:qmin,t,g.qmin)*u[i,t])
            @constraint(m,Qpoint[i,t,d]<=opinterval(c,g.name,:qmax,t,g.qmax)*u[i,t])
            @constraint(m,Ppoint[i,t,d]>=opinterval(c,g.name,:pmin,t,g.pmin)*u[i,t])
            @constraint(m,Ppoint[i,t,d]<=opinterval(c,g.name,:pmax,t,g.pmax)*u[i,t])
        end
    end
    for plant in c.system.plants,t in 1:T
        units=findall(g->g.plant==plant.name,c.system.generators)
        groups=findall(k->!isempty(intersect(units,members[k])),1:K)
        isempty(groups) && continue
        response=minimum(c.reserve_groups[k].response_hours for k in groups)
        a=(c.grid[t]+c.grid[t+1])/2
        for (d,direction) in ((1,:up),(2,:down))
            total_p=sum(Ppoint[i,t,d] for i in units);total_q=sum(Qpoint[i,t,d] for i in units)
            @constraint(m,total_p<=opinterval(c,plant.name,:pmax,t,plant.pmax))
            maxq=optional_operation(c,plant.name,:qmax,t;default=plant.qmax)
            maxq===nothing || @constraint(m,total_q<=maxq)
            for i in units
                @constraint(m,total_p>=opinterval(c,plant.name,:pmin,t,plant.pmin)*u[i,t])
                @constraint(m,total_q>=opinterval(c,plant.name,:qmin,t,plant.qmin)*u[i,t])
            end
            delta_p=d==1 ? total_p-sum(P[i,t] for i in units) : sum(P[i,t] for i in units)-total_p
            delta_q=d==1 ? total_q-sum(Q[i,t] for i in units) : sum(Q[i,t] for i in units)-total_q
            rate=ramp_allowance(c,plant,Symbol(:ramp_,direction),a,a+response)
            flow=ramp_allowance(c,plant,Symbol(:discharge_ramp_,direction),a,a+response)
            rate===nothing || @constraint(m,delta_p<=rate)
            flow===nothing || @constraint(m,delta_q<=flow)
        end
    end
    value=@expression(m,sum((c.grid[t+1]-c.grid[t])*(
        opinterval(c,group.name,:up_price,t,group.up_price)*sum(reserve_up[i,k,t] for i in 1:G)+
        opinterval(c,group.name,:down_price,t,group.down_price)*sum(reserve_down[i,k,t] for i in 1:G)-
        (group.shortfall_cost===nothing ? 0.0 : group.shortfall_cost)*
            (reserve_shortfall_up[k,t]+reserve_shortfall_down[k,t]))
        for (k,group) in enumerate(c.reserve_groups),t in 1:T))
    (;up=reserve_up,down=reserve_down,up_shortfall=reserve_shortfall_up,
        down_shortfall=reserve_shortfall_down,q=reserve_q,p=reserve_p,Qpoint,Ppoint,value)
end

function _local_reserve_points!(m,c,u,H,P,Q,turbine_ops,electrical_ops,tailwater_ops,
    outlet_ops,turbine_flow_ops)
    r=_reserve_variables!(m,c,u,P,Q)
    r===nothing && return nothing
    ix=nodeindex(c.system)
    for (i,g) in enumerate(c.system.generators),t in eachindex(c.prices),d in 1:2
        g.mode==:pump && continue
        plant=plantof(c.system,g)
        units=findall(z->z.plant==plant.name,c.system.generators)
        total=sum(r.Qpoint[j,t,d] for j in units)
        op=outlet_ops[plant.name]
        receiver=op===nothing ? H[ix[plant.target],t] : op isa Real ? op : op(H[ix[plant.target],t])
        tail=tailwater_ops[plant.name]===nothing ? 0.0 : tailwater_ops[plant.name](total)
        hd=H[ix[plant.source],t]-receiver-tail
        q,p=r.Qpoint[i,t,d],r.Ppoint[i,t,d]
        eta=turbine_ops[i]===nothing ? g.efficiency-g.qcurvature*((q-g.qbest)/g.qbest)^2-
            g.hcurvature*((hd-g.hbest)/g.hbest)^2 : turbine_ops[i](q,hd)
        electrical=electrical_ops[i]===nothing ? 1.0 : electrical_ops[i](p)
        if u[i,t]==0
            fix(r.q[i,t,d],0;force=true);fix(r.p[i,t,d],0;force=true)
        else
            @constraint(m,(p-0.00981*q*hd*eta*electrical)/40==0)
            @constraint(m,g.hmin<=hd<=g.hmax)
            @constraint(m,g.min_efficiency<=eta<=1.0)
            if turbine_flow_ops[i]!==nothing
                lo,hi=turbine_flow_ops[i]
                @constraint(m,q>=lo(hd));@constraint(m,q<=hi(hd))
            end
        end
    end
    (;r...,graphs=nothing)
end

function _global_reserve_points!(m,c,u,H,P,Q,node_bounds,plant_heads)
    r=_reserve_variables!(m,c,u,P,Q)
    r===nothing && return nothing
    G,T=length(c.system.generators),length(c.prices)
    graphs=Array{Union{Nothing,_GeneratorPowerGraph},3}(nothing,G,T,2)
    heads=Dict{Tuple{Symbol,Int,Int},_PlantHeadGraph}()
    for (i,g) in enumerate(c.system.generators),t in 1:T,d in 1:2
        g.mode==:pump && continue
        plant=plantof(c.system,g);key=(plant.name,t,d)
        if !haskey(heads,key)
            units=findall(z->z.plant==plant.name,c.system.generators)
            total=sum(r.Qpoint[j,t,d] for j in units)
            totalmax=sum(c.system.generators[j].qmax for j in units)
            tail=plant.tailwater_curve===nothing ? nothing :
                _global_tensor_table_graph!(m,plant.tailwater_curve,total,0.0,totalmax;
                    name=Symbol("reserve_tail_",i,"_",t,"_",d))
            nominal=plant_heads[(plant.name,t)]
            gross=nominal.head+(nominal.tailwater===nothing ? 0.0 : nominal.tailwater.y)
            src,dst=node_bounds[(plant.source,t)],node_bounds[(plant.target,t)]
            floor=plant.outlet_head_floor
            floor===nothing || (dst=(max(dst[1],floor),max(dst[2],floor)))
            tailbounds=plant.tailwater_curve===nothing ? (0.0,0.0) : extrema(plant.tailwater_curve.y)
            lo,hi=src[1]-dst[2]-tailbounds[2],src[2]-dst[1]-tailbounds[1]
            hd=@variable(m,lower_bound=lo,upper_bound=hi,base_name="reserve_head_$(i)_$(t)_$(d)")
            @constraint(m,hd==gross-(tail===nothing ? 0.0 : tail.y))
            heads[key]=_PlantHeadGraph(hd,nothing,tail)
        end
        headgraph=heads[key];hd=headgraph.head
        hlo,hhi=lower_bound(hd),upper_bound(hd)
        q,p=r.Qpoint[i,t,d],r.Ppoint[i,t,d]
        known=u[i,t] isa Number ? u[i,t] : is_fixed(u[i,t]) ? fix_value(u[i,t]) : nothing
        qmax=known==0 ? 0.0 : g.qmax
        known==0 && (fix(r.q[i,t,d],0;force=true);fix(r.p[i,t,d],0;force=true))
        turbine=nothing
        eta=if g.turbine_table===nothing
            lo,hi=_global_analytic_eta_bounds(g,0.0,qmax,hlo,hhi)
            z=@variable(m,lower_bound=lo,upper_bound=hi,base_name="reserve_eta_$(i)_$(t)_$(d)")
            @constraint(m,z==g.efficiency-g.qcurvature*((q-g.qbest)/g.qbest)^2-
                g.hcurvature*((hd-g.hbest)/g.hbest)^2)
            z
        else
            turbine=_global_tensor_turbine_graph!(m,g.turbine_table,q,hd,0.0,qmax,hlo,hhi;
                name=Symbol("reserve_turbine_",i,"_",t,"_",d))
            turbine.eta
        end
        electrical=g.generator_efficiency_curve===nothing ? nothing :
            _global_tensor_table_graph!(m,g.generator_efficiency_curve,p,0.0,g.pmax;
                name=Symbol("reserve_electrical_",i,"_",t,"_",d))
        @constraint(m,(p-0.00981*q*hd*eta*(electrical===nothing ? 1.0 : electrical.y))/40==0)
        @constraint(m,hd>=g.hmin-max(0.0,g.hmin-hlo)*(1-u[i,t]))
        @constraint(m,hd<=g.hmax+max(0.0,hhi-g.hmax)*(1-u[i,t]))
        @constraint(m,eta>=g.min_efficiency-max(0.0,g.min_efficiency-lower_bound(eta))*(1-u[i,t]))
        @constraint(m,eta<=1+max(0.0,upper_bound(eta)-1)*(1-u[i,t]))
        qlo=qhi=nothing
        if g.turbine_table!==nothing
            table=g.turbine_table
            qlo=_global_tensor_table_graph!(m,table.heads,
                h->turbine_qmin(table,h;extrapolation=:linear),hd,hlo,hhi;
                name=Symbol("reserve_qlo_",i,"_",t,"_",d))
            qhi=_global_tensor_table_graph!(m,table.heads,
                h->turbine_qmax(table,h;extrapolation=:linear),hd,hlo,hhi;
                name=Symbol("reserve_qhi_",i,"_",t,"_",d))
            @constraint(m,q>=qlo.y-max(0.0,upper_bound(qlo.y))*(1-u[i,t]))
            @constraint(m,q<=qhi.y+max(0.0,g.qmax-lower_bound(qhi.y))*(1-u[i,t]))
        end
        graphs[i,t,d]=_GeneratorPowerGraph(headgraph,eta,turbine,qlo,qhi,electrical)
    end
    (;r...,graphs)
end

function _reserve_values(c,r)
    r===nothing && return nothing
    up,down=value.(r.up),value.(r.down)
    su=[max(0.0,opinterval(c,g.name,:up_requirement,t,g.up_requirement)-sum(up[:,k,t]))
        for (k,g) in enumerate(c.reserve_groups),t in eachindex(c.prices)]
    sd=[max(0.0,opinterval(c,g.name,:down_requirement,t,g.down_requirement)-sum(down[:,k,t]))
        for (k,g) in enumerate(c.reserve_groups),t in eachindex(c.prices)]
    ReserveAllocation(up,down,value.(r.Qpoint[:,:,1]),value.(r.Qpoint[:,:,2]),su,sd)
end

function _start_reserves!(r,c,warm)
    r===nothing && return
    physical=warm===nothing ? nothing : _schedule_input(warm;case=c)
    allocation=physical===nothing ? nothing : physical.reserves
    for i in eachindex(c.system.generators),t in eachindex(c.prices),d in 1:2
        g=c.system.generators[i]
        q=allocation===nothing ? (physical===nothing ? 0.6g.qmax : physical.generator_q[i,t]) :
            (d==1 ? allocation.up_q[i,t] : allocation.down_q[i,t])
        set_start_value(r.q[i,t,d],q/50)
        set_start_value(r.p[i,t,d],physical===nothing ? power(g,q,g.hbest)/40 :
            (physical.power[i,t]+(allocation===nothing ? 0.0 :
                d==1 ? sum(allocation.up[i,:,t]) : -sum(allocation.down[i,:,t])))/40)
    end
    for i in axes(r.up,1),k in axes(r.up,2),t in axes(r.up,3)
        set_start_value(r.up[i,k,t],allocation===nothing ? 0.0 : allocation.up[i,k,t])
        set_start_value(r.down[i,k,t],allocation===nothing ? 0.0 : allocation.down[i,k,t])
    end
    for k in axes(r.up_shortfall,1),t in axes(r.up_shortfall,2)
        set_start_value(r.up_shortfall[k,t],allocation===nothing ? 0.0 : allocation.up_shortfall[k,t])
        set_start_value(r.down_shortfall[k,t],allocation===nothing ? 0.0 : allocation.down_shortfall[k,t])
    end
end

function _lift_reserves!(b,c,physical,put,table)
    r=b.reserves
    r===nothing && return
    allocation=physical.reserves
    allocation===nothing && error("reserve incumbent needs explicit allocation")
    for i in axes(r.up,1),k in axes(r.up,2),t in axes(r.up,3)
        put(r.up[i,k,t],allocation.up[i,k,t]);put(r.down[i,k,t],allocation.down[i,k,t])
    end
    for k in axes(r.up_shortfall,1),t in axes(r.up_shortfall,2)
        put(r.up_shortfall[k,t],allocation.up_shortfall[k,t])
        put(r.down_shortfall[k,t],allocation.down_shortfall[k,t])
    end
    ix=nodeindex(c.system)
    for (i,g) in enumerate(c.system.generators),t in eachindex(c.prices),d in 1:2
        q=d==1 ? allocation.up_q[i,t] : allocation.down_q[i,t]
        p=physical.power[i,t]+(d==1 ? sum(allocation.up[i,:,t]) : -sum(allocation.down[i,:,t]))
        plant=plantof(c.system,g)
        flow=d==1 ? allocation.up_q : allocation.down_q
        total=sum(flow[j,t] for (j,z) in enumerate(c.system.generators) if z.plant==plant.name)
        h=physical.head[ix[plant.source],t]-outlet_head(plant,physical.head[ix[plant.target],t])-tailwater(plant,total)
        graph=r.graphs[i,t,d]
        put(r.q[i,t,d],q/50);put(r.p[i,t,d],p/40)
        graph===nothing && continue
        put(graph.head.head,h)
        table(graph.head.tailwater,plant.tailwater_curve,total)
        if graph.turbine===nothing
            put(graph.efficiency,efficiency(g,q,h))
        else
            for (variable,val) in _global_tensor_turbine_values(graph.turbine,q,h)
                put(variable,val)
            end
            put(graph.efficiency,turbine_efficiency(g.turbine_table,q,h;extrapolation=:linear))
        end
        table(graph.flow_min,z->turbine_qmin(g.turbine_table,z;extrapolation=:linear),h)
        table(graph.flow_max,z->turbine_qmax(g.turbine_table,z;extrapolation=:linear),h)
        table(graph.electrical,g.generator_efficiency_curve,p)
    end
end

"""Independently check spinning capacity and shortfall, apart from nominal energy.

This checks unit/plant operating points at scheduled node heads and qualified
response/energy limits. It does not certify activated tunnel or river trajectories.
"""
function validate_reserves(c,x::Schedule;tolerance=1e-4)
    residuals=Dict{String,Float64}();errors=String[]
    record(name,value)=begin
        residuals[name]=max(0.0,value)
        (!isfinite(value) || value>tolerance) && push!(errors,"$name reserve residual $value")
    end
    finish()=Dict("valid"=>isempty(errors),"residuals"=>residuals,"errors"=>errors,
        "scope"=>"spinning capacity at scheduled hydraulic heads; qualified activation energy only")
    if isempty(c.reserve_groups)
        x.reserves===nothing || push!(errors,"reserve allocation requires reserve groups")
        return finish()
    end
    r=x.reserves
    if r===nothing
        push!(errors,"missing reserve allocation");return finish()
    end
    G,K,T=length(c.system.generators),length(c.reserve_groups),length(c.prices)
    for (A,shape) in ((r.up,(G,K,T)),(r.down,(G,K,T)),(r.up_q,(G,T)),
        (r.down_q,(G,T)),(r.up_shortfall,(K,T)),(r.down_shortfall,(K,T)))
        if size(A)!=shape || !all(isfinite,A)
            push!(errors,"malformed reserve allocation");return finish()
        end
    end
    members=reserve_members(c);ix=nodeindex(c.system)
    for (k,group) in enumerate(c.reserve_groups),t in 1:T
        for (direction,R,S) in ((:up,r.up,r.up_shortfall),(:down,r.down,r.down_shortfall))
            requirement=opinterval(c,group.name,Symbol(direction,:_requirement),t,getfield(group,Symbol(direction,:_requirement)))
            total=sum(R[:,k,t]);shortfall=max(0.0,requirement-total)
            cap=optional_operation(c,group.name,Symbol(direction,:_max),t;
                default=getfield(group,Symbol(direction,:_max)))
            cap===nothing && (cap=sum(c.system.generators[i].pmax for i in members[k]))
            record("$(group.name)_$(direction)_capacity",total-cap)
            record("$(group.name)_$(direction)_shortfall",abs(S[k,t]-shortfall))
            group.shortfall_cost===nothing && record("$(group.name)_$(direction)_obligation",shortfall)
            record("$(group.name)_$(direction)_slack_bound",max(-S[k,t],S[k,t]-requirement))
            for (i,g) in enumerate(c.system.generators)
                record("$(g.name)_$(group.name)_$(direction)_nonnegative",-R[i,k,t])
                available=opinterval(c,g.name,:reserve_available,t,g.reserve_available)
                ceiling=i in members[k] ? g.pmax*available*x.commitment[i,t] : 0.0
                record("$(g.name)_$(group.name)_$(direction)_eligibility",R[i,k,t]-ceiling)
            end
        end
    end
    for t in 1:T,d in 1:2
        direction=d==1 ? :up : :down
        R=d==1 ? r.up : r.down;Q=d==1 ? r.up_q : r.down_q
        actual=zeros(G)
        for (i,g) in enumerate(c.system.generators)
            if g.mode==:pump
                actual[i]=x.power[i,t]
                record("$(g.name)_$(direction)_witness",abs(Q[i,t]-x.generator_q[i,t]))
                continue
            end
            plant=plantof(c.system,g)
            units=findall(z->z.plant==plant.name,c.system.generators)
            totalq=sum(Q[j,t] for j in units)
            h=x.head[ix[plant.source],t]-outlet_head(plant,x.head[ix[plant.target],t])-tailwater(plant,totalq)
            q=Q[i,t];actual[i]=power(g,q,h)
            amount=sum(R[i,:,t]);u=x.commitment[i,t]
            qlo,qhi=opinterval(c,g.name,:qmin,t,g.qmin),opinterval(c,g.name,:qmax,t,g.qmax)
            plo,phi=opinterval(c,g.name,:pmin,t,g.pmin),opinterval(c,g.name,:pmax,t,g.pmax)
            record("$(g.name)_$(direction)_flow",max(qlo*u-q,q-qhi*u))
            record("$(g.name)_$(direction)_power",max(plo*u-actual[i],actual[i]-phi*u))
            record("$(g.name)_$(direction)_witness",abs(actual[i]-(x.power[i,t]+(d==1 ? amount : -amount))))
            record("$(g.name)_$(direction)_flow_direction",d==1 ? x.generator_q[i,t]-q : q-x.generator_q[i,t])
            if u==1
                eta=efficiency(g,q,h)
                record("$(g.name)_$(direction)_head",max(g.hmin-h,h-g.hmax))
                record("$(g.name)_$(direction)_efficiency",max(g.min_efficiency-eta,eta-1))
                g.turbine_table===nothing || record("$(g.name)_$(direction)_turbine",
                    max(turbine_qmin(g.turbine_table,h;extrapolation=:linear)-q,
                        q-turbine_qmax(g.turbine_table,h;extrapolation=:linear)))
            end
            groups=findall(k->!isempty(intersect(units,members[k])),1:K)
            if !isempty(groups)
                response=minimum(c.reserve_groups[k].response_hours for k in groups)
                a=(c.grid[t]+c.grid[t+1])/2
                rate=ramp_allowance(c,g,Symbol(:ramp_,direction),a,a+response)
                flow=ramp_allowance(c,g,Symbol(:discharge_ramp_,direction),a,a+response)
                rate===nothing || record("$(g.name)_$(direction)_response",amount-rate)
                flow===nothing || record("$(g.name)_$(direction)_flow_response",
                    (d==1 ? q-x.generator_q[i,t] : x.generator_q[i,t]-q)-flow)
                budget=getfield(g,Symbol(:reserve_,direction,:_energy))
                budget===nothing || record("$(g.name)_$(direction)_energy",sum(
                    c.reserve_groups[k].activation_hours*R[i,k,t] for k in groups)-budget)
            end
        end
        for plant in c.system.plants
            units=findall(g->g.plant==plant.name,c.system.generators)
            groups=findall(k->!isempty(intersect(units,members[k])),1:K)
            isempty(groups) && continue
            totalp=sum(actual[units]);totalq=sum(Q[units,t])
            pmin=opinterval(c,plant.name,:pmin,t,plant.pmin)
            qmin=opinterval(c,plant.name,:qmin,t,plant.qmin)
            on=any(x.commitment[i,t]==1 for i in units) ? 1 : 0
            record("$(plant.name)_$(direction)_plant_power",max(pmin*on-totalp,totalp-opinterval(c,plant.name,:pmax,t,plant.pmax)))
            maxq=optional_operation(c,plant.name,:qmax,t;default=plant.qmax)
            record("$(plant.name)_$(direction)_plant_flow",max(qmin*on-totalq,maxq===nothing ? 0.0 : totalq-maxq))
            response=minimum(c.reserve_groups[k].response_hours for k in groups)
            a=(c.grid[t]+c.grid[t+1])/2
            rate=ramp_allowance(c,plant,Symbol(:ramp_,direction),a,a+response)
            flow=ramp_allowance(c,plant,Symbol(:discharge_ramp_,direction),a,a+response)
            amount=sum(R[i,k,t] for i in units,k in 1:K)
            rate===nothing || record("$(plant.name)_$(direction)_plant_response",amount-rate)
            flow===nothing || record("$(plant.name)_$(direction)_plant_flow_response",
                (d==1 ? totalq-sum(x.generator_q[units,t]) : sum(x.generator_q[units,t])-totalq)-flow)
        end
    end
    finish()
end

"""Capacity controls with nominal-discharge guesses for their two operating points.
`up` and `down` are unit × group × interval MW arrays. Forward reconstruction
infers the physical witnesses and shortfalls; independent schedule validation
checks the stored witnesses without repairing them.
"""
function ReserveAllocation(;up,down,generator_q)
    size(up)==size(down) && size(up,1)==size(generator_q,1) &&
    size(up,3)==size(generator_q,2) || throw(DimensionMismatch("reserve capacity controls"))
    K,T=size(up,2),size(up,3)
    ReserveAllocation(_reserve_array(up),_reserve_array(down),
        Matrix{Float64}(generator_q),Matrix{Float64}(generator_q),zeros(K,T),zeros(K,T))
end

function _reserve_bisection(residual,lo,hi,tolerance)
    a,b=lo,hi;fa,fb=only(residual([a])),only(residual([b]))
    isfinite(fa) && isfinite(fb) || return nothing
    abs(fa)<=tolerance && return a
    abs(fb)<=tolerance && return b
    signbit(fa)==signbit(fb) && return nothing
    for iteration in 1:80
        mid=(a+b)/2;fm=only(residual([mid]))
        isfinite(fm) || return nothing
        abs(fm)<=tolerance && return mid
        if signbit(fm)==signbit(fa)
            a=mid;fa=fm
        else
            b=mid
        end
    end
    nothing
end

# Invert the plant power equations with bounded damped Newton. Scalar failures
# fall back to a sign bracket, without assuming power monotonicity. Bounds only
# restrict trial guesses; capacity controls are never clipped.
function reserve_witnesses(c,r::ReserveAllocation,u,trajectory;tolerance=1e-4)
    G,K,T=length(c.system.generators),length(c.reserve_groups),length(c.prices)
    size(r.up)==size(r.down)==(G,K,T) && size(r.up_q)==size(r.down_q)==(G,T) ||
        throw(DimensionMismatch("reserve capacity controls and witness guesses"))
    all(isfinite,r.up) && all(isfinite,r.down) && all(>=(0),r.up) && all(>=(0),r.down) ||
        throw(ArgumentError("finite nonnegative reserve capacity controls required"))
    ix=nodeindex(c.system)
    upq=copy(trajectory.generator_q);downq=copy(upq)
    for plant in c.system.plants,t in 1:T,d in 1:2
        members=findall(g->g.plant==plant.name,c.system.generators)
        all(c.system.generators[i].mode==:pump for i in members) && continue
        active=filter(i->u[i,t]==1,members)
        isempty(active) && continue
        capacities=d==1 ? r.up : r.down
        target=[trajectory.power[i,t]+(d==1 ? 1 : -1)*sum(capacities[i,:,t]) for i in active]
        guess=d==1 ? r.up_q : r.down_q
        lo=[d==1 ? trajectory.generator_q[i,t] : opinterval(c,c.system.generators[i].name,:qmin,t,c.system.generators[i].qmin) for i in active]
        hi=[d==1 ? opinterval(c,c.system.generators[i].name,:qmax,t,c.system.generators[i].qmax) : trajectory.generator_q[i,t] for i in active]
        q=clamp.([guess[i,t] for i in active],lo,hi)
        gross=trajectory.head[ix[plant.source],t]-outlet_head(plant,trajectory.head[ix[plant.target],t])
        residual(q)=begin
            h=gross-tailwater(plant,sum(q))
            [power(c.system.generators[i],q[j],h)-target[j] for (j,i) in enumerate(active)]
        end
        solve_tolerance=min(tolerance,1e-8)
        converged=maximum(abs,residual(q);init=0.0)<=solve_tolerance
        # Coupled nominal guesses can sit at stationary power maxima. Supplied
        # operating points retain their branch; only nominal guesses move inward.
        if !converged && length(q)>1
            for (j,i) in enumerate(active)
                q[j]==trajectory.generator_q[i,t] && target[j]!=trajectory.power[i,t] &&
                    (q[j]=(lo[j]+hi[j])/2)
            end
        end
        for iteration in 1:60
            converged && break
            f=residual(q)
            normf=maximum(abs,f;init=0.0)
            if isfinite(normf) && normf<=solve_tolerance
                converged=true;break
            end
            jacobian=zeros(length(q),length(q))
            for j in eachindex(q)
                step=1e-5*max(1.0,abs(q[j]))
                a,b=copy(q),copy(q);a[j]+=step;b[j]-=step
                jacobian[:,j]=(residual(a)-residual(b))/(2step)
            end
            change=try jacobian\f catch;break end
            all(isfinite,change) || break
            alpha=1.0
            for j in eachindex(q)
                change[j]>0 && (alpha=min(alpha,0.99*(q[j]-lo[j])/change[j]))
                change[j]<0 && (alpha=min(alpha,0.99*(hi[j]-q[j])/-change[j]))
            end
            accepted=false
            for backtrack in 1:40
                trial=q-alpha*change
                if all(lo .- 1e-9 .<= trial .<= hi .+ 1e-9) &&
                    maximum(abs,residual(trial);init=0.0)<normf
                    q=trial;accepted=true;break
                end
                alpha/=2
            end
            accepted || break
        end
        if !converged && length(q)==1
            root=_reserve_bisection(residual,only(lo),only(hi),solve_tolerance)
            if root!==nothing
                g=c.system.generators[only(active)]
                h=gross-tailwater(plant,root);eta=efficiency(g,root,h)
                admitted=g.hmin<=h<=g.hmax && g.min_efficiency<=eta<=1.0 &&
                    (g.turbine_table===nothing ||
                        turbine_qmin(g.turbine_table,h;extrapolation=:linear)<=root<=
                        turbine_qmax(g.turbine_table,h;extrapolation=:linear))
                admitted && (q[1]=root;converged=true)
            end
        end
        # A boundary point can miss an exact root by less than solver/audit
        # tolerance. Keep that physical point rather than altering its capacity.
        converged || (converged=maximum(abs,residual(q);init=0.0)<=tolerance)
        converged || throw(ArgumentError("reserve operating-point reconstruction failed at $(plant.name), interval $t"))
        destination=d==1 ? upq : downq
        destination[active,t]=q
    end
    su=[max(0.0,opinterval(c,g.name,:up_requirement,t,g.up_requirement)-sum(r.up[:,k,t]))
        for (k,g) in enumerate(c.reserve_groups),t in 1:T]
    sd=[max(0.0,opinterval(c,g.name,:down_requirement,t,g.down_requirement)-sum(r.down[:,k,t]))
        for (k,g) in enumerate(c.reserve_groups),t in 1:T]
    ReserveAllocation(copy(r.up),copy(r.down),upq,downq,su,sd)
end

Base.:(==)(a::ReserveAllocation,b::ReserveAllocation) =
    all(getfield(a,field)==getfield(b,field) for field in fieldnames(ReserveAllocation))
function Base.show(io::IO,::MIME"text/plain",r::ReserveAllocation)
    print(io,"ReserveAllocation(",size(r.up,1)," units, ",size(r.up,2),
        " groups, ",size(r.up,3)," intervals)")
end

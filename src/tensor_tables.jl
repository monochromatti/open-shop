# Exact tensor interpolation with independent one-dimensional SOS2 coordinates.
# Cubic corrections share discharge products across all head columns.
function _global_tensor_coordinate!(m, x, nodes; name)
    nn=Float64.(nodes)
    !isempty(nn) && all(isfinite,nn) && all(diff(nn).>0) ||
        throw(ArgumentError("tensor coordinates need finite increasing nodes"))
    cache=get!(m.ext,:global_tensor_coordinate_cache,Dict{Any,Any}())
    key=x isa VariableRef ? (x,Tuple(nn)) : nothing
    key!==nothing && haskey(cache,key) && return cache[key]
    n=length(nn)
    weights=if n==1
        @constraint(m,x==only(nn))
        [1.0]
    else
        ww=@variable(m,[1:n],lower_bound=0,upper_bound=1,base_name="$(name)_weight")
        @constraint(m,sum(ww)==1)
        # Centering avoids cancellation between large nearby coordinate values.
        @constraint(m,x==nn[1]+sum((nn[i]-nn[1])*ww[i] for i in 2:n))
        n>2 && @constraint(m,ww in SOS2(collect(1.0:n)))
        ww
    end
    data=(name=string(name),x=x,nodes=nn,weights=weights)
    push!(get!(m.ext,:global_tensor_coordinates,Any[]),data)
    key!==nothing && (cache[key]=data)
    data
end

function _global_tensor_nodes(grid,lo,hi)
    all(isfinite,(lo,hi)) && lo<=hi ||
        throw(ArgumentError("tensor coordinates need finite ordered bounds"))
    unique(vcat(Float64(lo),filter(v->lo<v<hi,grid),Float64(hi)))
end

function _global_tensor_coordinate_weights(nodes,x)
    isfinite(x) && first(nodes)<=x<=last(nodes) ||
        throw(DomainError(x,"tensor coordinate outside its graph domain"))
    n=length(nodes)
    n==1 && return [1.0]
    i=clamp(searchsortedlast(nodes,x),1,n-1)
    t=(x-nodes[i])/(nodes[i+1]-nodes[i])
    weights=zeros(n)
    weights[i]=1-t
    weights[i+1]=t
    weights
end

function _global_tensor_table!(m,c::TableCurve,x,lo,hi;name=gensym(:tensor_table))
    coord=_global_tensor_coordinate!(m,x,_global_tensor_nodes(c.x,lo,hi);name)
    values=[table_value(c,z;extrapolation=:linear) for z in coord.nodes]
    all(isfinite,values) || throw(ArgumentError("tensor table extension bounds overflow"))
    y=@variable(m,lower_bound=minimum(values),upper_bound=maximum(values),base_name=string(name))
    @constraint(m,y==sum(values[i]*coord.weights[i] for i in eachindex(values)))
    get!(m.ext,:global_tensor_tables,Dict{String,Any}())[string(name)]=(coordinate=coord,values=values,y=y)
    y
end

# Coefficients in the local [0,1] discharge coordinate at a supplied head.
# Each extension segment keeps its original secant, even for PCHIP tables.
function _global_tensor_polynomial(c,ql,qr,h)
    qm=ql+(qr-ql)/2
    i=_curve_segment(c.discharge,qm,:linear)
    j=_curve_segment(c.heads,h,:linear)
    cubic=c.interpolation==:pchip_discharge && first(c.discharge)<=qm<=last(c.discharge)
    a=_global_discharge_coefficients(c,i,j,cubic)
    b=_global_discharge_coefficients(c,i,j+1,cubic)
    u=(h-c.heads[j])/(c.heads[j+1]-c.heads[j])
    coefficients=ntuple(k->a[k]+u*(b[k]-a[k]),4)
    dq=c.discharge[i+1]-c.discharge[i]
    _global_shift_polynomial(coefficients,(ql-c.discharge[i])/dq,(qr-ql)/dq)
end

function _global_tensor_turbine!(m,c::TurbineTable,q,h,qlo,qhi,hlo,hhi;
    name=gensym(:tensor_turbine),quadratic=false)
    qc=_global_tensor_coordinate!(m,q,_global_tensor_nodes(c.discharge,qlo,qhi);name="$(name)_q")
    hc=_global_tensor_coordinate!(m,h,_global_tensor_nodes(c.heads,hlo,hhi);name="$(name)_h")
    nq,nh=length(qc.nodes),length(hc.nodes)
    nodal=[turbine_efficiency(c,x,y;extrapolation=:linear) for x in qc.nodes,y in hc.nodes]
    A=zeros(max(0,nq-1),nh)
    B=similar(A)
    ranges=Tuple{Float64,Float64}[]
    for j in 1:nh
        rr=Tuple{Float64,Float64}[]
        for i in 1:(nq-1)
            a=_global_tensor_polynomial(c,qc.nodes[i],qc.nodes[i+1],hc.nodes[j])
            # p(t) = p(0)(1-t)+p(1)t + A(1-t)^2*t + B(1-t)*t^2.
            A[i,j]=-a[3]-a[4]
            B[i,j]=-a[3]-2a[4]
            push!(rr,_global_polynomial_range(a,0.0,1.0))
        end
        lo,hi=nq==1 ? (nodal[1,j],nodal[1,j]) : (minimum(first,rr),maximum(last,rr))
        slack=1e-10*max(1.0,abs(lo),abs(hi))
        push!(ranges,(lo-slack,hi+slack))
    end
    all(isfinite,nodal) && all(isfinite,A) && all(isfinite,B) &&
        all(pair->all(isfinite,pair),ranges) ||
        throw(ArgumentError("tensor turbine extension bounds overflow"))
    active=[i for i in 1:(nq-1) if any(!iszero,view(A,i,:)) || any(!iszero,view(B,i,:))]
    w=Dict{Int,VariableRef}()
    r=Dict{Int,VariableRef}()
    s=Dict{Int,Any}()
    for i in active
        r[i]=@variable(m,lower_bound=0,upper_bound=4/27+1e-12,base_name="$(name)_r[$i]")
        if quadratic
            w[i]=@variable(m,lower_bound=0,upper_bound=1/4+1e-12,base_name="$(name)_w[$i]")
            @constraint(m,w[i]==qc.weights[i]*qc.weights[i+1])
            @constraint(m,r[i]==qc.weights[i]*w[i])
            # An active SOS2 pair sums to one; inactive pairs have w=r=0.
            # Thus the second cubic basis is exactly w-r, with no new variable.
            s[i]=w[i]-r[i]
            @constraint(m,0<=s[i]<=4/27+1e-12)
        else
            s[i]=@variable(m,lower_bound=0,upper_bound=4/27+1e-12,base_name="$(name)_s[$i]")
            @constraint(m,r[i]==qc.weights[i]^2*qc.weights[i+1])
            @constraint(m,s[i]==qc.weights[i]*qc.weights[i+1]^2)
        end
    end
    # At most one adjacent discharge pair contributes at an SOS2 feasible point.
    !isempty(w) && @constraint(m,sum(values(w))<=1/4+1e-12)
    columns=@variable(m,[1:nh],base_name="$(name)_column")
    for j in 1:nh
        set_lower_bound(columns[j],ranges[j][1])
        set_upper_bound(columns[j],ranges[j][2])
        @constraint(m,columns[j]==sum(nodal[i,j]*qc.weights[i] for i in 1:nq)+
            sum(A[i,j]*r[i]+B[i,j]*s[i] for i in active;init=0.0))
    end
    eta=@variable(m,lower_bound=minimum(first,ranges),upper_bound=maximum(last,ranges),base_name=string(name))
    @constraint(m,eta==sum(columns[j]*hc.weights[j] for j in 1:nh))
    get!(m.ext,:global_tensor_turbines,Dict{String,Any}())[string(name)]=
        (qcoordinate=qc,hcoordinate=hc,nodal=nodal,A=A,B=B,w=w,r=r,s=s,
         quadratic=quadratic,columns=columns,eta=eta)
    eta
end

# Values for every turbine-specific auxiliary. Shared coordinate weights are
# included so callers can merge these dictionaries with the physical start.
function _global_tensor_turbine_values(data,q,h)
    qweights=_global_tensor_coordinate_weights(data.qcoordinate.nodes,q)
    hweights=_global_tensor_coordinate_weights(data.hcoordinate.nodes,h)
    values=Dict{VariableRef,Float64}()
    for (coord,weights) in ((data.qcoordinate,qweights),(data.hcoordinate,hweights))
        for (v,x) in zip(coord.weights,weights)
            v isa VariableRef && (values[v]=x)
        end
    end
    for i in keys(data.r)
        values[data.r[i]]=qweights[i]^2*qweights[i+1]
        haskey(data.w,i) && (values[data.w[i]]=qweights[i]*qweights[i+1])
        data.s[i] isa VariableRef && (values[data.s[i]]=qweights[i]*qweights[i+1]^2)
    end
    eta=0.0
    for j in eachindex(data.columns)
        e=sum(data.nodal[i,j]*qweights[i] for i in eachindex(qweights))+
            sum(data.A[i,j]*values[data.r[i]]+data.B[i,j]*qweights[i]*qweights[i+1]^2
                for i in keys(data.r);init=0.0)
        values[data.columns[j]]=e
        eta+=e*hweights[j]
    end
    values[data.eta]=eta
    values
end

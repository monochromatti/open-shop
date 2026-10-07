module TableConsistency
using JuMP
export add_table_consistency!, lift_table_consistency!

_value(x,assigned)=x isa Number ? Float64(x) : JuMP.value(v->assigned[v],x)

"""Add transportation envelopes while retaining every exact table equation.

`mode=:bilinear` selects tables without cubic corrections; `:all` also links
the PCHIP basis products. `selected` can restrict the table names. Size limits
select tables in deterministic name order, without restricting physical domains.
"""
function add_table_consistency!(b,c;mode=:all,selected=nothing,
    max_tables=typemax(Int),max_added_variables=typemax(Int))
    mode in (:bilinear,:all) || throw(ArgumentError("table consistency mode must be :bilinear or :all"))
    max_tables>=0 && max_added_variables>=0 || throw(ArgumentError("table consistency size limits must be nonnegative"))
    m=b.m
    haskey(m.ext,:table_consistency) && throw(ArgumentError("table consistency has already been added"))
    records=Dict{String,Any}()
    m.ext[:table_consistency]=records
    tables=get(m.ext,:global_tensor_turbines,Dict())
    wanted=selected===nothing ? nothing : Set(string.(selected))
    before_variables=num_variables(m)
    before_constraints=num_constraints(m;count_variable_in_set_constraints=false)
    skipped_affine=0;skipped_mode=0;skipped_size=0
    for n in sort!(collect(keys(tables)))
        wanted!==nothing && !(n in wanted) && continue
        data=tables[n]
        nq,nh=size(data.nodal)
        if nq==1 || nh==1
            skipped_affine+=1
            continue
        end
        active=sort!(collect(keys(data.r)))
        if mode==:bilinear && !isempty(active)
            skipped_mode+=1
            continue
        end
        added=(nq+2length(active))*nh
        if length(records)>=max_tables || num_variables(m)-before_variables+added>max_added_variables
            skipped_size+=1
            continue
        end
        lambda=data.qcoordinate.weights
        mu=data.hcoordinate.weights
        gamma=@variable(m,[1:nq,1:nh],lower_bound=0,upper_bound=1,base_name="$(n)_joint")
        for i in 1:nq
            @constraint(m,sum(gamma[i,j] for j in 1:nh)==lambda[i])
        end
        for j in 1:nh
            @constraint(m,sum(gamma[i,j] for i in 1:nq)==mu[j])
        end
        R=Dict{Tuple{Int,Int},VariableRef}()
        S=Dict{Tuple{Int,Int},VariableRef}()
        for i in active
            # The source SOS2 graph bounds both cubic bases by4/27. Preserve
            # its outward guards; nodal/correction coefficients may be negative.
            rmax=upper_bound(data.r[i])
            smax=4/27+1e-12
            for j in 1:nh
                R[(i,j)]=@variable(m,lower_bound=0,upper_bound=rmax,base_name="$(n)_joint_r[$i,$j]")
                S[(i,j)]=@variable(m,lower_bound=0,upper_bound=smax,base_name="$(n)_joint_s[$i,$j]")
                @constraint(m,R[(i,j)]<=rmax*mu[j])
                @constraint(m,R[(i,j)]>=data.r[i]-rmax*(1-mu[j]))
                @constraint(m,S[(i,j)]<=smax*mu[j])
                @constraint(m,S[(i,j)]>=data.s[i]-smax*(1-mu[j]))
                # w=lambda_i*lambda_(i+1) is bounded by either endpoint weight.
                @constraint(m,R[(i,j)]+S[(i,j)]<=gamma[i,j])
                @constraint(m,R[(i,j)]+S[(i,j)]<=gamma[i+1,j])
            end
            @constraint(m,sum(R[(i,j)] for j in 1:nh)==data.r[i])
            @constraint(m,sum(S[(i,j)] for j in 1:nh)==data.s[i])
        end
        if !isempty(active)
            for j in 1:nh
                @constraint(m,sum(R[(i,j)]+S[(i,j)] for i in active)<=(1/4+1e-12)*mu[j])
            end
        end
        @constraint(m,data.eta==sum(data.nodal[i,j]*gamma[i,j] for i in 1:nq,j in 1:nh)+
            sum(data.A[i,j]*R[(i,j)]+data.B[i,j]*S[(i,j)] for i in active,j in 1:nh;init=0.0))
        records[n]=(data=data,gamma=gamma,R=R,S=S)
    end
    profile=Dict(
        "mode"=>string(mode),"tables_added"=>length(records),
        "bilinear_tables"=>count(r->isempty(r.data.r),values(records)),
        "pchip_tables"=>count(r->!isempty(r.data.r),values(records)),
        "variables_added"=>num_variables(m)-before_variables,
        "constraints_added"=>num_constraints(m;count_variable_in_set_constraints=false)-before_constraints,
        "skipped_affine"=>skipped_affine,"skipped_mode"=>skipped_mode,"skipped_size"=>skipped_size,
        "table_names"=>sort!(collect(keys(records))))
    m.ext[:table_consistency_profile]=profile
    !isempty(records) && push!(get!(m.ext,:experiment_start_lifters,Any[]),
        assigned->lift_table_consistency!(m,assigned))
    profile
end

"""Extend an assigned physical/table start with actual products, without clipping."""
function lift_table_consistency!(m,assigned)
    for record in values(get(m.ext,:table_consistency,Dict()))
        data=record.data
        lambda=[_value(v,assigned) for v in data.qcoordinate.weights]
        mu=[_value(v,assigned) for v in data.hcoordinate.weights]
        for i in eachindex(lambda),j in eachindex(mu)
            assigned[record.gamma[i,j]]=lambda[i]*mu[j]
        end
        for (i,j) in keys(record.R)
            assigned[record.R[(i,j)]]=_value(data.r[i],assigned)*mu[j]
            assigned[record.S[(i,j)]]=_value(data.s[i],assigned)*mu[j]
        end
    end
    assigned
end
end

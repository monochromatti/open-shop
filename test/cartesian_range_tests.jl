@testset "Cartesian turbine range strengthening" begin
    function cartesian_values(m,table,qv,hv)
        byname=Dict(name(v)=>v for v in all_variables(m))
        values=Dict(byname["q"]=>qv,byname["h"]=>hv,
            byname["eta"]=>OpenSHOP.turbine_efficiency(table,qv,hv;extrapolation=:linear))
        cells=m.ext[:global_turbine_cells]["eta"]
        selected=findfirst(pair->pair[1][1]<=qv<=pair[1][2] &&
            pair[2][1]<=hv<=pair[2][2],cells)
        for (k,((ql,qr),(hl,hr))) in enumerate(cells)
            i=OpenSHOP._curve_segment(table.discharge,ql+(qr-ql)/2,:linear)
            j=OpenSHOP._curve_segment(table.heads,hl+(hr-hl)/2,:linear)
            haskey(byname,"eta_cell[$k]") && (values[byname["eta_cell[$k]"]]=k==selected ? 1.0 : 0.0)
            values[byname["eta_q[$k]"]]=k==selected ?
                (qv-table.discharge[i])/(table.discharge[i+1]-table.discharge[i]) : 0.0
            values[byname["eta_h[$k]"]]=k==selected ?
                (hv-table.heads[j])/(table.heads[j+1]-table.heads[j]) : 0.0
        end
        values
    end
    function cartesian_residual(m,values)
        residual=0.0
        for ref in all_constraints(m;include_variable_in_set_constraints=true)
            obj=constraint_object(ref)
            x=obj.func isa Number ? obj.func : JuMP.value(v->values[v],obj.func)
            set=obj.set
            violation=set isa OpenSHOP.MOI.EqualTo ? abs(x-set.value) :
                set isa OpenSHOP.MOI.LessThan ? max(0.0,x-set.upper) :
                set isa OpenSHOP.MOI.GreaterThan ? max(0.0,set.lower-x) :
                set isa OpenSHOP.MOI.Interval ? max(0.0,set.lower-x,x-set.upper) :
                set isa OpenSHOP.MOI.ZeroOne ? max(abs(x-round(x)),max(0.0,-x,x-1)) :
                error("unsupported Cartesian test constraint")
            residual=max(residual,violation)
        end
        residual
    end
    for interpolation in (:bilinear,:pchip_discharge)
        table=TurbineTable([50.0,90.0,140.0],[2.0,5.0,10.0,16.0],
            [0.70 0.75 0.79;0.91 0.95 0.93;0.85 0.92 0.96;0.73 0.81 0.86],
            [2.0,2.0,2.0],[16.0,16.0,16.0];interpolation)
        for (qlo,qhi,hlo,hhi) in ((-1.0,19.0,30.0,170.0),(3.0,14.0,60.0,130.0),
                                  (5.0,5.0+1e-8,89.0,90.0),(7.0,7.0,100.0,100.0))
            models=[]
            etas=[]
            for (exact_bounds,range_cuts) in ((false,false),(true,false),(false,true))
                m=Model();@variable(m,qlo<=q<=qhi);@variable(m,hlo<=h<=hhi)
                eta=OpenSHOP._global_turbine!(m,table,q,h,qlo,qhi,hlo,hhi;
                    name=:eta,tightened=false,exact_bounds,range_cuts)
                push!(models,m);push!(etas,eta)
            end
            @test all(m->m.ext[:global_turbine_cells]==models[1].ext[:global_turbine_cells],models)
            @test all(m->!m.ext[:global_turbine_normalized]["eta"],models)
            @test all(m->num_variables(m)==num_variables(models[1]),models)
            counts=[length(all_constraints(m;include_variable_in_set_constraints=true)) for m in models]
            @test counts==[counts[1],counts[1],counts[1]+2]
            @test lower_bound(etas[2])>=lower_bound(etas[1])-1e-9
            @test upper_bound(etas[2])<=upper_bound(etas[1])+1e-9
            @test lower_bound(etas[3])==lower_bound(etas[2])
            @test upper_bound(etas[3])==upper_bound(etas[2])
            cells=models[1].ext[:global_turbine_cells]["eta"]
            qs=unique(vcat(qlo,qhi,vec([t==0 ? l : t==1 ? r : l+t*(r-l)
                for (l,r) in unique(first.(cells)),t in (0.,1/3,2/3,1.)])))
            hs=unique(vcat(hlo,hhi,vec([t==0 ? l : t==1 ? r : l+t*(r-l)
                for (l,r) in unique(last.(cells)),t in (0.,0.5,1.)])))
            for m in models
                max_residual=0.0
                complete=true
                for qv in qs,hv in hs
                    values=cartesian_values(m,table,qv,hv)
                    complete &= length(values)==num_variables(m)
                    max_residual=max(max_residual,cartesian_residual(m,values))
                end
                @test complete
                @test max_residual<1e-9
            end
        end
    end
    # The selector cuts contain information beyond the union efficiency bounds.
    table=TurbineTable([60.0,100.0],[2.0,4.0,8.0],
        [0.70 0.75;0.85 0.90;0.90 0.95],[2.0,2.0],[8.0,8.0])
    m=Model();@variable(m,0<=q<=8);@variable(m,60<=h<=100)
    eta=OpenSHOP._global_turbine!(m,table,q,h,0.0,8.0,60.0,100.0;
        name=:eta,tightened=false,range_cuts=true)
    cuts=[ref for ref in all_constraints(m;include_variable_in_set_constraints=false)
        if constraint_object(ref).func isa JuMP.GenericAffExpr &&
           coefficient(constraint_object(ref).func,eta)!=0 &&
           !(constraint_object(ref).set isa OpenSHOP.MOI.EqualTo)]
    @test length(cuts)==2
    values=cartesian_values(m,table,0.0,60.0)
    values[eta]=upper_bound(eta)
    @test any(cuts) do ref
        obj=constraint_object(ref)
        x=JuMP.value(v->values[v],obj.func)
        obj.set isa OpenSHOP.MOI.LessThan ? x>obj.set.upper+1e-3 : x<obj.set.lower-1e-3
    end
end

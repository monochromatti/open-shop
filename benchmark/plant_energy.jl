# At one shared plant head, electrical output cannot exceed hydraulic power
# times the greatest admissible turbine and electrical efficiencies. This is
# redundant for feasible schedules but avoids intermediate product variables
# in the aggregate upper bound. It uses the declared midpoint equations.
function add_plant_energy_bounds!(b,c)
    for plant in c.system.plants,t in eachindex(c.prices)
        units=findall(g->g.plant==plant.name,c.system.generators)
        isempty(units) && continue
        coefficients=Float64[]
        supported=true
        for i in units
            g=c.system.generators[i]
            eta=variable_by_name(b.m,g.turbine_table===nothing ? "eta_$(i)_$(t)" : "turbine_$(i)_$(t)")
            electrical=g.generator_efficiency_curve===nothing ? nothing :
                variable_by_name(b.m,"electrical_$(i)_$(t)")
            if electrical!==nothing && lower_bound(electrical)<0
                supported=false;break
            end
            emax=electrical===nothing ? 1. : upper_bound(electrical)
            push!(coefficients,max(0.,min(1.,upper_bound(eta)))*emax)
        end
        supported || continue
        head=b.shared_heads[(plant.name,t)]
        power=sum(b.P[i,t] for i in units)
        effective_flow=sum(a*b.GQ[i,t] for (a,i) in zip(coefficients,units))
        @constraint(b.m,(power-0.00981*head*effective_flow)/40<=0)
    end
end

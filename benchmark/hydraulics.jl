# Run with julia --project=. benchmark/hydraulics.jl [output.json].
using OpenSHOP, JuMP, Statistics, LinearAlgebra
include(joinpath(pkgdir(OpenSHOP),"examples","mixed_source_network.jl"))

function hydraulic_benchmark(output)
    rows=Any[]
    for distributed in (false,true)
        fixture=MixedSourceNetworkExample.mixed_source_fixture(;distributed)
        c=fixture.case
        x=dispatch_from_controls(c;commitment=fixture.u,
            generator_q=fixture.q,gate=fixture.gates)
        x["validation"]["valid"] || error("invalid benchmark controls")
        simulate(c;generator_q=fixture.q,gate=fixture.gates)
        samples=[@timed simulate(c;generator_q=fixture.q,gate=fixture.gates) for _ in 1:9]
        local_model=OpenSHOP._build_dispatch(c;u=fixture.u)
        global_model=OpenSHOP._build_global_dispatch(c;joint=true)
        counts(model)=Dict("variables"=>num_variables(model),
            "constraints"=>num_constraints(model;count_variable_in_set_constraints=false))
        push!(rows,Dict("case"=>c.name,"seconds"=>median(s.time for s in samples),
            "allocated_bytes"=>median(s.bytes for s in samples),
            "objective"=>x["objective"],"local_model"=>counts(local_model.m),
            "global_model"=>counts(global_model.m)))
    end
    writejson(output,Dict("julia"=>string(VERSION),"cpu"=>Sys.CPU_NAME,
        "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "warm_repetitions"=>9,"cases"=>rows))
end

if abspath(PROGRAM_FILE)==@__FILE__
    hydraulic_benchmark(isempty(ARGS) ? "hydraulics.json" : only(ARGS))
end

# Paired native SCIP diagnostics: one frozen case and one seed for every variant.
include("diagnose.jl")
input=isempty(ARGS) ? joinpath(@__DIR__,"..","examples","tokke_vinje","generated","operating") : abspath(ARGS[1])
output=length(ARGS)<2 ? joinpath(@__DIR__,"..","results","tokke_vinje_paired") : abspath(ARGS[2])
margin=length(ARGS)<3 ? .1 : parse(Float64,ARGS[3])
seconds=length(ARGS)<4 ? 60. : parse(Float64,ARGS[4])
repeats=length(ARGS)<5 ? 1 : parse(Int,ARGS[5])
# Accept one arbitrary JSON case or a directory of generated horizons, including
# seasonal variants. Mapping reports and seed/result files are not case inputs.
paths=if isfile(input)
    [input]
else
    filter(readdir(input;join=true)) do path
        endswith(path,".json") || return false
        try
            data=JSON3.read(read(path,String));haskey(data,"grid") && haskey(data,"prices") && haskey(data,"generators")
        catch
            false
        end
    end
end
isempty(paths) && error("no case JSON files found at $input")
paired_benchmark(sort(paths);output,time_limit=seconds,repeats,operational_margin=margin)

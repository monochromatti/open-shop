"""Run examples in separate namespaces while reusing compiled package code."""
function run_examples(files)
    for file in files
        path=abspath(file)
        println("Running ",file)
        flush(stdout)
        scope=Module(gensym(:Example))
        # Preserve script entry guards without sharing example globals.
        Core.eval(scope, :(const PROGRAM_FILE = $path))
        seconds=@elapsed Base.include(scope,path)
        println("Completed ",file," in ",round(seconds;digits=2)," s")
        flush(stdout)
    end
    nothing
end

if abspath(PROGRAM_FILE)==@__FILE__
    isempty(ARGS) && error("usage: run_examples.jl example.jl [example.jl ...]")
    run_examples(ARGS)
end

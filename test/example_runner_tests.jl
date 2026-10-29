include(joinpath(@__DIR__,"..","scripts","run_examples.jl"))

@testset "Example batches execute script entries in isolated namespaces" begin
    mktempdir() do dir
        files=[joinpath(dir,"example$i.jl") for i in 1:2]
        for (i,file) in enumerate(files)
            write(file,"""
                const example_value = $i
                if abspath(PROGRAM_FILE)==@__FILE__
                    write(joinpath(@__DIR__,"output$i"),string(example_value))
                end
                """)
        end
        run_examples(files)
        @test read(joinpath(dir,"output1"),String)=="1"
        @test read(joinpath(dir,"output2"),String)=="2"

        broken=joinpath(dir,"broken.jl")
        write(broken,"error(\"example failed\")")
        rm(joinpath(dir,"output2"))
        @test_throws LoadError run_examples([broken,files[2]])
        @test !isfile(joinpath(dir,"output2"))
    end
end

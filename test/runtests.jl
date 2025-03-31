using MapBField2Surface
using Test
using Aqua

@testset "MapBField2Surface.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(MapBField2Surface)
    end
    # Write your tests here.
end

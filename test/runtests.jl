using OffsetTables
using Test
using Aqua

@testset "OffsetTables.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(OffsetTables)
    end
    # Write your tests here.
end

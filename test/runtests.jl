using OffsetTables
using Test
using Aqua

@testset "OffsetTables.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(OffsetTables, deps_compat=false)
        Aqua.test_deps_compat(OffsetTables, check_extras=false)
    end

    @testset "Name" begin
        @test rand(OffsetTable([1])) == 1
        @test rand(OffsetTable([0])) == 1
        @test_throws ArgumentError OffsetTable(Int[])
        @test rand(OffsetTable([1, 1])) in [1, 2]
        @test_broken rand(OffsetTable([1, 0])) == 1 # TODO: support 100% probability
        @test_broken rand(OffsetTable([1, 1, 1])) in [1, 2, 3] # TODO: support non-power of 2 lengths
    end
end

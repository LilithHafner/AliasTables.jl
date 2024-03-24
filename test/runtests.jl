using OffsetTables
using Test
using Aqua

@testset "OffsetTables.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(OffsetTables, deps_compat=false)
        Aqua.test_deps_compat(OffsetTables, check_extras=false)
    end

    @testset "Basic" begin
        @test rand(OffsetTable([1])) == 1
        @test_throws ArgumentError OffsetTable([0])
        @test_throws ArgumentError OffsetTable(Int[])
        @test rand(OffsetTable([1, 1])) in [1, 2]
        @test rand(OffsetTable([1, 0])) == 1
        @test rand(OffsetTable([1, 1, 1])) in [1, 2, 3]
        @test_throws MethodError OffsetTable(OffsetTable([1]).probability_offset)
        @test rand(OffsetTable([1e-70, 1])) == 2
    end

    @testset "Invalid weight error messages" begin
        @test_throws ArgumentError("found negative weight -1 at index 2") OffsetTable([1, -1])
        @test_throws ArgumentError("all weights are zero") OffsetTable([0, 0])
        @test_throws ArgumentError("all weights are zero") OffsetTable([0])
        @test_throws ArgumentError("all weights are zero") OffsetTable(UInt[0, 0])
        @test_throws ArgumentError("all weights are zero") OffsetTable(UInt[0])
        @test_throws ArgumentError("weights must be non-empty") OffsetTable(Int[])
        @test_throws ArgumentError("weights must be non-empty") OffsetTable(UInt[])
    end

    @testset "Equality and hashing" begin
        a = OffsetTable([1, 2, 3])
        b = OffsetTable([1, 2, 3, 0, 0])
        @test a == b
    end
end

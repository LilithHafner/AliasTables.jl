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
        @test rand(OffsetTable([0, 1]), 3) == [2,2,2]
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

        data = [
            [
                OffsetTable([1, 2, 5]),
                OffsetTable([2, 4, 10]),
                OffsetTable([1, 2, 5] * 1729),
                OffsetTable([1, 2, 5, 0, 0]),
                OffsetTable{UInt16, Int32}([1,2,5]),
                OffsetTable{UInt16}([1,2,5]),
                OffsetTable{UInt, Int32}([1,2,5]),
                OffsetTable(UInt[unsigned(1)<<61, unsigned(2)<<61, unsigned(5)<<61]),
            ],[
                OffsetTable([1, 2.0001, 5]),
            ],[
                OffsetTable([0,0,0,0,1]),
                OffsetTable([0,0,0,0,1,0,0,0,0,0,0,0]),
                OffsetTable([1e-70,0,0,0,1]),
                OffsetTable([1e-70,0,0,0,1,0,0,0,1e-70]),
            ],[
                OffsetTable([1, 2, 3, 5]),
                OffsetTable{UInt64, Int8}([1, 2, 3, 5]),
            ],[
                OffsetTable{UInt16}([1, 2, 3, 5]),
                OffsetTable{UInt16, Int8}([1, 2, 3, 5]),
            ]
        ]

        for group1 in data, group2 in data, a in group1, b in group2
            if group1 === group2
                for a in [a, deepcopy(a)], b in [b, deepcopy(b)]
                    @test a == b
                    hash(a) == hash(b) || @show a, b
                    @test hash(a) == hash(b)
                    @test OffsetTables.probabilities(float, a) == OffsetTables.probabilities(float, b)
                    if eltype(OffsetTables.probabilities(a)) == eltype(OffsetTables.probabilities(b))
                        @test OffsetTables.probabilities(a) == OffsetTables.probabilities(b)
                    end
                end
            else
                @test a != b
                @test hash(a) != hash(b)
                @test OffsetTables.probabilities(float, a) != OffsetTables.probabilities(float, b)
                @test OffsetTables.probabilities(a) != OffsetTables.probabilities(b)
            end
        end
    end
end

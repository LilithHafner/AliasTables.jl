using AliasTables
using Test
using Aqua
using Random
using RegressionTests

@testset "AliasTables.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(AliasTables, deps_compat=false)
        Aqua.test_deps_compat(AliasTables, check_extras=false)
    end

    @testset "Basic" begin
        @test rand(AliasTable([1])) == 1
        @test_throws ArgumentError AliasTable([0])
        @test_throws ArgumentError AliasTable(Int[])
        @test rand(AliasTable([1, 1])) in [1, 2]
        @test rand(AliasTable([1, 0])) == 1
        @test rand(AliasTable([1, 1, 1])) in [1, 2, 3]
        @test_throws MethodError AliasTable(AliasTable([1]).probability_alias)
        @test rand(AliasTable([1e-70, 1])) == 2
        @test rand(AliasTable([0, 1]), 3)::Vector{Int} == [2,2,2]
        @test rand(AliasTable{UInt, Int8}([0, 1]), 3)::Vector{Int8} == [2,2,2]
        @test rand(AliasTable([typemax(Int)-10, 5, 5, 5])) == 1
    end

    @testset "Invalid weight error messages" begin
        @test_throws ArgumentError("found negative weight -1") AliasTable([1, -1])
        @test_throws ArgumentError("found negative weight -1") AliasTable([1, 1, -1])
        @test_throws ArgumentError("all weights are zero") AliasTable([0, 0])
        @test_throws ArgumentError("all weights are zero") AliasTable([0])
        @test_throws ArgumentError("all weights are zero") AliasTable(UInt[0, 0])
        @test_throws ArgumentError("all weights are zero") AliasTable(UInt[0])
        @test_throws ArgumentError("weights must be non-empty") AliasTable(Int[])
        @test_throws ArgumentError("weights must be non-empty") AliasTable(UInt[])
        @test_throws ArgumentError("sum(weights) is too low") AliasTable(UInt[123, 456], normalize=false)
        @test_throws ArgumentError("sum(weights) is too high") AliasTable(UInt[unsigned(3)<<62, unsigned(2)<<62, unsigned(3)<<62], normalize=false)
        @test_throws ArgumentError("sum(weights) overflows") AliasTable(UInt[unsigned(3)<<62, unsigned(2)<<62, unsigned(3)<<62])
        @test AliasTables.probabilities(float, AliasTable(UInt[unsigned(3)<<61, unsigned(2)<<61, unsigned(3)<<61])) == [3,2,3] ./ 8
        @test AliasTables.probabilities(float, AliasTable(UInt[unsigned(3)<<61, unsigned(2)<<61, unsigned(3)<<61], normalize=false)) == [3,2,3] ./ 8
    end

    @testset "probabilities()" begin
        @test AliasTables.probabilities(float, AliasTable([1, 2, 3])) == [1, 2, 3]/6
        @test AliasTables.probabilities(float, AliasTable([1, 2, 3, 0, 0])) == [1, 2, 3, 0, 0]/6
        @test AliasTables.probabilities(AliasTable([1, 2, 3, 0, 2])) == [1, 2, 3, 0, 2] .<< 61
    end

    @testset "sample()" begin
        @test Base.hasmethod(AliasTables.sample, Tuple{UInt, AliasTable{UInt, Int}})
        @test Base.hasmethod(AliasTables.sample, Tuple{Random.MersenneTwister, AliasTable{UInt, Int}})
        @test !Base.hasmethod(AliasTables.sample, Tuple{UInt32, AliasTable{UInt64, Int}})
    end

    @testset "Exact" begin
        for i in 1:100
            p = rand(i)
            at = AliasTable(p)
            @test maximum(abs, AliasTables.probabilities(at) ./ (big(typemax(UInt))+1) .- p ./ sum(big, p)) ≤ .5^64
            @test AliasTables.probabilities(float, at) ≈  p ./ sum(p)
            @test AliasTable(AliasTables.probabilities(at)) == at
            # @test AliasTable(AliasTables.probabilities(float, at)) == at

            if i == 1
                p2 = [typemax(UInt)]
            else
                p2 = floor.(UInt, (typemax(UInt)/sum(p)) .* p)
                p2[end] = typemax(UInt) - sum(p2[1:end-1]) + 1
            end
            at2 = AliasTable(p2)
            @test maximum(abs, AliasTables.probabilities(at2) ./ (big(typemax(UInt))+1) .- p2 ./ sum(big, p2)) ≤ .5^64
            @test AliasTables.probabilities(at2) == p2
            @test AliasTables.probabilities(float, at2) ≈  p ./ sum(p)
            # @test AliasTable(AliasTables.probabilities(float, at2)) == at2
            @test AliasTable(AliasTables.probabilities(at2)) == at2
        end

        function counts(x, levels)
            c = zeros(Int, levels)
            for v in x
                c[v] += 1
            end
            c
        end
        let at = AliasTable{UInt16}([10, 5, 1])
            @test counts(Iterators.map(x -> AliasTables.sample(x, at), typemin(UInt16):typemax(UInt16)), 3) ==
                2^16/16 * [10, 5, 1]
        end
    end

    @testset "Equality and hashing" begin
        a = AliasTable([1, 2, 3])
        b = AliasTable([1, 2, 3, 0, 0])
        @test a != b
        @test a.probability_alias == b.probability_alias

        data = [
            [
                AliasTable([1, 2, 5]),
                AliasTable([2, 4, 10]),
                AliasTable([1, 2, 5] * 1729),
                AliasTable{UInt16, Int32}([1,2,5]),
                AliasTable{UInt16}([1,2,5]),
                AliasTable{UInt, Int32}([1,2,5]),
                AliasTable(UInt[unsigned(1)<<61, unsigned(2)<<61, unsigned(5)<<61]),
            ],[
                AliasTable([1, 2, 5, 0, 0]),
            ],[
                AliasTable([1, 2.0001, 5]),
            ],[
                AliasTable([0,0,0,0,1]),
                AliasTable([1e-70,0,0,0,1]),
            ],[
                AliasTable([0,0,0,0,1,0,0,0,0,0,0,0]),
                AliasTable([1e-70,0,0,0,1,0,0,0,0,0,0,1e-70]),
            ],[
                AliasTable([1, 2, 3, 5]),
                AliasTable{UInt64, Int8}([1, 2, 3, 5]),
            ],[
                AliasTable{UInt16}([1, 2, 3, 5]),
                AliasTable{UInt16, Int8}([1, 2, 3, 5]),
            ]
        ]

        for group1 in data, group2 in data, a in group1, b in group2
            if group1 === group2
                for a in [a, deepcopy(a)], b in [b, deepcopy(b)]
                    @test a == b
                    hash(a) == hash(b) || @show a, b
                    @test hash(a) == hash(b)
                    @test AliasTables.probabilities(float, a) == AliasTables.probabilities(float, b)
                    if eltype(AliasTables.probabilities(a)) == eltype(AliasTables.probabilities(b))
                        @test AliasTables.probabilities(a) == AliasTables.probabilities(b)
                    end
                end
            else
                @test a != b
                @test hash(a) != hash(b)
                @test AliasTables.probabilities(float, a) != AliasTables.probabilities(float, b)
                @test AliasTables.probabilities(a) != AliasTables.probabilities(b)
            end
        end
    end

    @testset "RegressionTests" begin
        RegressionTests.test(skip_unsupported_platforms=true)
    end
end

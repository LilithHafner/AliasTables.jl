using AliasTables
using Test, Aqua, RegressionTests
using Random, OffsetArrays, StableRNGs

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
        @test AliasTable{UInt8}([0x80, 0x80]) ==
              AliasTable{UInt8}([0x81, 0x81]) ==
              AliasTable([1,1]) ==
              AliasTable{UInt8}(UInt128[typemax(UInt64), typemax(UInt64)] .<< 5) ==
              AliasTable{UInt8}([typemax(UInt32), typemax(UInt32)]) ==
              AliasTable(Float16[1, 1])
        @test rand(AliasTable{UInt8}(fill(0x80, 2^18))) in 1:2^18
        @test AliasTable{UInt8}(vcat(fill(0x00, 2^8), 0x80, 0x80)) == # Issue #34
                     AliasTable(vcat(fill(0x00, 2^8), 0x80, 0x80))
        @test length(AliasTable([1, 2, 3])) == 3
        @test rand(AliasTable{UInt8}(vcat(fill(0, 500), 1))) == 501 # PR #36
    end

    @testset "Invalid weight error messages" begin
        @test_throws ArgumentError("found negative weight -1") AliasTable([1, -1])
        @test_throws ArgumentError("found negative weight -1") AliasTable([1, 1, -1])
        @test_throws ArgumentError("found negative weight -1") AliasTable([3, typemax(Int), -1, typemax(Int)], _normalize=false)
        @test_throws ArgumentError("all weights are zero") AliasTable([0, 0])
        @test_throws ArgumentError("all weights are zero") AliasTable([0])
        @test_throws ArgumentError("all weights are zero") AliasTable(UInt[0, 0])
        @test_throws ArgumentError("all weights are zero") AliasTable(UInt[0])
        @test_throws ArgumentError("weights must be non-empty") AliasTable(Int[])
        @test_throws ArgumentError("weights must be non-empty") AliasTable(UInt[])
        @test_throws ArgumentError("sum(weights) is too low") AliasTable(UInt[123, 456], _normalize=false)
        @test_throws ArgumentError("sum(weights) is too high") AliasTable(UInt[unsigned(3)<<62, unsigned(2)<<62, unsigned(3)<<62], _normalize=false)
        @test_throws ArgumentError("sum(weights) overflows") AliasTable(UInt[unsigned(3)<<62, unsigned(2)<<62, unsigned(3)<<62])
        @test AliasTables.probabilities(float, AliasTable(UInt[unsigned(3)<<61, unsigned(2)<<61, unsigned(3)<<61])) == [3,2,3] ./ 8
        @test AliasTables.probabilities(float, AliasTable(UInt[unsigned(3)<<61, unsigned(2)<<61, unsigned(3)<<61], _normalize=false)) == [3,2,3] ./ 8
        @test_throws ArgumentError("offset arrays are not supported but got an array with index other than 1") AliasTable(OffsetVector([1,2], 1))
        @test AliasTable(OffsetVector([1,2], 0)) == AliasTable([1,2])
        @test_throws ArgumentError("sum(weights) is too high") AliasTable{UInt8}(vcat(fill(0x00, 2^8), 0x80, 0x81), _normalize=false) # _lookup_alias_table
        @test_throws ArgumentError("sum(weights) is too low") AliasTable{UInt8}(vcat(fill(0x00, 2^8), 0x80, 0x7f), _normalize=false) # _lookup_alias_table
    end

    @testset "probabilities()" begin
        @test AliasTables.probabilities(float, AliasTable([1, 2, 3])) == [1, 2, 3]/6
        @test AliasTables.probabilities(float, AliasTable([1, 2, 3, 0, 0])) == [1, 2, 3, 0, 0]/6
        @test AliasTables.probabilities(AliasTable([1, 2, 3, 0, 2])) == [1, 2, 3, 0, 2] .<< 61
        p = AliasTables.probabilities(AliasTable{UInt8}(fill(0x80, 2^18)))
        @test p isa Vector{UInt8}
        @test sum(p) == 256
        @test all(∈(0:1), p)
    end

    @testset "sample()" begin
        @test Base.hasmethod(AliasTables.sample, Tuple{UInt, AliasTable{UInt, Int}})
        @test !Base.hasmethod(AliasTables.sample, Tuple{Random.MersenneTwister, AliasTable{UInt, Int}})
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
        @testset "Sampling when there are more than typemax(T) weights" begin
            let at = AliasTable{UInt8}(vcat([1, 1], fill(0, 2^8)))
                @test counts(Iterators.map(x -> AliasTables.sample(x, at), 0x00:0xff), 258) ==
                    vcat([128, 128], zeros(2^8))
            end
            let at = AliasTable{UInt8}(vcat([1], fill(0, 2^8)))
                @test counts(Iterators.map(x -> AliasTables.sample(x, at), 0x00:0xff), 257) ==
                    vcat(256, zeros(2^8))
            end
        end
    end

    @testset "Equality and hashing" begin
        a = AliasTable([1, 2, 3])
        b = AliasTable([1, 2, 3, 0, 0])
        @test a != b
        @test a.probability_alias != b.probability_alias

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
            ],[
                AliasTable{UInt8}(vcat(fill(0x00, 2^8), 0x80, 0x80)), # Issue #34
                AliasTable(vcat(fill(0x00, 2^8), 0x80, 0x80))
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

    @testset "Display" begin
        function repr_test(x, str)
            @test repr(x) == str
            @test eval(Meta.parse(str)) == x
        end

        repr_test(AliasTable([1, 2, 3]), "AliasTable([0x2aaaaaaaaaaaaaab, 0x5555555555555555, 0x8000000000000000])")
        repr_test(AliasTable([1, 2, 5]), "AliasTable([0x2000000000000000, 0x4000000000000000, 0xa000000000000000])")
        repr_test(AliasTable{UInt8}([1, 2, 5, 0, 0]), "AliasTable{UInt8}([0x20, 0x40, 0xa0, 0x00, 0x00])")
        repr_test(AliasTable{UInt16, Int}([1, 2, 3]), "AliasTable{UInt16}([0x2aab, 0x5555, 0x8000])")
        repr_test(AliasTable{UInt64, Int16}([1, 2, 3]), "AliasTable{UInt64, Int16}([0x2aaaaaaaaaaaaaab, 0x5555555555555555, 0x8000000000000000])")
        repr_test(AliasTable{UInt16, Int16}([1, 2, 3]), "AliasTable{UInt16, Int16}([0x2aab, 0x5555, 0x8000])")

        if VERSION >= v"1.6" # Below this version, Base.typeinfo_implicit is not defined, so we display with redundant type info.
            repr_test([AliasTable([1, 2, 5]), AliasTable([0.0,1.0])], "[AliasTable([0x2000000000000000, 0x4000000000000000, 0xa000000000000000]), AliasTable([0x0000000000000000, 0xffffffffffffffff])]")
            repr_test([AliasTable{UInt16}([1, 2, 5]), AliasTable{UInt16}([0.0,1.0])]::Vector{AliasTable{UInt16, Int}}, "[AliasTable{UInt16}([0x2000, 0x4000, 0xa000]), AliasTable{UInt16}([0x0000, 0xffff])]")
            repr_test([AliasTable{UInt16}([1, 2, 5]), AliasTable([0.0,1.0])]::Vector{<:AliasTable}, "[AliasTable{UInt16}([0x2000, 0x4000, 0xa000]), AliasTable([0x0000000000000000, 0xffffffffffffffff])]")
        end
    end

    @testset "Consistency across Julia version and os" begin
        # This package does not guarantee reproducibility across minor versions of itself
        # so feel free to update the results here if the random stream changes.
        # However, we currentlty do prodcuce the same random stream on all Julia versions
        # and OSs. This test, running in CI, confirms that.

        at = AliasTable([0.4487, 0.46, 0.402, 0.062, 0.669, 0.1103, 0.9893, 0.1807, 0.7676])
        rng = StableRNGs.StableRNG(1729)
        @test rand(rng, at, 19) == [9, 3, 2, 2, 9, 7, 7, 9, 7, 5, 8, 5, 9, 7, 3, 2, 9, 9, 5]
    end

    @testset "Misc" begin
        probability_alias = AliasTables.Memory{Tuple{UInt8, Int}}(undef, 2)
        AliasTables._alias_table!(probability_alias, (0x01, 0xff)) == AliasTable([1,255])
    end

    @testset "RegressionTests" begin
        "CI" ∈ keys(ENV) && RegressionTests.test(skip_unsupported_platforms=true)
    end
end

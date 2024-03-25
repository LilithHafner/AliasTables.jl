# OffsetTables

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://LilithHafner.github.io/OffsetTables.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://offsettables.lilithhafner.com/dev)
[![Build Status](https://github.com/LilithHafner/OffsetTables.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/LilithHafner/OffsetTables.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/LilithHafner/OffsetTables.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/LilithHafner/OffsetTables.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/O/OffsetTables.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/O/OffsetTables.html)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

OffsetTables provides the `OffsetTable` type, which is an object that defines a probability
distribution over `1:n` for some `n`. They are efficient to construct and very efficient to
sample from.

An offset table can be combined with a dense vector of values to create a discrete
distribution over anything.

Internally, OffsetTables define a mapping from an unsigned integer type to the sampling
domain. To get a random sample according to the OffsetTable's distribution, one must provide
a random unsigned integer uniformly at random. One can also provide a `Random.AbstractRNG`
object instead and a random unsigned integer will be generated using that rng. When using
the random API, this latter approach is taken.

```julia
julia> using OffsetTables

julia> ot = OffsetTable([5,10,1])
OffsetTable([(0x2000000000000000, 1), (0x0000000000000000, 0), (0x3000000000000000, -2), (0x4000000000000000, -2)])

julia> rand(ot, 10)
10-element Vector{Int64}:
 2
 1
 2
 2
 2
 2
 1
 1
 3
 2

julia> using Chairmarks

julia> @b ot rand
2.990 ns

julia> @b rand(UInt)
2.734 ns

julia> @b rand(1000) OffsetTable
8.323 μs (5 allocs: 23.906 KiB)

julia> @b OffsetTable(rand(1000)) rand(_, 1000)
1.420 μs (3 allocs: 7.875 KiB)

julia> ot = OffsetTable{UInt16}([5,10,1])
OffsetTable{UInt16}([(0x2000, 1), (0x0000, 0), (0x3000, -2), (0x4000, -2)])

julia> countmap(OffsetTables.sample(x, ot) for x in typemin(UInt16):typemax(UInt16))
Dict{Any, Int64} with 3 entries:
  2 => 40960
  3 => 4096
  1 => 20480

julia> values(ans) .* 16/2^16
3-element Vector{Float64}:
 10.0
  1.0
  5.0
```
module AliasTables

using Random

export AliasTable
VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public sample, probabilities"))

const Memory = isdefined(Base, :Memory) ? Base.Memory : Vector # VERSION <= 1.10

if isdefined(Base, :top_set_bit)
    const top_set_bit = Base.top_set_bit
else
    top_set_bit(x::Integer) = 64 - leading_zeros(UInt64(x)) # VERSION <= 1.9
end

if isdefined(Base, :require_one_based_indexing) # VERSION == 1.0
    const require_one_based_indexing = Base.require_one_based_indexing
else
    require_one_based_indexing(A...) = !Base.has_offset_axes(A...) || throw(ArgumentError("offset arrays are not supported but got an array with index other than 1"))
end

"""
    AliasTable{T<:Unsigned=UInt, I<:Integer=Int}(weights::AbstractVector{<:Real})

An efficient data structure for sampling from a discrete distribution.

Maps every value representable by `T` to a value of type `I` in `eachindex(wights)` such
that the number of values maped to a given index of `weights` is proportional to the value
at that index.

The mapping can be accessed directly via
[`AliasTables.sample(x::T, at::AliasTable{T, I})`](@ref AliasTables.sample)
or indirectly via the `Random` API: `rand(at)`, `rand(rng, at)`, `rand(at, dims...)`, etc.

# Example

```jldoctest; filter=[r" [1-3]"]
julia> at = AliasTable([1, 3, 1])
AliasTable([0x3333333333333334, 0x9999999999999999, 0x3333333333333333])

julia> rand(at, 5)
5-element Vector{Int64}:
 2
 3
 2
 2
 3
```
"""
struct AliasTable{T <: Unsigned, I <: Integer}
    mask::T
    probability_alias::Memory{Tuple{T, I}}
    length::I
    """
        _AliasTable(probability_alias::Memory{Tuple{T, I}})

    Construct an `AliasTable` from a `Memory` of `(probability, alias)` pairs.

    Callers are responsible for ensuring that

        probability_alias == AliasTable(probabilities(_AliasTable(probability_alias))).probability_alias

    i.e. that wherever they got the `probability_alias` from, it is a form that would be
    produced by the internal default construcors (which are subject to change).

    If callers fail to do this, then equality and hashing may be broken.
    """
    function _AliasTable(probability_alias::Memory{Tuple{T, I}}, len) where {T, I}
        shift = 8sizeof(T) - top_set_bit(length(probability_alias)) + 1
        mask = (one(T) << shift) - one(T)
        new{T, I}(mask, probability_alias, len)
    end
    global _AliasTable
end

AliasTable(weights::AbstractVector{<:Real}; _normalize=true) = AliasTable{UInt64, Int}(weights; _normalize=_normalize)
AliasTable{T}(weights::AbstractVector{<:Real}; _normalize=true) where T <: Unsigned = AliasTable{T, Int}(weights; _normalize=_normalize)
function AliasTable{T, I}(weights; _normalize=true) where {T <: Unsigned, I <: Integer}
    require_one_based_indexing(weights)
    if _normalize
        (is_constant, sm) = checked_sum(weights)
        if is_constant
            _constant_alias_table(T, I, sm, length(weights))
        elseif sm-true == typemax(T) # pre-normalized
            _alias_table(T, I, weights)
        else
            norm = normalize_to_uint(T, weights, sm)
            _alias_table(T, I, norm)
        end
    else
        _alias_table(T, I, weights)
    end
end

function _constant_alias_table(::Type{T}, ::Type{I}, index, length) where {I, T}
    bitshift = top_set_bit(index - 1)
    len = 1 << bitshift
    points_per_cell = one(T) << (8*sizeof(T) - bitshift) # typemax(T)+1 / len

    probability_alias = Memory{Tuple{T, I}}(undef, len)
    for i in 1:len
        probability_alias[i] = (points_per_cell, index-i)
    end
    _AliasTable(probability_alias, length)
end

function throw_on_negatives(weights)
    for w in weights
        w < 0 && throw(ArgumentError("found negative weight $w"))
    end
end
function get_only_nonzero(weights)
    only_nonzero = -1
    for (i, w) in enumerate(weights)
        if w > 0
            if only_nonzero == -1
                only_nonzero = i
            else
                only_nonzero = -2
                break
            end
        end
    end
    only_nonzero == -1 && throw(ArgumentError("all weights are zero"))
    only_nonzero
end

"Like Iterators.Take, but faster"
struct HotTake{T<:Array}
    xs::T
    n::Int
    HotTake(xs::Array, n::Int) = new{typeof(xs)}(xs, min(n, length(xs)))
end
Base.iterate(A::HotTake, i=1) = ((i - 1)%UInt < A.n%UInt ? (@inbounds A.xs[i], i + 1) : nothing)
Base.length(A::HotTake) = A.n
hot_take(xs::Array, n) = HotTake(xs, n)
hot_take(xs, n) = Iterators.take(xs, n)

function _alias_table(::Type{T}, ::Type{I}, weights0) where {I, T}
    throw_on_negatives(weights0)
    onz = get_only_nonzero(weights0)
    onz == -2 || return _constant_alias_table(T, I, onz, length(weights0))

    weights = hot_take(weights0, findlast(!iszero, weights0))
    bitshift = top_set_bit(length(weights) - 1)
    len = 1 << bitshift # next_or_eq_power_of_two(length(weights))
    points_per_cell = one(T) << (8*sizeof(T) - bitshift) # typemax(T)+1 / len

    probability_alias = Memory{Tuple{T, I}}(undef, len)

    # @show sum(weights)
    # @show weights

    enum_weights = enumerate(weights)
    (surplus_i, surplus_desired), surplus_state = (thirsty_i, thirsty_desired), thirsty_state = iterate(enumerate(weights))
    current_i = surplus_i
    current_desired = surplus_desired

    while true
        # @show current_i, current_desired, points_per_cell
        if current_desired < points_per_cell # Surplus (strict)
            while true                                                            # Find the next thirsty cell
                ix = iterate(enum_weights, thirsty_state)
                ix === nothing && throw(ArgumentError("sum(weights) is too low")) # If there is no thirsty cell, there are more points than requsted by weights.
                (thirsty_i, thirsty_desired), thirsty_state = ix
                thirsty_desired >= points_per_cell && break
            end
            excess = points_per_cell - current_desired                            # Assign this many extra points
            probability_alias[current_i] = (excess, thirsty_i-current_i)          # To the targeted cell
            current_i = thirsty_i                                                 # Now we have to make sure that thristy cell gets exactly what it wants and no more
            current_desired = thirsty_desired - excess                            # It wants what it wants and hasn't already been transferred
        else                                 # Thirsty (loose)
            while true                                                            # Find the next surplus cell
                ix = iterate(enum_weights, surplus_state)
                ix === nothing && @goto break_outer                               # If there is no surplus cell, handle below
                (surplus_i, surplus_desired), surplus_state = ix
                surplus_desired < points_per_cell && break
            end
            excess = points_per_cell - surplus_desired                            # Assign this many extra points
            probability_alias[surplus_i] = (excess, current_i-surplus_i)          # From the cell with surplus to this cell
            current_desired -= excess                                             # We now don't want as many points (and may even no longer be thristy)
        end
    end
    @label break_outer
    # println()

    # There are no real surplus cells, but there may be synthetic surplus cells to round out
    # to a power of two. Those synthetic cells have structural weight of 0 so surplus value
    # of points_per_cell each. There may be remaining thristy cells, and the current cell is
    # thirsty.

    surplus_state_2 = length(weights)


    while true
        # @show current_i, current_desired, points_per_cell
        if current_desired < points_per_cell # Surplus (strict)
            while true                                                             # Find the next thirsty cell
                ix = iterate(enum_weights, thirsty_state)
                ix === nothing && throw(ArgumentError("sum(weights) is too low"))  # If there is no thirsty cell, there are more points than requsted by weights.
                (thirsty_i, thirsty_desired), thirsty_state = ix
                thirsty_desired >= points_per_cell && break
            end
            excess = points_per_cell - current_desired                             # Assign this many extra points
            probability_alias[current_i] = (excess, thirsty_i-current_i)          # To the targeted cell
            current_i = thirsty_i                                                  # Now we have to make sure that thristy cell gets exactly what it wants and no more
            current_desired = thirsty_desired - excess                             # It wants what it wants and hasn't already been transferred
        else                                 # Thirsty (loose)
            surplus_state_2 += true # Find the next surplus cell
            surplus_state_2 > len && break # If there is no surplus cell, handle below
            surplus_i = surplus_state_2
            excess = points_per_cell                                               # Assign all the points
            probability_alias[surplus_i] = (points_per_cell, current_i-surplus_i) # From the synthetic cell with surplus to this cell
            current_desired -= excess                                              # We now don't want as many points (and may even no longer be thristy)
        end
    end

    # @show probability_alias, points_per_cell, current_desired

    if points_per_cell < current_desired  # Strictly thirsty, and no surplus cells, so exceed the desired weight.
        throw(ArgumentError("sum(weights) is too high"))
    end
    probability_alias[current_i] = (0, 0)
    # Just right. There are no surplus cells left and no current surplus or thirst. All
    # that's left are future loosely thirsty cells, all of which should be a thirst of
    # exactly 0.
    while true         # Loop over all thirsty celss
        ix = iterate(enum_weights, thirsty_state)
        ix === nothing && break # Out of thirsty cells, yay!
        (thirsty_i, thirsty_desired), thirsty_state = ix
        points_per_cell < thirsty_desired && throw(ArgumentError("sum(weights) is too high")) # Strictly thirsty, with no surplus to draw from.
        points_per_cell == thirsty_desired && (probability_alias[thirsty_i] = (0, 0)) # loosely thirsty, but satisfied. Zero out the undef.
    end

    _AliasTable(probability_alias, length(weights0))
end

"""
    sample(x::T, at::AliasTable{T, I}) -> I

Sample from `at` using the seed `x`.

If `x` is chosen uniformly at random from the set of all values representable by `T` then
the output will be a random sample from the distribution represented by `at`. The mapping is
deterministic and not pseudo-random so for patterned input `x` the output will be patterned
as well.

See also [`AliasTable`](@ref), [`AliasTables.probabilities`](@ref)
"""
function sample(x::T, at::AliasTable{T, I}) where {T, I}
    shift = 8sizeof(T) - top_set_bit(length(at.probability_alias)) + 1
    cell = (x >> shift) + 1
    # @assert (one(T) << shift) - one(T) == at.mask
    val = x & at.mask
    @inbounds prob, alias = at.probability_alias[cell%Int]
    (((val < prob) * alias + cell)%I)::I
end

### Random API
Random.rand(rng::Random.AbstractRNG, at::Random.SamplerTrivial{<:AliasTable{T}}) where T = sample(rand(rng, T), at.self)
Random.gentype(::Type{AliasTable{T, I}}) where {T, I} = I

### Reconstruct probabilities
"""
    probabilities(at::AliasTable{T}) -> Vector{T}

Recover the exact sampling weights from a given `AliasTable`. The returned values will
sum to one more than `typemax(T)`, unless `at` is a constant distribution (e.g.
`AliasTable([0,1,0])`), in which case the weights will sum to `typemax(T)`.

See also [`AliasTable`](@ref), [`AliasTables.sample`](@ref)

# Examples

```jldoctest
julia> at = AliasTable([1, 3, 1])
AliasTable([0x3333333333333334, 0x9999999999999999, 0x3333333333333333])

julia> AliasTables.probabilities(at)
3-element Vector{UInt64}:
 0x3333333333333334
 0x9999999999999999
 0x3333333333333333

julia> AliasTables.probabilities(AliasTable([0, 1, 0]))
3-element Vector{UInt64}:
 0x0000000000000000
 0xffffffffffffffff
 0x0000000000000000
```
"""
function probabilities(at::AliasTable{T}) where T
    bitshift = top_set_bit(length(at.probability_alias) - 1)
    points_per_cell = one(T) << (8*sizeof(T) - bitshift)#typemax(T)+1 / len
    probs = zeros(T, at.length)
    for (i, (prob, alias)) in enumerate(at.probability_alias)
        probs[i + alias] += prob
        keep = points_per_cell - prob
        iszero(keep) || (probs[i] += keep)
    end
    if all(iszero, probs)
        probs[sample(zero(T), at)] = typemax(T)
    end
    probs
end

"""
    probabilities(float, at::AliasTable{T}) -> Vector{<:AbstractFloat}

Return the sampling probabilities of `at`. The returned vector will sum to 1.0, up to
rounding error.

# Example

```jldoctest
julia> AliasTables.probabilities(float, AliasTable([1, 3, 1]))
3-element Vector{Float64}:
 0.2
 0.6
 0.2
"""
probabilities(::typeof(float), at::AliasTable{T}) where T =
    probabilities(at) ./ (float(typemax(T))+1)

### Length accessor
"""
    length(at::AliasTable)

Get the number of weights that `at` was constructed with, including trailing zeros.

# Example

```jldoctest
julia> length(AliasTable([1, 3, 0]))
3
```
"""
Base.length(at::AliasTable) = at.length

### Show
function Base.show(io::IO, at::AliasTable{T, I}) where {T, I}
    print(io, AliasTable)
    if I != Int
        print(io, "{", T, ", ", I, "}")
    elseif T != UInt64
        print(io, "{", T, "}")
    end
    print(io, "(")
    print(IOContext(io, :typeinfo=>Vector{T}), probabilities(at))
    print(io, ")")
end
isdefined(Base, :typeinfo_implicit) && (Base.typeinfo_implicit(::Type{<:AliasTable}) = true)

### Equality and hashing

# These naive implementations are equivalent to computing equality
# based on probabilities because the constrors are deterministic w.r.t
# the input weights exclusind trailing zeros and the length is tracked.
Base.:(==)(at1::AliasTable{T}, at2::AliasTable{T}) where T = at1.length == at2.length && at1.probability_alias == at2.probability_alias
function Base.:(==)(at1::AliasTable{T1}, at2::AliasTable{T2}) where {T1, T2}
    at1.length == at2.length || return false
    length(at1.probability_alias) == length(at2.probability_alias) || return false
    bitshift = 8(sizeof(T1) - sizeof(T2))
    for (po1, po2) in zip(at1.probability_alias, at2.probability_alias)
        po1[2] == po2[2] &&
        if bitshift > 0
            po1[1] == T1(po2[1]) << bitshift
        else
            T2(po1[1]) << -bitshift == po2[1]
        end || return false
    end
    true
end
struct MapVector{T, F, P} <: AbstractVector{T}
    parent::P
    f::F
end
Base.size(mv::MapVector) = size(mv.parent)
Base.getindex(mv::MapVector{T, F, P}, i) where {T, F, P} = mv.f(mv.parent[i])
function Base.hash(at::AliasTable, h::UInt)
    h ⊻= Sys.WORD_SIZE == 32 ? 0x7719cd5e : 0x0a0c5cfeeb10f090
    h = hash(at.length, h)
    # isempty(at.probability_alias) && return hash(0, h) # This should never happen, but it makes first not throw.
    po1 = first(at.probability_alias)
    norm(x) = (ldexp(float(x[1]), -8sizeof(x[1])), x[2])
    hash(MapVector{typeof(norm(po1)), typeof(norm), typeof(at.probability_alias)}(at.probability_alias, norm), h)
end

## Normalization

maybe_unsigned(x) = x # this is necessary to avoid calling unsigned on BigInt, Flaot64, etc.
maybe_unsigned(x::Base.BitSigned) = unsigned(x)

maybe_add_with_overflow(x::Base.BitInteger, y::Base.BitInteger) = Base.Checked.add_with_overflow(x, convert(typeof(x), y))
maybe_add_with_overflow(x, y) = x+y, false

####

# 2-4 passes (skip first two if nomralize = false)
# Initial sum, check for overflow, negatives, allzero & exactness
# If not exact, compute the mapped sum and the error there
# Dual pass for construction

# First pass
widen_to_word(x) = x
widen_to_word(x::Base.BitSignedSmall) = Int(x)
widen_to_word(x::Base.BitUnsignedSmall) = Int(x)
sum_prepare(x) = maybe_unsigned(widen_to_word(x))

function checked_sum(weights)
    xi = iterate(weights)
    xi === nothing && throw(ArgumentError("weights must be non-empty"))
    x, i = xi
    nonzero_index = 1
    while iszero(x)
        nonzero_index += 1
        xi = iterate(weights, i)
        xi === nothing && throw(ArgumentError("all weights are zero"))
        x, i = xi
    end
    x < 0 && throw(ArgumentError("found negative weight $x"))
    x0 = x
    while true
        xi = iterate(weights, i)
        xi === nothing && return (true, nonzero_index)
        x, i = xi
        iszero(x) || break
    end
    # There are two nonzero elements
    x < 0 && throw(ArgumentError("found negative weight $x"))
    sm, overflow = maybe_add_with_overflow(sum_prepare(x0), sum_prepare(x))
    overflow && !iszero(sm) && throw(ArgumentError("sum(weights) overflows"))
    while true
        xi = iterate(weights, i)
        xi === nothing && break
        x, i = xi
        x < 0 && throw(ArgumentError("found negative weight $x"))
        sm, o = maybe_add_with_overflow(sm, sum_prepare(x))
        overflow |= o
        overflow && !iszero(sm) && throw(ArgumentError("sum(weights) overflows"))
    end
    isfinite(sm) || throw(ArgumentError("sum(weights) == $sm which is not finite"))
    (false, sm)
end

# Second phase

# Slower than allocating:
# function normalize_to_uint_lazy_frac_div(::Type{T}, weights, sm) where T <: Unsigned
#     sm2 = zero(T)
#     for x in weights
#         sm2 += frac_div(T(x), sm)
#     end
#     sm2_copy = sm2 # lolz https://github.com/JuliaLang/julia/issues/15276

#     bonus = typemax(sm2)-sm2+1
#     (frac_div(T(x), sm) + (i <= bonus) for (i,x) in enumerate(weights))
# end

widen_float(T, x) = typemax(T) < floatmax(x) ? x : widen_float(T, widen(x))
function normalize_to_uint(::Type{T}, v, sm) where {T <: Unsigned}
    if sm isa AbstractFloat
        shift = 8sizeof(T)-exponent(sm + sqrt(eps(sm)))-1
        v2 = res = [floor(T, ldexp(widen_float(T, x), shift)) for x in v]
        onz = get_only_nonzero(v2)
        onz != -2 && return res
        sm2 = sum(res)
    else
        v2 = v
        res = Vector{T}(undef, length(v))
        sm2 = sm
    end

    sm3 = zero(T)

    T2 = promote_type(widen(T), typeof(sm2))
    for (i,x) in enumerate(v2)
        # @assert x < sm2
        # @assert sm2 != 0
        val = div(T2(maybe_unsigned(x)) << 8sizeof(T), sm2) % T
        sm3 += val
        res[i] = val
    end

    sm3 == 0 && any(!iszero(res)) && return res

    for i in sm3:typemax(sm3)
        res[typemax(sm3)-i+1] += true
    end

    res
end

end

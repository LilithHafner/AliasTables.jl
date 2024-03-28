module AliasTables

using Random

export AliasTable
VERSION >= v"1.11.0-DEV.469" && eval(Meta.parse("public sample"))

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
    AliasTable{T<:Unsigned=UInt, I<:Integer=Int}(weights::AbstractVector{<:Real}; normalize=true)

An efficient data structure for sampling from a discrete distribution.

Maps every value representable by `T` to a value of type `I` in `eachindex(wights)` such
that the number of values maped to a given index of `weights` is proportional to the value
at that index.

The mapping can be accessed directly via
[`AliasTables.sample(x::T, ot::AliasTable{T, I})`](@ref AliasTables.sample(::T, ::AliasTable{T, I}) where {T, I})
or indirectly by passing a random number generator which will be used to generate a
random input of type `T` for you via
[`AliasTables.sample(rng::Random.AbstractRNG, ot::AliasTable{T, I})`](@ref AliasTables.sample(::Random.AbstractRNG, ::AliasTable))
or simply via the `Random` API: `rand(ot)`, `rand(rng, ot)`, `rand(ot, dims...)`, etc.

Set `normalize = false` for incrased performance when the weights are already normalized to
sum to exactly the number of values representable by `T` (i.e. `typemax(T)+1`). A different
sum will result in an error unless exactly one weight is non-zero, in which case the sum is
not checked and the `AliasTable` represents a constant distribution which always produces
the index of the nonzero weight.
"""
struct AliasTable{T <: Unsigned, I <: Integer}
    probability_alias::Memory{Tuple{T, I}}
    """
        _AliasTable(probability_alias::Memory{Tuple{T, I}})

    Construct an `AliasTable` from a `Memory` of `(probability, alias)` pairs.

    Callers are responsible for ensuring that

        probability_alias == AliasTable(probabilities(_AliasTable(probability_alias))).probability_alias

    i.e. that wherever they got the `probability_alias` from, it is a form that would be
    produced by the internal default construcors (which are subject to change).

    If callers fail to do this, then equality and hashing may be broken.
    """
    _AliasTable(probability_alias::Memory{Tuple{T, I}}) where {T, I} = new{T, I}(probability_alias)
    global _AliasTable
end

AliasTable(weights::AbstractVector{<:Real}; normalize=true) = AliasTable{UInt, Int}(weights; normalize=normalize)
AliasTable{T}(weights::AbstractVector{<:Real}; normalize=true) where T <: Unsigned = AliasTable{T, Int}(weights; normalize=normalize)
function AliasTable{T, I}(weights; normalize=true) where {T <: Unsigned, I <: Integer}
    # function _AliasTable(::Type{T}, ::Type{I}, weights; normalize=true) where {T <: Unsigned, I <: Integer}
    require_one_based_indexing(weights)
    if normalize
        (is_constant, sm) = checked_sum(weights)
        if is_constant
            _constant_alias_table(T, I, sm)
        elseif sm == 0 # pre-normalized
            _alias_table(T, I, weights)
        else
            # norm = normalize_to_uint_lazy_frac_div(T, weights, sm)
            norm = normalize_to_uint_frac_div(T, weights, sm)
            _alias_table(T, I, norm)
        end
    else
        _alias_table(T, I, weights)
    end
end

function _constant_alias_table(::Type{T}, ::Type{I}, index) where {I, T}
    bitshift = top_set_bit(index - 1)
    len = 1 << bitshift
    points_per_cell = one(T) << (8*sizeof(T) - bitshift)#typemax(T)+1 / len

    probability_alias = Memory{Tuple{T, I}}(undef, len)
    for i in 1:len
        probability_alias[i] = (points_per_cell, index-i)
    end
    _AliasTable(probability_alias)
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
        elseif w < 0
            throw(ArgumentError("found negative weight $w at index $i"))
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
    onz = get_only_nonzero(weights0)
    onz == -2 || return _constant_alias_table(T, I, onz)

    # weights = Iterators.take(weights0, findlast(!iszero, weights0))
    weights = hot_take(weights0, findlast(!iszero, weights0))
    # weights = weights0
    # while iszero(last(weights))
        # pop!(weights)
    # end
    bitshift = top_set_bit(length(weights) - 1)
    len = 1 << bitshift#next_or_eq_power_of_two(length(weights))
    points_per_cell = one(T) << (8*sizeof(T) - bitshift)#typemax(T)+1 / len

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
            while true                                                             # Find the next surplus cell
                ix = iterate(enum_weights, surplus_state)
                ix === nothing && @goto break_outer                                # If there is no surplus cell, handle below
                (surplus_i, surplus_desired), surplus_state = ix
                surplus_desired < points_per_cell && break
            end
            excess = points_per_cell - surplus_desired                             # Assign this many extra points
            probability_alias[surplus_i] = (excess, current_i-surplus_i)          # From the cell with surplus to this cell
            current_desired -= excess                                              # We now don't want as many points (and may even no longer be thristy)
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

    _AliasTable(probability_alias)
end

"""
    sample(x::T, ot::AliasTable{T, I}) -> I

Sample from `ot` using the seed `x`.

If `x` is chosen uniformly at random from the set of all values representable by `T` then
the output will be a random sample from the distribution represented by `ot`. The mapping is
deterministic and not pseudo-random so for patterned input `x` the output will be patterned
as well.

See also [`AliasTable`](@ref)
"""
function sample(x::T, ot::AliasTable{T, I}) where {T, I}
    count_ones(length(ot.probability_alias)) == 1 || return zero(I) # This should never happen, but it makes the @inbounds safe
    shift = 8sizeof(T) - top_set_bit(length(ot.probability_alias)) + 1
    cell = (x >> shift) + 1
    val = x & ((one(T) << shift) - one(T))
    @inbounds prob, alias = ot.probability_alias[cell%Int]
    # (val < prob ? (alias+cell) : I(cell))::I
    # I((val < prob) * alias + cell)::I
    (((val < prob) * alias + cell)%I)::I
end

"""
    sample(rng::Random.AbstractRNG, ot::AliasTable{T, I}) -> I

Sample from `ot` using randomness drawn from `rng`.

Produces a random sample from the distribution represented by `ot`.

See also [`AliasTable`](@ref), `Random.rand`
"""
function sample(rng::Random.AbstractRNG, ot::AliasTable{T, I}) where {T, I}
    sample(rand(rng, T), ot)
end

### Random API
Random.rand(rng::Random.AbstractRNG, ot::Random.SamplerTrivial{<:AliasTable}) = sample(rng, ot.self)
Random.gentype(::Type{AliasTable{T, I}}) where {T, I} = I

### Reconstruct probabilities
function probabilities(ot::AliasTable{T}) where T
    bitshift = top_set_bit(length(ot.probability_alias) - 1)
    points_per_cell = one(T) << (8*sizeof(T) - bitshift)#typemax(T)+1 / len
    probs = zeros(T, length(ot.probability_alias))
    for (i, (prob, alias)) in enumerate(ot.probability_alias)
        probs[i + alias] += prob
        probs[i] += points_per_cell - prob
    end
    li = findlast(!iszero, probs)
    if li != nothing
        resize!(probs, li)
    else # overflow
        resize!(probs, sample(zero(T), ot))
        probs[end] = typemax(T)
    end
    probs
end

probabilities(::typeof(float), ot::AliasTable{T}) where T =
    probabilities(ot) ./ (float(typemax(T))+1)


### Show
function Base.show(io::IO, ot::AliasTable{T, I}) where {T, I}
    print(io, AliasTable)
    if get(io, :typeinfo, nothing) != AliasTable{T, I} && (T != UInt || I != Int)
        if I == Int
            print(io, "{", T, "}")
        else
            print(io, "{", T, ", ", I, "}")
        end
    end
    print(io, "(")
    # print(IOContext(io, :typeinfo=>Vector{T}), probabilities(ot))
    print(IOContext(io, :typeinfo=>Memory{Tuple{T, I}}), ot.probability_alias)
    print(io, ")")
end

### Equality and hashing

# These naive implementations are equivalent to computing equality based on
# `WithZeros(probabilities, Inf)` because the constrors are deterministic w.r.t
# the input weights exclusind trailing zeros.
Base.:(==)(ot1::AliasTable{T}, ot2::AliasTable{T}) where T = ot1.probability_alias == ot2.probability_alias
function Base.:(==)(ot1::AliasTable{T1}, ot2::AliasTable{T2}) where {T1, T2}
    bitshift = 8(sizeof(T1) - sizeof(T2))
    length(ot1.probability_alias) == length(ot2.probability_alias) || return false
    for (po1, po2) in zip(ot1.probability_alias, ot2.probability_alias)
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
function Base.hash(ot::AliasTable, h::UInt)
    h ⊻= Sys.WORD_SIZE == 32 ? 0x7719cd5e : 0x0a0c5cfeeb10f090
    # isempty(ot.probability_alias) && return hash(0, h) # This should never happen, but it makes first not throw.
    po1 = first(ot.probability_alias)
    norm(x) = (ldexp(float(x[1]), -8sizeof(x[1])), x[2])
    hash(MapVector{typeof(norm(po1)), typeof(norm), typeof(ot.probability_alias)}(ot.probability_alias, norm), h)
end

## Normalization

maybe_unsigned(x) = x # this is necessary to avoid calling unsigned on BigInt, Flaot64, etc.
maybe_unsigned(x::Base.BitSigned) = unsigned(x)

maybe_add_with_overflow(x::Base.BitInteger, y::Base.BitInteger) = Base.Checked.add_with_overflow(x, convert(typeof(x), y))
maybe_add_with_overflow(x, y) = x+y, false

function normalize_to_uint(::Type{T}, v::AbstractVector{<:Real}) where {T <: Unsigned}
    only_nonzero = -1
    sm = maybe_unsigned(Base.add_sum(zero(eltype(v)), zero(eltype(v))))
    overflow = false
    for (i, w) in enumerate(v)
        if only_nonzero != 0 && w > 0
            if only_nonzero == -1
                only_nonzero = i
            else
                only_nonzero = 0
            end
        elseif w < 0
            throw(ArgumentError("found negative weight $w at index $i"))
        end
        sm, o = maybe_add_with_overflow(sm, w)
        sm != 0 && o && throw(ArgumentError("sum(weights) overflows"))
        overflow |= o
    end
    only_nonzero == -1 && throw(ArgumentError("all weights are zero"))
    if only_nonzero != 0
        res = zeros(T, length(v))
        res[only_nonzero] = 1
        return res
    end
    overflow && typemax(T) == typemax(eltype(v)) && return v
    overflow && throw(ArgumentError("sum(weights) overflows"))

    length(v) <= 1 && return ones(T, length(v))
    sm2 = sm # https://github.com/JuliaLang/julia/issues/15276
    res = [x <= typemax(T) ? floor(T, x) : typemax(T) for x in (ldexp(x/sm2, 8sizeof(T)) for x in v)] # TODO: make this lazy & non-allocating

    count_nonzero = 0
    for x in res
        if x != 0
            count_nonzero += 1
            count_nonzero > 1 && break
        end
    end
    @assert count_nonzero != 0
    count_nonzero == 1 && return res # No need to normalize

    leftover = signed(sum(widen, res) - typemax(T)-1)
    argmx = argmax(res)
    if (res[argmx] -= leftover) < 0
        throw(ArgumnetError("normalization failed")) # TODO: eliminate this
    end
    res
end

function frac_div(x::T, y::T) where T <: Unsigned
    # @assert x < y
    # @assert y != 0
    div(widen(x) << 8sizeof(T), y) % T
end

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

function normalize_to_uint_frac_div(::Type{T}, v, sm=sum(T,v)) where {T <: Unsigned}
    if sm isa AbstractFloat
        shift = 8sizeof(T)-exponent(sm + sqrt(eps(sm)))-1
        v2 = res = [floor(T, ldexp(x, shift)) for x in v]
        onz = get_only_nonzero(v2)
        onz != -2 && return res
        sm2 = sum(res)
    else
        v2 = v
        res = Vector{T}(undef, length(v))
        sm2 = sm
    end

    sm3 = zero(T)

    for (i,x) in enumerate(v2)
        val = frac_div(T(x), T(sm2))
        sm3 += val
        res[i] = val
    end

    sm3 == 0 && return res

    for i in sm3:typemax(sm3)
        res[typemax(sm3)-i+1] += true
    end

    res
end

end

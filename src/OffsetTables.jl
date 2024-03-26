module OffsetTables

using Random

export OffsetTable
public sample

"""
    OffsetTable{T<:Unsigned=UInt, I<:Integer=Int}(weights::AbstractVector{<:Real}; normalize=true)

An efficient data structure for sampling from a discrete distribution.

Maps every value representable by `T` to a value of type `I` in `eachindex(wights)` such
that the number of values maped to a given index of `weights` is proportional to the value
at that index.

The mapping can be accessed directly via
[`OffsetTables.sample(x::T, ot::OffsetTable{T, I})`](@ref OffsetTables.sample(::T, ::OffsetTable{T, I}) where {T, I})
or indirectly by passing a random number generator which will be used to generate a
random input of type `T` for you via
[`OffsetTables.sample(rng::Random.AbstractRNG, ot::OffsetTable{T, I})`](@ref OffsetTables.sample(::Random.AbstractRNG, ::OffsetTable))
or simply via the `Random` API: `rand(ot)`, `rand(rng, ot)`, `rand(ot, dims...)`, etc.

Set `normalize = false` for incrased performance when the weights are already normalized to
sum to exactly the number of values representable by `T` (i.e. `typemax(T)+1`). A different
sum will result in an error unless exactly one weight is non-zero, in which case the sum is
not checked and the `OffsetTable` represents a constant distribution which always produces
the index of the nonzero weight.
"""
struct OffsetTable{T <: Unsigned, I <: Integer}
    probability_offset::Memory{Tuple{T, I}}
    """
        _OffsetTable(probability_offset::Memory{Tuple{T, I}})

    Construct an `OffsetTable` from a `Memory` of `(probability, offset)` pairs.

    Callers are responsible for ensuring that

        probability_offset == OffsetTable(probabilities(_OffsetTable(probability_offset))).probability_offset

    i.e. that wherever they got the `probability_offset` from, it is a form that would be
    produced by the internal default construcors (which are subject to change).

    If callers fail to do this, then equality and hashing may be broken.
    """
    _OffsetTable(probability_offset::Memory{Tuple{T, I}}) where {T, I} = new{T, I}(probability_offset)
    global _OffsetTable
end

OffsetTable(weights::AbstractVector{<:Real}; normalize=true) = OffsetTable{UInt, Int}(weights; normalize)
OffsetTable{T}(weights::AbstractVector{<:Real}; normalize=true) where T <: Unsigned = OffsetTable{T, Int}(weights; normalize)
function OffsetTable{T, I}(weights::AbstractVector{<:Real}; normalize=true) where {T <: Unsigned, I <: Integer}
    Base.require_one_based_indexing(weights)
    isempty(weights) && throw(ArgumentError("weights must be non-empty"))
    if normalize
        _offset_table(T, I, normalize_to_uint(T, weights))
    else
        _offset_table(T, I, weights)
    end
end

struct WithZeros{T, P <: AbstractVector{T}} <: AbstractVector{T}
    parent::P
    length::Int
end
Base.size(wz::WithZeros) = (wz.length,)
Base.getindex(wz::WithZeros{T}, i) where T = i <= length(wz.parent) ? wz.parent[i] : zero(T)

function _constant_offset_table(::Type{I}, ::Type{T}, index) where {I, T}
    bitshift = Base.top_set_bit(index - 1)
    len = 1 << bitshift#next_or_eq_power_of_two(length(weights))
    points_per_cell = one(T) << (8*sizeof(T) - bitshift)#typemax(T)+1 / len

    probability_offset = Memory{Tuple{T, I}}(undef, len)
    for i in 1:len
        probability_offset[i] = (points_per_cell, index-i)
    end
    _OffsetTable(probability_offset)
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

function _offset_table(::Type{T}, ::Type{I}, weights::AbstractVector{<:Unsigned}) where {I, T}
    onz = get_only_nonzero(weights)
    onz == -2 || return _constant_offset_table(I, T, onz)

    bitshift = Base.top_set_bit(findlast(!iszero, weights) - 1)
    len = 1 << bitshift#next_or_eq_power_of_two(length(weights))
    points_per_cell = one(T) << (8*sizeof(T) - bitshift)#typemax(T)+1 / len

    probability_offset = Memory{Tuple{T, I}}(undef, len)

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
            probability_offset[current_i] = (excess, thirsty_i-current_i)          # To the targeted cell
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
            probability_offset[surplus_i] = (excess, current_i-surplus_i)          # From the cell with surplus to this cell
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
            probability_offset[current_i] = (excess, thirsty_i-current_i)          # To the targeted cell
            current_i = thirsty_i                                                  # Now we have to make sure that thristy cell gets exactly what it wants and no more
            current_desired = thirsty_desired - excess                             # It wants what it wants and hasn't already been transferred
        else                                 # Thirsty (loose)
            surplus_state_2 += true # Find the next surplus cell
            surplus_state_2 > len && break # If there is no surplus cell, handle below
            surplus_i = surplus_state_2
            excess = points_per_cell                                               # Assign all the points
            probability_offset[surplus_i] = (points_per_cell, current_i-surplus_i) # From the synthetic cell with surplus to this cell
            current_desired -= excess                                              # We now don't want as many points (and may even no longer be thristy)
        end
    end

    # @show probability_offset, points_per_cell, current_desired

    if points_per_cell < current_desired  # Strictly thirsty, and no surplus cells, so exceed the desired weight.
        throw(ArgumentError("sum(weights) is too high"))
    end
    probability_offset[current_i] = (0, 0)
    # Just right. There are no surplus cells left and no current surplus or thirst. All
    # that's left are future loosely thirsty cells, all of which should be a thirst of
    # exactly 0.
    while true         # Loop over all thirsty celss
        ix = iterate(enum_weights, thirsty_state)
        ix === nothing && break # Out of thirsty cells, yay!
        (thirsty_i, thirsty_desired), thirsty_state = ix
        points_per_cell < thirsty_desired && throw(ArgumentError("sum(weights) is too high")) # Strictly thirsty, with no surplus to draw from.
        probability_offset[thirsty_i] = (0, 0)
    end

    _OffsetTable(probability_offset)
end

"""
    sample(x::T, ot::OffsetTable{T, I}) -> I

Sample from `ot` using the seed `x`.

If `x` is chosen uniformly at random from the set of all values representable by `T` then
the output will be a random sample from the distribution represented by `ot`. The mapping is
deterministic and not pseudo-random so for patterned input `x` the output will be patterned
as well.

See also [`OffsetTable`](@ref)
"""
function sample(x::T, ot::OffsetTable{T, I}) where {T, I}
    count_ones(length(ot.probability_offset)) == 1 || return zero(I) # This should never happen, but it makes the @inbounds safe
    shift = 8sizeof(T) - Base.top_set_bit(length(ot.probability_offset)) + 1
    cell = (x >> shift) + 1
    val = x & ((one(T) << shift) - one(T))
    @inbounds prob, offset = ot.probability_offset[cell%Int]
    # (val < prob ? (offset+cell) : I(cell))::I
    # I((val < prob) * offset + cell)::I
    (((val < prob) * offset + cell)%I)::I
end

"""
    sample(rng::Random.AbstractRNG, ot::OffsetTable{T, I}) -> I

Sample from `ot` using randomness drawn from `rng`.

Produces a random sample from the distribution represented by `ot`.

See also [`OffsetTable`](@ref), `Random.rand`
"""
function sample(rng::Random.AbstractRNG, ot::OffsetTable{T, I}) where {T, I}
    sample(rand(rng, T), ot)
end

### Random API
Random.rand(rng::Random.AbstractRNG, ot::Random.SamplerTrivial{<:OffsetTable}) = sample(rng, ot.self)
Random.gentype(::Type{OffsetTable{T, I}}) where {T, I} = I

### Reconstruct probabilities
function probabilities(ot::OffsetTable{T}) where T
    bitshift = Base.top_set_bit(length(ot.probability_offset) - 1)
    points_per_cell = one(T) << (8*sizeof(T) - bitshift)#typemax(T)+1 / len
    probs = zeros(T, length(ot.probability_offset))
    for (i, (prob, offset)) in enumerate(ot.probability_offset)
        probs[i + offset] += prob
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

probabilities(::typeof(float), ot::OffsetTable{T}) where T =
    probabilities(ot) ./ (float(typemax(T))+1)


### Show
function Base.show(io::IO, ot::OffsetTable{T, I}) where {T, I}
    print(io, OffsetTable)
    if get(io, :typeinfo, nothing) != OffsetTable{T, I} && (T != UInt || I != Int)
        if I == Int
            print(io, "{", T, "}")
        else
            print(io, "{", T, ", ", I, "}")
        end
    end
    print(io, "(")
    # print(IOContext(io, :typeinfo=>Vector{T}), probabilities(ot))
    print(IOContext(io, :typeinfo=>Memory{Tuple{T, I}}), ot.probability_offset)
    print(io, ")")
end

### Equality and hashing

# These naive implementations are equivalent to computing equality based on
# `WithZeros(probabilities, Inf)` because the constrors are deterministic w.r.t
# the input weights exclusind trailing zeros.
Base.:(==)(ot1::OffsetTable{T}, ot2::OffsetTable{T}) where T = ot1.probability_offset == ot2.probability_offset
function Base.:(==)(ot1::OffsetTable{T1}, ot2::OffsetTable{T2}) where {T1, T2}
    bitshift = 8(sizeof(T1) - sizeof(T2))
    length(ot1.probability_offset) == length(ot2.probability_offset) || return false
    for (po1, po2) in zip(ot1.probability_offset, ot2.probability_offset)
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
function Base.hash(ot::OffsetTable, h::UInt)
    h ⊻= Sys.WORD_SIZE == 32 ? 0x7719cd5e : 0x0a0c5cfeeb10f090
    # isempty(ot.probability_offset) && return hash(0, h) # This should never happen, but it makes first not throw.
    po1 = first(ot.probability_offset)
    norm(x) = (ldexp(float(x[1]), -8sizeof(x[1])), x[2])
    hash(MapVector{typeof(norm(po1)), typeof(norm), typeof(ot.probability_offset)}(ot.probability_offset, norm), h)
end

## Mediocre float handling

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

function normalize_to_uint_widen(::Type{T}, v::AbstractVector{<:Integer}) where {T <: Unsigned}
    T2 = widen(T)
    den = sum(T2, v)
    num = T2(typemax(T)) + one(T2)
    res = Vector{T}(undef, length(v))
    for (i,x) in enumerate(v)
        chosen = Base.udiv_int(num*x, den)
        num -= chosen
        den -= x
        res[i] = chosen
    end
    res
end

function normalize_to_uint_dual_overflow(::Type{T}, v::AbstractVector{<:Integer}) where {T <: Unsigned}
    sm = sum(T, v)
    q, r = divrem(typemax(T), sm)
    if r == sm-true
        q += true
        r = zero(r)
    else
        r += true
    end
    # q, r = divrem(typemax(T)+1, sm)

    error = zero(T)
    res = Vector{T}(undef, length(v))
    for (i,x) in enumerate(v)
        tx = T(x)
        # @show q, r, tx, r*tx
        error += r*tx # TODO: check overflow
        carry, error = divrem(error, sm)
        res[i] = q*tx + carry
    end
    res
end

function mul_hi_lo(x::T, y::T) where T <: Unsigned
    hbits = 4sizeof(T)
    lo_mask = T(1) << hbits - 1
    lo_x = x & lo_mask
    lo_y = y & lo_mask
    hi_x = x >> hbits
    hi_y = y >> hbits

    lo = lo_x * lo_y
    hi = hi_x * hi_y

    mid = lo_x * hi_y
    lo += mid << hbits
    hi += mid >> hbits

    mid = hi_x * lo_y
    lo += mid << hbits
    hi += mid >> hbits

    hi, lo
end

function mul_hi_lo_2(x::T, y::T) where T <: Unsigned
    hbits = 4sizeof(T)
    lo_mask = T(1) << hbits - 1
    lo_x = x & lo_mask
    lo_y = y & lo_mask
    hi_x = x >> hbits
    hi_y = y >> hbits

    hi = hi_x * hi_y

    mid = lo_x * hi_y
    hi += mid >> hbits

    mid = hi_x * lo_y
    hi += mid >> hbits

    hi, x*y
end

function mul_hi_lo_3(x::T, y::T) where T <: Unsigned
    xy = widemul(x, y)
    (x >> 8sizeof(T)) % T, x % T
end

function normalize_to_uint_dual(::Type{T}, v::AbstractVector{<:Integer}) where {T <: Unsigned}
    sm = sum(T, v)
    q, r = divrem(typemax(T), sm)
    if r == sm-true
        q += true
        r = zero(r)
    else
        r += true
    end
    # q, r = divrem(typemax(T)+1, sm)

    error = zero(T)
    res = Vector{T}(undef, length(v))
    for (i,x) in enumerate(v)
        tx = T(x)
        # @show q, r, tx, r*tx
        hi, lo = mul_hi_lo_3(r, tx)
        error += r*tx # TODO: check overflow
        carry, error = divrem(error, sm)
        res[i] = q*tx + carry
    end
    res
end


function frac_div(x::T, y::T) where T <: Unsigned
    # @assert x < y
    # @assert y != 0
    div(widen(x) << 8sizeof(T), y) % T
end
function frac_div_2(x::T, y::T) where T <: Unsigned
    # @assert x < y
    # @assert y != 0
    Base.udiv_int(promote(widen(x) << 8sizeof(T), y)...) % T
end

function normalize_to_uint_frac_div(::Type{T}, v::AbstractVector{<:Integer}) where {T <: Unsigned}
    sm = sum(T, v)

    sm2 = zero(T)

    res = Vector{T}(undef, length(v))
    for (i,x) in enumerate(v)
        val = frac_div(T(x), sm)
        sm2 += val
        res[i] = val
    end

    for i in sm2:typemax(sm2)
        res[typemax(sm2)-i+1] += true
    end

    res
end

####

# 2-4 passes (skip first two if nomralize = false)
# Initial sum, check for overflow, negatives, allzero & exactness
# If not exact, compute the mapped sum and the error there
# Dual pass for construction

## Strict (fragile) & efficient float handling

# struct RescaleView{T, P} <: AbstractVector{T}
#     parent::P
# end
# Base.getindex(rv::RescaleView{T}, i) where T = T(ldexp(rv.parent[i], 8*sizeof(T)))
# Base.size(rv::RescaleView) = size(rv.parent)
# OffsetTable{T, I}(weights::AbstractVector{<:Real}) where {T, I} = _offset_table(I, RescaleView{T, typeof(weights)}(weights))

## Better float handling (wip)

# function offset_table_float(::Type{T}, ::Type{I}, weights::AbstractVector{<:Real}) where {T <: Unsigned, I <: Integer}
#     Base.require_one_based_indexing(weights)

#     bitshift = Base.top_set_bit(length(weights) - 1)
#     len = 1 << bitshift#next_or_eq_power_of_two(length(weights))
#     points_per_cell = one(T) << (8*sizeof(T) - bitshift)#typemax(T)+1 / len

#     thirsty_i = surplus_i = current_i = firstindex(weights)
#     current_desired = weights[current_i]

#     exp = exponent(sum(weights))
#     real_to_uint(x) = floor(T, ldexp(x, 8sizeof(T)-exp))
#     # sum(real_to_uint, weights) could be as low as 0 if weights = fill(1, 2^33) and T=UInt32, though in that case we can assign every bin either 0 or 1 point at will
#     # it could be equal to typemax(T)+1 if weights = [1] or there is otherwise no rounding down.
#     sm = sum(widen ∘ real_to_uint, weights)



#     probability_offset = Memory{Tuple{T, I}}(undef, len)


#     while true
#         # @show current_i, current_desired, points_per_cell
#         if current_desired <= points_per_cell # Surplus
#             excess = points_per_cell - current_desired                     # Assign this many extra points
#             thirsty_i = findnext(>(points_per_cell), weights, thirsty_i+1) # Which is the next available thirsty cell
#             thirsty_i === nothing && break                                 # If there is no thirsty cell, handle below
#             probability_offset[current_i] = (excess, thirsty_i-current_i)  # To the targeted cell
#             current_i = thirsty_i                                          # Now we have to make sure that thristy cell gets exactly what it wants and no more
#             current_desired = weights[current_i] - excess                  # It wants what it wants nad hasn't already been transferred
#         else                                  # Thirsty (strictly)
#             surplus_i = findnext(<=(points_per_cell), weights, surplus_i+1) # Find the next surplus cell
#             surplus_i === nothing && throw(ArgumentError("sum(weights) is too high")) # Lacking points, so unnable to reach the desired weight
#             excess = weights[surplus_i] - points_per_cell                   # Assign this many extra points
#             probability_offset[surplus_i] = (excess, current_i-surplus_i)   # From the cell with surplus to this cell
#             current_desired -= excess                                       # We now don't want as many points (and may even no longer be thristy)
#         end
#     end
#     if current_desired < points_per_cell # Surplus points, so exceed the desired weight
#         throw(ArgumentError("sum(weights) is too low"))
#     end
#     probability_offset[current_i] = (0, 0)
#     # Just right. There are no thristy cells left and no current surplus. All that's left
#     # are future surplus cells, all of which should be a surplus of exactly 0.
#     while true
#         surplus_i = findnext(<=(points_per_cell), weights, surplus_i+1) # Find the next surplus cell
#         surplus_i === nothing && break
#         weights[surplus_i] == points_per_cell || throw(ArgumentError("sum(weights) is too low"))
#         probability_offset[surplus_i] = (0, 0)
#     end

#     OffsetTable(probability_offset)
# end



end

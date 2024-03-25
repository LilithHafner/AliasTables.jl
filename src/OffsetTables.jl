module OffsetTables

export OffsetTable

struct OffsetTable{T, I}
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

OffsetTable(weights::AbstractVector{<:Real}) = OffsetTable{UInt, Int}(weights)
OffsetTable(weights::AbstractVector{T}) where T <: Unsigned = OffsetTable{T, Int}(weights)

OffsetTable{T}(weights::AbstractVector{<:Real}) where T = OffsetTable{T, Int}(weights)

OffsetTable{T, I}(weights::AbstractVector{<:Real}) where {T, I} = _offset_table(I, normalize_to_uint(T, weights))
OffsetTable{T, I}(weights::AbstractVector{T}) where {T<:Real, I} = _offset_table(I, weights)

struct WithZeros{T} <: AbstractVector{T}
    parent::AbstractVector{T}
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

function _offset_table(::Type{I}, weights::AbstractVector{<:Unsigned}) where I
    T = eltype(weights)
    Base.require_one_based_indexing(weights)

    isempty(weights) && throw(ArgumentError("weights must be non-empty"))

    onz = get_only_nonzero(weights)
    onz == -2 || return _constant_offset_table(I, T, onz)

    bitshift = Base.top_set_bit(findlast(!iszero, weights) - 1)
    len = 1 << bitshift#next_or_eq_power_of_two(length(weights))
    points_per_cell = one(T) << (8*sizeof(T) - bitshift)#typemax(T)+1 / len

    probability_offset = Memory{Tuple{T, I}}(undef, len)

    weights_extended = WithZeros(weights, len)
    thirsty_i = surplus_i = current_i = firstindex(weights_extended)
    current_desired = weights_extended[current_i]

    while true
        # @show current_i, current_desired, points_per_cell
        if current_desired <= points_per_cell # Surplus
            excess = points_per_cell - current_desired                     # Assign this many extra points
            thirsty_i = findnext(>(points_per_cell), weights_extended, thirsty_i+1) # Which is the next available thirsty cell
            thirsty_i === nothing && break                                 # If there is no thirsty cell, handle below
            probability_offset[current_i] = (excess, thirsty_i-current_i)  # To the targeted cell
            current_i = thirsty_i                                          # Now we have to make sure that thristy cell gets exactly what it wants and no more
            current_desired = weights_extended[current_i] - excess                  # It wants what it wants nad hasn't already been transferred
        else                                  # Thirsty (strictly)
            surplus_i = findnext(<=(points_per_cell), weights_extended, surplus_i+1) # Find the next surplus cell
            surplus_i === nothing && throw(ArgumentError("sum(weights) is too high")) # Lacking points, so unnable to reach the desired weight
            excess = points_per_cell - weights_extended[surplus_i]                   # Assign this many extra points
            probability_offset[surplus_i] = (excess, current_i-surplus_i)   # From the cell with surplus to this cell
            current_desired -= excess                                       # We now don't want as many points (and may even no longer be thristy)
        end
    end
    if current_desired < points_per_cell # Surplus points, so exceed the desired weight
        throw(ArgumentError("sum(weights) is too low"))
    end
    probability_offset[current_i] = (0, 0)
    # Just right. There are no thristy cells left and no current surplus. All that's left
    # are future surplus cells, all of which should be a surplus of exactly 0.
    while true
        surplus_i = findnext(<=(points_per_cell), weights_extended, surplus_i+1) # Find the next surplus cell
        surplus_i === nothing && break
        weights_extended[surplus_i] == points_per_cell || throw(ArgumentError("sum(weights) is too low"))
        probability_offset[surplus_i] = (0, 0)
    end

    _OffsetTable(probability_offset)
end

function sample(rng, ot::OffsetTable{T, I}) where {T, I}
    sample(rand(rng, T), ot)
end

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

### Random API
using Random
Random.rand(rng::Random.AbstractRNG, ot::Random.SamplerTrivial{<:OffsetTable}) = sample(rng, ot.self)

### Reconstruct probabilities
function probabilities(ot::OffsetTable{T}) where T
    bitshift = Base.top_set_bit(length(ot.probability_offset) - 1)
    points_per_cell = one(T) << (8*sizeof(T) - bitshift)#typemax(T)+1 / len
    probs = zeros(T, length(ot.probability_offset))
    for (i, (prob, offset)) in enumerate(ot.probability_offset)
        probs[i + offset] += prob
        probs[i] += points_per_cell - prob
    end
    resize!(probs, findlast(!iszero, probs))
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

function normalize_to_uint(::Type{T}, v::AbstractVector{<:Real}) where {T <: Unsigned}
    isempty(v) && throw(ArgumentError("weights must be non-empty"))

    onz = get_only_nonzero(v)
    if onz != -2
        res = zeros(T, length(v))
        res[onz] = 1
        return res
    end

    length(v) <= 1 && return ones(T, length(v))
    sm = sum(v)
    res = [x <= typemax(T) ? floor(T, x) : typemax(T) for x in (ldexp(x/sm, 8sizeof(T)) for x in v)] # TODO: make this lazy & non-allocating

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

module OffsetTables

export OffsetTable

struct OffsetTable{T, I}
    probability_offset::Memory{Tuple{T, I}}
end
struct RescaleView{T, P} <: AbstractVector{T}
    parent::P
end
Base.getindex(rv::RescaleView{T}, i) where T = T(ldexp(rv.parent[i], 8*sizeof(T)))
Base.size(rv::RescaleView) = size(rv.parent)

OffsetTable(weights::AbstractVector{<:Real}) = OffsetTable{UInt, Int}(weights)
OffsetTable(weights::AbstractVector{T}) where T <: Unsigned = OffsetTable{T, Int}(weights)

OffsetTable{T}(weights::AbstractVector{<:Real}) where T = OffsetTable{T, Int}(weights)

OffsetTable{T, I}(weights::AbstractVector{<:Real}) where {T, I} = _offset_table(I, RescaleView{T, typeof(weights)}(weights))
OffsetTable{T, I}(weights::AbstractVector{T}) where {T<:Real, I} = _offset_table(I, weights)

function _offset_table(::Type{I}, weights::AbstractVector{<:Unsigned}) where I
    T = eltype(weights)
    Base.require_one_based_indexing(weights)

    bitshift = Base.top_set_bit(length(weights) - 1)
    len = 1 << bitshift#next_or_eq_power_of_two(length(weights))
    points_per_cell = one(T) << (8*sizeof(T) - bitshift)#typemax(T)+1 / len

    thirsty_i = surplus_i = current_i = firstindex(weights)
    current_desired = weights[current_i]

    probability_offset = Memory{Tuple{T, I}}(undef, len)

    while true
        # @show current_i, current_desired, points_per_cell
        if current_desired <= points_per_cell # Surplus
            excess = points_per_cell - current_desired                     # Assign this many extra points
            thirsty_i = findnext(>(points_per_cell), weights, thirsty_i+1) # Which is the next available thirsty cell
            thirsty_i === nothing && break                                 # If there is no thirsty cell, handle below
            probability_offset[current_i] = (excess, thirsty_i-current_i)  # To the targeted cell
            current_i = thirsty_i                                          # Now we have to make sure that thristy cell gets exactly what it wants and no more
            current_desired = weights[current_i] - excess                  # It wants what it wants nad hasn't already been transferred
        else                                  # Thirsty (strictly)
            surplus_i = findnext(<=(points_per_cell), weights, surplus_i+1) # Find the next surplus cell
            surplus_i === nothing && throw(ArgumentError("sum(weights) is too high")) # Lacking points, so unnable to reach the desired weight
            excess = weights[surplus_i] - points_per_cell                   # Assign this many extra points
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
        surplus_i = findnext(<=(points_per_cell), weights, surplus_i+1) # Find the next surplus cell
        surplus_i === nothing && break
        weights[surplus_i] == points_per_cell || throw(ArgumentError("sum(weights) is too low"))
        probability_offset[surplus_i] = (0, 0)
    end

    OffsetTable(probability_offset)
end

function sample(rng, ot::OffsetTable{T, I}) where {T, I}
    sample(rand(rng, T), ot)
end

function sample(x::T, ot::OffsetTable{T, I}) where {T, I}
    count_ones(x) == 1 || return zero(I) # This should never happen, but it makes the @inbounds safe
    shift = 8sizeof(T) - Base.top_set_bit(length(ot.probability_offset)) + 1
    cell = (x >> shift) + 1
    val = x & ((one(T) << shift) - one(T))
    @inbounds prob, offset = ot.probability_offset[cell%Int]
    # (val < prob ? (offset+cell) : I(cell))::I
    # I((val < prob) * offset + cell)::I
    (((val < prob) * offset + cell)%I)::I
end

### Reconstruct probabilities
function probabilities(ot::OffsetTable{T}) where T
    bitshift = Base.top_set_bit(length(ot.probability_offset) - 1)
    points_per_cell = one(T) << (8*sizeof(T) - bitshift)#typemax(T)+1 / len
    probs = zeros(T, length(ot.probability_offset))
    for (i, (prob, offset)) in enumerate(ot.probability_offset)
        probs[i + offset] += prob
        probs[i] += points_per_cell - prob
    end
    probs
end

probabilities(::typeof(float), ot::OffsetTable{T}) where T =
    probabilities(ot) ./ (float(typemax(T))+1)


### Show
function Base.show(io::IO, ot::OffsetTable{T, I}) where {T, I}
    print(io, typeof(ot), "(")
    print(IOContext(io, :typeinfo=>Memory{Tuple{T, I}}), ot.probability_offset)
    print(io, ")")
end

end

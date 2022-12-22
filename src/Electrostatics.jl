module Electrostatics

using ..LastHomework: DiscreteLaplacian

export ReshapeVector,
    getbcindices,
    getsquareindices,
    getchargeindices,
    checkbc,
    checksquare,
    checkcharges,
    setbc!,
    setsquare!,
    setcharges!

struct ReshapeVector{T} <: AbstractVector{T}
    data::Vector{T}
    size::NTuple{2,Int64}
    function ReshapeVector{T}(data, size) where {T}
        if length(data) != prod(size)
            throw(
                DimensionMismatch(
                    "dimensions $size must be consistent with array size $(length(data))!"
                ),
            )
        end
        return new(data, size)
    end
end
ReshapeVector(data::AbstractVector{T}, dims) where {T} = ReshapeVector{T}(data, dims)
ReshapeVector(data::AbstractVector{T}, dims...) where {T} = ReshapeVector{T}(data, dims)

abstract type FixedRegion end
struct Boundary <: FixedRegion end
struct InternalSquare <: FixedRegion end
struct PointCharges <: FixedRegion end

abstract type PartiallyFixedVector{T} <: AbstractVector{T} end
struct SolutionVector{T} <: PartiallyFixedVector{T}
    parent::Vector{T}
end
struct ResidualVector{T} <: PartiallyFixedVector{T}
    parent::Vector{T}
end

function getindices(ϕ::AbstractMatrix, ::Boundary)
    cartesian_indices = CartesianIndices(ϕ)
    # Note the geometry of the region and the matrix rows/columns ordering are the same!
    # See https://discourse.julialang.org/t/how-to-get-the-cartesian-indices-of-a-row-column-in-a-matrix/91940/2
    return vcat(
        cartesian_indices[begin, :],  # Bottom
        cartesian_indices[end, :],  # Top
        cartesian_indices[:, begin],  # Left
        cartesian_indices[:, end],  # Right
    )
end
function getindices(ϕ::AbstractMatrix, ::InternalSquare)
    M, N = size(ϕ)
    xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ = map(Int64, (M / 2, M * 3//4, N * 5//8, N * 7//8))
    return map(Iterators.product(xₘᵢₙ:xₘₐₓ, yₘᵢₙ:yₘₐₓ)) do (i, j)
        CartesianIndex(j, i)  # Note y -> row, x -> column
    end
end
function getindices(ρ::AbstractMatrix, ::PointCharges)
    M, N = size(ρ)
    x₁, x₂, y = map(Int64, (M / 4, M * 3//4, N / 8))
    return map(CartesianIndex, ((y, x₁), (y, x₂)))  # Note y -> row, x -> column
end
function getindices(ρ::AbstractMatrix, ::PointCharges)
    M, N = size(ρ)
    x₁, x₂, y = map(Int64, (M / 4, M * 3//4, N / 8))
    return map(CartesianIndex, ((y, x₁), (y, x₂)))  # Note y -> row, x -> column
end
# See See https://discourse.julialang.org/t/how-to-convert-cartesianindex-n-values-to-int64/15074/4
# and http://docs.julialang.org/en/v1/base/arrays/#Base.LinearIndices
function getindices(vec::ReshapeVector, region::FixedRegion)
    mat = reshape(vec)
    linear_indices = LinearIndices(mat)
    cartesian_indices = collect(getindices(vec, region))  # `getindex` only accepts vector indices
    return linear_indices[cartesian_indices]
end

function _checkequal(data::AbstractVecOrMat, value, region::FixedRegion)
    indices = getindices(data, region)
    for index in indices
        @assert data[index] == value
    end
    return nothing
end

function _setconst!(data::AbstractVecOrMat, value, region::FixedRegion)
    indices = getindices(data, region)
    for index in indices
        data[index] = value
    end
    return data
end

Base.parent(vec::ReshapeVector) = vec.data

Base.size(vec::ReshapeVector) = size(parent(vec))

Base.IndexStyle(::Type{ReshapeVector{T}}) where {T} = IndexLinear()

Base.getindex(vec::ReshapeVector, i) = getindex(parent(vec), i)

Base.setindex!(vec::ReshapeVector, v, i) = setindex!(parent(vec), v, i)

# Base.similar(::ReshapeVector, ::Type{T}, dims::Dims) where {T} =
#     ReshapeVector(Vector{T}(undef, dims), dims)

Base.reshape(vec::ReshapeVector) = reshape(vec.data, vec.size)

Base.parent(vec::PartiallyFixedVector) = vec.data

Base.size(vec::PartiallyFixedVector) = size(parent(vec))

Base.IndexStyle(::Type{PartiallyFixedVector{T}}) where {T} = IndexLinear()

Base.getindex(vec::PartiallyFixedVector, i) = getindex(parent(vec), i)

Base.setindex!(vec::PartiallyFixedVector, v, i) = setindex!(parent(vec), v, i)

Base.similar(::PartiallyFixedVector, ::Type{T}, dims::Dims) where {T} =
    PartiallyFixedVector(Vector{T}(undef, dims))

function Base.:*(A::DiscreteLaplacian, 𝐯::PartiallyFixedVector)
    𝐯′ = A * 𝐯
    _setconst!(f, 𝐯, 1)
    return 𝐯′
end

end

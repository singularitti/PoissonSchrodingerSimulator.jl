module Electrostatics

using ..LastHomework: DiscreteLaplacian

export Boundary,
    InternalSquare,
    PointCharges,
    SolutionVector,
    ResidualVector,
    getindices,
    checkequal,
    set!

struct Region
    dims::NTuple{2,Int64}
end
Region(m, n) = Region((m, n))

const BOX = Region(128, 128)

abstract type FixedValueRegion{T} end
struct Boundary{T} <: FixedValueRegion{T}
    value::T
end
struct InternalSquare{T} <: FixedValueRegion{T}
    value::T
end
struct PointCharges{T} <: FixedValueRegion{T}
    value::T
end

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
# See See https://discourse.julialang.org/t/how-to-convert-cartesianindex-n-values-to-int64/15074/4
# and http://docs.julialang.org/en/v1/base/arrays/#Base.LinearIndices
function getindices(vec::PartiallyFixedVector, region::FixedValueRegion)
    mat = reshape(vec, BOX.dims)
    linear_indices = LinearIndices(mat)
    cartesian_indices = collect(getindices(vec, region))  # `getindex` only accepts vector indices
    return linear_indices[cartesian_indices]
end

function checkequal(data::AbstractVecOrMat, region::FixedValueRegion)
    indices = getindices(data, region)
    for index in indices
        @assert data[index] == region.value
    end
    return nothing
end

function set!(data::AbstractVecOrMat, region::FixedValueRegion)
    indices = getindices(data, region)
    for index in indices
        data[index] = region.value
    end
    return data
end

Base.parent(vec::PartiallyFixedVector) = vec.parent

Base.size(vec::PartiallyFixedVector) = size(parent(vec))

Base.IndexStyle(::Type{<:PartiallyFixedVector}) = IndexLinear()

Base.getindex(vec::PartiallyFixedVector, i) = getindex(parent(vec), i)

Base.setindex!(vec::PartiallyFixedVector, v, i) = setindex!(parent(vec), v, i)

Base.similar(::PartiallyFixedVector, ::Type{T}, dims::Dims) where {T} =
    PartiallyFixedVector(Vector{T}(undef, dims))

function Base.:*(A::DiscreteLaplacian, 𝐯::SolutionVector)
    𝐯′ = A * 𝐯
    set!(𝐯, Boundary(0))
    set!(𝐯, InternalSquare(5))
    return 𝐯′
end
function Base.:*(A::DiscreteLaplacian, 𝐯::ResidualVector)
    𝐯′ = A * 𝐯
    set!(𝐯, PointCharges(-20))
    return 𝐯′
end

end

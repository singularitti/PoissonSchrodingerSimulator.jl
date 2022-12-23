module Electrostatics

using ..LastHomework: DiscreteLaplacian

export Boundary, InternalSquare, PointCharges, getindices, checkequal, set

abstract type FixedValueRegion{T} end
struct Boundary{T} <: FixedValueRegion{T}
    boxsize::NTuple{2,Int}
    value::T
end
struct InternalSquare{T} <: FixedValueRegion{T}
    boxsize::NTuple{2,Int}
    value::T
end
struct PointCharge{T} <: FixedValueRegion{T}
    boxsize::NTuple{2,Int}
    value::T
end
struct PointCharges{T} <: FixedValueRegion{T}
    boxsize::NTuple{2,Int}
    value::T
end

N::Int64 = 33

BOUNDARY = Boundary((N, N), 0)
SQUARE = InternalSquare((N, N), 5)

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
    M, N = size(ϕ) .- 1
    xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ = map(Int64, (M / 2, M * 3//4, N * 5//8, N * 7//8))
    return map(Iterators.product(xₘᵢₙ:xₘₐₓ, yₘᵢₙ:yₘₐₓ)) do (i, j)
        CartesianIndex(j + 1, i + 1)   # Note y -> row, x -> column
    end
end
function getindices(ρ::AbstractMatrix, ::PointCharges)
    M, N = size(ρ) .- 1
    x₁, x₂, y = map(Int64, (M / 4, M * 3//4, N / 8)) .+ 1
    return map(CartesianIndex, ((y, x₁), (y, x₂)))  # Note y -> row, x -> column
end
# See See https://discourse.julialang.org/t/how-to-convert-cartesianindex-n-values-to-int64/15074/4
# and http://docs.julialang.org/en/v1/base/arrays/#Base.LinearIndices
function getindices(vec::AbstractVector, region::FixedValueRegion)
    mat = reshape(vec, region.boxsize)
    linear_indices = LinearIndices(mat)
    cartesian_indices = collect(getindices(mat, region))  # `getindex` only accepts vector indices
    return linear_indices[cartesian_indices]
end

function checkequal(data, region::FixedValueRegion)
    indices = getindices(data, region)
    for index in indices
        @assert data[index] == region.value
    end
    return nothing
end

function set(data, region::FixedValueRegion)
    indices = getindices(data, region)
    for index in indices
        data[index] = region.value
    end
    return vec(data)
end

function Base.:*(A::DiscreteLaplacian, 𝐩ₙ::AbstractVector)
    𝐩ₙ = set(𝐩ₙ, BOUNDARY)
    𝐩ₙ = set(𝐩ₙ, SQUARE)
    return parent(A) * parent(𝐩ₙ)
end

end

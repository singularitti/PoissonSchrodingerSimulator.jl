module Electrostatics

using ..LastHomework: DiscreteLaplacian

export Boundary, InternalSquare, PointCharges, getindices, checkequal, set!

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

function checkequal(data::PartiallyFixedMatrix, region::FixedValueRegion)
    indices = getindices(data, region)
    for index in indices
        @assert data[index] == region.value
    end
    return nothing
end

function set!(data::PartiallyFixedMatrix, region::FixedValueRegion)
    indices = getindices(data, region)
    for index in indices
        data[index] = region.value
    end
    return data
end

end

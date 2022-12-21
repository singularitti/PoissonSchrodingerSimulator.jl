module Electrostatics

using ..LastHomework: DiscreteLaplacianPBCs

export ReshapeVector, checkbc, checksquare, checkcharges, setbc!, setsquare!, setcharges!

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
ReshapeVector(data::AbstractVector{T}, size) where {T} = ReshapeVector{T}(data, size)

function getbcindices(ϕ::AbstractMatrix)
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
getbcindices(𝛟::ReshapeVector) = _getindices(getbcindices, 𝛟)

checkbc(ϕ, ϕ₀) = _checkequal(getbcindices, ϕ, ϕ₀)

setbc!(ϕ, ϕ₀) = _setconst!(getbcindices, ϕ, ϕ₀)

function getsquareindices(ϕ::AbstractMatrix)
    M, N = size(ϕ)
    xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ = map(Int64, (M / 2, M * 3//4, N * 5//8, N * 7//8))
    return map(Iterators.product(xₘᵢₙ:xₘₐₓ, yₘᵢₙ:yₘₐₓ)) do (i, j)
        CartesianIndex(j, i)  # Note y -> row, x -> column
    end
end
getsquareindices(𝛟::ReshapeVector) = _getindices(getsquareindices, 𝛟)

checksquare(ϕ, ϕ₀) = _checkequal(getsquareindices, ϕ, ϕ₀)

setsquare!(ϕ, ϕ₀) = _setconst!(getsquareindices, ϕ, ϕ₀)

function getchargeindices(ρ::AbstractMatrix)
    M, N = size(ρ)
    x₁, x₂, y = map(Int64, (M / 4, M * 3//4, N / 8))
    return map(CartesianIndex, ((y, x₁), (y, x₂)))  # Note y -> row, x -> column
end
getchargeindices(𝛒::ReshapeVector) = _getindices(getchargeindices, 𝛒)

checkcharges(ρ, ρ₀) = _checkequal(getchargeindices, ρ, ρ₀)

setcharges!(ρ, ρ₀) = _setconst!(getchargeindices, ρ, ρ₀)

# See See https://discourse.julialang.org/t/how-to-convert-cartesianindex-n-values-to-int64/15074/4
# and http://docs.julialang.org/en/v1/base/arrays/#Base.LinearIndices
function _getindices(f::Function, vec::ReshapeVector)
    vec = reshape(vec)
    linear_indices = LinearIndices(vec)
    cartesian_indices = f(vec)
    return linear_indices[cartesian_indices]
end

function _checkequal(f::Function, data::AbstractVecOrMat, value)
    indices = f(data)
    for index in indices
        @assert data[index] == value
    end
    return nothing
end

function _setconst!(f, data::AbstractVecOrMat, value)
    indices = f(data)
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

Base.similar(::ReshapeVector, ::Type{T}, dims::Dims) where {T} =
    ReshapeVector(Vector{T}(undef, dims), dims)

Base.reshape(vec::ReshapeVector) = reshape(vec.data, vec.size)

end

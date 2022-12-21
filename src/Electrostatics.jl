module Electrostatics

using ..LastHomework: DiscreteLaplacianPBCs

export checkbc, checksquare, checkcharges, setbc!, setsquare!, setcharges!

struct ReshapeVector{T} <: AbstractVector{T}
    data::Vector{T}
    size::NTuple{2,Int64}
end

function checkbc(ϕ::AbstractMatrix, ϕ₀)
    @assert ϕ[begin, :] == ϕ₀  # Top
    @assert ϕ[end, :] == ϕ₀  # Bottom
    @assert ϕ[:, begin] == ϕ₀  # Left
    @assert ϕ[:, end] == ϕ₀  # Right
    return nothing
end

function setbc!(ϕ::AbstractMatrix, ϕ₀)
    ϕ[begin, :] = ϕ₀  # Top
    ϕ[end, :] = ϕ₀  # Bottom
    ϕ[:, begin] = ϕ₀  # Left
    ϕ[:, end] = ϕ₀  # Right
    return ϕ
end

function getsquareindices(ϕ::AbstractMatrix)
    M, N = size(ϕ)
    xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ = map(Int64, (M / 2, M * 3//4, N * 5//8, N * 7//8))
    return map(xₘᵢₙ:xₘₐₓ, yₘᵢₙ:yₘₐₓ) do i, j
        CartesianIndex(i, j)
    end
end
getsquareindices(𝛟::ReshapeVector) = _getindices(getsquareindices, 𝛟)

checksquare(ϕ, ϕ₀) = _checkequal(getsquareindices, ϕ, ϕ₀)

setsquare!(ϕ, ϕ₀) = _setconst!(setsquare!, ϕ, ϕ₀)

function getchargeindices(ρ::AbstractMatrix)
    M, N = size(ρ)
    x₁, x₂, y = map(Int64, (M / 4, M * 3//4, N / 8))
    return map(CartesianIndex, ((x₁, y), (x₂, y)))
end

checkcharges(ρ, ρ₀) = _checkequal(getchargeindices, ρ, ρ₀)

setcharges!(ρ, ρ₀) = _setconst!(setcharges!, ρ, ρ₀)

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

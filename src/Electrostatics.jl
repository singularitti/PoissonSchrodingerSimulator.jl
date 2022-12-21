module Electrostatics

using ..LastHomework: DiscreteLaplacianPBCs

export checkbc, checksquare, checkcharges, setbc!, setsquare!, setcharges!

function checkbc(ϕ::AbstractMatrix, ϕ₀)
    @assert ϕ[begin, :] == ϕ₀  # Top
    @assert ϕ[end, :] == ϕ₀  # Bottom
    @assert ϕ[:, begin] == ϕ₀  # Left
    @assert ϕ[:, end] == ϕ₀  # Right
    return nothing
end
checkbc(𝛟::AbstractVector, M, N, ϕ₀) = _checkvec(checkbc, 𝛟, M, N, ϕ₀)

function setbc!(ϕ::AbstractMatrix, ϕ₀)
    ϕ[begin, :] = ϕ₀  # Top
    ϕ[end, :] = ϕ₀  # Bottom
    ϕ[:, begin] = ϕ₀  # Left
    ϕ[:, end] = ϕ₀  # Right
    return ϕ
end
setbc!(𝛟::AbstractVector, M, N, ϕ₀) = _setvec!(setbc!, 𝛟, M, N, ϕ₀)

function getsquareindices(ϕ::AbstractMatrix)
    M, N = size(ϕ)
    xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ = map(Int64, (M / 2, M * 3//4, N * 5//8, N * 7//8))
    return map(xₘᵢₙ:xₘₐₓ, yₘᵢₙ:yₘₐₓ) do i, j
        CartesianIndex(i, j)
    end
end
getsquareindices(𝛟::AbstractVector, M, N) =
    _getindices(getsquareindices, 𝛟::AbstractVector, M, N)

function checksquare(ϕ::AbstractMatrix, ϕ₀)
    indices = getsquareindices(ϕ)
    for index in indices
        @assert ϕ[index] == ϕ₀
    end
    return nothing
end
checksquare(𝛟::AbstractVector, M, N, ϕ₀) = _checkvec(checksquare, 𝛟, M, N, ϕ₀)

function setsquare!(ϕ::AbstractMatrix, ϕ₀)
    indices = getsquareindices(ϕ)
    for index in indices
        ϕ[index] = ϕ₀
    end
    return ϕ
end
setsquare!(𝐯::AbstractVector, M, N, ϕ₀) = _setvec!(setsquare!, 𝐯, M, N, ϕ₀)

function getchargeindices(ρ::AbstractMatrix)
    M, N = size(ρ)
    x₁, x₂, y = map(Int64, (M / 4, M * 3//4, N / 8))
    return map(CartesianIndex, ((x₁, y), (x₂, y)))
end

function checkcharges(ρ::AbstractMatrix, ρ₀)
    indices = getchargeindices(ρ)
    for index in indices
        @assert ρ[index] == ρ₀
    end
    return nothing
end
checkcharges(𝛒::AbstractVector, M, N, ρ₀) = _checkvec(checkcharges, 𝛒, M, N, ρ₀)

function setcharges!(ρ::AbstractMatrix, ρ₀)
    M, N = size(ρ)
    x₁, x₂, y = map(Int64, (M / 4, M * 3//4, N / 8))
    ρ[x₁, y] = ρ₀
    ρ[x₂, y] = ρ₀
    return ρ
end
setcharges!(𝛒::AbstractVector, M, N, ρ₀) = _setvec!(setcharges!, 𝛒, M, N, ρ₀)

_checkvec(f::Function, 𝐯::AbstractVector, M, N, value) = f(reshape(𝐯, M, N), value)

# See See https://discourse.julialang.org/t/how-to-convert-cartesianindex-n-values-to-int64/15074/4
# and http://docs.julialang.org/en/v1/base/arrays/#Base.LinearIndices
function _getindices(f::Function, 𝐯::AbstractVector, M, N)
    v = reshape(𝐯, M, N)
    linear_indices = LinearIndices(v)
    cartesian_indices = f(v)
    return linear_indices[cartesian_indices]
end

function _setvec!(f::Function, 𝐯::AbstractVector, M, N, value)
    indices = f(𝐯, M, N)
    for index in indices
        𝐯[index] = value
    end
    return 𝐯
end

end

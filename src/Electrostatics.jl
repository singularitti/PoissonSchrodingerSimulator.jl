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
checkbc(𝛟::AbstractVector, M, N, ϕ₀) = _checkvector(checkbc, 𝛟, M, N, ϕ₀)

function setbc!(ϕ::AbstractMatrix, ϕ₀)
    ϕ[begin, :] = ϕ₀  # Top
    ϕ[end, :] = ϕ₀  # Bottom
    ϕ[:, begin] = ϕ₀  # Left
    ϕ[:, end] = ϕ₀  # Right
    return ϕ
end
function setbc!(𝛟::AbstractVector, M, N, ϕ₀)
    ϕ = reshape(𝛟, M, N)
    ϕ = setbc!(ϕ, ϕ₀)
    return reshape(ϕ, length(ϕ))
end

function checksquare(ϕ::AbstractMatrix, ϕ₀)
    M, N = size(ϕ)
    xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ = map(Int64, (M / 2, M * 3//4, N * 5//8, N * 7//8))
    for i in xₘᵢₙ:xₘₐₓ
        for j in yₘᵢₙ:yₘₐₓ
            @assert ϕ[i, j] == ϕ₀
        end
    end
    return nothing
end
checksquare(𝛟::AbstractVector, M, N, ϕ₀) = _checkvector(checksquare, 𝛟, M, N, ϕ₀)

function setsquare!(ϕ::AbstractMatrix, ϕ₀)
    M, N = size(ϕ)
    xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ = map(Int64, (M / 2, M * 3//4, N * 5//8, N * 7//8))
    for i in xₘᵢₙ:xₘₐₓ
        for j in yₘᵢₙ:yₘₐₓ
            ϕ[i, j] = ϕ₀
        end
    end
    return ϕ
end
function setsquare!(𝛟::AbstractVector, M, N, ϕ₀)
    ϕ = reshape(𝛟, M, N)
    ϕ = setsquare!(ϕ, ϕ₀)
    return reshape(ϕ, length(ϕ))
end

function checkcharges(ρ::AbstractMatrix, ρ₀)
    M, N = size(ρ)
    x₁, x₂, y = map(Int64, (M / 4, M * 3//4, N / 8))
    @assert ρ[x₁, y] == ρ₀
    @assert ρ[x₂, y] == ρ₀
    return nothing
end
checkcharges(𝛒::AbstractVector, M, N, ρ₀) = _checkvector(checkcharges, 𝛒, M, N, ρ₀)

function setcharges!(ρ::AbstractMatrix, ρ₀)
    M, N = size(ρ)
    x₁, x₂, y = map(Int64, (M / 4, M * 3//4, N / 8))
    ρ[x₁, y] = ρ₀
    ρ[x₂, y] = ρ₀
    return ρ
end
function setcharges!(𝛒::AbstractVector, M, N, ρ₀)
    ρ = reshape(𝛒, M, N)
    ρ = setsquare!(ρ, ρ₀)
    return reshape(ρ, length(ρ))
end

_checkvector(f::Function, 𝐯::AbstractVector, M, N, value) = f(reshape(𝐯, M, N), value)

end

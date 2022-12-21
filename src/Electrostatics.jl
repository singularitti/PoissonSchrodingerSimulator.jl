module Electrostatics

using ..LastHomework: DiscreteLaplacianPBCs

export setbc!, setsquare!, setcharges!

function setbc!(ϕ::AbstractMatrix, v=zero(eltype(ϕ)))
    ϕ[begin, :] = v  # Top
    ϕ[end, :] = v  # Bottom
    ϕ[:, begin] = v  # Left
    ϕ[:, end] = v  # Right
    return ϕ
end
function setbc!(𝛟::AbstractVector, M, N, v=zero(eltype(𝛟)))
    ϕ = reshape(𝛟, M, N)
    ϕ = setbc!(ϕ, v)
    return reshape(ϕ, length(ϕ))
end

function setsquare!(ϕ::AbstractMatrix, v=oneunit(eltype(ϕ)))
    M, N = size(ϕ)
    xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ = map(Int64, (M / 2, M * 3//4, N * 5//8, N * 7//8))
    for i in xₘᵢₙ:xₘₐₓ
        for j in yₘᵢₙ:yₘₐₓ
            ϕ[i, j] = v
        end
    end
    return ϕ
end
function setsquare!(𝛟::AbstractVector, M, N, v=oneunit(eltype(𝛟)))
    ϕ = reshape(𝛟, M, N)
    ϕ = setsquare!(ϕ, v)
    return reshape(ϕ, length(ϕ))
end

function setcharges!(ρ::AbstractMatrix, v=oneunit(eltype(ρ)))
    M, N = size(ρ)
    x₁, x₂, y = map(Int64, (M / 4, M * 3//4, N / 8))
    ρ[x₁, y] = v
    ρ[x₂, y] = v
    return ρ
end
function setcharges!(𝛒::AbstractVector, M, N, v=oneunit(eltype(𝛒)))
    ρ = reshape(𝛒, M, N)
    ρ = setsquare!(ρ, v)
    return reshape(ρ, length(ρ))
end

end

module Electrostatics

using ..LastHomework: DiscreteLaplacianPBCs

export checkbc, checksquare, checkcharges, setbc!, setsquare!, setcharges!

function checkbc(ϕ::AbstractMatrix, v)
    @assert ϕ[begin, :] == v  # Top
    @assert ϕ[end, :] == v  # Bottom
    @assert ϕ[:, begin] == v  # Left
    @assert ϕ[:, end] == v  # Right
    return nothing
end
checkbc(𝛟::AbstractVector, M, N, v) = _checkvector(checkbc, 𝛟, M, N, v)

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

function checksquare(ϕ::AbstractMatrix, v)
    M, N = size(ϕ)
    xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ = map(Int64, (M / 2, M * 3//4, N * 5//8, N * 7//8))
    for i in xₘᵢₙ:xₘₐₓ
        for j in yₘᵢₙ:yₘₐₓ
            @assert ϕ[i, j] == v
        end
    end
    return nothing
end
checksquare(𝛟::AbstractVector, M, N, v) = _checkvector(checksquare, 𝛟, M, N, v)

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

function checkcharges(ρ::AbstractMatrix, v)
    M, N = size(ρ)
    x₁, x₂, y = map(Int64, (M / 4, M * 3//4, N / 8))
    @assert ρ[x₁, y] == v
    @assert ρ[x₂, y] == v
    return nothing
end
checkcharges(𝛒::AbstractVector, M, N, v) = _checkvector(checkcharges, 𝛒, M, N, v)

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

_checkvector(f::Function, 𝐯::AbstractVector, M, N, value) = f(reshape(𝐯, M, N), value)

end

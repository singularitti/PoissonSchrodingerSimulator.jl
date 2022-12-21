module ConjugateGradient

using LinearAlgebra: dot, norm
using OffsetArrays: OffsetVector, Origin

export solve, isconverged, eachstep

mutable struct IterationStep
    n::UInt64
    alpha::Float64
    beta::Float64
    x::Vector{Float64}
    r::Vector{Float64}
    p::Vector{Float64}
end

mutable struct ConvergenceHistory
    maxiter::UInt64
    isconverged::Bool
    data::OffsetVector{IterationStep}
end

function solve(A, 𝐛, 𝐱₀=zeros(length(𝐛)); atol=eps(), maxiter=2000)
    history = ConvergenceHistory(maxiter, false, OffsetVector([], Origin(0)))
    𝐱ₙ = 𝐱₀
    𝐫ₙ = 𝐛 - A * 𝐱ₙ  # Initial residual, 𝐫₀
    𝐩ₙ = 𝐫ₙ  # Initial momentum, 𝐩₀
    for n in 0:maxiter
        if norm(𝐫ₙ) < atol
            history.isconverged = true
            break
        end
        αₙ = compute_alpha(A, 𝐫ₙ, 𝐩ₙ)
        𝐱ₙ₊₁ = 𝐱ₙ + αₙ * 𝐩ₙ
        𝐫ₙ₊₁ = 𝐫ₙ - αₙ * A * 𝐩ₙ
        βₙ = compute_beta(𝐫ₙ₊₁, 𝐫ₙ)
        𝐩ₙ₊₁ = 𝐫ₙ₊₁ + βₙ * 𝐩ₙ
        push!(history.data, IterationStep(n, αₙ, βₙ, 𝐱ₙ, 𝐫ₙ, 𝐩ₙ))
        # Prepare for a new iteration
        𝐱ₙ, 𝐫ₙ, 𝐩ₙ = 𝐱ₙ₊₁, 𝐫ₙ₊₁, 𝐩ₙ₊₁
    end
    return 𝐱ₙ, history
end

compute_alpha(A, 𝐫, 𝐩) = dot(𝐫, 𝐫) / dot(𝐩, A, 𝐩)

compute_beta(A, 𝐫, 𝐩) = -dot(𝐩, A, 𝐫) / dot(𝐩, A, 𝐩)
compute_beta(𝐫ₙ₊₁, 𝐫ₙ) = dot(𝐫ₙ₊₁, 𝐫ₙ₊₁) / dot(𝐫ₙ, 𝐫ₙ)

isconverged(ch::ConvergenceHistory) = ch.isconverged

struct EachStep
    history::ConvergenceHistory
end

eachstep(ch::ConvergenceHistory) = EachStep(ch)

Base.iterate(iter::EachStep) = iterate(iter.history.data)
Base.iterate(iter::EachStep, state) = iterate(iter.history.data, state)

Base.eltype(::EachStep) = IterationStep

Base.length(iter::EachStep) = length(iter.history.data)

Base.size(iter::EachStep, dim...) = size(iter.history.data, dim...)

Base.getindex(iter::EachStep, i) = getindex(iter.history.data, i)

Base.firstindex(iter::EachStep) = firstindex(iter.history.data)

Base.lastindex(iter::EachStep) = lastindex(iter.history.data)

function Base.show(io::IO, step::IterationStep)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(step)
        Base.show_default(IOContext(io, :limit => true), step)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, summary(step))
        println(io, " n = ", Int(step.n))
        println(io, " α = ", step.alpha)
        println(io, " β = ", step.beta)
        println(io, " 𝐱 = ", step.x)
        println(io, " 𝐫 = ", step.r)
        println(io, " 𝐩 = ", step.p)
    end
end

end
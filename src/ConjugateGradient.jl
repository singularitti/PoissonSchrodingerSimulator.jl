module ConjugateGradient


struct IterationStep
    n::UInt64
    alpha::Float64
    beta::Float64
    x::Vector{Float64}
    r::Vector{Float64}
    p::Vector{Float64}
end

mutable struct ConvergenceHistory
    isconverged::Bool
    history::Vector{IterationStep}
end

end

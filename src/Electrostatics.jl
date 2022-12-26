module Electrostatics

using LinearAlgebra: norm, ⋅

using ..LastHomework:
    DiscreteLaplacian, Boundary, InternalSquare, PointCharges, validate, setvalues!
using ..ConjugateGradient: IterationStep, setconverged!, log!

import ..ConjugateGradient: solve!

function solve!(
    logger,
    A::DiscreteLaplacian,
    𝐛,
    𝐱₀;
    atol=eps(),
    maxiter=2000,
    charge=-20,
    bc=0,
    ext_pot=5,
)
    N = Int(sqrt(length(𝐛)))
    BOUNDARY = Boundary((N, N), bc)
    SQUARE = InternalSquare((N, N), ext_pot)
    SQUARE_RESIDUAL = InternalSquare((N, N), 0)
    setvalues!(𝐱₀, BOUNDARY)
    setvalues!(𝐱₀, SQUARE)
    setvalues!(𝐛, PointCharges((N, N), charge))
    𝐱ₙ = copy(𝐱₀)
    𝐫ₙ = 𝐛 - A * 𝐱ₙ  # Initial residual, 𝐫₀
    𝐩ₙ = 𝐫ₙ  # Initial momentum, 𝐩₀, notice that if you cahnge 𝐩ₙ, 
    for n in 0:maxiter
        if norm(𝐫ₙ) < atol
            setconverged!(logger)
            break
        end
        setvalues!(𝐩ₙ, BOUNDARY)
        setvalues!(𝐩ₙ, SQUARE_RESIDUAL)  # Set 𝐩ₙ and 𝐫₀ to zero in the square
        A𝐩ₙ = A * 𝐩ₙ  # Avoid running it multiple times
        αₙ = 𝐫ₙ ⋅ 𝐫ₙ / 𝐩ₙ ⋅ A𝐩ₙ  # `⋅` means dot product between two vectors
        𝐱ₙ₊₁ = 𝐱ₙ + αₙ * 𝐩ₙ
        αₙA𝐩ₙ = αₙ * A𝐩ₙ
        setvalues!(αₙA𝐩ₙ, BOUNDARY)
        setvalues!(αₙA𝐩ₙ, SQUARE_RESIDUAL)  # Set αₙA𝐩ₙ to 0 in the square
        𝐫ₙ₊₁ = 𝐫ₙ - αₙA𝐩ₙ
        βₙ = 𝐫ₙ₊₁ ⋅ 𝐫ₙ₊₁ / 𝐫ₙ ⋅ 𝐫ₙ
        𝐩ₙ₊₁ = 𝐫ₙ₊₁ + βₙ * 𝐩ₙ
        log!(logger, IterationStep(n, αₙ, βₙ, 𝐱ₙ, 𝐫ₙ, 𝐩ₙ))
        𝐱ₙ, 𝐫ₙ, 𝐩ₙ = 𝐱ₙ₊₁, 𝐫ₙ₊₁, 𝐩ₙ₊₁  # Prepare for a new iteration
        # validate(𝐫ₙ, SQUARE_RESIDUAL)
        # validate(𝐱ₙ, SQUARE)  # Check if 𝐱ₙ is 5 in the square
    end
    return 𝐱ₙ
end

end

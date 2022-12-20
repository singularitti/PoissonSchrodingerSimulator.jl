module Electrostatics

using LinearAlgebra: Symmetric, diagm
using SparseArrays: sparse

export discretize_laplacian

function discretize_laplacian(N)
    A = diagm(
        0 => fill(-4, N^2),
        1 => fill(1, N^2 - 1),
        N - 1 => fill(1, N^2 - N + 1),
        N => fill(1, N^2 - N),
        N^2 - N => fill(1, N),
    )  # An upper triangular matrix
    return sparse(Symmetric(A))
end

end

module Electrostatics

using LinearAlgebra: Symmetric
using SparseArrays: SparseMatrixCSC, spdiagm

export DiscreteLaplacian

struct DiscreteLaplacian <: AbstractMatrix{Int64}
    data::SparseMatrixCSC{Int64,Int64}
end
function DiscreteLaplacian(N::Integer)
    A = spdiagm(
        0 => fill(-4, N^2),
        1 => fill(1, N^2 - 1),
        N - 1 => fill(1, N^2 - N + 1),
        N => fill(1, N^2 - N),
        N^2 - N => fill(1, N),
    )  # An upper triangular matrix
    return DiscreteLaplacian(Symmetric(A))
end

end

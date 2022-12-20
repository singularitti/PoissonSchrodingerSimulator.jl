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

Base.parent(S::DiscreteLaplacian) = S.data

Base.size(S::DiscreteLaplacian) = size(parent(S))

Base.IndexStyle(::Type{DiscreteLaplacian}) = IndexLinear()

Base.getindex(S::DiscreteLaplacian, i) = getindex(parent(S), i)

Base.setindex!(S::DiscreteLaplacian, v, i) = setindex!(parent(S), v, i)

end

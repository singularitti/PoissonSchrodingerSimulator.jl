using LinearAlgebra: Symmetric, issymmetric, isposdef, det
using SparseArrays: AbstractSparseMatrix, SparseMatrixCSC, spdiagm

using ToggleableAsserts: @toggled_assert

export DiscreteLaplacian

struct DiscreteLaplacian <: AbstractSparseMatrix{Int64,Int64}
    parent::SparseMatrixCSC{Int64,Int64}
    function DiscreteLaplacian(parent::AbstractMatrix)
        @toggled_assert issymmetric(parent)
        @toggled_assert all(iszero(sum(parent; dims=2)))
        @toggled_assert isposdef(parent)
        return new(parent)
    end
end
function DiscreteLaplacian(N::Integer)
    A = spdiagm(
        0 => fill(4, N^2),
        1 => [mod(i, N) == N - 1 ? 0 : -1 for i in 0:(N^2 - 2)],
        N - 1 => [mod(i, N) == 0 ? -1 : 0 for i in 0:(N^2 - N)],
        N => fill(-1, N^2 - N),
        N^2 - N => fill(-1, N),
    )  # An upper triangular matrix
    return DiscreteLaplacian(Symmetric(A))
end

leading_principal_minor(A::AbstractMatrix, k) = A[begin:k, begin:k]

Base.parent(S::DiscreteLaplacian) = S.parent

Base.size(S::DiscreteLaplacian) = size(parent(S))

Base.IndexStyle(::Type{DiscreteLaplacian}) = IndexLinear()

Base.getindex(S::DiscreteLaplacian, i) = getindex(parent(S), i)

Base.setindex!(S::DiscreteLaplacian, v, i) = setindex!(parent(S), v, i)

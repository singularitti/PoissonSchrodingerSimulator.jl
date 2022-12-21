using LinearAlgebra: Symmetric
using SparseArrays: AbstractSparseMatrix, SparseMatrixCSC, spdiagm

export DiscreteLaplacianPBCs

struct DiscreteLaplacianPBCs <: AbstractSparseMatrix{Int64,Int64}
    data::SparseMatrixCSC{Int64,Int64}
end
function DiscreteLaplacianPBCs(N::Integer)
    A = spdiagm(
        0 => fill(-4, N^2),
        1 => fill(1, N^2 - 1),
        N - 1 => fill(1, N^2 - N + 1),
        N => fill(1, N^2 - N),
        N^2 - N => fill(1, N),
    )  # An upper triangular matrix
    return DiscreteLaplacianPBCs(Symmetric(A))
end

Base.parent(S::DiscreteLaplacianPBCs) = S.data

Base.size(S::DiscreteLaplacianPBCs) = size(parent(S))

Base.IndexStyle(::Type{DiscreteLaplacianPBCs}) = IndexLinear()

Base.getindex(S::DiscreteLaplacianPBCs, i) = getindex(parent(S), i)

Base.setindex!(S::DiscreteLaplacianPBCs, v, i) = setindex!(parent(S), v, i)

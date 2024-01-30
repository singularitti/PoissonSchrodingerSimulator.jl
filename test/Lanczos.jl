module Lanczos

using PoissonSchrodingerSimulator.Lanczos: lanczos
using LinearAlgebra: ishermitian
using Test: @testset, @test

@testset "" begin
    N = 100
    A = rand(N, N) .- 1 / 2
    A = (A + A') / 2
    @test ishermitian(A)
    T, Q = lanczos(A)
    @test T == Q' * A * Q
    @test A == Q * T * Q'
end

end

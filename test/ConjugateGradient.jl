module ConjugateGradient

using LastHomework: solve

# Example is from https://en.wikipedia.org/wiki/Conjugate_gradient_method#Numerical_example
@testset "Test Wikipedia example" begin
    A = [
        4 1
        1 3
    ]
    𝐛 = [1, 2]
    𝐱₀ = [2, 1]
    𝐱, ch = solve(A, 𝐛, 𝐱₀, 1e-24)
    @test 𝐱 ≈ [1 / 11, 7 / 11]  # Compare with the exact solution
    @test norm(A * 𝐱 - 𝐛) / norm(𝐛) ≤ 1e-12
    @test isconverged(ch) == true
    steps = eachstep(ch)
    @test steps[0].r == steps[0].p == -[8, 3]
    @test steps[0].alpha == 73 / 331
    @test steps[0].beta ≈ 0.008771369374138607
    @test steps[1].x == [2, 1] - 73 / 331 * [8, 3]
    @test steps[1].r == -[8, 3] + 73 / 331 * [4 1; 1 3] * [8, 3]
    @test steps[1].p ≈ [-0.3511377223647101, 0.7229306048685207]
    @test steps[1].alpha ≈ 0.4122042341220423
    @test steps[2].x ≈ [0.09090909090909094, 0.6363636363636365]
end

end

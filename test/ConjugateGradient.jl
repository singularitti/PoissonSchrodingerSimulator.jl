module ConjugateGradient

using LastHomework.ConjugateGradient: Logger, solve, solve!, isconverged, eachstep
using LinearAlgebra: norm
using Test: @testset, @test

# Example is from https://en.wikipedia.org/wiki/Conjugate_gradient_method#Numerical_example
@testset "Test Wikipedia example" begin
    A = [
        4 1
        1 3
    ]
    𝐛 = [1, 2]
    𝐱₀ = [2, 1]
    logger = Logger(2000)
    𝐱 = solve!(logger, A, 𝐛, 𝐱₀; atol=1e-24)
    @test 𝐱 ≈ [1 / 11, 7 / 11]  # Compare with the exact solution
    @test norm(A * 𝐱 - 𝐛) / norm(𝐛) ≤ 1e-12
    @test isconverged(logger) == true
    steps = eachstep(logger)
    @test steps[0].r == steps[0].p == -[8, 3]
    @test steps[0].alpha == 73 / 331
    @test steps[0].beta ≈ 0.008771369374138607
    @test steps[1].x == [2, 1] - 73 / 331 * [8, 3]
    @test steps[1].r == -[8, 3] + 73 / 331 * [4 1; 1 3] * [8, 3]
    @test steps[1].p ≈ [-0.3511377223647101, 0.7229306048685207]
    @test steps[1].alpha ≈ 0.4122042341220423
    @test steps[2].x ≈ [0.09090909090909094, 0.6363636363636365]
end

# Example is from https://optimization.mccormick.northwestern.edu/index.php/Conjugate_gradient_methods#Numerical_Example_of_the_method
@testset "Test Erik Zuehlke's example" begin
    A = [
        5 1
        1 2
    ]
    𝐛 = [2, 2]
    𝐱₀ = [1, 2]
    logger = Logger(2000)
    𝐱 = solve!(logger, A, 𝐛, 𝐱₀; atol=1e-24)
    @test 𝐱 ≈ [0.2222222222222221, 0.8888888888888891]  # Compare with other's result
    @test norm(A * 𝐱 - 𝐛) / norm(𝐛) == 0
    @test isconverged(logger) == true
    steps = eachstep(logger)
    @test steps[0].r == steps[0].p == -[5, 3]
    @test steps[0].alpha == 34 / 173
    @test steps[0].beta ≈ 0.028099836279194094  # The example's result is wrong
    @test steps[1].x == [1, 2] - 34 / 173 * [5, 3]
    @test steps[1].r == -[5, 3] + 34 / 173 * [5 1; 1 2] * [5, 3]
    @test steps[1].p ≈ [0.3623909920144345, -0.9224497978549232]
end

# See https://towardsdatascience.com/complete-step-by-step-conjugate-gradient-algorithm-from-scratch-202c07fb52a8
@testset "Test Albers Uzila's examples" begin
    @testset "Problem 1" begin
        A = [
            2.5409 -0.0113
            -0.0113 0.5287
        ]
        𝐛 = [1.3864, 0.3719]
        logger = Logger(2000)
        𝐱 = solve!(logger, A, 𝐛, -[3, 4])
        @test 𝐱 ≈ [0.5488138979502294, 0.7151533895344008]
        @test norm(A * 𝐱 - 𝐛) / norm(𝐛) < 1e-15
        @test isconverged(logger) == true
        steps = eachstep(logger)
        @test steps[1].x ≈ [0.742786502583181, -2.975857971024216]
        @test norm(steps[1].r) ≈ 2.025447442457243
        @test steps[2].x ≈ [0.5488138979502315, 0.7151533895344007]
        @test norm(steps[2].r) < 1e-14
    end
    @testset "Problem 2" begin
        A = [
            0.7444 -0.5055 -0.0851
            -0.5055 3.4858 0.0572
            -0.0851 0.0572 0.4738
        ]
        𝐛 = [-0.0043, 2.2501, 0.2798]
        logger = Logger(2000)
        𝐱 = solve!(logger, A, 𝐛, [3, 1, -7])
        @test 𝐱 ≈ [0.5488032997143618, 0.7151992261015149, 0.6027728262403653]
        @test norm(A * 𝐱 - 𝐛) / norm(𝐛) < 1e-15
        @test isconverged(logger) == true
    end
    @testset "Problem 3" begin
        A = [
            3.4430 -0.3963 2.5012 0.9525 0.6084 -1.2728
            -0.3963 0.6015 -0.4108 -0.1359 -0.0295 0.2630
            2.5012 -0.4108 2.5927 0.7072 0.5587 -1.0613
            0.9525 -0.1359 0.7072 1.1634 0.1920 -0.4344
            0.6084 -0.0295 0.5587 0.1920 0.7636 -0.3261
            -1.2728 0.2630 -1.0613 -0.4344 -0.3261 1.0869
        ]
        𝐛 = [3.0685, 0.0484, 2.5783, 1.2865, 0.8671, -0.8230]
        𝐱 = solve(A, 𝐛, [9, 0, -2, 3, -2, 5])
        @test 𝐱 ≈ [
            0.5488252073566335,
            0.7152045853108671,
            0.6027868107757592,
            0.5448522879867495,
            0.4236962375185705,
            0.6459055453301382,
        ]
    end
end

end

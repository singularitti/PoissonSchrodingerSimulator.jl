using LinearAlgebra: norm
using RecipesBase: @userplot, @recipe, @series

using .ConjugateGradient: ConvergenceHistory

export residualplot

@userplot ResidualPlot
@recipe function f(plot::ResidualPlot)
    # See http://juliaplots.org/RecipesBase.jl/stable/types/#User-Recipes-2
    history = plot.args[end]
    residuals = [norm(step.r) for step in history.data]
    # If we are passed two args, we use the first as labels
    steps = length(plot.args) == 2 ? plot.args[1] : eachindex(residuals)
    size --> (700, 400)
    markersize --> 2
    markerstrokecolor --> :auto
    markerstrokewidth --> 0
    xlims --> extrema(steps)
    ylims --> extrema(residuals)
    xguide --> raw"iteration step ($n$)"
    yguide --> raw"residual $\Vert \mathbf{b} - \mathrm{A} \mathbf{x} \Vert$"
    guidefontsize --> 10
    tickfontsize --> 8
    legendfontsize --> 8
    legend_foreground_color --> nothing
    legend_position --> :topright
    frame --> :box
    palette --> :tab20
    grid --> nothing
    @series begin
        seriestype --> :scatter
        steps, residuals
    end
    @series begin
        seriestype --> :path
        z_order --> :back
        label := ""
        steps, residuals
    end
end

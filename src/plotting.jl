using LinearAlgebra: norm
using RecipesBase: @userplot, @recipe, @series

using .ConjugateGradient: ConvergenceHistory

export regionheatmap, residualplot, surfaceplot

@recipe function f(laplacian::DiscreteLaplacian)
    size --> (800, 800)
    seriestype --> :heatmap
    yflip --> true  # Set the origin to the upper left corner, see https://github.com/MakieOrg/Makie.jl/issues/46
    xlims --> extrema(axes(laplacian, 1)) .+ (-0.5, 0.5)  # See https://discourse.julialang.org/t/can-plots-jl-heatmap-coordinates-start-at-1-instead-of-0-5/90385/3
    ylims --> extrema(axes(laplacian, 2)) .+ (-0.5, 0.5)
    xguide --> raw"columns"  # Note that it is reversed!
    xguideposition --> :top  # Place xguide along top axis
    xmirror --> true  # Place xticks along top axis, see https://github.com/JuliaPlots/Plots.jl/issues/337
    yguide --> raw"rows"  # Note that it is reversed!
    tick_direction --> :out
    guidefontsize --> 12
    tickfontsize --> 10
    color --> :RdGy_3  # Tricolor colormap
    color_limits --> extrema(laplacian)
    colorbar --> :none
    frame --> :box
    aspect_ratio --> :equal
    margins --> (0, :mm)  # See https://github.com/JuliaPlots/Plots.jl/issues/4522#issuecomment-1318511879
    right_margin --> (1.5, :mm)
    return axes(laplacian)..., laplacian
end

@userplot RegionHeatMap
@recipe function f(plot::RegionHeatMap)
    data = plot.args[end]
    size --> (800, 800)
    seriestype --> :heatmap
    yflip --> true  # Set the origin to the upper left corner, see https://github.com/MakieOrg/Makie.jl/issues/46
    xlims --> extrema(axes(data, 1)) .+ (-0.5, 0.5)  # See https://discourse.julialang.org/t/can-plots-jl-heatmap-coordinates-start-at-1-instead-of-0-5/90385/3
    ylims --> extrema(axes(data, 2)) .+ (-0.5, 0.5)
    xguide --> raw"$x$"
    xguideposition --> :top  # Place xguide along top axis
    xmirror --> true  # Place xticks along top axis, see https://github.com/JuliaPlots/Plots.jl/issues/337
    yguide --> raw"$y$"
    tick_direction --> :out
    guidefontsize --> 12
    tickfontsize --> 10
    color --> :thermometer
    frame --> :box
    aspect_ratio --> :equal
    margins --> (0, :mm)  # See https://github.com/JuliaPlots/Plots.jl/issues/4522#issuecomment-1318511879
    return axes(data)..., data
end

@userplot SurfacePlot
@recipe function f(plot::SurfacePlot)
    data = plot.args[end]
    size --> (800, 800)
    seriestype --> :surface
    yflip --> true  # Set the origin to the upper left corner, see https://github.com/MakieOrg/Makie.jl/issues/46
    xlims --> extrema(axes(data, 1)) .+ (-0.5, 0.5)  # See https://discourse.julialang.org/t/can-plots-jl-heatmap-coordinates-start-at-1-instead-of-0-5/90385/3
    ylims --> extrema(axes(data, 2)) .+ (-0.5, 0.5)
    xguide --> raw"$x$"
    yguide --> raw"$y$"
    zguide --> raw"$z$"
    tick_direction --> :out
    guidefontsize --> 10
    tickfontsize --> 8
    color --> :thermometer
    frame --> :box
    aspect_ratio --> :equal  # See https://docs.juliaplots.org/latest/gallery/gr/generated/gr-ref060/
    margins --> (0, :mm)  # See https://github.com/JuliaPlots/Plots.jl/issues/4522#issuecomment-1318511879
    return axes(data)..., data
end

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

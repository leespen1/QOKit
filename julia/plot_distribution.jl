using GLMakie, NPZ
"""
Plot 2D slices of N(c', d, c) as slices

If two parameters are varied, then a 3D surface is plotted.

If all three parameters are varied, then the surfaces will be plotted on top of each other
transparently.

WARNING - Technically these results are off by one, since N[1,1,1] corresponds to N(0,0,0)

"""
function plot_dist_surface(N, c_prime::Union{Colon, Integer}=:, d::Union{Colon, Integer}=:, c::Union{Colon, Integer}=:,
        ; shading=true, colormap=:viridis, num_cols=4, common_zlims=true, common_colorrange=true)
    if isa(N, String) # Convert from numpy file to array
        N_ary = npzread(N)
    else
        N_ary = N
    end
    # Colormaps: deep, darkterrain, viridis, inferno, etc. Can check available ones with Makie.available_gradients()
    params = (c_prime, d, c)
    param_labels = ("c'", "d", "c")
    num_plot_params = count(isequal(:), (c_prime, d, c))
    @assert num_plot_params >= 2 "At least two of the three parameters must be ':'. " 

    first_param_i = findfirst(isequal(:), params)
    second_param_i = isnothing(first_param_i) ? nothing : findnext(isequal(:), params, first_param_i+1)
    third_param_i = isnothing(first_param_i) ? nothing : findnext(isequal(:), params, second_param_i+1)

    xlabel = param_labels[first_param_i]
    ylabel = param_labels[second_param_i]
    zlabel = "N(c'; d, c)"


    if num_plot_params == 2
        title = ""
        fig, ax, surf = surface(
            N_ary[c_prime, d, c],
            axis=(type=Axis3, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel),
            shading=shading, colormap=colormap
        )
        Colorbar(fig[1,2], surf)
        
    else
        fig = Figure()
        max_N = maximum(N_ary)
        colorrange = common_colorrange ? (0, max_N) : Makie.automatic
        for (i, n) in enumerate(axes(N_ary, 1))
            fig_row, fig_col = divrem(i-1, num_cols) .+ 1
            ax = Axis3(fig[1+fig_row, 1+fig_col], limits=(nothing, nothing, (0, max_N)), xlabel="c", ylabel="d", zlabel="N(c'; d, c)")
            #surface!(ax, N_ary[n, :, :], colorrange=(0, max_N))
            surface!(ax, N_ary[n, :, :], colorrange=Makie.automatic)
            ax.title = "c'=$n"
        end
    end


    return fig
end

function plot_dist_volume(N; colormap=:viridis, transparency=false)
    if isa(N, String) # Convert from numpy file to array
        N_ary = npzread(N)
    else
        N_ary = N
    end
    xlabel="c'"
    ylabel="d"
    zlabel="c"
    return volume(
        N_ary,
        axis=(type=Axis3, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel),
        transparency=false, colormap=colormap
    )
end


using GLMakie, NPZ, JuliaQAOA, NativeFileDialog
"""
Plot 2D slices of N(c', d, c) as slices

If two parameters are varied, then a 3D surface is plotted.

If all three parameters are varied, then the surfaces will be plotted on top of each other
transparently.

WARNING - Technically these results are off by one, since N[1,1,1] corresponds to N(0,0,0)

"""
function plot_dist_surface(N, c_prime::Union{Colon, Integer}=:, d::Union{Colon, Integer}=:, c::Union{Colon, Integer}=:,
        ; shading=true, colormap=:viridis, num_cols=4, common_zlims=true, common_colorrange=true)
    if N isa String # Convert from numpy file to array
        N_ary = npzread(N)
    elseif N isa JuliaQAOA.TriangleProxy
        N_ary = JuliaQAOA.cpu_compute_homodist(N)
    elseif N isa Array
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
        title = "TODO display fixed parameter here"
        fig, ax, surf = surface(
            N_ary[c_prime, d, c],
            axis=(type=Axis3, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel),
            shading=shading, colormap=colormap
        )
        Colorbar(fig[1,2], surf)
        return fig        
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
    if N isa String # Convert from numpy file to array
        N_ary = npzread(N)
    elseif N isa JuliaQAOA.TriangleProxy
        N_ary = JuliaQAOA.cpu_compute_homodist(N)
    elseif N isa Array
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

function interactive_plot(;shading=true, colormap=:viridis)
    # Create observables for parameters
    p1 = Observable(10) # num constraints
    p2 = Observable(5) # num qubits
    p3 = Observable(1.0) # height adjustment
    p4 = Observable(0.5) # centeradjustment
    p5 = Observable(0.25) # left angle
    p6 = Observable(0.25) # right angle
    c_prime_obs = Observable(7) 


    proxy_obs = lift(p1, p2, p3, p4, p5, p6) do p1_, p2_, p3_, p4_, p5_, p6_
        IntuitiveTriangleProxy(p1_, p2_, p3_, p4_, p5_, p6_) 
    end

    N_ary_obs = lift(
        (proxy, c_prime) -> JuliaQAOA.cpu_compute_homodist(proxy)[c_prime,:,:],
        proxy_obs, c_prime_obs
    )

    title = lift(c_prime -> "c' = $(c_prime)", c_prime_obs)
    fig, ax, surf = surface(
        N_ary_obs,
        axis=(type=Axis3, xlabel="d", ylabel="c", zlabel="N(c'; d, c)", title=title),
        shading=shading, colormap=colormap
    )
    Colorbar(fig[1,1][1,2], surf)

    on(proxy_obs) do _
        autolimits!(ax)
    end

    # UI fields for input
    fig[2, 1] = grid = GridLayout()

    int_fields = [
        Textbox(fig, placeholder="Num Constraints", validator=Int64, tellwidth=false),
        Textbox(fig, placeholder="Num Nodes",       validator=Int64, tellwidth=false),
        Textbox(fig, placeholder="Fixed c'",        validator=Int64, tellwidth=false),
    ]

    float_fields = [
        Textbox(fig, placeholder="Height Adjustment", validator=Float64, tellwidth=false),
        Textbox(fig, placeholder="Center Adjustment", validator=Float64, tellwidth=false),
        Textbox(fig, placeholder="Left Angle",        validator=Float64, tellwidth=false),
        Textbox(fig, placeholder="Right Angle",       validator=Float64, tellwidth=false)
    ]

    for (i, box) in enumerate(int_fields)
        grid[1, i] = box
    end

    for (i, box) in enumerate(float_fields)
        grid[2, i] = box
    end

    # Connect text boxes to observables
    on(int_fields[1].stored_string) do str
        try p1[] = parse(Int, str) catch end
    end
    on(int_fields[2].stored_string) do str
        try p2[] = parse(Int, str) catch end
    end
    on(int_fields[3].stored_string) do str
        try c_prime_obs[] = parse(Int, str) catch end
    end

    on(float_fields[1].stored_string) do str
        try p3[] = parse(Float64, str) catch end
    end
    on(float_fields[2].stored_string) do str
        try p4[] = parse(Float64, str) catch end
    end
    on(float_fields[3].stored_string) do str
        try p5[] = parse(Float64, str) catch end
    end
    on(float_fields[4].stored_string) do str
        try p6[] = parse(Float64, str) catch end
    end

    return fig
end

"""
Parse string as either an integer, a colon, or a unitrange (n1:n2).
"""
function parse_token(s::AbstractString)
    s = strip(s)

    # literal colon
    if s == ":"
        return Colon()  # or just :
    end

    # try integer
    if occursin(r"^-?\d+$", s)
        return parse(Int, s)
    end

    # try unitrange a:b
    m = match(r"^(-?\d+)\s*:\s*(-?\d+)$", s)
    if m !== nothing
        lo = parse(Int, m.captures[1])
        hi = parse(Int, m.captures[2])
        return lo:hi
    end

    error("Unrecognized token: $s")
end



function interactive_sample(;shading=true, colormap=:viridis,  num_cols=4,
    common_zlims=true, common_colorrange=true
)

    cprime_range_obs = Observable{Any}(1)

    fig = Figure()
    axfig = fig[1,1]
    fig[2, 1] = grid = GridLayout()

    cprime_range_box = Textbox(fig, placeholder="c' range", tellwidth=false)
    grid[1,1] = cprime_range_box

    on(cprime_range_box.stored_string) do str
        cprime_range_obs[] = parse_token(str)
    end

    filename = ""

    function update()
        if isempty(filename) # Do nothing if empty filename
            return
        end

        # Hacky way to clear old figures. Once layout is expanded, can't be retracted
        axfig_contents = contents(axfig)
        if !isempty(axfig_contents)
            for gridlayout in axfig_contents
                foreach(delete!, contents(gridlayout))
            end
        end

        N_ary = npzread(filename)
        max_N = maximum(N_ary)
        colorrange = common_colorrange ? (0, max_N) : Makie.automatic

        cprime_range = cprime_range_obs[]
        indices = isa(cprime_range, Colon) ? (1:size(N_ary, 1)) : cprime_range

        for (i, n) in enumerate(indices)
            fig_row, fig_col = divrem(i-1, num_cols) .+ 1
            #empty!(contents(axfig[1+fig_row, 1+fig_col]))
            ax = Axis3(axfig[1+fig_row, 1+fig_col], limits=(nothing, nothing, (0, max_N)), xlabel="d", ylabel="c", zlabel="N(c'; d, c)")
            ax.title = "c'=$n"
            #push!(axis3_vec, ax)
            surface!(ax, N_ary[n, :, :], colorrange=Makie.automatic)
        end
    end

    btn_RUN = Button(grid[1,2], label = " Open file... ")
    on(btn_RUN.clicks) do c
        #@async begin
        @sync begin
            filename = fetch(Threads.@spawn pick_file(""))
        end
        update()
    end
    on(cprime_range_obs) do _
        update()
    end
    return fig
end

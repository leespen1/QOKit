#!/usr/bin/env julia
#
# plot_Nc.jl — Read N(c) distributions from HDF5 files produced by the Python
# Nc_*.py scripts and plot the cost distributions.
#
# Usage:
#   julia plot_Nc.jl                     # scan current directory
#   julia plot_Nc.jl path/to/datadir     # scan a specific directory
#   julia plot_Nc.jl file1.h5 file2.h5   # plot specific files
#
# Requires: HDF5, CairoMakie

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "julia"))

using HDF5
using CairoMakie
using Statistics
using Printf

"""
    find_nc_h5_files(root)

Recursively find all `.h5` files under `root` that contain N(c) distributions
(i.e. exclude files ending in `_costs.h5` and `_homodist.h5`).
"""
function find_nc_h5_files(root::AbstractString)
    files = String[]
    for (dirpath, _, filenames) in walkdir(root)
        for f in filenames
            if endswith(f, ".h5") && !endswith(f, "_costs.h5") && !endswith(f, "_homodist.h5")
                push!(files, joinpath(dirpath, f))
            end
        end
    end
    sort!(files)
    return files
end

"""
    read_nc_data(filepath) -> (Nc, attrs)

Read an N(c) HDF5 file. Returns the Nc matrix (num_graphs × num_costs) and a
Dict of attributes.
"""
function read_nc_data(filepath::AbstractString)
    h5open(filepath, "r") do f
        if !haskey(f, "Nc")
            error("File $filepath does not contain an 'Nc' dataset")
        end
        Nc = read(f, "Nc")
        attrs = Dict{String,Any}()
        for key in keys(HDF5.attributes(f))
            attrs[key] = read(HDF5.attributes(f), key)
        end
        return Nc, attrs
    end
end

"""
    trim_padding(row) -> trimmed

Remove trailing zeros that are just padding (keep at most one trailing zero
to preserve the shape of distributions that genuinely reach zero).
"""
function trim_padding(row::AbstractVector)
    last_nonzero = findlast(!iszero, row)
    if last_nonzero === nothing
        return row[1:1]
    end
    # keep one zero past the last nonzero if it exists, to show the tail
    endpoint = min(last_nonzero + 1, length(row))
    return row[1:endpoint]
end

"""
    make_label(attrs) -> String

Build a human-readable label from the HDF5 attributes.
"""
function make_label(attrs::Dict)
    parts = String[]
    gt = get(attrs, "graphType", "Unknown")
    push!(parts, gt)
    if haskey(attrs, "numNodes")
        push!(parts, "n=$(attrs["numNodes"])")
    end
    if haskey(attrs, "edgesPerNode")
        push!(parts, "m=$(attrs["edgesPerNode"])")
    end
    if haskey(attrs, "edgeProbability")
        push!(parts, "p=$(attrs["edgeProbability"])")
    end
    if haskey(attrs, "nearestNeighbors")
        push!(parts, "k=$(attrs["nearestNeighbors"])")
    end
    if haskey(attrs, "rewiringProbability")
        push!(parts, "r=$(attrs["rewiringProbability"])")
    end
    return join(parts, ", ")
end

"""
    make_filename(attrs) -> String

Build a filename-safe string from the HDF5 attributes.
"""
function make_filename(attrs::Dict)
    parts = String[]
    gt = get(attrs, "graphType", "Unknown")
    push!(parts, gt)
    if haskey(attrs, "numNodes")
        push!(parts, "n=$(attrs["numNodes"])")
    end
    if haskey(attrs, "edgesPerNode")
        push!(parts, "m=$(attrs["edgesPerNode"])")
    end
    if haskey(attrs, "edgeProbability")
        push!(parts, "p=$(attrs["edgeProbability"])")
    end
    if haskey(attrs, "nearestNeighbors")
        push!(parts, "k=$(attrs["nearestNeighbors"])")
    end
    if haskey(attrs, "rewiringProbability")
        push!(parts, "r=$(attrs["rewiringProbability"])")
    end
    return join(parts, "_")
end

"""
    plot_nc!(ax, Nc, num_graphs)

Plot N(c) data on the given axis. Chooses bar plot, individual lines + mean,
or mean ± std band depending on the number of graphs.
"""
function plot_nc!(ax, Nc, num_graphs)
    if num_graphs == 1
        # Single graph: simple bar plot
        nc_vec = trim_padding(vec(Nc[:, 1]))
        costs = 0:(length(nc_vec) - 1)
        barplot!(ax, collect(costs), Float64.(nc_vec); color = :steelblue)
    elseif num_graphs <= 10
        # Few graphs: show individual lines + mean
        for g in 1:num_graphs
            nc_vec = trim_padding(vec(Nc[:, g]))
            costs = 0:(length(nc_vec) - 1)
            lines!(ax, collect(costs), Float64.(nc_vec); color = (:gray60, 0.5), linewidth = 1)
        end
        # Mean
        mean_nc = Float64.(mean(Nc, dims = 2)[:, 1])
        mean_trimmed = trim_padding(mean_nc)
        costs = 0:(length(mean_trimmed) - 1)
        lines!(ax, collect(costs), mean_trimmed; color = :steelblue, linewidth = 2.5, label = "Mean")
        axislegend(ax; position = :rt)
    else
        # Many graphs: show mean ± std as a band
        mean_nc = Float64.(mean(Nc, dims = 2)[:, 1])
        std_nc = Float64.(std(Nc, dims = 2)[:, 1])
        mean_trimmed = trim_padding(mean_nc)
        n = length(mean_trimmed)
        std_trimmed = std_nc[1:n]
        costs = collect(0:(n - 1))
        band!(ax, costs, mean_trimmed .- std_trimmed, mean_trimmed .+ std_trimmed;
              color = (:steelblue, 0.25))
        lines!(ax, costs, mean_trimmed; color = :steelblue, linewidth = 2.5, label = "Mean ± σ")
        axislegend(ax; position = :rt)
    end
end

function main()
    # Determine which files to plot
    args = ARGS
    h5files = String[]

    if isempty(args)
        # Default: scan subdirectories of the current directory
        for entry in readdir(".")
            if isdir(entry) && startswith(entry, "Data_")
                append!(h5files, find_nc_h5_files(entry))
            end
        end
    else
        for arg in args
            if isdir(arg)
                append!(h5files, find_nc_h5_files(arg))
            elseif isfile(arg) && endswith(arg, ".h5")
                push!(h5files, arg)
            else
                @warn "Skipping $arg (not a directory or .h5 file)"
            end
        end
    end

    if isempty(h5files)
        println("No N(c) HDF5 files found.")
        return
    end

    println("Found $(length(h5files)) N(c) file(s):")
    for f in h5files
        println("  $f")
    end

    # Read all datasets
    datasets = []
    for filepath in h5files
        try
            Nc, attrs = read_nc_data(filepath)
            push!(datasets, (; filepath, Nc, attrs))
        catch e
            @warn "Skipping $filepath" exception = e
        end
    end

    if isempty(datasets)
        println("No valid datasets found.")
        return
    end

    # Save individual plots
    individual_dir = "Nc_individual_plots"
    mkpath(individual_dir)
    for ds in datasets
        Nc = ds.Nc
        attrs = ds.attrs
        label = make_label(attrs)
        num_graphs = size(Nc, 2)

        ifig = Figure(size = (900, 350))
        ax = Axis(
            ifig[1, 1];
            xlabel = "Cost c",
            ylabel = "N(c)",
            title = "$label  ($num_graphs graph$(num_graphs > 1 ? "s" : ""))",
        )
        plot_nc!(ax, Nc, num_graphs)

        fname = "Nc_$(make_filename(attrs)).png"
        ipath = joinpath(individual_dir, fname)
        save(ipath, ifig)
        println("  Saved $ipath")
    end

    # Save a combined plot for a group of datasets
    function save_group_plot(group, outpath)
        n = length(group)
        fig = Figure(size = (900, 300 * n + 50))
        for (i, ds) in enumerate(group)
            Nc = ds.Nc
            attrs = ds.attrs
            label = make_label(attrs)
            num_graphs = size(Nc, 2)
            ax = Axis(
                fig[i, 1];
                xlabel = "Cost c",
                ylabel = "N(c)",
                title = "$label  ($num_graphs graph$(num_graphs > 1 ? "s" : ""))",
            )
            plot_nc!(ax, Nc, num_graphs)
        end
        save(outpath, fig)
        println("  Saved $outpath")
    end

    # Helper: get non-numNodes parameter keys for a dataset
    other_param_keys = ["edgesPerNode", "edgeProbability", "nearestNeighbors", "rewiringProbability"]

    combined_dir = "Nc_combined_plots"
    mkpath(combined_dir)

    # Group by graph type
    by_type = Dict{String,Vector}()
    for ds in datasets
        gt = get(ds.attrs, "graphType", "Unknown")
        push!(get!(by_type, gt, []), ds)
    end
    for (gt, group) in sort(collect(by_type))
        save_group_plot(group, joinpath(combined_dir, "Nc_$(gt).png"))
    end

    # Group by (graph type, numNodes)
    by_nodes = Dict{String,Vector}()
    for ds in datasets
        gt = get(ds.attrs, "graphType", "Unknown")
        nn = get(ds.attrs, "numNodes", "?")
        key = "$(gt)_n=$(nn)"
        push!(get!(by_nodes, key, []), ds)
    end
    for (key, group) in sort(collect(by_nodes))
        save_group_plot(group, joinpath(combined_dir, "Nc_$(key).png"))
    end

    # Group by (graph type, other params) — varying numNodes
    by_other = Dict{String,Vector}()
    for ds in datasets
        gt = get(ds.attrs, "graphType", "Unknown")
        parts = [gt]
        for k in other_param_keys
            if haskey(ds.attrs, k)
                push!(parts, "$(k)=$(ds.attrs[k])")
            end
        end
        key = join(parts, "_")
        push!(get!(by_other, key, []), ds)
    end
    for (key, group) in sort(collect(by_other))
        save_group_plot(group, joinpath(combined_dir, "Nc_$(key).png"))
    end

    # Grid plots: rows = numNodes, columns = other parameters, one per graph type
    for (gt, group) in sort(collect(by_type))
        # Build a lookup: (numNodes, other_params_string) -> dataset
        lookup = Dict{Tuple{Any,String},Any}()
        all_nodes = Set()
        all_params = Set{String}()
        for ds in group
            nn = get(ds.attrs, "numNodes", "?")
            parts = String[]
            for k in other_param_keys
                if haskey(ds.attrs, k)
                    push!(parts, "$(k)=$(ds.attrs[k])")
                end
            end
            param_key = isempty(parts) ? "" : join(parts, ", ")
            push!(all_nodes, nn)
            push!(all_params, param_key)
            lookup[(nn, param_key)] = ds
        end

        sorted_nodes = sort(collect(all_nodes))
        sorted_params = sort(collect(all_params))
        nrows = length(sorted_nodes)
        ncols = length(sorted_params)

        fig = Figure(size = (max(400, 350 * ncols), 250 * nrows + 80))
        Label(fig[0, 1:ncols], gt; fontsize = 20, font = :bold)

        for (ci, param) in enumerate(sorted_params)
            # Column header
            col_label = isempty(param) ? "" : param
            if !isempty(col_label)
                Label(fig[1, ci, Top()], col_label; fontsize = 14, padding = (0, 0, 5, 0))
            end

            for (ri, nn) in enumerate(sorted_nodes)
                ds = get(lookup, (nn, param), nothing)
                ax = Axis(
                    fig[ri, ci];
                    xlabel = "Cost c",
                    ylabel = "N(c)",
                    title = "n=$nn",
                    titlesize = 12,
                )
                if ds !== nothing
                    plot_nc!(ax, ds.Nc, size(ds.Nc, 2))
                end
            end
        end

        outpath = joinpath(combined_dir, "Nc_grid_$(gt).png")
        save(outpath, fig)
        println("  Saved $outpath")
    end

    println("\nAll plots saved.")
end

main()

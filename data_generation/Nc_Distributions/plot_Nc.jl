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

using HDF5
using GLMakie
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

    # Create figure — one row per dataset
    n_datasets = length(datasets)
    fig = Figure(size = (900, 300 * n_datasets + 50))

    for (i, ds) in enumerate(datasets)
        Nc = ds.Nc  # shape: (num_costs, num_graphs) — HDF5.jl reads column-major
        attrs = ds.attrs
        label = make_label(attrs)
        num_graphs = size(Nc, 2)

        ax = Axis(
            fig[i, 1];
            xlabel = "Cost c",
            ylabel = "N(c)",
            title = "$label  ($num_graphs graph$(num_graphs > 1 ? "s" : ""))",
        )

        if num_graphs == 1
            # Single graph: simple bar plot
            nc_vec = trim_padding(vec(Nc[:, 1]))
            costs = 0:(length(nc_vec) - 1)
            barplot!(ax, collect(costs), Float64.(nc_vec); color = :steelblue)
        elseif num_graphs <= 10
            # Few graphs: show individual lines + mean
            max_len = 0
            for g in 1:num_graphs
                nc_vec = trim_padding(vec(Nc[:, g]))
                max_len = max(max_len, length(nc_vec))
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

    outpath = "Nc_distributions.png"
    save(outpath, fig; px_per_unit = 2)
    println("\nPlot saved to $outpath")
end

main()

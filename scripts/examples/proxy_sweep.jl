using DrWatson
@quickactivate "JuliaQAOA"
using NPZ: npzwrite, npzread
using Dates: now
using DelimitedFiles, JuliaQAOA, Logging, Dates, ArgParse
using Printf: @sprintf

"""
Given a filename of the form X1=Y1_X2=Y2_[...].ext,
return a dictionary Dict(X1 => Y1, X2 => Y2, [...])
"""
function filename_dict(in_filename)
    in_filename_no_ext = splitext(basename(in_filename))[1]
    keyval_pair_strs = split(in_filename_no_ext, '_')
    keyval_dict = Dict((kv[1], kv[2]) for kv in split.(keyval_pair_strs, '='))
    return keyval_dict
end

function parse_filename_regex(filename::String)
    name = splitext(basename(filename))[1]

    # Regular expressions for each field
    backend     = match(r"backend=([^_]+)", name).captures[1]
    graphs      = parse(Int, match(r"graphs=(\d+)", name).captures[1])
    graphtype   = match(r"graphtype=([^_]+)", name).captures[1]
    neighbors   = parse(Int, match(r"neighbors=(\d+)", name).captures[1])
    numnodes    = parse(Int, match(r"numnodes=(\d+)", name).captures[1])
    probability = parse(Float64, match(r"probability=([\d.]+)", name).captures[1])
    seedstart   = parse(Int, match(r"seedstart=(\d+)", name).captures[1])

    # Return a named tuple
    return (
        backend     = backend,
        graphs      = graphs,
        graphtype   = graphtype,
        neighbors   = neighbors,
        numnodes    = numnodes,
        probability = probability,
        seedstart   = seedstart,
    )
end

function sweep_parameters(in_filename, N_gridpoints=100, use_gpu=false, streamed_version=false, batch_size=1_000, graph_type="Any")

    sampled_homodist = npzread(in_filename)

    num_constraints = size(sampled_homodist, 1) - 1

    in_filename_no_ext = splitext(basename(in_filename))[1]
    keyval_dict = filename_dict(in_filename_no_ext)

    num_nodes = parse(Int, keyval_dict["numnodes"]) 
    probability = parse(Float64, keyval_dict["probability"])
    filename_graphtype = keyval_dict["graphtype"]

    if (graph_type != "Any") && (graph_type != filename_graphtype)
        println("Graph type $filename_graphtype does not match required graph type $graph_type. Skipping.")
        return
    end

    @assert ndims(sampled_homodist) == 3 "Homogeneous distribution must be 3D."
    @assert size(sampled_homodist, 1) == size(sampled_homodist,3) "1st and 3rd dimensions of homogeneous distribution must be same length."
    @assert size(sampled_homodist, 2) == num_nodes + 1 "num_nodes from filename does not agree with number of hamming distances in the homogeneous distribution."


    # Ideas for scale to make gridpoints more sensible:
    # - l/right_angle can be changed to an angle determining the slope, so range is 0:pi
    # - ceneter_adjustment_add_range can be done as a fraction of the total number of possible costs, so range is -0.5:0.5
    # - height_adjustment is the hardest

    height_adjustment_range = LinRange(0, 20, N_gridpoints)
    ceneter_adjustment_add_range = LinRange(0, 1, N_gridpoints)
    left_angle_range = LinRange(0, 0.5, N_gridpoints)
    right_angle_range = LinRange(0, 0.5, N_gridpoints)

    rangestr(r) = "$(minimum(r)):$(maximum(r))"
    out_filename_base = in_filename_no_ext * "_htweaksub=$(rangestr(height_adjustment_range))_hctweakadd=$(rangestr(ceneter_adjustment_add_range))_ltweakmul=$(rangestr(left_angle_range))_rtweakmul=$(rangestr(right_angle_range))_ngridpoints=$(N_gridpoints)"
    out_filename_npy = "results/" * "parametersweep_" * out_filename_base * ".npy"
    out_filename_csv = "results/" * "parametersweep_" * out_filename_base * ".csv"
    out_filename_opt_csv = "results/" * "optimalparams_" * out_filename_base * ".csv"
    mkpath("results")

    params = Iterators.product(height_adjustment_range, ceneter_adjustment_add_range, left_angle_range, right_angle_range)


    # This may be unnecessary, just make big file but handle writes in chunks
    if streamed_version
        mean_mse = 0.0
        min_mse = Inf
        min_mse_params = (NaN, NaN, NaN, NaN)

        for params_chunk in Iterators.partition(params, 1_000)
            proxies = [JuliaQAOA.IntuitiveTriangleProxy(num_constraints, num_nodes, params...) for params in params_chunk]

            if use_gpu
                mses_vec = JuliaQAOA.gpu_multi_proxy_mse(proxies, sampled_homodist, batch_size=batch_size) |> Array |> vec
            else
                mses_vec = JuliaQAOA.cpu_multi_proxy_mse(proxies, sampled_homodist, batch_size=batch_size)
            end

            mean_mse += sum(x -> isnan(x) ? 0 : x, mses_vec) # Technically, should also reduce number I divide mean_MSE by later on
            local_min_mse, min_mse_index = findmin(mses_vec)
            local_min_mse, min_mse_index = findmin(x -> isnan(x) ? Inf : x, mses_vec)
            if local_min_mse < min_mse
                min_mse = local_min_mse
                min_mse_params = params_chunk[min_mse_index]
            end
        end
        mean_mse /= N_gridpoints^4
        println("Min MSE = $min_mse\nMean MSE = $mean_mse\nOptimal parameters = $min_mse_params")

        optimal_params_mat = hcat(min_mse_params..., min_mse, mean_mse)
        writedlm(out_filename_opt_csv, optimal_params_mat)
    else
        params_vec = params |> collect |> vec

        proxies = [JuliaQAOA.IntuitiveTriangleProxy(num_constraints, num_nodes, params...) for params in params_vec]

        # End supression of warnings

        #proxies_mat = [
        #    TriangleProxy(num_constraints, num_nodes, height_adjustment, ceneter_adjustment_add, left_angle, right_angle)
        #    for height_adjustment  in height_adjustment_range,
        #        ceneter_adjustment_add in ceneter_adjustment_add_range,
        #        left_angle  in left_angle_range,
        #        right_angle  in right_angle_range
        #]
        
        if use_gpu
            mses_vec = JuliaQAOA.gpu_multi_proxy_mse(proxies, sampled_homodist, batch_size=batch_size) |> Array |> vec
        else
            mses_vec = JuliaQAOA.cpu_multi_proxy_mse(proxies, sampled_homodist, batch_size=batch_size)
        end

        mean_mse = sum(x -> isnan(x) ? 0 : x, mses_vec) / length(mses_vec)
        min_mse, min_mse_index = findmin(x -> isnan(x) ? Inf : x, mses_vec)
        min_mse_params = params_vec[min_mse_index]
        println("Min MSE = $min_mse\nMean MSE = $mean_mse\nOptimal parameters = $min_mse_params")


        mses_4D_array = reshape(mses_vec, N_gridpoints, N_gridpoints, N_gridpoints, N_gridpoints)
        params_mat = [params[i] for params in params_vec, i in 1:4]
        mses_csv_mat = hcat(params_mat, mses_vec)

        # Save results
        npzwrite(out_filename_npy, mses_4D_array)
        # output will be csv, fist four columns are paremeters, last column is mse.
        writedlm(out_filename_csv, mses_csv_mat, ',')

        optimal_params_mat = hcat(min_mse_params..., min_mse, mean_mse)
        writedlm(out_filename_opt_csv, optimal_params_mat)
    end

    return optimal_params_mat
end

function collect_parameter_sweeps(directory, N_gridpoints=3, use_gpu=false, streamed_version=false, batch_size=1_000, graph_type="Any")
    files = readdir(directory, join=true)

    header = hcat("backend", "graphs", "graphtype", "neighbors", "numnodes", "probability", "seedstart", "p1", "p2", "p3", "p4", "minmse", "meanmse")
    header = hcat(
        @sprintf("%-10s", "backend"),
        @sprintf("%-6s",  "graphs"),
        @sprintf("%-15s", "graphtype"),
        @sprintf("%-9s",  "neighbors"),
        @sprintf("%-9s",  "numnodes"),
        @sprintf("%-12s", "probability"),
        @sprintf("%-9s",  "seedstart"),
        @sprintf("%-17s", "param1"),
        @sprintf("%-17s", "param2"),
        @sprintf("%-17s", "param3"),
        @sprintf("%-17s", "param4"),
        @sprintf("%-17s", "minMSE"),
        @sprintf("%-17s", "meanMSE"),
    )
    full_data_mat = Matrix{Any}(undef, 0, length(header))
    full_data_mat = vcat(full_data_mat, header)
    for numpy_file in filter(x -> isfile(x) && endswith(x, ".npy"), files)
        println("Performing $N_gridpoints-point parameter sweep for file ", numpy_file, ".")
        println("Starting at time ", now())
        flush(stdout)
        start_sec = time()

        optimal_params_mat = sweep_parameters(numpy_file, N_gridpoints, use_gpu, streamed_version, batch_size, graph_type)

        println("Finished at time ", now(), ".")
        end_sec = time()


        data_params = parse_filename_regex(numpy_file)
        formatted_data_params = [
            @sprintf("%-10s",   data_params.backend),
            @sprintf("%-6d",    data_params.graphs),
            @sprintf("%-15s",   data_params.graphtype),
            @sprintf("%-9d",    data_params.neighbors),
            @sprintf("%-9d",    data_params.numnodes),
            @sprintf("%-12.3f", data_params.probability),
            @sprintf("%-9d",    data_params.seedstart),
        ]
        formatted_row = hcat(
            @sprintf("%-10s",   data_params.backend),
            @sprintf("%-6d",    data_params.graphs),
            @sprintf("%-15s",   data_params.graphtype),
            @sprintf("%-9d",    data_params.neighbors),
            @sprintf("%-9d",    data_params.numnodes),
            @sprintf("%-12.3f", data_params.probability),
            @sprintf("%-9d",    data_params.seedstart),
            @sprintf("%-17.17g", optimal_params_mat[1]),
            @sprintf("%-17.17g", optimal_params_mat[2]),
            @sprintf("%-17.17g", optimal_params_mat[3]),
            @sprintf("%-17.17g", optimal_params_mat[4]),
            @sprintf("%-17.17g", optimal_params_mat[5]),
            @sprintf("%-17.17g", optimal_params_mat[6]),
        )
        row = hcat(formatted_data_params..., optimal_params_mat...)
        #full_data_mat = vcat(full_data_mat, row)
        full_data_mat = vcat(full_data_mat, formatted_row)

        println("Parameter sweep took $(end_sec - start_sec) seconds.\n", "-"^20, "\n\n")
        flush(stdout)
    end

    writedlm(stdout, full_data_mat, '\t')
    open("results/mega_optimal_params.csv", "a") do fileio
        writedlm(fileio, full_data_mat, '\t')
    end
end

if abspath(PROGRAM_FILE) == @__FILE__ # Only run this if file is being run from the command line
    s = ArgParseSettings()
    @add_arg_table s begin
        "directory"
            help = "Directory of input npy files containing sampled homogeneous distributions."
            required = true
        "ngridpoints"
            help = "Number of gridpoints to use along each axis."
            arg_type = Int
            required = true
        "--gpu", "-g"
            help = "Flag for using GPU"
            action = :store_true
        "--stream", "-s"
            help = "Flag to use 'streamed' version, which will find optimal parameters but not store the results of the entire parameter sweep."
            action = :store_true
        "--batch_size", "-b"
            help = "Number of proxies to do at once on GPU"
            arg_type = Int
            default = 1_000
        "--graph_type", "-t"
            help = "Only process graphs of a certain type."
            arg_type = String
            default = "Any"
    end

    parsed_args = parse_args(ARGS, s)
    collect_parameter_sweeps(parsed_args["directory"], parsed_args["ngridpoints"], parsed_args["gpu"], parsed_args["stream"], parsed_args["batch_size"], parsed_args["graph_type"])
end



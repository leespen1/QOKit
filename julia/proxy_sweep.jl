using NPZ: npzwrite, npzread
using Dates: now
using DelimitedFiles, JuliaQAOA, Logging, Dates, ArgParse


function sweep_parameters(in_filename, N_gridpoints=100, use_gpu=false, streamed_version=false)

    sampled_homodist = npzread(in_filename)

    num_constraints = size(sampled_homodist, 1) - 1

    in_filename_no_ext = split(basename(in_filename),'.')[1]
    keyval_pair_strs = split(in_filename_no_ext, '_')
    keyval_dict = Dict((kv[1], kv[2]) for kv in split.(keyval_pair_strs, '='))

    num_nodes = parse(Int, keyval_dict["numnodes"]) 
    probability = parse(Float64, keyval_dict["probability"])
    graphtype = keyval_dict["graphtype"]

    @assert ndims(sampled_homodist) == 3 "Homogeneous distribution must be 3D."
    @assert size(sampled_homodist, 1) == size(sampled_homodist,3) "1st and 3rd dimensions of homogeneous distribution must be same length."
    @assert size(sampled_homodist, 2) == num_nodes + 1 "num_nodes from filename does not agree with number of hamming distances in the homogeneous distribution."

    bounds = [
        (0, (num_nodes^2) / 3),   # h_tweak_sub >= 0
        (-10, 10),   # hc_tweak_add can be small positive/negative
        (0.005, 2),   # l_tweak_mul > 0
        (0.05, 2),   # r_tweak_mul > 0
    ]

    # Ideas for scale to make gridpoints more sensible:
    # - l/r_tweak_mul can be changed to an angle determining the slope, so range is 0:pi
    # - hc_tweak_add_range can be done as a fraction of the total number of possible costs, so range is -0.5:0.5
    # - h_tweak_sub is the hardest

    h_tweak_sub_range = LinRange(0, num_nodes^2 / 3, N_gridpoints)
    hc_tweak_add_range = LinRange(-10, 10, N_gridpoints)
    l_tweak_mul_range = LinRange(0.005, 2, N_gridpoints)
    r_tweak_mul_range = LinRange(0.005, 2, N_gridpoints)

    rangestr(r) = "$(minimum(r)):$(maximum(r))"
    out_filename_base = in_filename_no_ext * "_htweaksub=$(rangestr(h_tweak_sub_range))_hctweakadd=$(rangestr(hc_tweak_add_range))_ltweakmul=$(rangestr(l_tweak_mul_range))_rtweakmul=$(rangestr(r_tweak_mul_range))_ngridpoints=$(N_gridpoints)"
    out_filename_npy = "results/" * "parametersweep_" * out_filename_base * ".npy"
    out_filename_csv = "results/" * "parametersweep_" * out_filename_base * ".csv"
    out_filename_opt_csv = "results/" * "optimalparams_" * out_filename_base * ".csv"

    params = Iterators.product(h_tweak_sub_range, hc_tweak_add_range, l_tweak_mul_range, r_tweak_mul_range)


    # This may be unnecessary, just make big file but handle writes in chunks
    if streamed_version
        mean_mse = 0.0
        min_mse = Inf
        min_mse_params = (NaN, NaN, NaN, NaN)

        for params_chunk in Iterators.partition(params, 10_000)
            old_logger = global_logger()
            global_logger(ConsoleLogger(stderr, Logging.Error))
            proxies = [TriangleProxy(num_constraints, num_nodes, params...) for params in params_chunk]
            global_logger(old_logger)

            if use_gpu
                mses_vec = JuliaQAOA.gpu_multi_proxy_mse(proxies, sampled_homodist) |> Array |> vec
            else
                mses_vec = JuliaQAOA.cpu_multi_proxy_mse(proxies, sampled_homodist, batch_size=1000)
            end

            mean_mse += sum(mses_vec)
            local_min_mse, min_mse_index = findmin(mses_vec)
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

        # Temporarily supress warnings (which happen when h_tweak_sub makes the peak negative)
        old_logger = global_logger()
        global_logger(ConsoleLogger(stderr, Logging.Error))

        proxies = [TriangleProxy(num_constraints, num_nodes, params...) for params in params_vec]

        global_logger(old_logger)
        # End supression of warnings

        #proxies_mat = [
        #    TriangleProxy(num_constraints, num_nodes, h_tweak_sub, hc_tweak_add, l_tweak_mul, r_tweak_mul)
        #    for h_tweak_sub  in h_tweak_sub_range,
        #        hc_tweak_add in hc_tweak_add_range,
        #        l_tweak_mul  in l_tweak_mul_range,
        #        r_tweak_mul  in r_tweak_mul_range
        #]
        
        if use_gpu
            mses_vec = JuliaQAOA.gpu_multi_proxy_mse(proxies, sampled_homodist) |> Array |> vec
        else
            mses_vec = JuliaQAOA.cpu_multi_proxy_mse(proxies, sampled_homodist, batch_size=1000)
        end

        mean_mse = sum(mses_vec) / length(mses_vec)
        min_mse, min_mse_index = findmin(mses_vec)
        min_mse_params = params_vec[min_mse_index]
        println("Min MSE = $min_mse\nMean MSE = $mean_mse\nOptimal parameters = $min_mse_params")


        mses_4D_array = reshape(mses_vec, N_gridpoints, N_gridpoints, N_gridpoints, N_gridpoints)
        params_mat = [params[i] for params in params_vec, i in 1:4]
        mses_csv_mat = hcat(params_mat, mses_vec)

        # Save results
        mkpath("results")
        npzwrite(out_filename_npy, mses_4D_array)
        # output will be csv, fist four columns are paremeters, last column is mse.
        writedlm(out_filename_csv, mses_csv_mat, ',')

        optimal_params_mat = hcat(min_mse_params..., min_mse, mean_mse)
        writedlm(out_filename_opt_csv, optimal_params_mat)
    end

end

function collect_parameter_sweeps(directory, N_gridpoints=3, use_gpu=false, streamed_version=false)
    files = readdir(directory, join=true)
    for numpy_file in filter(x -> isfile(x) && endswith(x, ".npy"), files)
        println("Performing parameter sweep for file ", numpy_file, ".")
        println("Starting at time ", now())
        start_sec = time()

        sweep_parameters(numpy_file, N_gridpoints, use_gpu, streamed_version)

        println("Finished at time ", now(), ".")
        end_sec = time()
        println("Parameter sweep took $(end_sec - start_sec) seconds.\n", "-"^20, "\n\n")
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
    end

    parsed_args = parse_args(ARGS, s)
    collect_parameter_sweeps(parsed_args["directory"], parsed_args["ngridpoints"], parsed_args["gpu"], parsed_args["stream"])
end



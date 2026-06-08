using DrWatson
@quickactivate "JuliaQAOA"
using DelimitedFiles, JuliaQAOA
using Printf: @sprintf
using NPZ: npzread

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

function main(directory)
    files = readdir(directory, join=true)

    header = hcat(
        @sprintf("%-17s", "NormalizedMSE"),
        @sprintf("%-15s", "graphtype"),
        @sprintf("%-9s",  "numnodes"),
        @sprintf("%-12s", "probability"),
        @sprintf("%-9s",  "neighbors"),
        @sprintf("%-17s", "SampledDistSum"),
        @sprintf("%-17s", "ProxyDistSum"),
        @sprintf("%-17s", "MeanSquaredError"),
        @sprintf("%-17s", "SumSquaredError"),
        @sprintf("%-17s", "NormalizedSSE"),
        @sprintf("%-6s",  "graphs"),
        @sprintf("%-9s",  "seedstart"),
        @sprintf("%-10s", "backend"),
    )
    full_data_mat = header
    for numpy_file in filter(x -> isfile(x) && endswith(x, ".npy"), files)
        data_params = parse_filename_regex(numpy_file)
        sampled_homodist = npzread(numpy_file)
        num_constraints = size(sampled_homodist, 1) - 1

        # Question how to get the 
        paper_proxy = PaperProxy(
            num_constraints,
            data_params.numnodes,
            data_params.probability
        )

        n_entries = length(sampled_homodist)

        proxy_homodist = JuliaQAOA.cpu_compute_homodist(paper_proxy)

        sampled_homodist_sum = sum(sampled_homodist)
        proxy_homodist_sum = sum(proxy_homodist)

        mean_sampled_value = sampled_homodist_sum / n_entries
        mean_proxy_value = proxy_homodist_sum / n_entries
         
        sum_squared_error = mapreduce((x, y) -> (x-y)^2, +, sampled_homodist, proxy_homodist)
        mean_squared_error = sum_squared_error / n_entries

        normal_sampled_homodist = sampled_homodist ./ sampled_homodist_sum
        normal_proxy_homodist = proxy_homodist ./ proxy_homodist_sum

        n_sum_squared_error = mapreduce((x, y) -> (x-y)^2, +, normal_sampled_homodist, normal_proxy_homodist)
        n_mean_squared_error = n_sum_squared_error / n_entries
        
        row = hcat(
            @sprintf("%-17.17g", n_mean_squared_error),
            @sprintf("%-15s", data_params.graphtype),
            @sprintf("%-9d", data_params.numnodes),
            @sprintf("%-12.3f", data_params.probability),
            @sprintf("%-9d", data_params.neighbors),
            @sprintf("%-17.17g", sampled_homodist_sum),
            @sprintf("%-17.17g", proxy_homodist_sum),
            @sprintf("%-17.17g", mean_squared_error),
            @sprintf("%-17.17g", n_sum_squared_error),
            @sprintf("%-17.17g", sum_squared_error),
            @sprintf("%-6d", data_params.graphs),
            @sprintf("%-9d", data_params.seedstart),
            @sprintf("%-10s", data_params.backend),
        )

        full_data_mat = vcat(full_data_mat, row)
    end
    writedlm(stdout, full_data_mat, '\t')
    open("results/paperproxy_mse_data.csv", "a") do fileio
        writedlm(fileio, full_data_mat, '\t')
    end

    return full_data_mat
end

if abspath(PROGRAM_FILE) == @__FILE__ # Only run this if file is being run from the command line
    @assert length(ARGS) == 1 "Must provide exactly one argument: the directory from which to read the sampled homogeneous distribution files."
    main(ARGS[1])
end

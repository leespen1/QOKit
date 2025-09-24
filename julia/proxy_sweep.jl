using NPZ: npzwrite

function sweep_parameters(in_filename)

    in_filename_no_ext = split(in_filename,'.')[1]
    
    keyval_pair_strs = split(in_filename_no_ext, '_')
    keyval_dict = Dict((kv[1], kv[2]) for kv in split(keyval_pair_strs, '='))
    num_nodes = parse(Int, keyval_dict["numnodes"]) 

    bounds = [
        (0, (num_nodes^2) / 3),   # h_tweak_sub >= 0
        (-10, 10),   # hc_tweak_add can be small positive/negative
        (0.005, 2),   # l_tweak_mul > 0
        (0.05, 2),   # r_tweak_mul > 0
    ]

    # Ideas for scale to make gridpoints easier:
    # - l/r_tweak_mul can be changed to an angle determining the slope
    # - hc_tweak_add_range can be changed so that the 

    N_gridpoints = 101

    h_tweak_sub_range = LinRange(0, num_nodes^2 / 3, N_gridpoints)
    hc_tweak_add_range = LinRange(-10, 10, N_gridpoints)
    l_tweak_mul_range = LinRange(0.005, 2, N_gridpoints)
    r_tweak_mul_range = LinRange(0.005, 2, N_gridpoints)
    proxies = [
        TriangleProxy(num_constraints, num_qubits, h_tweak_sub, hc_tweak_add, l_tweak_mul, r_tweak_mul)
        for h_tweak_sub  in h_tweak_sub_range
            hc_tweak_add in hc_tweak_add_range
            l_tweak_mul  in l_tweak_mul_range
            r_tweak_mul  in r_tweak_mul_range
    ]

    rangestr(r) = "$(minimum(r)):$(maximum(r))"
    out_filename =  in_filename_no_ext * "_htweaksub=$(rangestr(h_tweak_sub_range))_hctweakadd=$(rangestr(hc_tweak_add_range))_ltweakmul=$(rangestr(l_tweak_mul_range))_rtweakmul=$(rangestr(r_tweak_mul_range)).npy"
    println(filename)

end

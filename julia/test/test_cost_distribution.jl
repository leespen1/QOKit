using JuliaQAOA, Test, CUDA
import Random: MersenneTwister

@testset showtiming=true "hamming_distance" begin
    @testset "Basic examples" begin
        # Same bitstrings have distance 0
        @test hamming_distance(0b000, 0b000) == 0
        @test hamming_distance(0b111, 0b111) == 0

        # Complementary bitstrings
        @test hamming_distance(0b000, 0b111) == 3
        @test hamming_distance(0b0000, 0b1111) == 4

        # Single bit differences
        @test hamming_distance(0b100, 0b000) == 1
        @test hamming_distance(0b010, 0b000) == 1
        @test hamming_distance(0b001, 0b000) == 1

        # Two bit differences
        @test hamming_distance(0b101, 0b000) == 2
        @test hamming_distance(0b110, 0b000) == 2
        @test hamming_distance(0b011, 0b000) == 2

        # Mixed examples
        @test hamming_distance(0b010, 0b101) == 3
        @test hamming_distance(0b1010, 0b0101) == 4
        @test hamming_distance(0b1100, 0b0011) == 4
    end

    @testset "Symmetry" begin
        for _ in 1:10
            a = rand(MersenneTwister(0), 0:1023)
            b = rand(MersenneTwister(1), 0:1023)
            @test hamming_distance(a, b) == hamming_distance(b, a)
        end
    end

    @testset "Triangle inequality" begin
        for _ in 1:10
            a = rand(MersenneTwister(0), 0:255)
            b = rand(MersenneTwister(1), 0:255)
            c = rand(MersenneTwister(2), 0:255)
            @test hamming_distance(a, c) <= hamming_distance(a, b) + hamming_distance(b, c)
        end
    end
end


@testset showtiming=true "get_real_distribution_from_costs" begin
    @testset "Simple 2-vertex example" begin
        # 2 vertices, 1 edge (complete graph K2)
        # Bitstrings: 00, 01, 10, 11
        # For MaxCut on K2: costs are [0, 1, 1, 0] (00 and 11 have cost 0, 01 and 10 have cost 1)
        costs = Float64[0, 1, 1, 0]
        num_edges = 1
        num_vertices = 2

        n_dist = get_real_distribution_from_costs(costs, num_edges, num_vertices)

        # Shape should be (4 bitstrings, 3 distances, 2 costs)
        @test size(n_dist) == (4, 3, 2)

        # For bitstring 00 (index 1, cost 0):
        # - distance 0: only 00 itself, cost 0 → n[1,1,1] = 1
        # - distance 1: 01 (cost 1) and 10 (cost 1) → n[1,2,2] = 2
        # - distance 2: only 11, cost 0 → n[1,3,1] = 1
        @test n_dist[1, 1, 1] == 1  # d=0, c=0
        @test n_dist[1, 1, 2] == 0  # d=0, c=1
        @test n_dist[1, 2, 1] == 0  # d=1, c=0
        @test n_dist[1, 2, 2] == 2  # d=1, c=1
        @test n_dist[1, 3, 1] == 1  # d=2, c=0
        @test n_dist[1, 3, 2] == 0  # d=2, c=1

        # For bitstring 01 (index 2, cost 1):
        # - distance 0: only 01 itself, cost 1 → n[2,1,2] = 1
        # - distance 1: 00 (cost 0) and 11 (cost 0) → n[2,2,1] = 2
        # - distance 2: only 10, cost 1 → n[2,3,2] = 1
        @test n_dist[2, 1, 1] == 0  # d=0, c=0
        @test n_dist[2, 1, 2] == 1  # d=0, c=1
        @test n_dist[2, 2, 1] == 2  # d=1, c=0
        @test n_dist[2, 2, 2] == 0  # d=1, c=1
        @test n_dist[2, 3, 1] == 0  # d=2, c=0
        @test n_dist[2, 3, 2] == 1  # d=2, c=1

        # For bitstring 10 (index 3, cost 1):
        # - distance 0: only 10 itself, cost 1 → n[3,1,2] = 1
        # - distance 1: 00 (cost 0) and 11 (cost 0) → n[3,2,1] = 2
        # - distance 2: only 01, cost 1 → n[3,3,2] = 1
        @test n_dist[3, 1, 1] == 0  # d=0, c=0
        @test n_dist[3, 1, 2] == 1  # d=0, c=1
        @test n_dist[3, 2, 1] == 2  # d=1, c=0
        @test n_dist[3, 2, 2] == 0  # d=1, c=1
        @test n_dist[3, 3, 1] == 0  # d=2, c=0
        @test n_dist[3, 3, 2] == 1  # d=2, c=1

        # For bitstring 11 (index 4, cost 0):
        # - distance 0: only 11 itself, cost 0 → n[4,1,1] = 1
        # - distance 1: 01 (cost 1) and 10 (cost 1) → n[4,2,2] = 2
        # - distance 2: only 00, cost 0 → n[4,3,1] = 1
        @test n_dist[4, 1, 1] == 1  # d=0, c=0
        @test n_dist[4, 1, 2] == 0  # d=0, c=1
        @test n_dist[4, 2, 1] == 0  # d=1, c=0
        @test n_dist[4, 2, 2] == 2  # d=1, c=1
        @test n_dist[4, 3, 1] == 1  # d=2, c=0
        @test n_dist[4, 3, 2] == 0  # d=2, c=1

        # Total count for any bitstring should be 2^n = 4
        for x in 1:4
            @test sum(n_dist[x, :, :]) == 4
        end
    end

    @testset "Row sums equal binomial coefficients" begin
        # For any bitstring x, the number of bitstrings at distance d is C(n, d)
        num_vertices = 4
        num_edges = 3  # arbitrary
        costs = rand(MersenneTwister(42), 0:num_edges, 2^num_vertices) .|> Float64

        n_dist = get_real_distribution_from_costs(costs, num_edges, num_vertices)

        for x in 1:2^num_vertices
            for d in 0:num_vertices
                # Sum over all costs at distance d should equal C(n, d)
                @test sum(n_dist[x, d+1, :]) == binomial(num_vertices, d)
            end
        end
    end

    @testset "max_num_edges padding" begin
        costs = Float64[0, 1, 1, 0]
        num_edges = 1
        num_vertices = 2

        n_dist = get_real_distribution_from_costs(costs, num_edges, num_vertices; max_num_edges=5)

        # Should have 6 costs (0 to 5) instead of 2
        @test size(n_dist) == (4, 3, 6)

        # Extra columns should be zeros
        @test all(n_dist[:, :, 3:6] .== 0)
    end
end


@testset showtiming=true "get_homogeneous_distribution_from_costs" begin
    @testset "Simple 2-vertex example" begin
        costs = Float64[0, 1, 1, 0]
        num_edges = 1
        num_vertices = 2

        n_dist = get_real_distribution_from_costs(costs, num_edges, num_vertices)
        N_dist = get_homogeneous_distribution_from_costs(costs, n_dist)

        # Shape should be (2 costs, 3 distances, 2 costs)
        @test size(N_dist) == (2, 3, 2)

        # For c'=0 (bitstrings 00 and 11):
        # Average of n(00; d, c) and n(11; d, c)
        # 00: at d=0 has itself (c=0), at d=1 has 01,10 (c=1), at d=2 has 11 (c=0)
        # 11: at d=0 has itself (c=0), at d=1 has 01,10 (c=1), at d=2 has 00 (c=0)
        # They should be the same by symmetry
        @test N_dist[1, 1, 1] == 1.0  # d=0, c=0: average of [1, 1] = 1
        @test N_dist[1, 2, 2] == 2.0  # d=1, c=1: average of [2, 2] = 2
        @test N_dist[1, 3, 1] == 1.0  # d=2, c=0: average of [1, 1] = 1
    end

    @testset "Consistency: sum over c equals binomial" begin
        num_vertices = 4
        num_edges = 4
        costs = rand(MersenneTwister(123), 0:num_edges, 2^num_vertices) .|> Float64

        n_dist = get_real_distribution_from_costs(costs, num_edges, num_vertices)
        N_dist = get_homogeneous_distribution_from_costs(costs, n_dist)

        # For each c', sum over c at distance d should still equal C(n, d)
        for c_prime in 0:num_edges
            # Check if any bitstrings have this cost
            if any(costs .== c_prime)
                for d in 0:num_vertices
                    @test sum(N_dist[c_prime+1, d+1, :]) ≈ binomial(num_vertices, d) atol=1e-10
                end
            end
        end
    end
end


@testset showtiming=true "get_homogeneous_distribution_from_costs_direct" begin
    @testset "Matches two-step computation" begin
        num_vertices = 5
        num_edges = 6
        costs = rand(MersenneTwister(456), 0:num_edges, 2^num_vertices) .|> Float64

        # Two-step method
        n_dist = get_real_distribution_from_costs(costs, num_edges, num_vertices)
        N_dist_twostep = get_homogeneous_distribution_from_costs(costs, n_dist)

        # Direct method
        N_dist_direct = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)

        @test N_dist_twostep ≈ N_dist_direct atol=1e-10
    end

    @testset "Matches two-step with padding" begin
        num_vertices = 4
        num_edges = 3
        max_num_edges = 8
        costs = rand(MersenneTwister(789), 0:num_edges, 2^num_vertices) .|> Float64

        n_dist = get_real_distribution_from_costs(costs, num_edges, num_vertices; max_num_edges=max_num_edges)
        N_dist_twostep = get_homogeneous_distribution_from_costs(costs, n_dist; max_num_edges=max_num_edges)
        N_dist_direct = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices; max_num_edges=max_num_edges)

        @test size(N_dist_twostep) == size(N_dist_direct)
        @test N_dist_twostep ≈ N_dist_direct atol=1e-10
    end
end


@testset showtiming=true "GPU real distribution functions" begin
    if CUDA.has_cuda_gpu()
        @testset "gpu_get_real_distribution_from_costs matches CPU" begin
            num_vertices = 6
            num_edges = 8
            costs = rand(MersenneTwister(111), 0:num_edges, 2^num_vertices) .|> Float64

            cpu_result = get_real_distribution_from_costs(costs, num_edges, num_vertices)
            gpu_result = gpu_get_real_distribution_from_costs(costs, num_edges, num_vertices) |> Array

            @test size(cpu_result) == size(gpu_result)
            @test cpu_result ≈ gpu_result atol=1e-10
        end

        @testset "gpu_get_homogeneous_distribution_from_costs_direct matches CPU" begin
            num_vertices = 6
            num_edges = 8
            costs = rand(MersenneTwister(222), 0:num_edges, 2^num_vertices) .|> Float64

            cpu_result = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)
            gpu_result = gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices) |> Array

            @test size(cpu_result) == size(gpu_result)
            @test cpu_result ≈ gpu_result atol=1e-10
        end

        @testset "GPU with max_num_edges padding" begin
            num_vertices = 5
            num_edges = 4
            max_num_edges = 10
            costs = rand(MersenneTwister(333), 0:num_edges, 2^num_vertices) .|> Float64

            cpu_result = get_homogeneous_distribution_from_costs_direct(
                costs, num_edges, num_vertices; max_num_edges=max_num_edges
            )
            gpu_result = gpu_get_homogeneous_distribution_from_costs_direct(
                costs, num_edges, num_vertices; max_num_edges=max_num_edges
            ) |> Array

            @test size(cpu_result) == size(gpu_result)
            @test cpu_result ≈ gpu_result atol=1e-10
        end
    else
        @warn "Skipping GPU tests because no GPU detected"
    end
end


@testset showtiming=true "Utility functions" begin
    @testset "pad_to_shape" begin
        arr = [1 2; 3 4; 5 6]  # 3x2
        padded = pad_to_shape(arr, (5, 4))

        @test size(padded) == (5, 4)
        @test padded[1:3, 1:2] == arr
        @test all(padded[4:5, :] .== 0)
        @test all(padded[:, 3:4] .== 0)
    end

    @testset "pad_to_match" begin
        a = ones(2, 3, 4)
        b = ones(3, 4, 5)

        a_padded, b_result = pad_to_match(a, b)

        @test size(a_padded) == (3, 4, 5)
        @test b_result === b  # b unchanged
        @test a_padded[1:2, 1:3, 1:4] == a

        # Test reverse order: when b is first arg and a is second,
        # a (the smaller one) gets padded
        b_result2, a_padded2 = pad_to_match(b, a)
        @test b_result2 === b  # b unchanged since it's larger
        @test size(a_padded2) == size(b)  # a was padded to match b
        @test a_padded2[1:2, 1:3, 1:4] == a
    end

    @testset "pad_and_stack" begin
        arr1 = ones(2, 3, 4)
        arr2 = 2 .* ones(3, 4, 5)
        arr3 = 3 .* ones(2, 2, 3)

        stacked = pad_and_stack([arr1, arr2, arr3])

        # Result shape: (num_arrays, max_dim1, max_dim2, max_dim3)
        @test size(stacked) == (3, 3, 4, 5)
        # Check that original values are preserved in the non-padded region
        @test all(stacked[1, 1:2, 1:3, 1:4] .== 1.0)
        @test all(stacked[2, 1:3, 1:4, 1:5] .== 2.0)
        @test all(stacked[3, 1:2, 1:2, 1:3] .== 3.0)
        # Check that padded regions are zeros
        @test all(stacked[1, 3:3, :, :] .== 0.0)
        @test all(stacked[3, :, 3:4, :] .== 0.0)
    end

    @testset "average_distributions" begin
        arr1 = ones(2, 3, 4) .* 2.0
        arr2 = ones(2, 3, 4) .* 4.0

        avg = average_distributions([arr1, arr2])

        @test size(avg) == (2, 3, 4)
        @test all(avg .≈ 3.0)
    end

    @testset "stddev_distributions" begin
        arr1 = ones(2, 3, 4) .* 2.0
        arr2 = ones(2, 3, 4) .* 4.0

        sd = stddev_distributions([arr1, arr2])

        @test size(sd) == (2, 3, 4)
        # std of [2, 4] with N normalization = 1.0
        @test all(sd .≈ 1.0)
    end

    @testset "distributions_mean_and_stddev" begin
        arr1 = ones(2, 3, 4) .* 2.0
        arr2 = ones(2, 3, 4) .* 4.0

        m, s = distributions_mean_and_stddev([arr1, arr2])

        @test all(m .≈ 3.0)
        @test all(s .≈ 1.0)
    end

    @testset "distribution_array_to_dict" begin
        arr = zeros(3, 2, 2)
        arr[1, 1, 1] = 5.0
        arr[2, 2, 1] = 3.0
        arr[3, 1, 2] = 7.0

        d = distribution_array_to_dict(arr)

        @test length(d) == 3
        # Keys are 0-indexed for Python compatibility
        @test d[(0, 0, 0)] == 5.0
        @test d[(1, 1, 0)] == 3.0
        @test d[(2, 0, 1)] == 7.0
    end

    @testset "get_pearson_correlation_coefficients" begin
        # Identical distributions should have correlation 1
        arr1 = rand(MersenneTwister(0), 3, 4, 5)
        correlations = get_pearson_correlation_coefficients(arr1, arr1)

        @test length(correlations) == 3
        @test all(correlations .≈ 1.0)

        # Negatively correlated
        arr2 = -arr1
        correlations_neg = get_pearson_correlation_coefficients(arr1, arr2)
        @test all(correlations_neg .≈ -1.0)
    end
end


#==============================================================================#
#              Additional N(c'; d, c) Tests                                     #
#==============================================================================#

@testset showtiming=true "N(c'; d, c) properties" begin
    @testset "Uniform costs - all bitstrings same cost" begin
        # When all bitstrings have the same cost, N(c'; d, c) should have a specific structure
        num_vertices = 4
        num_edges = 5
        uniform_cost = 3
        costs = fill(Float64(uniform_cost), 2^num_vertices)

        N_dist = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)

        # Only one c' value should have non-zero entries (the uniform cost)
        for c_prime in 0:num_edges
            if c_prime != uniform_cost
                @test all(N_dist[c_prime + 1, :, :] .== 0.0)
            end
        end

        # For the uniform cost c', all neighbors at any distance also have that cost
        # So N(uniform_cost; d, c) should be binomial(n, d) when c == uniform_cost, else 0
        for d in 0:num_vertices
            @test N_dist[uniform_cost + 1, d + 1, uniform_cost + 1] == binomial(num_vertices, d)
            for c in 0:num_edges
                if c != uniform_cost
                    @test N_dist[uniform_cost + 1, d + 1, c + 1] == 0.0
                end
            end
        end
    end

    @testset "Symmetry: N is symmetric in exchange of bitstring and complement" begin
        # For MaxCut-like costs where c(x) = num_edges - c(~x), the distribution
        # should exhibit certain symmetries. Here we test a simpler property:
        # Sum over all c' of N(c'; d, c) weighted by count(c') should equal
        # binomial(n, d) * 2^n / (num_edges + 1) approximately for well-distributed costs
        num_vertices = 5
        num_edges = 6
        costs = rand(MersenneTwister(555), 0:num_edges, 2^num_vertices) .|> Float64

        N_dist = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)

        # For each d, the weighted average of N(c'; d, c) over c' should sum to binomial(n, d)
        cost_counts = [count(==(c), costs) for c in 0:num_edges]

        for d in 0:num_vertices
            # Total count at distance d from all bitstrings = 2^n * binomial(n, d)
            total_from_N = 0.0
            for c_prime in 0:num_edges
                if cost_counts[c_prime + 1] > 0
                    total_from_N += cost_counts[c_prime + 1] * sum(N_dist[c_prime + 1, d + 1, :])
                end
            end
            expected = 2^num_vertices * binomial(num_vertices, d)
            @test total_from_N ≈ expected atol=1e-10
        end
    end

    @testset "Binary costs (0 or 1 only)" begin
        num_vertices = 4
        num_edges = 1  # Only costs 0 or 1
        costs = rand(MersenneTwister(666), 0:1, 2^num_vertices) .|> Float64

        N_dist = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)

        # Shape should be (2, 5, 2)
        @test size(N_dist) == (2, num_vertices + 1, 2)

        # For each c', sum over c at each d should be binomial(n, d)
        for c_prime in 0:1
            if any(costs .== c_prime)
                for d in 0:num_vertices
                    @test sum(N_dist[c_prime + 1, d + 1, :]) ≈ binomial(num_vertices, d) atol=1e-10
                end
            end
        end
    end

    @testset "Larger system correctness" begin
        # Test with larger system to catch indexing issues
        num_vertices = 8
        num_edges = 12
        costs = rand(MersenneTwister(777), 0:num_edges, 2^num_vertices) .|> Float64

        n_dist = get_real_distribution_from_costs(costs, num_edges, num_vertices)
        N_dist_twostep = get_homogeneous_distribution_from_costs(costs, n_dist)
        N_dist_direct = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)

        # Both methods should agree
        @test N_dist_twostep ≈ N_dist_direct atol=1e-10

        # Basic properties should hold
        @test size(N_dist_direct) == (num_edges + 1, num_vertices + 1, num_edges + 1)

        # Check row sums for each c' that exists
        for c_prime in 0:num_edges
            if any(costs .== c_prime)
                for d in 0:num_vertices
                    @test sum(N_dist_direct[c_prime + 1, d + 1, :]) ≈ binomial(num_vertices, d) atol=1e-10
                end
            end
        end
    end

    @testset "Distance 0 diagonal property" begin
        # At distance 0, we're only counting the bitstring itself
        # So N(c'; d=0, c) = 1 if c == c', else 0
        num_vertices = 5
        num_edges = 7
        costs = rand(MersenneTwister(888), 0:num_edges, 2^num_vertices) .|> Float64

        N_dist = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)

        for c_prime in 0:num_edges
            if any(costs .== c_prime)
                # At d=0, only c == c' should have count 1
                @test N_dist[c_prime + 1, 1, c_prime + 1] ≈ 1.0 atol=1e-10
                for c in 0:num_edges
                    if c != c_prime
                        @test N_dist[c_prime + 1, 1, c + 1] ≈ 0.0 atol=1e-10
                    end
                end
            end
        end
    end

    @testset "Distance n (maximum) property" begin
        # At distance n, we're only counting the complement bitstring
        # The complement of x has some cost c_complement
        num_vertices = 4
        num_edges = 6
        costs = rand(MersenneTwister(999), 0:num_edges, 2^num_vertices) .|> Float64

        N_dist = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)

        # For each c', sum at d=n should be 1 (only one bitstring at max distance)
        for c_prime in 0:num_edges
            if any(costs .== c_prime)
                @test sum(N_dist[c_prime + 1, num_vertices + 1, :]) ≈ 1.0 atol=1e-10
            end
        end
    end
end


@testset showtiming=true "N(c'; d, c) edge cases" begin
    @testset "Single bitstring (n=1)" begin
        # n=1: bitstrings are 0 and 1
        costs = Float64[0, 1]  # cost 0 for bitstring 0, cost 1 for bitstring 1
        num_edges = 1
        num_vertices = 1

        N_dist = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)

        # Shape: (2 costs, 2 distances, 2 costs)
        @test size(N_dist) == (2, 2, 2)

        # N(0; 0, 0) = 1 (bitstring 0 at distance 0 from itself)
        @test N_dist[1, 1, 1] == 1.0
        # N(0; 1, 1) = 1 (bitstring 1 at distance 1 from bitstring 0)
        @test N_dist[1, 2, 2] == 1.0
        # N(1; 0, 1) = 1 (bitstring 1 at distance 0 from itself)
        @test N_dist[2, 1, 2] == 1.0
        # N(1; 1, 0) = 1 (bitstring 0 at distance 1 from bitstring 1)
        @test N_dist[2, 2, 1] == 1.0
    end

    @testset "All zeros costs" begin
        num_vertices = 3
        num_edges = 0
        costs = zeros(Float64, 2^num_vertices)

        N_dist = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)

        @test size(N_dist) == (1, num_vertices + 1, 1)

        # All entries at c'=0, c=0 should be binomial coefficients
        for d in 0:num_vertices
            @test N_dist[1, d + 1, 1] == binomial(num_vertices, d)
        end
    end

    @testset "Maximum cost everywhere" begin
        num_vertices = 3
        num_edges = 5
        costs = fill(Float64(num_edges), 2^num_vertices)

        N_dist = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)

        # Only max cost slice should have non-zeros
        for c_prime in 0:(num_edges - 1)
            @test all(N_dist[c_prime + 1, :, :] .== 0.0)
        end

        # Max cost slice should have binomial structure
        for d in 0:num_vertices
            @test N_dist[num_edges + 1, d + 1, num_edges + 1] == binomial(num_vertices, d)
        end
    end

    @testset "Sparse costs (only a few values used)" begin
        num_vertices = 4
        num_edges = 10
        # Only use costs 0, 5, and 10
        sparse_costs = [0, 5, 10]
        costs = [Float64(sparse_costs[rand(MersenneTwister(i), 1:3)]) for i in 1:2^num_vertices]

        N_dist = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)

        # Unused cost indices should have all zeros
        for c_prime in 0:num_edges
            if !(c_prime in sparse_costs)
                @test all(N_dist[c_prime + 1, :, :] .== 0.0)
            end
        end

        # Used cost indices should satisfy binomial property
        for c_prime in sparse_costs
            for d in 0:num_vertices
                @test sum(N_dist[c_prime + 1, d + 1, :]) ≈ binomial(num_vertices, d) atol=1e-10
            end
        end
    end
end


#==============================================================================#
#              GPU vs CPU Agreement Tests                                       #
#==============================================================================#

@testset showtiming=true "GPU vs CPU agreement - comprehensive" begin
    if CUDA.has_cuda_gpu()
        @testset "Real distribution: various sizes" begin
            for num_vertices in [4, 5, 6, 7]
                num_edges = num_vertices + 2
                costs = rand(MersenneTwister(num_vertices * 100), 0:num_edges, 2^num_vertices) .|> Float64

                cpu_result = get_real_distribution_from_costs(costs, num_edges, num_vertices)
                gpu_result = gpu_get_real_distribution_from_costs(costs, num_edges, num_vertices) |> Array

                @test size(cpu_result) == size(gpu_result)
                @test cpu_result ≈ gpu_result atol=1e-10
            end
        end

        @testset "Homogeneous distribution: various sizes" begin
            for num_vertices in [4, 5, 6, 7, 8]
                num_edges = num_vertices + 3
                costs = rand(MersenneTwister(num_vertices * 200), 0:num_edges, 2^num_vertices) .|> Float64

                cpu_result = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)
                gpu_result = gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices) |> Array

                @test size(cpu_result) == size(gpu_result)
                @test cpu_result ≈ gpu_result atol=1e-10
            end
        end

        @testset "GPU with uniform costs" begin
            num_vertices = 6
            num_edges = 4
            uniform_cost = 2
            costs = fill(Float64(uniform_cost), 2^num_vertices)

            cpu_result = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)
            gpu_result = gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices) |> Array

            @test cpu_result ≈ gpu_result atol=1e-10

            # Verify the expected structure
            for d in 0:num_vertices
                @test gpu_result[uniform_cost + 1, d + 1, uniform_cost + 1] ≈ binomial(num_vertices, d) atol=1e-10
            end
        end

        @testset "GPU with binary costs" begin
            num_vertices = 7
            num_edges = 1
            costs = rand(MersenneTwister(1234), 0:1, 2^num_vertices) .|> Float64

            cpu_result = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)
            gpu_result = gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices) |> Array

            @test cpu_result ≈ gpu_result atol=1e-10
        end

        @testset "GPU with sparse costs" begin
            num_vertices = 6
            num_edges = 15
            # Only use costs 0, 7, 15
            sparse_costs = [0, 7, 15]
            costs = [Float64(sparse_costs[rand(MersenneTwister(i + 5000), 1:3)]) for i in 1:2^num_vertices]

            cpu_result = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)
            gpu_result = gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices) |> Array

            @test cpu_result ≈ gpu_result atol=1e-10
        end

        @testset "GPU with max_num_edges padding - various sizes" begin
            for num_vertices in [4, 5, 6]
                num_edges = num_vertices
                max_num_edges = num_vertices * 3
                costs = rand(MersenneTwister(num_vertices * 300), 0:num_edges, 2^num_vertices) .|> Float64

                cpu_result = get_homogeneous_distribution_from_costs_direct(
                    costs, num_edges, num_vertices; max_num_edges=max_num_edges
                )
                gpu_result = gpu_get_homogeneous_distribution_from_costs_direct(
                    costs, num_edges, num_vertices; max_num_edges=max_num_edges
                ) |> Array

                @test size(cpu_result) == size(gpu_result)
                @test size(cpu_result, 1) == max_num_edges + 1
                @test cpu_result ≈ gpu_result atol=1e-10

                # Padded region should be zeros
                @test all(cpu_result[num_edges + 2:end, :, :] .== 0.0)
                @test all(gpu_result[num_edges + 2:end, :, :] .≈ 0.0)
            end
        end

        @testset "GPU real distribution with padding" begin
            num_vertices = 5
            num_edges = 4
            max_num_edges = 12
            costs = rand(MersenneTwister(4321), 0:num_edges, 2^num_vertices) .|> Float64

            cpu_result = get_real_distribution_from_costs(costs, num_edges, num_vertices; max_num_edges=max_num_edges)
            gpu_result = gpu_get_real_distribution_from_costs(costs, num_edges, num_vertices; max_num_edges=max_num_edges) |> Array

            @test size(cpu_result) == size(gpu_result)
            @test size(cpu_result, 3) == max_num_edges + 1
            @test cpu_result ≈ gpu_result atol=1e-10
        end

        @testset "GPU consistency across multiple runs" begin
            # Ensure GPU results are deterministic
            num_vertices = 6
            num_edges = 8
            costs = rand(MersenneTwister(9999), 0:num_edges, 2^num_vertices) .|> Float64

            gpu_result1 = gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices) |> Array
            gpu_result2 = gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices) |> Array
            gpu_result3 = gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices) |> Array

            @test gpu_result1 ≈ gpu_result2 atol=1e-14
            @test gpu_result2 ≈ gpu_result3 atol=1e-14
        end
    else
        @warn "Skipping comprehensive GPU tests because no GPU detected"
    end
end


@testset showtiming=true "GPU vs CPU - special cost patterns" begin
    if CUDA.has_cuda_gpu()
        @testset "Alternating costs" begin
            # Costs alternate: 0, 1, 0, 1, ...
            num_vertices = 6
            num_edges = 1
            costs = Float64[i % 2 for i in 0:(2^num_vertices - 1)]

            cpu_result = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)
            gpu_result = gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices) |> Array

            @test cpu_result ≈ gpu_result atol=1e-10
        end

        @testset "Linearly increasing costs" begin
            num_vertices = 5
            num_edges = 2^num_vertices - 1
            costs = Float64.(0:(2^num_vertices - 1))

            cpu_result = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)
            gpu_result = gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices) |> Array

            @test cpu_result ≈ gpu_result atol=1e-10

            # Each bitstring has unique cost, so N(c'; d, c) = n(x; d, c) for the unique x with cost c'
            # This means each c' slice should have exactly one count per distance
            for c_prime in 0:num_edges
                for d in 0:num_vertices
                    total = sum(cpu_result[c_prime + 1, d + 1, :])
                    @test total ≈ binomial(num_vertices, d) atol=1e-10
                end
            end
        end

        @testset "Popcount-based costs" begin
            # Cost = number of 1s in the bitstring (popcount)
            num_vertices = 6
            num_edges = num_vertices  # max cost is n
            costs = Float64[count_ones(i) for i in 0:(2^num_vertices - 1)]

            cpu_result = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)
            gpu_result = gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices) |> Array

            @test cpu_result ≈ gpu_result atol=1e-10

            # Verify: number of bitstrings with cost c is binomial(n, c)
            for c in 0:num_vertices
                @test count(==(c), costs) == binomial(num_vertices, c)
            end
        end

        @testset "MaxCut-like cost structure" begin
            # Simulate MaxCut costs: c(x) + c(~x) = num_edges for all x
            # We create costs where this property roughly holds
            num_vertices = 5
            num_edges = 10
            num_bitstrings = 2^num_vertices

            # Generate costs respecting MaxCut symmetry
            costs = zeros(Float64, num_bitstrings)
            for x in 0:(num_bitstrings ÷ 2 - 1)
                c = rand(MersenneTwister(x + 10000), 0:num_edges)
                complement_x = (2^num_vertices - 1) - x
                costs[x + 1] = Float64(c)
                costs[complement_x + 1] = Float64(num_edges - c)
            end

            cpu_result = get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices)
            gpu_result = gpu_get_homogeneous_distribution_from_costs_direct(costs, num_edges, num_vertices) |> Array

            @test cpu_result ≈ gpu_result atol=1e-10
        end
    else
        @warn "Skipping special pattern GPU tests because no GPU detected"
    end
end

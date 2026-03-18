using JuliaQAOA, Test, CUDA
import Random: MersenneTwister



@testset showtiming=true "_expand(X)" begin
    @test _expand(1) == 1
    @test _expand([1,2]) == [1 2]
    @test _expand([1; 2;; 3; 4]) == [1;; 2;;; 3;; 4]
end

@testset showtiming=true "get_β_factors" begin
    @testset showtiming=true "Agrees with manufactured solution" begin
        n = 4
        @test get_β_factors(0,       n, pi_units=true) == [1,0,0,0,0]
        @test get_β_factors(0.25,    n, pi_units=true) ≈ fill(0.25, 5) .* [1, -1im, -1, 1im, 1] atol=1e-15
        @test get_β_factors(0.5,     n, pi_units=true) == [0,0,0,0,1]
        @test get_β_factors(0,       n, pi_units=false) == [1,0,0,0,0]
        @test get_β_factors(0.25*pi, n, pi_units=false) ≈ fill(0.25, 5) .* [1, -1im, -1, 1im, 1] atol=1e-15
        @test get_β_factors(0.5*pi,  n, pi_units=false) ≈ [0,0,0,0,1] atol=1e-15
    end

    @testset showtiming=true "Multidimensional β broadcasting consistent with scalar β." begin
        n = 4
        β_vec = [0, pi/6, pi/4]
        @test get_β_factors(β_vec, n) ≈ reduce(hcat, [get_β_factors(β, n) for β in β_vec]) atol=1e-15
        β_mat = [0    pi/6 pi/4
                 pi/3 pi/2 pi]
        @test get_β_factors(β_mat, n) ≈ reduce((x,y) -> cat(x, y; dims=3), [get_β_factors(β_vec, n) for β_vec in eachcol(β_mat)]) atol=1e-15
    end

    @testset showtiming=true "GPU result agrees with CPU result" begin
        if CUDA.has_cuda_gpu()
            n = 4
            β_Array = rand(MersenneTwister(0), 4, 5, 6) 
            @test get_β_factors(β_Array, n, pi_units=true) ≈ Array(get_β_factors(CuArray(β_Array), n, pi_units=true)) atol=1-14
            βpi_Array = β_Array .* pi
            @test get_β_factors(βpi_Array, n, pi_units=false) ≈ Array(get_β_factors(CuArray(βpi_Array), n, pi_units=false)) atol=1-14
        else
            @warn "Skipping GPU test because no GPU detected"
        end 

    end
end

@testset "get_γ_factors" begin
    @testset showtiming=true "Agrees with manufactured solution" begin
        # Convention: get_γ_factors computes exp(-iγc/2), so γ values are 2× the
        # phase angle per unit cost. This matches QOKit's exp(-0.5j*gamma*hc_diag).
        n = 4
        @test get_γ_factors(0,       n, pi_units=true) == ones(n+1)
        @test get_γ_factors(1.0,    n, pi_units=true) == [1, -1im, -1, 1im, 1]
        @test get_γ_factors(2.0,     n, pi_units=true) == [1, -1, 1, -1, 1]
        @test get_γ_factors(0,       n, pi_units=false) == ones(n+1)
        @test get_γ_factors(1.0*pi, n, pi_units=false) ≈ [1, -1im, -1, 1im, 1] atol=1e-15
        @test get_γ_factors(2*pi,  n, pi_units=false) ≈ [1, -1, 1, -1, 1] atol=1e-15
    end

    @testset showtiming=true "Multidimensional γ broadcasting consistent with scalar γ." begin
        n = 4
        γ_vec = [0, pi/6, pi/4]
        @test get_γ_factors(γ_vec, n) ≈ reduce(hcat, [get_γ_factors(γ, n) for γ in γ_vec]) atol=1e-15
        γ_mat = [0    pi/6 pi/4
                 pi/3 pi/2 pi]
        @test get_γ_factors(γ_mat, n) ≈ reduce((x,y) -> cat(x, y; dims=3), [get_γ_factors(γ_vec, n) for γ_vec in eachcol(γ_mat)]) atol=1e-15
    end

    @testset showtiming=true "CUDA result agrees with CPU result" begin
        if CUDA.has_cuda_gpu()
            n = 4
            γ_Array = rand(MersenneTwister(0), 4, 5, 6) 
            @test get_γ_factors(γ_Array, n, pi_units=true) ≈ Array(get_γ_factors(CuArray(γ_Array), n, pi_units=true)) atol=1e-14
            γpi_Array = γ_Array .* pi
            @test get_γ_factors(γpi_Array, n, pi_units=false) ≈ Array(get_γ_factors(CuArray(γpi_Array), n, pi_units=false)) atol=1e-14
        else
            @warn "Skipping GPU test because no GPU detected"
        end 
    end
end

@testset showtiming=true "QAOA algorithm" begin
    @testset showtiming=true "Agrees with manufactured solution." begin
        N = zeros(2,2,2)
        N[1,:,:] .= [1 3
                     2 4]
        N[2,:,:] .= [5 7
                     6 8]

        # in units of pi (γ uses QOKit convention: phase gate is exp(-iγC/2))
        γs = [2]
        βs = [3/4]

        manufactured_solution = [1+1im, 1+1im]

        @test QAOA_proxy_basic(N, γs, βs, pi_units=true)[end] ≈ manufactured_solution
        @test QAOA_proxy_single(N, γs, βs, pi_units=true, blas=true)[end] ≈ manufactured_solution
        @test QAOA_proxy_single(N, γs, βs, pi_units=true, blas=false)[end] ≈ manufactured_solution
        @test vec(QAOA_proxy_multi(N, _expand(γs), _expand(βs), pi_units=true, blas=true)[end]) ≈ manufactured_solution
    end
    @testset showtiming=true "Implementations agree for larger, random example." begin
        m = 50
        n = 10
        p = 5
        N = rand(MersenneTwister(0), 1+m, 1+n, 1+m)
        γs = rand(MersenneTwister(1), p)
        βs = rand(MersenneTwister(2), p)
        basic_result = reduce(hcat, QAOA_proxy_basic(N, γs, βs))
        single_result = reduce(hcat, QAOA_proxy_single(N, γs, βs))
        multi_result = reduce(hcat, QAOA_proxy_multi(N, _expand(γs), _expand(βs)))

        @test basic_result ≈ single_result
        @test basic_result ≈ multi_result
        @test single_result ≈ multi_result
    end
    @testset showtiming=true "QAOA_proxy_multi on multiple parameter sets agrees with" begin
        m = 50
        n = 10
        p = 5
        N = rand(MersenneTwister(0), 1+m, 1+n, 1+m)
        num_param_sets = 4
        γs = rand(MersenneTwister(1), num_param_sets, p)
        βs = rand(MersenneTwister(2), num_param_sets, p)
        Qs_single_collection = [QAOA_proxy_single(N, γs[j,:], βs[j,:]) for j in 1:num_param_sets]
        # Put in same format as result of QAOA_proxy_mult
        Qs_single = [hcat(getindex.(Qs_single_collection, 1+j)...) for j in 0:p]
        Qs_multi = QAOA_proxy_multi(N, γs, βs)

        @test cat(Qs_single..., dims=3) ≈ cat(Qs_multi..., dims=3)
    end
    @testset showtiming=true "GPU result agrees with CPU result" begin
        @testset "QAOA_proxy_single"  begin
            if CUDA.has_cuda_gpu()
                m = 50
                n = 10
                p = 5
                N = rand(MersenneTwister(0), Float32, 1+m, 1+n, 1+m)
                γs = rand(MersenneTwister(1), Float32, p)
                βs = rand(MersenneTwister(2), Float32, p)
                Qs_CPU = reduce(hcat, QAOA_proxy_single(N, γs, βs))
                Qs_GPU = reduce(hcat, QAOA_proxy_single(cu(N), cu(γs), cu(βs))) |> Array
                @test Qs_CPU ≈ Qs_GPU
            else
                @warn "Skipping GPU test because no GPU detected"
            end
        end
        @testset showtiming=true "QAOA_proxy_multi"  begin
            if CUDA.has_cuda_gpu()
                m = 50
                n = 10
                p = 5
                num_param_sets = 4
                N = rand(MersenneTwister(0), Float32, 1+m, 1+n, 1+m)
                γs = rand(MersenneTwister(1), Float32, num_param_sets, p)
                βs = rand(MersenneTwister(2), Float32, num_param_sets, p)
                Qs_CPU = cat(QAOA_proxy_multi(N, γs, βs)..., dims=3)
                Qs_GPU = cat(QAOA_proxy_multi(cu(N), cu(γs), cu(βs))..., dims=3) |> Array
                @test Qs_CPU ≈ Qs_GPU
            else
                @warn "Skipping GPU test because no GPU detected"
            end
        end

    end
end

@testset showtiming=true "QAOA expectation value calculation" begin
    @testset showtiming=true "Agrees with manufactured solution" begin
        Q = [1, 2im, 3+4im]
        Q_mat = reshape(Q, :, 1)
        P = [4, 5, 6]
        n = 2
        @test expectation(Q, P, n) isa Real
        @test size(expectation(Q_mat, P, n)) == (1,1)
        @test expectation(Q, P, n) == 4*(0*1*4 + 1*4*5 + 2*25*6)
        @test expectation(Q_mat, P, n) == [4*(0*1*4 + 1*4*5 + 2*25*6);;]
    end
    @testset showtiming=true "Matrix Q version agrees with vector Q version" begin
        m = 5
        n = 3
        num_states = 4
        Q = rand(MersenneTwister(0), ComplexF64, 1+m, num_states)
        P = rand(MersenneTwister(1), 1+m)

        vector_expectations = reduce(hcat, expectation(Qcol, P, n) for Qcol in eachcol(Q))
        matrix_expectations = expectation(Q, P, n)
        @test vector_expectations ≈ matrix_expectations atol=1e-14
    end
    @testset showtiming=true "GPU result agrees with CPU result"  begin
        if CUDA.has_cuda_gpu()
            m = 5
            n = 3
            num_states = 4
            Qvec = rand(MersenneTwister(0), ComplexF32, 1+m)
            Qmat = rand(MersenneTwister(0), ComplexF32, 1+m, num_states)
            P = rand(MersenneTwister(1), Float32, 1+m)

            @test expectation(cu(Qvec), cu(P), n) ≈ expectation(Qvec, P, n) rtol=1e-7
            @test Array(expectation(cu(Qmat), cu(P), n)) ≈ expectation(Qmat, P, n) rtol=1e-7
        else
            @warn "Skipping GPU test because no GPU detected"
        end
    end

end



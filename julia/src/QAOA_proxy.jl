"""
Expand X into a shape suitable for broadcasting.

- Number       -> leave as scalar
- Vector       -> make it a 1×N row vector
- Matrix       -> make it a 1×M×N array
- Higher ranks -> prepend a singleton dimension
"""
_expand(X::Number) = X
_expand(X::AbstractVector) = reshape(X, 1, :)
_expand(X::AbstractMatrix) = reshape(X, 1, size(X, 1), size(X, 2))
_expand(X::AbstractArray) = reshape(X, 1, size(X)...)

"""
Compute cos(β)^(n-d)(-i sinβ)^d  for d=0:n.

For N-dimensional β, return N+1-dimensional array, where first dimension
corresponds to varying d, rest correspond to varying β.

By default, β is given in units of pi.

n = problem size / number of nodes / number of qubits
"""
function get_β_factors(β, n::Integer; pi_units=true)
    @assert eltype(β) <: Real "β values must be real!"
    expanded_β = _expand(β)
    sinβ = pi_units ? sinpi.(expanded_β) : sin.(expanded_β)
    cosβ = pi_units ? cospi.(expanded_β) : cos.(expanded_β)

    # cu(range) converts to CUDA version *only if needed.*
    return @. cosβ ^ (n:-1:0) * ((-im*sinβ) ^ (0:n))
end


"""
Compute exp(-iγc) for c=0:m.

For N-dimensional γ, return N+1-dimensional array, where first dimension
corresponds to varying c, rest correspond to varying γ.

By default, γ is given in units of pi.

m = number of constraints / number of edges / number of distinct costs minus one
"""
function get_γ_factors(γ, m::Integer; pi_units=true)
    @assert eltype(γ) <: Real "γ values must be real!"
    expanded_γ = _expand(γ)
    # cu(range) converts to CUDA version *only if needed.*
    # cis(x) equals exp(im*x)
    if pi_units
        return @. cispi((-1 * expanded_γ) * (0:m))
    else
        return @. cis((-1 * expanded_γ) * (0:m))
    end
end

"""
Basic version of the QAOA proxy, which follows the algorithm from the paper
most closely. Should be less efficient than other versions, but most likely to
be implemented correctly.
"""
function QAOA_proxy_basic(
        N::AbstractArray{<: Real, 3},
        γs::AbstractVector{<: Real},
        βs::AbstractVector{<: Real}, 
        ; pi_units=true,
    )
    @assert length(βs) == length(γs) "βs and γs must both have the same length (p)." 
    @assert size(N, 1) == size(N, 3) "1st and 3rd dimensions of homogeneous distribution must be the same length (1+m)."
    m = size(N, 1) - 1
    n = size(N, 2) - 1
    p = length(γs)

    Q0 = fill(ComplexF64(1/sqrt(2^n)), 1+m)
    Qs = [Q0] # Store each Q in a vector
    for ℓ in 1:p
        γ = γs[ℓ]
        β = βs[ℓ]

        Qprev = Qs[end]
        Qℓ = zeros(ComplexF64, 1+m)
        for c_prime in 0:m
            if pi_units
                Qℓ[1+c_prime] = sum(
                    cospi(β)^(n-d) * (-1im*sinpi(β))^d * cispi(-γ*c) * Qprev[1+c] *N[1+c_prime, 1+d, 1+c]
                    for d in 0:n, c in 0:m
                )
            else
                Qℓ[1+c_prime] = sum(
                    cos(β)^(n-d) * (-1im*sin(β))^d * cis(-γ*c) * Qprev[1+c] *N[1+c_prime, 1+d, 1+c]
                    for d in 0:n, c in 0:m
                )
            end
        end
        push!(Qs, Qℓ)
    end
    return Qs
end

"""
Version using a fused mat-vec mult approach. Should be faster due to use of
BLAS and multithreading.
"""
function QAOA_proxy_single(
        N::AbstractArray{<: Real, 3},
        γs::AbstractVector{<: Real},
        βs::AbstractVector{<: Real}, 
        ; pi_units::Bool=true, blas::Bool=true
    )
    @assert length(βs) == length(γs) "βs and γs must both have the same length (p)." 
    @assert size(N, 1) == size(N, 3) "1st and 3rd dimensions of homogeneous distribution must be the same length (1+m)."
    m = size(N, 1) - 1
    n = size(N, 2) - 1
    p = length(γs)

    γ_factors = get_γ_factors(γs, m, pi_units=pi_units)
    β_factors = get_β_factors(βs, n, pi_units=pi_units)
    # Reshape N into matrix, second dim is multi-index of (d,c)
    M = reshape(N, 1+m, :) 

    # Q₀
    Q0 = similar(N, complex(eltype(N)), 1+m)
    Q0 .= 1/sqrt(2^n)
    Qs = [Q0] # Store each Q in a vector

    for ℓ in 1:p
        γ_factors_ℓ = γ_factors[:,ℓ]
        β_factors_ℓ = β_factors[:,ℓ]
        
        Qprev = Qs[end]
        v_2D = β_factors_ℓ .* _expand(γ_factors_ℓ .* Qprev)
        v = reshape(v_2D, :) # length (1+m)*(1+d)
        Qℓ = similar(Q0)
        if blas
            mul!(Qℓ, M, v)
        else
            Qℓ .= M*v
        end
        push!(Qs, Qℓ)
    end
    return Qs
end

"""
Version using a fused mat-mat mult approach. Should be faster due to use of
BLAS and multithreading.

For γs and βs, note that the *column index* increases with p.
Contrast with QAOA_proxy_single, where the *row index* increases with p.
"""
function QAOA_proxy_multi(
        N::AbstractArray{<: Real, 3},
        γs::AbstractVecOrMat{<: Real},
        βs::AbstractVecOrMat{<: Real}, 
        ; pi_units::Bool=true, blas::Bool=true
    )
    @assert size(βs) == size(γs) "βs and γs must both have the same dimensions." 
    @assert size(N, 1) == size(N, 3) "1st and 3rd dimensions of homogeneous distribution must be the same length (1+m)."
    m = size(N, 1) - 1
    n = size(N, 2) - 1
    num_param_sets = size(γs, 1)
    p = size(γs, 2)

    γ_factors = get_γ_factors(γs, m, pi_units=pi_units)
    β_factors = get_β_factors(βs, n, pi_units=pi_units)
    # Reshape N into matrix, second dim is multi-index of (d,c)
    M = reshape(N, 1+m, (1+m)*(1+n)) 

    # Q₀
    Q0 = similar(N, complex(eltype(N)), 1+m, num_param_sets)
    Q0 .= 1/sqrt(2^n)
    Qs = [Q0] # Store each Q in a vector

    for ℓ in 1:p
        γ_factors_ℓ = γ_factors[:,:,ℓ]
        β_factors_ℓ = reshape(β_factors[:,:,ℓ], 1+n, 1, num_param_sets)
        
        Qprev = Qs[end]
        v_3D = β_factors_ℓ .* _expand(γ_factors_ℓ .* Qprev)
        v = reshape(v_3D, (1+m)*(1+n), num_param_sets) # size (1+m)*(1+n) × num_param_sets
        Qℓ = similar(Q0)
        if blas
            mul!(Qℓ, M, v)
        else
            Qℓ .= M*v
        end
        push!(Qs, Qℓ)
    end
    return Qs
end


#=

"""
Given a set of proxy amplitudes (i.e. the homomgeneous proxy of the state
vector.), compute the expectation value associated with those amplitudes.
"""
function QAOA_proxy_expectation(
        proxy,
        proxy_amplitudes::AbstractVector{<: Number},
        drop_first_N_costs::Integer=0
    )::Float64

    proxy_expectation_value::Float64 = 0.0
    num_unique_costs = length(proxy_amplitudes)

    for cost in drop_first_N_costs:num_unique_costs-1
        proxy_expectation_value += cost * N_cost_distribution(proxy, cost) * abs2(proxy_amplitudes[1+cost])
    end

    return proxy_expectation_value
end



"""
Inverse here is the additive inverse, not multiplicative inverse.

I.e. we are taking -x, not 1/x.

Scipy allows us to minimize functions, but we want to maximize our cost
(the number of edges cut), so we have to flip the +/- sign to turn our
problem into a minimization problem.
"""
function inverse_proxy_objective_function(
        proxy,
        num_QAOA_layers::Integer,
        expectations::Union{Nothing, AbstractVector}=nothing
    )

    function inverse_objective(args...)::Float64
        # Note that args[0] in Python is args[1] in Julia 
        gammas_betas_combined_vec = pyconvert(Vector, args[1]) # Conversion may not be necessary. Remove if impacts performance significantly.

        @assert length(gammas_betas_combined_vec) == 2*num_QAOA_layers
        gammas = @view gammas_betas_combined_vec[1:num_QAOA_layers]
        betas = @view gammas_betas_combined_vec[num_QAOA_layers+1:end]

        proxy_amplitudes = QAOA_proxy(
            proxy,
            gammas,
            betas,
        )

        final_proxy_amplitudes = @view proxy_amplitudes[:,end]

        expectation = QAOA_proxy_expectation(
            proxy,
            final_proxy_amplitudes
        )
        current_time = time()

        if !isnothing(expectations)
            push!(expectations, (current_time, expectation))
        end

        return -expectation
    end

    return inverse_objective
end
=#

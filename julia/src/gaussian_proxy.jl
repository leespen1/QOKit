#=
gaussian_proxy.jl — Gaussian blob proxy for N(c';d,c).

For each cost class c', N(c';d,c) is modeled as a sum of contributions from
a Gaussian-like shape in (d,c) space. More expressive than TriangleProxy
(which uses piecewise linear), but still fully parameterized.

Key design choices:
- For each c' and distance d, the distribution over c is Gaussian with:
  - Mean: linearly interpolated from c' (at d=0) to a global center (at d=n/2)
  - Stddev: scales with sqrt(d) (like a diffusion process)
  - Height: scales with distance (peaks at d=n/2, minimal at d=0)
- The binomial weight C(n,d) is folded into the height
- Consistent P(c') is derived from the marginal of N

Parameters (7 total, more expressive than Triangle's 4):
- center_target: where the mean converges at d=n/2 (fraction of m)
- sigma_base: base standard deviation at d=1
- sigma_scale: how fast sigma grows with sqrt(d)
- height_base: height at d=0 (should be ~1 for the delta-like d=0 contribution)
- height_scale: how fast height grows with d
- height_power: exponent for height growth (1=linear, 0.5=sqrt, etc.)
- center_bias: asymmetric shift of center toward low or high cost
=#

using Distributions: Binomial, pdf

struct GaussianProxy <: AbstractProxy
    num_constraints::Int64
    num_qubits::Int64
    center_target::Float64    # Where mean cost converges at d=n/2 (as fraction of m)
    sigma_base::Float64       # Base sigma at d=1
    sigma_scale::Float64      # Sigma growth rate with sqrt(d)
    height_base::Float64      # Height at d=0
    height_scale::Float64     # Height scaling factor
    height_power::Float64     # Height growth exponent
    center_bias::Float64      # Asymmetric center shift
end

"""Default constructor with reasonable initial parameters."""
function GaussianProxy(num_constraints::Integer, num_qubits::Integer;
    center_target::Real=0.5,
    sigma_base::Real=1.0,
    sigma_scale::Real=1.0,
    height_base::Real=1.0,
    height_scale::Real=1.0,
    height_power::Real=0.5,
    center_bias::Real=0.0,
)
    return GaussianProxy(num_constraints, num_qubits,
        Float64(center_target), Float64(sigma_base), Float64(sigma_scale),
        Float64(height_base), Float64(height_scale), Float64(height_power),
        Float64(center_bias))
end

"""
P(c') — cost distribution. For GaussianProxy, we use a simple approximation:
a Gaussian centered at m*center_target with appropriate width.
This ensures N/P consistency.
"""
function P_cost_distribution(proxy::GaussianProxy, cost::Integer)::Float64
    m = proxy.num_constraints
    n = proxy.num_qubits
    # Use binomial as base (like PaperProxy) but with adjusted center
    # For consistency, P should be related to the proxy's view of cost distribution
    # Use a Gaussian approximation centered at m*center_target
    μ = m * proxy.center_target
    # Variance from binomial-like: m * p * (1-p) where p = center_target
    σ² = m * proxy.center_target * (1 - proxy.center_target)
    σ = sqrt(max(σ², 0.1))
    # Gaussian PDF, discretized
    val = exp(-0.5 * ((cost - μ) / σ)^2) / (σ * sqrt(2π))
    return val
end

"""
N(c') — expected number of bitstrings with cost c'.
"""
function N_cost_distribution(proxy::GaussianProxy, cost::Integer)::Float64
    return P_cost_distribution(proxy, cost) * (1 << proxy.num_qubits)
end

"""
N(c'; d, c) — the Gaussian blob model.

For fixed c' and d:
- The distribution over c is Gaussian centered at μ_c(c', d) with stddev σ(d)
- μ_c slides from c' at d=0 toward m*center_target at d=n/2
- σ grows as sigma_base + sigma_scale * sqrt(d)
- Height scales as height_base + height_scale * d^height_power,
  weighted by binomial(n, d)
"""
function N_cost_distance_distribution(proxy::GaussianProxy,
    cost_1::Integer, distance::Integer, cost_2::Integer)::Float64

    n = proxy.num_qubits
    m = proxy.num_constraints

    # Reflect distance to [0, n/2] (symmetric around n/2)
    d = distance > div(n, 2) ? n - distance : distance

    # Special case: d=0 means cost_2 must equal cost_1
    if d == 0
        return cost_2 == cost_1 ? proxy.height_base : 0.0
    end

    half_n = n / 2.0
    t = d / half_n  # t ∈ [0, 1] as d goes from 0 to n/2

    # Mean cost at this distance: interpolate from cost_1 to center
    center = m * proxy.center_target + proxy.center_bias * (cost_1 - m * proxy.center_target)
    μ_c = cost_1 * (1 - t) + center * t

    # Standard deviation grows with sqrt(distance)
    σ = proxy.sigma_base + proxy.sigma_scale * sqrt(d)
    σ = max(σ, 0.1)  # prevent degenerate cases

    # Height: how many bitstrings at distance d
    # Base it on binomial(n, d) shape but with tunable scaling
    # C(n,d) / 2^n using log to avoid overflow
    log_binom = sum(log(n - i + 1) - log(i) for i in 1:d; init=0.0) - n * log(2)
    binom_weight = exp(log_binom)
    height = (proxy.height_base + proxy.height_scale * d^proxy.height_power) * binom_weight * (1 << n)

    # Gaussian value
    z = (cost_2 - μ_c) / σ
    gaussian = exp(-0.5 * z^2) / (σ * sqrt(2π))

    return height * gaussian
end

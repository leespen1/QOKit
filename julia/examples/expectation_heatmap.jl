using JuliaQAOA, GLMakie, LaTeXStrings
import Random: MersenneTwister

# Set up problem
n = 30
m = 450
prob_edge = 0.5
n_gridpoints = 100
println("n=$n, m=$m, n_gridpoints=$n_gridpoints")

#proxy = PaperProxy(m, n, prob_edge)
proxy = IntuitiveTriangleProxy(m, n)
t1 = time()
P = P_cost_distribution.(proxy, 0:m)
t2 = time()
println("Time taken to compute P: $(t2-t1)")

γ_range = range(0, 2, length=n_gridpoints)
β_range = range(0, 2, length=n_gridpoints)
γs = repeat(γ_range, inner=length(β_range))
βs = repeat(β_range, outer=length(γ_range))

t1 = time()
N = JuliaQAOA.cpu_compute_homodist(proxy)
t2 = time()
println("Time taken to compute N: $(t2-t1)")

# Collect heatmap data
t1 = time()
Qf = QAOA_proxy_multi(N, γs, βs; pi_units=true) |> last
t2 = time()
println("Time taken to run QAOA proxy: $(t2-t1)")

t1 = time()
expectations = reshape(expectation(Qf, P, n), length(β_range), length(γ_range))
t2 = time()
println("Time taken to compute expectations: $(t2-t1)")


# Plot heatmap
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel=L"\gamma \times \pi",
    ylabel=L"\beta \times \pi",
)


hm = heatmap!(ax, γ_range, β_range, expectations)
Colorbar(fig[1, end+1], hm)

fig

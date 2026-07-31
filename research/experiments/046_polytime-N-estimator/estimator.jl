# E046 prototype: polynomial-time (in n) Monte Carlo estimator of the
# homogeneous distribution N(c'; d, c), the class counts M_c, and the attained
# cost classes, using only the edge list (never the 2^n cost array).
#
# Pipeline:
#   1. Wang-Landau over cost classes c in 0..m  -> ln g(c) (density of states)
#   2. Flat-histogram production chain (pi(x) ∝ 1/g(c(x))) -> up to S distinct
#      members per attained cost class
#   3. Per member x and distance d: sample K random d-subsets of flip
#      positions, evaluate costs from scratch (O(m) each), and scale by
#      C(n,d): n_hat(x; d, c) = C(n,d) * freq(c). Distances with
#      C(n,d) <= K are enumerated exactly (Gosper's hack).
#   4. counts_hat[c] = 2^n * softmax(ln g)[c]
#
# Total cost is polynomial in n and m; nothing enumerates 2^n states.

using Random: AbstractRNG, Xoshiro
using Base.Threads: @threads

"Cut cost of state `s` (bits = partition) for 0-based `edges`. O(m)."
@inline function cut_cost(s::Integer, edges::Vector{Tuple{Int,Int}})
    c = 0
    @inbounds for (i, j) in edges
        c += ((s >> i) & 1) != ((s >> j) & 1)
    end
    return c
end

"Adjacency lists (0-based vertices) from the edge list."
function adjacency(n::Int, edges::Vector{Tuple{Int,Int}})
    adj = [Int[] for _ in 1:n]
    for (i, j) in edges
        push!(adj[i + 1], j)
        push!(adj[j + 1], i)
    end
    return adj
end

"Cost change from flipping bit `v` (0-based) of state `s`. O(deg(v))."
@inline function flip_delta(s::Integer, v::Int, adj::Vector{Vector{Int}})
    δ = 0
    bv = (s >> v) & 1
    @inbounds for u in adj[v + 1]
        δ += ((s >> u) & 1) == bv ? 1 : -1
    end
    return δ
end

"""
    wang_landau(n, edges; rng, flat=0.8, lnf_final=1e-3,
                check_every=10_000*n, stage_cap=2_000_000) -> lng

Wang-Landau estimate of ln g(c) for c in 0..m (Vector of length m+1;
`-Inf` marks classes never visited). Flatness is checked over visited
classes only; a stage that hits `stage_cap` steps advances anyway (the
estimate keeps refining in later stages).
"""
function wang_landau(n::Int, edges::Vector{Tuple{Int,Int}}; rng::AbstractRNG,
                     flat::Float64=0.8, lnf_final::Float64=1e-3,
                     check_every::Int=10_000 * n, stage_cap::Int=2_000_000)
    m = length(edges)
    adj = adjacency(n, edges)
    lng = fill(-Inf, m + 1)
    hist = zeros(Int, m + 1)
    s = rand(rng, 0:(1 << n) - 1)
    c = cut_cost(s, edges)
    lng[c + 1] = 0.0
    lnf = 1.0
    while lnf > lnf_final
        steps = 0
        while steps < stage_cap
            for _ in 1:check_every
                v = rand(rng, 0:n-1)
                c2 = c + flip_delta(s, v, adj)
                lg2 = lng[c2 + 1]
                if isinf(lg2) || log(rand(rng)) < lng[c + 1] - lg2
                    s ⊻= 1 << v
                    c = c2
                    if isinf(lng[c + 1]); lng[c + 1] = 0.0; end
                end
                lng[c + 1] += lnf
                hist[c + 1] += 1
            end
            steps += check_every
            vis = findall(isfinite, lng)
            hv = view(hist, vis)
            if minimum(hv) >= flat * (sum(hv) / length(vis))
                break
            end
        end
        fill!(hist, 0)
        lnf /= 2
    end
    return lng
end

"""
    sample_class_members(n, edges, lng; rng, samples_per_class,
                         thin=n, sweep_cap=300*samples_per_class) -> members

Flat-histogram Metropolis chain with stationary pi(x) ∝ 1/g(c(x)); records
(after `thin`-flip thinning) up to `samples_per_class` distinct states per
attained class. Runs until every visited class is full or the step budget
(`sweep_cap` thinned records per class) is exhausted; tiny classes simply
yield fewer members.
"""
function sample_class_members(n::Int, edges::Vector{Tuple{Int,Int}},
                              lng::Vector{Float64}; rng::AbstractRNG,
                              samples_per_class::Int, thin::Int=n,
                              sweep_cap::Int=300 * samples_per_class)
    m = length(edges)
    adj = adjacency(n, edges)
    vis = findall(isfinite, lng)
    members = [Set{Int}() for _ in 1:(m + 1)]
    s = rand(rng, 0:(1 << n) - 1)
    c = cut_cost(s, edges)
    @assert isfinite(lng[c + 1]) "production chain started in a class Wang-Landau never visited"
    budget = sweep_cap * length(vis)
    for _ in 1:budget
        for _ in 1:thin
            v = rand(rng, 0:n-1)
            c2 = c + flip_delta(s, v, adj)
            lg2 = lng[c2 + 1]
            if isfinite(lg2) && log(rand(rng)) < lng[c + 1] - lg2
                s ⊻= 1 << v
                c = c2
            end
        end
        mem = members[c + 1]
        length(mem) < samples_per_class && push!(mem, s)
        all(ci -> length(members[ci]) >= samples_per_class, vis) && break
    end
    return members
end

"Next integer with the same popcount (Gosper's hack)."
@inline function next_popcount(x::Int)
    u = x & -x
    v = x + u
    return v + (((v ⊻ x) ÷ u) >> 2)
end

"""
    profile_row!(row, x, n, edges, d, K, rng)

Accumulate the estimate of n(x; d, c) over c into `row` (length m+1).
Exact enumeration when C(n,d) <= K, else K sampled d-subsets scaled by
C(n,d). Each row sums to C(n,d) exactly by construction.
"""
function profile_row!(row::AbstractVector{Float64}, x::Int, n::Int,
                      edges::Vector{Tuple{Int,Int}}, d::Int, K::Int,
                      rng::AbstractRNG)
    total = binomial(n, d)
    if total <= K
        mask = (1 << d) - 1
        last = mask << (n - d)
        while true
            row[cut_cost(x ⊻ mask, edges) + 1] += 1.0
            mask == last && break
            mask = next_popcount(mask)
        end
    else
        w = total / K
        pos = collect(0:n-1)
        for _ in 1:K
            mask = 0
            for i in 1:d   # partial Fisher-Yates -> d distinct positions
                j = rand(rng, i:n)
                pos[i], pos[j] = pos[j], pos[i]
                mask |= 1 << pos[i]
            end
            row[cut_cost(x ⊻ mask, edges) + 1] += w
        end
    end
    return row
end

"""
    polytime_homogeneous_distribution(n, edges; samples_per_class=10,
        K_subsets=200, rng, wl_kwargs...) -> (N, counts, members)

Polynomial-time estimate of the homogeneous distribution N(c'; d, c)
(shape (m+1, n+1, m+1), zero rows for classes Wang-Landau never reached),
the class counts M_c (Float64, summing to 2^n over visited classes), and
the per-class member sets used. Uses only `edges`; never touches the 2^n
state space.
"""
function polytime_homogeneous_distribution(n::Int, edges::Vector{Tuple{Int,Int}};
        samples_per_class::Int=10, K_subsets::Int=200, rng::AbstractRNG,
        wl_kwargs...)
    m = length(edges)
    lng = wang_landau(n, edges; rng, wl_kwargs...)
    members = sample_class_members(n, edges, lng; rng, samples_per_class)
    # rng is consumed serially per (class, member, d) job before threading
    jobs = Tuple{Int,Int}[]           # (ci, x)
    for ci in 1:(m + 1), x in members[ci]
        push!(jobs, (ci, x))
    end
    rngs = [Xoshiro(rand(rng, UInt64)) for _ in jobs]
    N = zeros(Float64, m + 1, n + 1, m + 1)
    profs = Vector{Matrix{Float64}}(undef, length(jobs))
    @threads for k in eachindex(jobs)
        (_, x) = jobs[k]
        P = zeros(Float64, n + 1, m + 1)
        for d in 0:n
            profile_row!(view(P, d + 1, :), x, n, edges, d, K_subsets, rngs[k])
        end
        @assert isapprox(sum(P), 2.0^n; rtol=1e-10) "profile of x=$x does not sum to 2^n"
        profs[k] = P
    end
    for (k, (ci, _)) in enumerate(jobs)
        @views N[ci, :, :] .+= profs[k]
    end
    for ci in 1:(m + 1)
        nm = length(members[ci])
        nm > 0 && (@views N[ci, :, :] ./= nm)
    end
    # counts: 2^n * softmax(ln g)
    lmax = maximum(filter(isfinite, lng))
    w = [isfinite(l) ? exp(l - lmax) : 0.0 for l in lng]
    counts = (2.0^n / sum(w)) .* w
    @assert isapprox(sum(counts), 2.0^n; rtol=1e-10)
    return N, counts, members
end

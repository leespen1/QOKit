# Batched FUR X-mixer using KernelAbstractions shared memory.
#
# Instead of one kernel launch per qubit, processes groups of NQ qubits in a
# single launch using a butterfly network in shared memory.  This matches the
# algorithm in QOKit's furx.cu (shared-memory path, NQ > 6).
#
# Included by JuliaQAOAKernelAbstractionsExt.jl
#
# KA CPU codegen note: the CPU backend splits the kernel body at each
# @synchronize into separate WorkgroupLoop blocks, each with its own @index
# scope.  So @index must be called fresh in each block — values from a
# previous block are NOT available after @synchronize.

using KernelAbstractions: @kernel, @index, @Const, @localmem, @synchronize,
    @groupsize, @uniform, get_backend, synchronize, KernelAbstractions


# ─── Shared-memory batched X-mixer kernel ─────────────────────────────────────

@kernel function _batched_x_mixer_kernel!(state, cosβ, sinβ,
                                           q_offset, ::Val{NQ}) where {NQ}
    # half_tile = 2^(NQ-1) = number of threads per workgroup
    # Each workgroup handles 2*half_tile = 2^NQ amplitudes.
    # Computed from NQ (a type parameter) so it's a compile-time constant for GPU.

    @uniform half_tile = 1 << (NQ - 1)
    shmem = @localmem eltype(state) (1 << NQ,)

    # ── Block 0: Load from global → shared ────────────────────────────
    i = @index(Local, Linear)   # 1..half_tile
    g = @index(Group, Linear)   # 1..n_workgroups

    stride = 1 << q_offset
    index_mask = stride - 1
    stride_mask = ~index_mask
    base = ((stride_mask & (g - 1)) << NQ) | (index_mask & (g - 1))

    @inbounds begin
        shmem[i] = state[base + (i - 1) * stride + 1]
        shmem[i + half_tile] = state[base + (i - 1 + half_tile) * stride + 1]
    end

    @synchronize()

    # ── Blocks 1..NQ: Butterfly rotations ─────────────────────────────
    # Each iteration needs a fresh @index call because @synchronize starts a
    # new WorkgroupLoop on the CPU backend.

    for q in 0:(NQ - 1)
        i2 = @index(Local, Linear)
        tid0 = i2 - 1  # 0-based thread index

        mask1 = (1 << q) - 1
        mask2 = (half_tile - 1) - mask1

        ia = ((tid0 & mask1) | ((tid0 & mask2) << 1)) + 1  # 1-based
        ib = ia + (1 << q)

        @inbounds begin
            va = shmem[ia]
            vb = shmem[ib]
            shmem[ia] = cosβ * va - im * sinβ * vb
            shmem[ib] = -im * sinβ * va + cosβ * vb
        end

        @synchronize()
    end

    # ── Final block: Store shared → global ────────────────────────────
    i3 = @index(Local, Linear)
    g3 = @index(Group, Linear)

    stride3 = 1 << q_offset
    index_mask3 = stride3 - 1
    stride_mask3 = ~index_mask3
    base3 = ((stride_mask3 & (g3 - 1)) << NQ) | (index_mask3 & (g3 - 1))

    @inbounds begin
        state[base3 + (i3 - 1) * stride3 + 1] = shmem[i3]
        state[base3 + (i3 - 1 + half_tile) * stride3 + 1] = shmem[i3 + half_tile]
    end
end


# ─── Dispatch function ────────────────────────────────────────────────────────

"""
    gpu_apply_x_mixer_batched!(state, β, n; group_size=10)

Apply the X-mixer e^{-iβΣ_j X_j} using batched shared-memory FUR on GPU.

Groups `group_size` consecutive qubits per kernel launch (max 11 due to shared
memory limits).  Each launch processes 2^group_size amplitudes per workgroup
using a butterfly network in shared memory, reducing kernel launches from n to
ceil(n / group_size) and eliminating redundant global memory traffic.

Falls back to the single-qubit kernel for remainder groups of 1 qubit.
"""
function JuliaQAOA.gpu_apply_x_mixer_batched!(
    state::AbstractVector{<:Complex}, β::Real, n::Int; group_size::Int=10
)
    @assert 1 <= group_size <= 11 "group_size must be between 1 and 11"

    backend = get_backend(state)
    T = real(eltype(state))
    cosβ = T(cos(β))
    sinβ = T(sin(β))
    n_states = length(state)

    last_group = n % group_size
    full_groups = n - last_group

    # Process full groups
    for q_offset in 0:group_size:(full_groups - 1)
        _launch_batched_kernel!(backend, state, cosβ, sinβ, group_size, q_offset, n_states)
    end

    # Process remainder group
    if last_group > 0
        q_offset = full_groups
        if last_group == 1
            bit_mask = 1 << q_offset
            _x_mixer_qubit_kernel!(backend)(state, cosβ, sinβ, bit_mask,
                                            ndrange=n_states ÷ 2)
            synchronize(backend)
        else
            _launch_batched_kernel!(backend, state, cosβ, sinβ, last_group, q_offset, n_states)
        end
    end
end


function _launch_batched_kernel!(backend, state, cosβ, sinβ, nq, q_offset, n_states)
    half_tile = 1 << (nq - 1)
    n_workgroups = n_states >> nq
    ndrange = n_workgroups * half_tile

    _batched_x_mixer_kernel!(backend, half_tile)(
        state, cosβ, sinβ, q_offset, Val(nq), ndrange=ndrange
    )
    synchronize(backend)
end

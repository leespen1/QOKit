# Warp-shuffle FUR X-mixer using CUDA.jl directly.
#
# For small qubit groups (NQ <= 6, so <= 32 threads per sub-warp), uses warp
# shuffles instead of shared memory for intra-warp communication.  This matches
# the algorithm in QOKit's furx.cu (warp_furx_kernel path, NQ <= 6).
#
# Not portable — CUDA only.  Included by JuliaQAOACUDAExt.jl

using CUDA


# ─── Single-qubit fallback kernel (CUDA) ──────────────────────────────────────

function _single_qubit_x_kernel!(state, cosβ, sinβ, bit_mask)
    tid = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x - Int32(1)
    n_pairs = length(state) ÷ Int32(2)
    if tid < n_pairs
        t = tid
        x = (t & (bit_mask - Int32(1))) | ((t & ~(bit_mask - Int32(1))) << Int32(1))
        y = x | bit_mask
        @inbounds begin
            ax = state[x + 1]
            ay = state[y + 1]
            state[x + 1] = cosβ * ax - im * sinβ * ay
            state[y + 1] = -im * sinβ * ax + cosβ * ay
        end
    end
    return nothing
end


# ─── Warp-shuffle X-mixer kernel ─────────────────────────────────────────────

# We need separate kernels per NQ because the sub-warp size is a compile-time
# property of the shuffle width.  Julia's @generated or manual dispatch via
# Val{NQ} handles this.

function _warp_furx_kernel!(state, cosβ, sinβ, ::Val{NQ}, q_offset) where {NQ}
    HALF_TILE = 1 << (NQ - 1)   # threads per sub-warp (sub-warp size)

    # Multiple sub-warps can share a CUDA block.  Compute which sub-warp this
    # thread belongs to and which tile (block_idx) it processes.
    global_tid = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x - Int32(1)
    block_idx = global_tid ÷ Int32(HALF_TILE)
    tid = global_tid % Int32(HALF_TILE)

    # Compute global base offset (same as furx.cu get_offset)
    stride = Int32(1) << q_offset
    index_mask = stride - Int32(1)
    stride_mask = ~index_mask
    base = ((stride_mask & block_idx) << NQ) | (index_mask & block_idx)

    # Load 2 elements into registers
    load_offset = base + tid * Int32(2) * stride
    @inbounds v0 = state[load_offset + 1]
    @inbounds v1 = state[load_offset + stride + 1]

    # First rotation (on the pair loaded by this thread)
    new_v0 = cosβ * v0 - im * sinβ * v1
    new_v1 = -im * sinβ * v0 + cosβ * v1
    v0 = new_v0
    v1 = new_v1

    # Butterfly via warp shuffles for remaining NQ-1 qubits
    for q in Int32(0):Int32(NQ - 2)
        warp_stride = Int32(1) << q
        positive = (tid & warp_stride) == Int32(0)

        # Exchange the "outgoing" value with partner thread
        lane_idx = positive ? (tid + warp_stride) : (tid - warp_stride)
        # shfl_sync with width = HALF_TILE (sub-warp shuffle)
        if positive
            v1_re = CUDA.shfl_sync(0xFFFFFFFF, real(v1), lane_idx + Int32(1), Int32(HALF_TILE))
            v1_im = CUDA.shfl_sync(0xFFFFFFFF, imag(v1), lane_idx + Int32(1), Int32(HALF_TILE))
            v1 = Complex(v1_re, v1_im)
        else
            v0_re = CUDA.shfl_sync(0xFFFFFFFF, real(v0), lane_idx + Int32(1), Int32(HALF_TILE))
            v0_im = CUDA.shfl_sync(0xFFFFFFFF, imag(v0), lane_idx + Int32(1), Int32(HALF_TILE))
            v0 = Complex(v0_re, v0_im)
        end

        new_v0 = cosβ * v0 - im * sinβ * v1
        new_v1 = -im * sinβ * v0 + cosβ * v1
        v0 = new_v0
        v1 = new_v1
    end

    # Store results back.  After the butterfly, thread tid owns elements at
    # positions tid and tid + HALF_TILE within the tile.
    @inbounds state[base + tid * stride + 1] = v0
    @inbounds state[base + (tid + Int32(HALF_TILE)) * stride + 1] = v1

    return nothing
end


# ─── Dispatch function ────────────────────────────────────────────────────────

"""
    gpu_apply_x_mixer_warp!(state::CuArray, β, n; group_size=6)

Apply the X-mixer using warp-shuffle FUR on GPU (CUDA only).

For groups of up to 6 qubits, uses warp shuffles for intra-thread communication
instead of shared memory.  Each sub-warp of 2^(NQ-1) threads processes one tile
of 2^NQ amplitudes entirely in registers.

For remainder groups of 1 qubit, falls back to the single-qubit kernel.
"""
function JuliaQAOA.gpu_apply_x_mixer_warp!(
    state::CuArray{<:Complex}, β::Real, n::Int; group_size::Int=6
)
    @assert 1 <= group_size <= 6 "warp group_size must be between 1 and 6"

    T = real(eltype(state))
    cosβ = T(cos(β))
    sinβ = T(sin(β))
    n_states = length(state)

    last_group = n % group_size
    full_groups = n - last_group

    for q_offset in 0:group_size:(full_groups - 1)
        _launch_warp_kernel!(state, cosβ, sinβ, Val(group_size), q_offset, n_states)
    end

    if last_group > 0
        q_offset = full_groups
        if last_group == 1
            # Single qubit: use a simple CUDA kernel directly
            n_pairs = n_states ÷ 2
            block_size = min(256, n_pairs)
            n_blocks = cld(n_pairs, block_size)
            @cuda threads=block_size blocks=n_blocks _single_qubit_x_kernel!(
                state, cosβ, sinβ, Int32(1 << q_offset)
            )
            CUDA.synchronize()
        else
            _launch_warp_kernel!(state, cosβ, sinβ, Val(last_group), q_offset, n_states)
        end
    end
end


function _launch_warp_kernel!(state, cosβ, sinβ, ::Val{NQ}, q_offset, n_states) where {NQ}
    half_tile = 1 << (NQ - 1)
    n_tiles = n_states >> NQ
    total_threads = n_tiles * half_tile
    block_size = min(256, total_threads)
    n_blocks = cld(total_threads, block_size)

    @cuda threads=block_size blocks=n_blocks _warp_furx_kernel!(
        state, cosβ, sinβ, Val(NQ), q_offset
    )
    CUDA.synchronize()
end

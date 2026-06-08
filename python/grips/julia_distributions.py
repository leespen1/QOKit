"""
Helper module for computing N(c', d, c) distributions using Julia.

This module provides Python wrappers around the Julia distribution functions
in JuliaQAOA, handling the index ordering differences between Julia (column-major)
and Python (row-major).
"""

import numpy as np
import os

# Lazy initialization of Julia - only initialize when actually needed
_jl = None
_julia_initialized = False


def _init_julia():
    """Initialize Julia and load JuliaQAOA module."""
    global _jl, _julia_initialized
    if _julia_initialized:
        return _jl

    # Set thread count before Julia starts (must be set before first import)
    if "JULIA_NUM_THREADS" not in os.environ:
        os.environ["JULIA_NUM_THREADS"] = "auto"

    from juliacall import Main as jl

    _jl = jl

    grips_dir = os.path.dirname(os.path.abspath(__file__))
    # grips/ now lives at python/grips/; the JuliaQAOA project is the repo root,
    # two levels up (python/grips -> python -> repo root).
    julia_project_dir = os.path.normpath(os.path.join(grips_dir, "..", ".."))

    jl.seval(
        f"""
    using Pkg
    Pkg.activate("{julia_project_dir}")
    try
        using JuliaQAOA
        using CUDA
    catch e
        println("Encountered error during 'using JuliaQAOA', instantiating environment at {julia_project_dir}")
        Pkg.instantiate()
        using JuliaQAOA
        using CUDA
    end
    """
    )

    _julia_initialized = True
    return _jl


def has_cuda_gpu() -> bool:
    """Check if a CUDA GPU is available."""
    jl = _init_julia()
    return bool(jl.seval("CUDA.has_cuda_gpu()"))


def get_homogeneous_distribution(costs: np.ndarray, num_edges: int, num_vertices: int, use_gpu: bool = False, max_num_edges: int = 0) -> np.ndarray:
    """
    Compute the homogeneous distribution N(c'; d, c) from costs using Julia.

    Args:
        costs: 1D array of costs for each bitstring (length 2^num_vertices)
        num_edges: Number of edges in the graph
        num_vertices: Number of vertices in the graph
        use_gpu: Whether to use GPU acceleration
        max_num_edges: Optional, allocate extra space for compatibility

    Returns:
        3D numpy array of shape (num_costs, num_distances, num_costs)
        representing N(c'; d, c) where indices are [c', d, c]

    Note:
        Julia uses 1-based indexing internally but returns 0-based indexed results.
        The returned array is in Python's row-major order.
    """
    jl = _init_julia()

    # Convert costs to Float64 for Julia
    costs_float = costs.astype(np.float64)

    if use_gpu and has_cuda_gpu():
        # Use GPU version
        N_dist_julia = jl.gpu_get_homogeneous_distribution_from_costs_direct(costs_float, num_edges, num_vertices, max_num_edges=max_num_edges)
        # Convert CuArray to Array, then to numpy
        N_dist = np.array(jl.Array(N_dist_julia))
    else:
        # Use CPU version
        N_dist_julia = jl.get_homogeneous_distribution_from_costs_direct(costs_float, num_edges, num_vertices, max_num_edges=max_num_edges)
        N_dist = np.array(N_dist_julia)

    # Julia arrays are column-major, numpy is row-major
    # Julia returns shape (num_costs, num_distances, num_costs)
    # We need to handle the memory layout difference
    # The to_numpy() or np.array() already handles this, but we should verify shape
    # Expected shape: (num_costs, num_distances, num_costs)

    return N_dist


def get_real_distribution(costs: np.ndarray, num_edges: int, num_vertices: int, use_gpu: bool = False, max_num_edges: int = 0) -> np.ndarray:
    """
    Compute the real distribution n(x; d, c) from costs using Julia.

    Args:
        costs: 1D array of costs for each bitstring (length 2^num_vertices)
        num_edges: Number of edges in the graph
        num_vertices: Number of vertices in the graph
        use_gpu: Whether to use GPU acceleration
        max_num_edges: Optional, allocate extra space for compatibility

    Returns:
        3D numpy array of shape (num_bitstrings, num_distances, num_costs)
        representing n(x; d, c) where indices are [x, d, c]

    Warning:
        This is memory intensive for large num_vertices (O(2^n * n * m) storage).
    """
    jl = _init_julia()

    # Convert costs to Float64 for Julia
    costs_float = costs.astype(np.float64)

    if use_gpu and has_cuda_gpu():
        # Use GPU version
        n_dist_julia = jl.gpu_get_real_distribution_from_costs(costs_float, num_edges, num_vertices, max_num_edges=max_num_edges)
        # Convert CuArray to Array, then to numpy
        n_dist = np.array(jl.Array(n_dist_julia))
    else:
        # Use CPU version
        n_dist_julia = jl.get_real_distribution_from_costs(costs_float, num_edges, num_vertices, max_num_edges=max_num_edges)
        n_dist = np.array(n_dist_julia)

    return n_dist


def pad_distribution(dist: np.ndarray, target_shape: tuple) -> np.ndarray:
    """Pad a distribution array with zeros to match target shape."""
    if dist.shape == target_shape:
        return dist

    result = np.zeros(target_shape, dtype=dist.dtype)
    slices = tuple(slice(0, s) for s in dist.shape)
    result[slices] = dist
    return result


class HomogeneousDistributionAccumulator:
    """
    Accumulates N(c'; d, c) distributions across multiple graphs and computes
    their average. Memory efficient - only stores the running sum and count.
    """

    def __init__(self, max_num_edges: int = 0):
        """
        Args:
            max_num_edges: Pre-allocate for this many edges to avoid resizing.
                           If 0, will grow as needed.
        """
        self.max_num_edges = max_num_edges
        self._sum = None  # Running sum of distributions
        self._count = 0  # Number of distributions added
        self._shape = None

    def add(self, costs: np.ndarray, num_edges: int, num_vertices: int, use_gpu: bool = False):
        """
        Compute N(c'; d, c) for the given costs and add to accumulator.

        Args:
            costs: 1D array of costs for each bitstring
            num_edges: Number of edges in the graph
            num_vertices: Number of vertices in the graph
            use_gpu: Whether to use GPU acceleration
        """
        # Compute distribution for this graph
        N_dist = get_homogeneous_distribution(costs, num_edges, num_vertices, use_gpu=use_gpu, max_num_edges=self.max_num_edges)

        if self._sum is None:
            # First distribution - initialize
            self._sum = N_dist.astype(np.float64)
            self._shape = N_dist.shape
        else:
            # Pad if needed to match shapes
            if N_dist.shape != self._shape:
                # Determine new shape (max of each dimension)
                new_shape = tuple(max(s1, s2) for s1, s2 in zip(self._shape, N_dist.shape))
                self._sum = pad_distribution(self._sum, new_shape)
                N_dist = pad_distribution(N_dist, new_shape)
                self._shape = new_shape

            self._sum += N_dist

        self._count += 1

    @property
    def count(self) -> int:
        """Number of distributions accumulated."""
        return self._count

    @property
    def mean(self) -> np.ndarray:
        """Compute the average distribution."""
        if self._count == 0:
            raise ValueError("No distributions have been added yet")
        return self._sum / self._count

    @property
    def sum(self) -> np.ndarray:
        """Return the sum of all distributions."""
        if self._count == 0:
            raise ValueError("No distributions have been added yet")
        return self._sum.copy()

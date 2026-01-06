# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

QOKit (Quantum Optimization Toolkit) is a Python library for simulating and
benchmarking the Quantum Approximate Optimization Algorithm (QAOA). The `grips/`
directory contains G-RIPS 2024 Sendai research on improving QAOA
parameter-setting heuristics based on the paper "Parameter-setting heuristic for
the quantum alternating operator ansatz." The users of Claude Code in this
project are the Sendai researchers. 

**Branches**: Main branch is `grips` (stable). Active development on `gpu` branch.

**Virtual Environment**: At least on Spencer's computer, a virtual environment
(for the python portion of the project) is provided in `../qokitvenv`, which has
the necessary python packages installed.

## Common Commands

```bash
# Install (development mode)
pip install -e .

# Install with GPU support (CUDA 12.x)
pip install -e .[GPU-CUDA12]

# Python-only install (if C compilation fails)
QOKIT_PYTHON_ONLY=1 pip install -e .

# Run all tests with coverage
pytest --cov=qokit --cov-fail-under=75 -rs tests

# Run a single test file
pytest tests/test_qaoa_objective_maxcut.py

# Run GRIPS-specific tests
pytest grips_tests/

# Check formatting
black --check .

# Format code
black .
```

## Architecture

### Simulation Pipeline

The core QAOA simulation uses Fast Unitary Rotation (FUR) algorithm with multiple backends:
- **GPU (CUDA + Numba)**: Fastest, requires `cupy` - `qokit/fur/nbcuda/`
- **C-compiled**: Fast CPU - `qokit/fur/c/`
- **Python fallback**: Reference implementation - `qokit/fur/python/`

`qokit.fur.choose_simulator()` auto-selects the fastest available backend.

### GRIPS Proxy System

The proxy system approximates QAOA state evolution without full quantum simulation:

1. **Real Distribution** (`grips/real_distribution.py`): Computes `n(x; d, c)` - count of bitstrings at Hamming distance `d` with cost `c` from bitstring `x`

2. **Homogeneous Distribution**: Cost-averaged `N(c'; d, c)` enables parameter prediction across graph instances

3. **Proxy Classes** (all in `grips/`):
   - `PaperProxy` (`paper_proxy.py`): Original paper's binomial/multinomial approach
   - `NormalProxy` (`normal_proxy.py`): Multivariate normal approximation
   - `TriangleProxy` (`triangle_proxy.py`): GRIPS contribution - simplified triangle distribution

4. **Interface** (`grips/QAOA_proxy_interface.py`): Unified proxy API with Python (Numba JIT) and Julia backends

### Key Module Relationships

```
grips/QAOA_simulator.py     - Main simulation interface (QAOA_run, get_simulator)
grips/QAOA_proxy_interface.py - Proxy algorithm entry point
  └── paper_proxy.py / normal_proxy.py / triangle_proxy.py
grips/real_distribution.py  - Statistical distribution computation
grips/sendai_opt.py        - Parameter optimization (fit_proxy_to_real)
qokit/fur/__init__.py      - Simulator backend selection
```

### Julia Backend

Julia implementations in `julia/src/` provide high-performance alternatives. Configure with `USE_JULIA=True` in proxy files. Setup via `grips/setup_juliacall.py`.

## Code Style

- **Formatter**: Black with `line-length=160`
- **Python**: 3.10, 3.11
- **Type checking**: pyright
- **License headers**: Apache 2.0, SPDX headers required (checked by `addheader`)

## Key Concepts

- **QAOA**: Quantum Approximate Optimization Algorithm - hybrid quantum-classical algorithm
- **MaxCut**: Graph partitioning problem, primary benchmark in this codebase
- **LABS**: Low Autocorrelation Binary Sequences problem
- **Hamming Distance**: Number of differing bits between two bitstrings
- **FUR**: Fast Unitary Rotation - efficient state evolution algorithm

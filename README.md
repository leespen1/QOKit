# QOKit — Julia-first QAOA toolkit (G-RIPS Sendai)

![Julia tests](https://github.com/jpmorganchase/QOKit/actions/workflows/julia-test.yml/badge.svg)
![Python tests](https://github.com/jpmorganchase/QOKit/actions/workflows/qokit-package.yml/badge.svg)
[![arXiv](https://img.shields.io/badge/arXiv-2309.04841-b31b1b.svg?style=plastic)](https://arxiv.org/abs/2309.04841)

This repository began as the Python **Quantum Optimization Toolkit** (QOKit) and
has been restructured to be **Julia-first**. The repo root is now the
[`JuliaQAOA`](src/JuliaQAOA.jl) package, laid out as a
[DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/stable/) scientific
project. The original Python QOKit library and the G-RIPS research code now live
under [`python/`](python/).

## Layout

```
.                      JuliaQAOA package + DrWatson project
├── src/  ext/  test/  Julia package (module, GPU extensions, tests)
├── scripts/           Julia scripts: paper_figures/, benchmark/, examples/
├── plots/             Figure output (DrWatson plotsdir())
├── data/              Datasets and data-generation scripts
├── papers/            Paper sources (OverleafPaper/, References/)
├── notebooks/         Julia notebooks
├── research_log/      Experiment log (IEEE Quantum Week paper)
└── python/            All Python code (QOKit + grips), with its own README
    ├── qokit/  grips/  tests/  grips_tests/  grips_examples/
    ├── pyproject.toml  setup.py
    └── qokitvenv/      Python virtualenv lives here (gitignored)
```

## Julia quickstart

This is a DrWatson project. From the repo root:

```julia
using DrWatson
@quickactivate "JuliaQAOA"   # activates this project; sets projectdir(), datadir(), plotsdir(), ...
using JuliaQAOA
```

Run the tests, benchmarks, and paper figures (all from the repo root):

```bash
julia --project -e 'using Pkg; Pkg.test()'                 # full test suite
julia --project scripts/benchmark/benchmark_QAOA.jl        # a benchmark
bash  scripts/paper_figures/run_all.sh                     # all paper figures -> plots/
```

GPU acceleration loads automatically when a GPU package is available in your
environment: `using CUDA, JuliaQAOA` activates the CUDA/KernelAbstractions
extensions. See [`CLAUDE.md`](CLAUDE.md) for the full command reference and
architecture notes.

## Python package

All Python code is under [`python/`](python/) and retains its own
[README](python/README.md). Create the virtualenv **inside `python/`** and
install in development mode from the repo root:

```bash
python -m venv python/qokitvenv
source python/qokitvenv/bin/activate
pip install -U pip
pip install -e python              # add 'python[GPU-CUDA12]' for GPU, or
                                   # QOKIT_PYTHON_ONLY=1 to skip the C build
pytest python/tests                # run the QOKit test suite
```

The Python ↔ Julia bridge (`USE_JULIA=True` in `python/grips/`) activates the
repo-root `JuliaQAOA` project automatically via `juliacall`.

## Cite

For the simulators and other software tools, please cite:

```
@inproceedings{Lykov2023,
  series = {SC-W 2023},
  title = {Fast Simulation of High-Depth QAOA Circuits},
  url = {http://dx.doi.org/10.1145/3624062.3624216},
  DOI = {10.1145/3624062.3624216},
  booktitle = {Proceedings of the SC ’23 Workshops of The International Conference on High Performance Computing,  Network,  Storage,  and Analysis},
  publisher = {ACM},
  author = {Lykov,  Danylo and Shaydulin,  Ruslan and Sun,  Yue and Alexeev,  Yuri and Pistoia,  Marco},
  year = {2023},
  month = nov,
  collection = {SC-W 2023}
}
```

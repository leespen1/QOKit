#=
find_python.jl — Auto-detect the qokitvenv Python for PythonCall.

Include this file BEFORE `using JuliaQAOA` or `using PythonCall` to
automatically configure PythonCall to use the project's virtual environment.

    include("find_python.jl")  # or adjust path as needed
    using JuliaQAOA

The venv is expected at ../qokitvenv/ relative to the QOKit repo root.
If JULIA_PYTHONCALL_EXE is already set, this is a no-op.
=#

if !haskey(ENV, "JULIA_PYTHONCALL_EXE")
    let
        # @__DIR__ is julia/ (where this file and Project.toml live)
        # The venv is at ../../qokitvenv/ relative to julia/
        venv_python = normpath(joinpath(@__DIR__, "..", "..", "qokitvenv", "bin", "python"))
        if isfile(venv_python)
            ENV["JULIA_PYTHONCALL_EXE"] = venv_python
            ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
        end
    end
end

"""
GRIPS: A library for running QAOA and QAOA proxies
This file imports functions from all src files in the grips directory

With this in place, the user should be able to do "from grips import function_name", instead
of needing to do "from grips.file_name import function_name".
"""
import os

# Make the functions written in Julia available (call with `jl.function_name`)
# (by having this in __init__.py, should be able to simply do from grips import jl)
from juliacall import Main as jl
grips_dir = os.path.dirname(os.path.abspath(__file__))
julia_project_dir = os.path.normpath(os.path.join(grips_dir, "../julia"))
jl.seval(f'''
using Pkg
Pkg.activate("{julia_project_dir}")
try
    using JuliaQAOA
catch e
    if isa(e, ArgumentError)
        println("Encountered error during 'using JuliaQAOA', instantiating ...")
        Pkg.instantiate()
        using JuliaQAOA
    else
        rethrow(e)
    end
end
''')

from .normal_proxy import *
from .paper_proxy import *
from .plot_utils import *
from .QAOA_proxy_interface import *
from .QAOA_simulator import *
from .real_distribution import *
from .scipy_additional_optimizers import *
from .triangle_proxy import *
from .solve_maxcut_exact import *
from .sendai_opt import *

#from .stddev_mean_heatmaps import *

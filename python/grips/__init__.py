"""
GRIPS: A library for running QAOA and QAOA proxies
This file imports functions from all src files in the grips directory

With this in place, the user should be able to do "from grips import function_name", instead
of needing to do "from grips.file_name import function_name".
"""
import os

from .normal_proxy import *
from .paper_proxy import *
from .plot_utils import *
#from .QAOA_proxy_interface import * # Removing this so that julicall only needs to be loaded when using the proxy interface
from .QAOA_simulator import *
from .real_distribution import *
from .scipy_additional_optimizers import *
from .triangle_proxy import *
from .solve_maxcut_exact import *
from .sendai_opt import *
from .utils import *

#from .stddev_mean_heatmaps import *

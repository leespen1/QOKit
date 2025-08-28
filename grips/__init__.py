"""
GRIPS: A library for running QAOA and QAOA proxies
This file imports functions from all src files in the grips directory

With this in place, the user should be able to do "from grips import function_name", instead
of needing to do "from grips.file_name import function_name".
"""

from .normal_proxy import *
from .paper_proxy import *
from .plot_utils import *
from .QAOA_proxy_interface import *
from .QAOA_simulator import *
from .real_distribution import *
from .scipy_additional_optimizers import *
from .triangle_proxy import *

#from .stddev_mean_heatmaps import *
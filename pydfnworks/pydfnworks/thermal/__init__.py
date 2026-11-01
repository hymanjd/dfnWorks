"""
pydfnworks.thermal
==================

Heat extraction from fractured rock for hydraulic-fracture spacing
optimization. Two backends share one interface,

    heat_extraction(spacing, m_dot, T_res, T_inj, rock_props, fluid_props,
                    n_fractures, fracture_height, t_eval) -> Q(t) [W]

mpf_gringarten.heat_extraction  analytic Gringarten (1975) surrogate, ~100 us
dfn_heat.heat_extraction        DFN graph-transport backend (not implemented)

optimize_spacing() evaluates candidate spacings against any backend and
objective; cumulative_heat() and time_to_drawdown() are provided.
"""

from pydfnworks.thermal import mpf_gringarten, dfn_heat
from pydfnworks.thermal.common import GRANITE_PROPS, WATER_PROPS, SECONDS_PER_YEAR
from pydfnworks.thermal.mpf_gringarten import heat_extraction, dimensionless_parameters
from pydfnworks.thermal.dfn_heat import thermal_to_tdrw_parameters
from pydfnworks.thermal.optimize_spacing import (optimize_spacing, cumulative_heat,
                                                 time_to_drawdown, plot_optimization)

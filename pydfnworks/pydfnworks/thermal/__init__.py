"""
pydfnworks.thermal
==================

DFN-based heat extraction. Graph particle tracking supplies each path's
advective time and fracture-rock exchange integral; an analytic slab
kernel turns those into the produced temperature T_out(t) and heat rate
Q(t) for a reservoir with parallel hydraulic fractures at spacing D in a
natural fracture network. dataset.py builds such networks and writes
per-sample HDF5 files for downstream (ML) use.
"""

from pydfnworks.thermal.dfn_heat import (heat_extraction, thermal_breakthrough,
                                         slab_drawdown, thermal_to_tdrw_parameters,
                                         GRANITE_PROPS, WATER_PROPS, SECONDS_PER_YEAR)
from pydfnworks.thermal.dataset import (hydraulic_fracture_positions, add_hydraulic_fractures,
                                        graph_flow_for_rate, write_sample, append_manifest)

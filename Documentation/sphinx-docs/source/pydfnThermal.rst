.. _dfnWorks-python-chapter-thermal:

pydfnworks: thermal
========================================

DFN-based heat extraction. Graph particle tracking supplies each path's
advective time and fracture-rock exchange integral; an analytic slab kernel
turns those into the produced temperature and heat extraction rate for
parallel hydraulic fractures at a given spacing in a natural fracture
network.

Heat extraction from graph transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: pydfnworks.thermal.dfn_heat
    :members: heat_extraction, thermal_breakthrough, slab_drawdown, thermal_to_tdrw_parameters

Dataset helpers
^^^^^^^^^^^^^^^^
.. automodule:: pydfnworks.thermal.dataset
    :members: add_hydraulic_fractures, hydraulic_fracture_positions, write_sample, append_manifest

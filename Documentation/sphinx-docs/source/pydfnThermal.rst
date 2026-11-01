.. _dfnWorks-python-chapter-thermal:

pydfnworks: thermal
========================================

Heat extraction from fractured rock for hydraulic-fracture spacing
optimization. Two backends share one ``heat_extraction()`` interface so an
optimizer can call either interchangeably: the analytic Gringarten,
Witherspoon and Ohnishi (1975) multiple parallel fractures surrogate and
the DFN graph-transport backend (interface only in this version).

Gringarten multiple parallel fractures surrogate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: pydfnworks.thermal.mpf_gringarten
    :members: heat_extraction, dimensionless_parameters, dimensionless_drawdown, single_fracture_drawdown, laplace_drawdown

DFN-based heat transport
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: pydfnworks.thermal.dfn_heat
    :members: heat_extraction, thermal_to_tdrw_parameters

Spacing optimizer
^^^^^^^^^^^^^^^^^^
.. automodule:: pydfnworks.thermal.optimize_spacing
    :members: optimize_spacing, cumulative_heat, time_to_drawdown, plot_optimization

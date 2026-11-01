"""
.. module:: common.py
   :synopsis: Property constants and argument validation shared by the
              pydfnworks.thermal heat-extraction backends
.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
"""

import numpy as np

# Representative properties: Utah FORGE granitoid and liquid water
GRANITE_PROPS = {"rho": 2750.0, "c": 790.0, "k": 3.05}
WATER_PROPS = {"rho": 1000.0, "c": 4200.0}
SECONDS_PER_YEAR = 365.25 * 24 * 3600.0


def check_inputs(spacing, m_dot, T_res, T_inj, rock_props, fluid_props,
                 n_fractures, fracture_height, t_eval):
    """ Validate the positional arguments of heat_extraction()

    Parameters
    ----------
        See pydfnworks.thermal.mpf_gringarten.heat_extraction

    Returns
    -------
        t : numpy array
            t_eval as a 1D float array [s]

    Notes
    -----
        Raises ValueError (not print_log 'error', which exits) so the
        backends can sit inside an optimizer loop. Both backends call
        this so they reject the same inputs.
    """
    for name, value in (("spacing", spacing), ("m_dot", m_dot),
                        ("fracture_height", fracture_height)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive, got {value}")
    if not (np.isfinite(T_res) and np.isfinite(T_inj)) or T_inj >= T_res:
        raise ValueError(f"need finite T_inj < T_res, got {T_inj}, {T_res}")
    for name, props, keys in (("rock_props", rock_props, ("rho", "c", "k")),
                              ("fluid_props", fluid_props, ("rho", "c"))):
        for key in keys:
            if key not in props or not props[key] > 0:
                raise ValueError(f"{name}['{key}'] must be a positive number")
    if int(n_fractures) != n_fractures or n_fractures < 1:
        raise ValueError(f"n_fractures must be a positive integer, got {n_fractures}")
    t = np.atleast_1d(np.asarray(t_eval, dtype=float)).ravel()
    if t.size == 0 or not np.all(np.isfinite(t)) or np.any(t < 0):
        raise ValueError("t_eval must be non-empty, finite and non-negative [s]")
    return t

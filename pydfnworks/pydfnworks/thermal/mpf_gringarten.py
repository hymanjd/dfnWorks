"""
.. module:: mpf_gringarten.py
   :synopsis: Analytic heat-extraction surrogate: the Gringarten,
              Witherspoon and Ohnishi (1975) multiple parallel fractures
              (MPF) model, inverted numerically from the Laplace domain
.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>

Notes
-----
    n identical planar fractures of height H (flow direction) and width w
    are spaced a distance D apart in an infinite rock mass. Water enters
    every fracture at T_inj at the same rate; rock heat reaches the
    fracture by conduction normal to its plane, so each fracture draws on
    a slab of half-thickness D/2. The Laplace transform of the outlet
    drawdown T_wD = (T_res - T_out) / (T_res - T_inj) is (Eq. A17)

        T_wD*(s) = (1/s) exp[-sqrt(s) tanh(beta sqrt(s))],   t_D = t / t_scale

        Q_f     = m_dot / (rho_f n)                          per-fracture flow
        t_scale = 4 k (rho c)_r (w H)^2 / ((rho c)_f Q_f)^2  time scale
        beta    = (rho c)_f Q_f D / (4 k w H)                spacing parameter

    beta * t_scale is the piston time of one slab (heat stored in D w H of
    rock divided by the heat-capacity flow of the water). beta -> inf
    gives the single-fracture semi-infinite solution erfc(1/(2 sqrt t_D))
    (Lauwerier 1955; Bodvarsson and Tsang 1982); beta -> 0 gives a step
    at t_D = beta.

    The Laplace expression and parameter grouping follow GEOPHIRES-X
    (Beckers and McCabe 2019, MPFReservoir.py, MIT). Inversion is the
    double-precision Gaver-Stehfest routine shared with the TDRW models,
    18 terms (chosen against a 40-digit de Hoog reference: error in T_wD
    ~3e-5 for beta >= 1, ~6e-3 at beta = 0.1; more terms lose accuracy in
    double precision). Below beta ~ 0.05 the step at t_D = beta shows
    Stehfest ringing of a few percent; integrated heat stays within 1e-3.
    A 361-point curve costs about 100 microseconds.

    References
        Gringarten, Witherspoon, Ohnishi (1975). J. Geophys. Res. 80(8), 1120.
        Beckers, McCabe (2019). Geothermal Energy 7, 5.
        Bodvarsson, Tsang (1982). J. Geophys. Res. 87(B2), 1031.
"""

import numpy as np
from scipy.special import erfc

from pydfnworks.dfnGraph.transport.tdrw.stehfest import stehfest_coefficients
from pydfnworks.thermal.common import check_inputs

STEHFEST_TERMS = 18


def laplace_drawdown(s, beta):
    """ Laplace transform of the MPF outlet drawdown, (1/s) exp(-sqrt(s) tanh(beta sqrt(s)))

    Parameters
    ----------
        s : float or numpy array
            Laplace variable conjugate to dimensionless time

        beta : float
            dimensionless spacing parameter

    Returns
    -------
        F : float or numpy array
    """
    q = np.sqrt(s)
    return np.exp(-q * np.tanh(beta * q)) / s


def single_fracture_drawdown(t_D):
    """ Wide-spacing (beta -> inf) limit, erfc(1 / (2 sqrt(t_D)))

    Parameters
    ----------
        t_D : array_like
            dimensionless time >= 0

    Returns
    -------
        T_wD : numpy array
            drawdown in [0, 1]
    """
    t_D = np.atleast_1d(np.asarray(t_D, dtype=float))
    out = np.zeros_like(t_D)
    pos = t_D > 0
    out[pos] = erfc(0.5 / np.sqrt(t_D[pos]))
    return out


def dimensionless_parameters(spacing, m_dot, rock_props, fluid_props,
                             n_fractures, fracture_height, fracture_width=None):
    """ Dimensionless groups of the MPF model

    Parameters
    ----------
        spacing : float
            fracture spacing D [m]

        m_dot : float
            total mass flow rate [kg/s], shared equally by the fractures

        rock_props : dict
            {"rho": kg/m^3, "c": J/kg/K, "k": W/m/K}

        fluid_props : dict
            {"rho": kg/m^3, "c": J/kg/K}

        n_fractures : int
            number of fractures

        fracture_height : float
            fracture extent in the flow direction H [m]

        fracture_width : float or None
            extent perpendicular to flow w [m]; None uses H (square)

    Returns
    -------
        params : dict
            "beta" [-], "t_scale" [s], "t_piston" = beta * t_scale [s],
            "heat_in_place" = (rho c)_r n D w H [J/K]
    """
    w = fracture_height if fracture_width is None else fracture_width
    H = fracture_height
    rho_c_f = fluid_props["rho"] * fluid_props["c"]
    rho_c_r = rock_props["rho"] * rock_props["c"]
    k = rock_props["k"]
    Q_f = m_dot / (fluid_props["rho"] * n_fractures)
    t_scale = 4.0 * k * rho_c_r * (w * H)**2 / (rho_c_f * Q_f)**2
    beta = rho_c_f * Q_f * spacing / (4.0 * k * w * H)
    return {
        "beta": beta,
        "t_scale": t_scale,
        "t_piston": beta * t_scale,
        "heat_in_place": rho_c_r * n_fractures * spacing * w * H,
    }


def dimensionless_drawdown(t_D, beta, n_terms=STEHFEST_TERMS):
    """ Outlet drawdown T_wD(t_D; beta) by Gaver-Stehfest inversion

    f(t) ~ (ln 2 / t) sum_i V_i F(i ln 2 / t)

    Parameters
    ----------
        t_D : array_like
            dimensionless times >= 0

        beta : float
            dimensionless spacing parameter

        n_terms : int
            number of Stehfest coefficients (even). Default 18.

    Returns
    -------
        T_wD : numpy array
            drawdown clipped to [0, 1], zero at t_D = 0
    """
    t_D = np.atleast_1d(np.asarray(t_D, dtype=float))
    T_wD = np.zeros_like(t_D)
    pos = t_D > 0
    if np.any(pos):
        V = stehfest_coefficients(n_terms)
        a = np.log(2.0) / t_D[pos]
        s = np.outer(a, np.arange(1, n_terms + 1))
        T_wD[pos] = a * (laplace_drawdown(s, beta) @ V)
    return np.clip(T_wD, 0.0, 1.0)


def heat_extraction(spacing, m_dot, T_res, T_inj, rock_props, fluid_props,
                    n_fractures, fracture_height, t_eval, *,
                    fracture_width=None):
    """ Heat extraction rate Q(t) of the Gringarten MPF surrogate

    Parameters
    ----------
        spacing : float
            uniform fracture spacing [m]

        m_dot : float
            total mass flow rate [kg/s]

        T_res : float
            initial reservoir temperature [deg C]

        T_inj : float
            injection temperature [deg C]

        rock_props : dict
            {"rho": kg/m^3, "c": J/kg/K, "k": W/m/K}

        fluid_props : dict
            {"rho": kg/m^3, "c": J/kg/K}

        n_fractures : int
            number of fractures sharing the flow

        fracture_height : float
            fracture extent in the flow direction [m]

        t_eval : array_like
            times at which to evaluate Q [s]

        fracture_width : float or None
            keyword-only; extent perpendicular to flow [m], default = height

    Returns
    -------
        Q : numpy array
            m_dot c_f (T_out - T_inj) at each time in t_eval [W]

    Notes
    -----
        The nine positional arguments are the interface shared with
        pydfnworks.thermal.dfn_heat.heat_extraction(). The produced
        temperature is T_inj + Q / (m_dot c_f).
    """
    t = check_inputs(spacing, m_dot, T_res, T_inj, rock_props, fluid_props,
                     n_fractures, fracture_height, t_eval)
    p = dimensionless_parameters(spacing, m_dot, rock_props, fluid_props,
                                 n_fractures, fracture_height, fracture_width)
    T_wD = dimensionless_drawdown(t / p["t_scale"], p["beta"])
    return m_dot * fluid_props["c"] * (T_res - T_inj) * (1.0 - T_wD)

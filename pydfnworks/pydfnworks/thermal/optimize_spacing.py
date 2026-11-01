"""
.. module:: optimize_spacing.py
   :synopsis: Fracture spacing optimizer over a heat_extraction() backend
.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>

Notes
-----
    The optimizer evaluates every candidate spacing and keeps the best,
    which is all a one-dimensional, cheap (100 microsecond) and, for the
    surrogate, monotone objective needs; it also leaves the whole
    objective curve available. The backend is an argument, so the
    surrogate and the DFN backend are interchangeable here.
"""

import numpy as np

from pydfnworks.thermal import mpf_gringarten

_trapezoid = getattr(np, "trapezoid", None) or np.trapz


def cumulative_heat(t, Q):
    """ Objective: heat extracted over the window, integral of Q dt [J] """
    return float(_trapezoid(Q, t))


def time_to_drawdown(t, Q, fraction=0.9):
    """ Objective: first time Q falls to fraction * Q[0] [s]

    Linearly interpolated between samples; returns t[-1] if the threshold
    is not reached. Use functools.partial to change fraction.
    """
    threshold = fraction * Q[0]
    below = np.nonzero(Q <= threshold)[0]
    if below.size == 0:
        return float(t[-1])
    i = below[0]
    if i == 0 or Q[i - 1] == Q[i]:
        return float(t[i])
    return float(t[i - 1] + (Q[i - 1] - threshold) / (Q[i - 1] - Q[i]) * (t[i] - t[i - 1]))


def optimize_spacing(spacings, m_dot, T_res, T_inj, rock_props, fluid_props,
                     n_fractures, fracture_height, t_eval,
                     heat_extraction=mpf_gringarten.heat_extraction,
                     objective=cumulative_heat, **backend_kwargs):
    """ Find the candidate spacing that maximizes an objective of Q(t)

    Parameters
    ----------
        spacings : array_like
            candidate spacings [m], e.g. np.geomspace(2, 200, 40) or
            lateral_length / np.arange(2, 60) for integer stage counts

        m_dot, T_res, T_inj, rock_props, fluid_props, fracture_height, t_eval :
            passed to heat_extraction()

        n_fractures : int or callable
            fixed count, or n_fractures(spacing) (e.g. lateral / spacing)

        heat_extraction : callable
            backend with the pydfnworks.thermal signature. Default: the
            Gringarten surrogate.

        objective : callable
            objective(t, Q) -> float, maximized. Default cumulative_heat.

        backend_kwargs :
            forwarded to heat_extraction()

    Returns
    -------
        result : dict
            "spacing", "n_fractures", "objective", "Q" (at the best
            spacing), "t_eval", "spacings" and "objectives" (whole curve)
    """
    n_of = n_fractures if callable(n_fractures) else (lambda s: n_fractures)
    t = np.atleast_1d(np.asarray(t_eval, dtype=float))
    spacings = np.sort(np.asarray(spacings, dtype=float).ravel())
    counts = [int(n_of(s)) for s in spacings]
    curves = [heat_extraction(s, m_dot, T_res, T_inj, rock_props, fluid_props,
                              n, fracture_height, t, **backend_kwargs)
              for s, n in zip(spacings, counts)]
    objectives = np.array([objective(t, Q) for Q in curves])
    i = int(np.argmax(objectives))
    return {"spacing": float(spacings[i]), "n_fractures": counts[i],
            "objective": float(objectives[i]), "Q": curves[i], "t_eval": t,
            "spacings": spacings, "objectives": objectives}


def plot_optimization(result, filename, reference=None, seconds_per_year=None):
    """ Two-panel figure: objective vs spacing, and Q(t) at the best spacing

    Parameters
    ----------
        result : dict
            output of optimize_spacing()

        filename : str
            figure file name

        reference : tuple (label, Q) or None
            an extra Q(t) curve on the same t_eval, e.g. the DFN-derived spacing

        seconds_per_year : float or None
            if given, the time axis is in years
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = result["t_eval"] / (seconds_per_year or 1.0)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(result["spacings"], result["objectives"] / result["objectives"].max(),
             color="#2a78d6", lw=2)
    ax1.plot(result["spacing"], 1.0, "o", color="#2a78d6", ms=8, mec="white")
    ax1.set_xscale("log")
    ax1.set_xlabel("fracture spacing [m]")
    ax1.set_ylabel("objective / max")
    ax2.plot(t, result["Q"] / 1e6, color="#2a78d6", lw=2,
             label=f"best: {result['spacing']:0.1f} m, n = {result['n_fractures']}")
    if reference is not None:
        ax2.plot(t, reference[1] / 1e6, color="#eb6834", lw=2, ls="--", label=reference[0])
    ax2.set_xlabel("time [yr]" if seconds_per_year else "time [s]")
    ax2.set_ylabel("heat extraction rate [MW]")
    ax2.legend(frameon=False)
    for ax in (ax1, ax2):
        ax.grid(alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    plt.close(fig)

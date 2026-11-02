"""
.. module:: dfn_heat.py
   :synopsis: DFN-based heat extraction: produced temperature and heat
              rate from graph particle tracking with analytic fracture-
              rock heat exchange
.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>

Notes
-----
    Heat obeys the same fracture / rock-matrix equations as a solute with
    matrix porosity phi_m -> (rho c)_r / (rho c)_f and matrix diffusivity
    D_m -> alpha = k / (rho c)_r (thermal_to_tdrw_parameters). Graph
    particle tracking records, for every particle, its advective travel
    time tau and the exchange integral

        beta = sum over edges of 2 tau_i / b_i        [s/m]

    (b_i the edge aperture; Particle.delta_beta). For rock slabs of a
    single half-width B everywhere, the thermal retention time of a whole
    path has the Laplace transform

        exp[ -phi_m sqrt(D_m s) beta tanh( B sqrt(s / D_m) ) ]

    (Tang, Frind and Sudicky 1981; Lauwerier 1955 for B -> inf; the same
    kernel as the Gringarten et al. 1975 parallel-fracture model), so the
    thermal arrival distribution of each particle follows from (tau, beta)
    without event-by-event sampling. With t_s = (phi_m sqrt(D_m) beta)^2
    and beta_G = B / (phi_m D_m beta) it is the slab drawdown
    T_wD((t - tau) / t_s; beta_G) inverted here by Gaver-Stehfest. The
    fraction of produced water that has cooled is the flux-weighted sum of
    those curves over the particles (the tracer-to-thermal mapping of
    Shook 2001), giving

        T_out(t) = T_res - (T_res - T_inj) F(t),
        Q(t)     = m_dot c_f (T_res - T_inj) [1 - F(t)],

    with m_dot = rho_f times the graph-flow rate. B is the hydraulic
    fracture half-spacing D/2 and applies to every edge, natural fractures
    included; B = None gives the semi-infinite rock limit.

    References
        Lauwerier (1955). Appl. Sci. Res. A5, 145.
        Tang, Frind, Sudicky (1981). Water Resour. Res. 17(3), 555.
        Gringarten, Witherspoon, Ohnishi (1975). J. Geophys. Res. 80(8), 1120.
        Shook (2001). Geothermics 30(6), 573.
"""

import numpy as np
from scipy.special import erfc

from pydfnworks.dfnGraph.transport.tdrw.stehfest import stehfest_coefficients

# Representative properties: Utah FORGE granitoid and liquid water
GRANITE_PROPS = {"rho": 2750.0, "c": 790.0, "k": 3.05}
WATER_PROPS = {"rho": 1000.0, "c": 4200.0}
SECONDS_PER_YEAR = 365.25 * 24 * 3600.0
STEHFEST_TERMS = 18


def thermal_to_tdrw_parameters(rock_props, fluid_props):
    """ Map thermal properties onto the TDRW matrix parameters

    Parameters
    ----------
        rock_props : dict
            {"rho": kg/m^3, "c": J/kg/K, "k": W/m/K}

        fluid_props : dict
            {"rho": kg/m^3, "c": J/kg/K}

    Returns
    -------
        tdrw : dict
            "matrix_porosity" = (rho c)_r / (rho c)_f [-] and
            "matrix_diffusivity" = k / (rho c)_r [m^2/s]
    """
    rho_c_r = rock_props["rho"] * rock_props["c"]
    rho_c_f = fluid_props["rho"] * fluid_props["c"]
    return {"matrix_porosity": rho_c_r / rho_c_f,
            "matrix_diffusivity": rock_props["k"] / rho_c_r}


def slab_drawdown(t_D, beta_G, n_terms=STEHFEST_TERMS):
    """ Drawdown of a fracture between finite rock slabs, inverted from
    (1/u) exp[-sqrt(u) tanh(beta_G sqrt(u))] by Gaver-Stehfest

    Parameters
    ----------
        t_D : array_like
            dimensionless time t / t_s, >= 0

        beta_G : array_like, broadcastable with t_D
            slab half-width in conduction lengths, B / sqrt(D_m t_s).
            np.inf gives the semi-infinite closed form erfc(1/(2 sqrt t_D)).

        n_terms : int
            Stehfest coefficients (even). Default 18.

    Returns
    -------
        T_wD : numpy array
            drawdown clipped to [0, 1]; zero for t_D <= 0

    Notes
    -----
        Against a 40-digit reference the error is ~3e-5 for beta_G >= 1
        and ~6e-3 at beta_G = 0.1; below ~0.05 the response is a step at
        t_D = beta_G with Stehfest ringing of a few percent.
    """
    t_D, beta_G = np.broadcast_arrays(np.asarray(t_D, float), np.asarray(beta_G, float))
    out = np.zeros(t_D.shape)
    pos = t_D > 0
    inf = pos & ~np.isfinite(beta_G)
    fin = pos & np.isfinite(beta_G)
    if np.any(inf):
        out[inf] = erfc(0.5 / np.sqrt(t_D[inf]))
    if np.any(fin):
        V = stehfest_coefficients(n_terms)
        a = np.log(2.0) / t_D[fin]
        bg = beta_G[fin]
        acc = np.zeros(a.shape)
        for i, Vi in enumerate(V, 1):
            u = i * a
            q = np.sqrt(u)
            acc += Vi * np.exp(-q * np.tanh(bg * q)) / u
        out[fin] = a * acc
    return np.clip(out, 0.0, 1.0)


def thermal_breakthrough(advect_time, beta, t_eval, rock_props=GRANITE_PROPS,
                         fluid_props=WATER_PROPS, slab_half_width=None, weights=None):
    """ Fraction F(t) of produced water that has cooled, from particle paths

    Parameters
    ----------
        advect_time : array_like
            advective travel time of each exited particle [s]

        beta : array_like
            exchange integral sum(2 tau_i / b_i) of each particle [s/m]

        t_eval : array_like
            times at which to evaluate F [s]

        rock_props, fluid_props : dict
            thermal properties, see thermal_to_tdrw_parameters()

        slab_half_width : float or None
            rock slab half-width B [m], hydraulic fracture spacing / 2.
            None: semi-infinite rock.

        weights : array_like or None
            particle weights (flux); None = equal

    Returns
    -------
        F : numpy array
            cooled fraction in [0, 1] at each time in t_eval
    """
    tau = np.asarray(advect_time, float).ravel()
    beta = np.asarray(beta, float).ravel()
    t = np.atleast_1d(np.asarray(t_eval, float)).ravel()
    if tau.size == 0:
        raise ValueError("no exited particles")
    w = np.ones(tau.size) if weights is None else np.asarray(weights, float).ravel()
    w = w / w.sum()
    p = thermal_to_tdrw_parameters(rock_props, fluid_props)
    phi, D = p["matrix_porosity"], p["matrix_diffusivity"]
    t_s = (phi * np.sqrt(D) * beta)**2
    if slab_half_width is None:
        beta_G = np.full(tau.size, np.inf)
    else:
        beta_G = slab_half_width / (phi * D * beta)
    t_D = (t[:, None] - tau[None, :]) / t_s[None, :]
    return slab_drawdown(t_D, beta_G[None, :]) @ w


def heat_extraction(dfn, G, T_res, T_inj, t_eval, rock_props=GRANITE_PROPS,
                    fluid_props=WATER_PROPS, slab_half_width=None, nparticles=10**4,
                    partime_file="thermal_partime", **transport_kwargs):
    """ Produced temperature and heat extraction rate of a DFN

    Parameters
    ----------
        dfn : DFNWORKS object

        G : NetworkX graph
            flow graph from run_graph_flow() (inflow to outflow boundary)

        T_res : float
            initial reservoir temperature [deg C]

        T_inj : float
            injection temperature [deg C]

        t_eval : array_like
            times at which to evaluate the curves [s]

        rock_props, fluid_props : dict
            thermal properties

        slab_half_width : float or None
            rock slab half-width B [m]; hydraulic spacing / 2, None = infinite

        nparticles : int
            particles for graph transport (flux-weighted release). Default 10**4.

        partime_file : str
            graph transport output file stem

        transport_kwargs :
            forwarded to run_graph_transport()

    Returns
    -------
        result : dict
            "t_eval" [s], "F" cooled fraction, "T_out" [deg C], "Q" [W],
            "m_dot" [kg/s] from the graph flow, "flow_rate" [m^3/s],
            "advect_time" and "beta" of the exited particles,
            "n_exit", "n_stuck"

    Notes
    -----
        Particle tracking is purely advective (tdrw_flag=False); the
        thermal retention is applied analytically per particle, see the
        module notes. Particles that do not exit are dropped.
    """
    if T_inj >= T_res:
        raise ValueError(f"need T_inj < T_res, got {T_inj}, {T_res}")
    t = np.atleast_1d(np.asarray(t_eval, float)).ravel()
    particles = dfn.run_graph_transport(G, nparticles, partime_file, None,
                                        initial_positions="flux", tdrw_flag=False,
                                        **transport_kwargs)
    exited = [p for p in particles if p.exit_flag]
    tau = np.array([p.advect_time for p in exited])
    beta = np.array([p.beta for p in exited])
    inlet = [v for v in G.nodes if G.nodes[v]["inletflag"]]
    flow_rate = float(sum(G.edges[u, v]["vol_flow_rate"] for u in inlet for v in G.successors(u)))
    m_dot = fluid_props["rho"] * flow_rate
    F = thermal_breakthrough(tau, beta, t, rock_props, fluid_props, slab_half_width)
    dT = T_res - T_inj
    dfn.print_log(f"--> Thermal breakthrough from {tau.size} particles "
                  f"({len(particles) - tau.size} stuck); m_dot = {m_dot:0.3f} kg/s, "
                  f"cooled fraction at t_end = {F[-1]:0.3f}")
    return {"t_eval": t, "F": F, "T_out": T_res - dT * F, "Q": m_dot * fluid_props["c"] * dT * (1 - F),
            "m_dot": m_dot, "flow_rate": flow_rate, "advect_time": tau, "beta": beta,
            "n_exit": int(tau.size), "n_stuck": int(len(particles) - tau.size)}

"""
.. module:: dfn_heat.py
   :synopsis: DFN graph-transport heat-extraction backend (Phase 2).
              Interface and parameter mapping only; the kernel raises
              NotImplementedError.
.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>

Notes
-----
    Heat obeys the same fracture/matrix equations as a solute with
    matrix_porosity -> (rho c)_r / (rho c)_f and matrix_diffusivity ->
    k / (rho c)_r (thermal_to_tdrw_parameters), so the backend is a call
    to run_graph_transport() with those parameters followed by the
    tracer-to-thermal conversion of the arrival times (Shook 2001). Two
    things are left: that conversion, and a per-edge draw from the
    finite-slab kernel. The existing finite-slab TDRW models sample
    excursions event by event at a rate proportional to phi_m D_m, about
    1e7 times the solute value under this mapping, so one particle takes
    minutes; only tdrw_model='infinite' is practical with thermal
    parameters today.

    Reference: Shook (2001). Geothermics 30(6), 573-589.
"""

from pydfnworks.thermal.common import check_inputs


def thermal_to_tdrw_parameters(rock_props, fluid_props):
    """ Map thermal properties onto the TDRW parameters of run_graph_transport()

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


def heat_extraction(spacing, m_dot, T_res, T_inj, rock_props, fluid_props,
                    n_fractures, fracture_height, t_eval, *,
                    dfn=None, G=None, **transport_kwargs):
    """ Heat extraction rate Q(t) from DFN graph transport (not implemented)

    Parameters
    ----------
        spacing, ..., t_eval :
            as in pydfnworks.thermal.mpf_gringarten.heat_extraction()

        dfn : DFNWORKS object
            keyword-only; DFN whose flow graph is used

        G : NetworkX graph
            keyword-only; flow graph from run_graph_flow()

        transport_kwargs :
            forwarded to run_graph_transport()

    Returns
    -------
        Q : numpy array
            heat extraction rate at each time in t_eval [W]
    """
    check_inputs(spacing, m_dot, T_res, T_inj, rock_props, fluid_props,
                 n_fractures, fracture_height, t_eval)
    tdrw = thermal_to_tdrw_parameters(rock_props, fluid_props)
    raise NotImplementedError(
        "pydfnworks.thermal.dfn_heat: the DFN heat kernel is not implemented; "
        "use mpf_gringarten.heat_extraction. Planned: run_graph_transport(G, "
        f"..., tdrw_flag=True, matrix_porosity={tdrw['matrix_porosity']:.3f}, "
        f"matrix_diffusivity={tdrw['matrix_diffusivity']:.2e}, "
        f"fracture_spacing={spacing:.3g}) then arrival times -> Q(t).")

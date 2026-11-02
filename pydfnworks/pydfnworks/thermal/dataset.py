"""
.. module:: dataset.py
   :synopsis: Helpers to build a DFN heat-extraction dataset: parallel
              hydraulic fracture placement and per-sample HDF5 output
.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
"""

import os
import csv
import json
import numpy as np
import h5py


def hydraulic_fracture_positions(length, spacing):
    """ Centers of parallel hydraulic fractures filling a length at uniform spacing

    Parameters
    ----------
        length : float
            domain extent along the fracture normal [m]

        spacing : float
            fracture spacing D [m]

    Returns
    -------
        positions : numpy array
            n = floor(length / spacing) coordinates, centered on zero,
            spacing apart, with a half-spacing margin to each boundary
    """
    n = int(np.floor(length / spacing))
    if n < 1:
        raise ValueError(f"spacing {spacing} m exceeds the domain length {length} m")
    return (np.arange(n) - 0.5 * (n - 1)) * spacing


def add_hydraulic_fractures(dfn, spacing, length, radius, permeability, normal=(0, 1, 0)):
    """ Add parallel rectangular hydraulic fractures to a DFN

    Parameters
    ----------
        dfn : DFNWORKS object
            before check_input() / create_network()

        spacing : float
            fracture spacing [m]

        length : float
            domain extent along the normal [m]

        radius : float
            rectangle half-size [m]; larger than the domain cross-section
            so each fracture spans it

        permeability : float
            fracture permeability [m^2] (aperture by the cubic law)

        normal : tuple
            unit normal of the fracture planes. Default (0, 1, 0).

    Returns
    -------
        positions : numpy array
            fracture center coordinates along the normal [m]
    """
    positions = hydraulic_fracture_positions(length, spacing)
    n = np.asarray(normal, float)
    n = n / np.linalg.norm(n)
    # plain Python floats: the values are written verbatim into dfnGen's input
    normal_vector = [float(v) for v in n]
    for s in positions:
        dfn.add_user_fract(shape="rect", radii=float(radius),
                           translation=[float(v) for v in s * n],
                           normal_vector=normal_vector, permeability=permeability)
    return positions


def graph_flow_for_rate(dfn, inflow, outflow, m_dot, fluid_props, pressure_out=1.0e6):
    """ Solve graph flow with the pressure drop that delivers a target mass flow rate

    Parameters
    ----------
        dfn : DFNWORKS object

        inflow, outflow : str
            boundary names for run_graph_flow()

        m_dot : float
            target mass flow rate through the network [kg/s]

        fluid_props : dict
            {"rho": kg/m^3, "c": J/kg/K}

        pressure_out : float
            outlet pressure [Pa]. Default 1e6.

    Returns
    -------
        G : NetworkX graph
            flow graph at the rescaled pressure drop

        pressure_drop : float
            inlet minus outlet pressure that gives m_dot [Pa]

    Notes
    -----
        Graph flow is linear in the pressure drop, so a unit-drop solve
        fixes the network conductance and a second solve at the rescaled
        drop gives the target rate. Imposing the rate rather than the drop
        makes samples with different fracture counts comparable.
    """
    G = dfn.run_graph_flow(inflow, outflow, pressure_out + 1.0, pressure_out)
    inlet = [v for v in G.nodes if G.nodes[v]["inletflag"]]
    rate_unit = sum(G.edges[u, v]["vol_flow_rate"] for u in inlet for v in G.successors(u))
    pressure_drop = m_dot / (fluid_props["rho"] * rate_unit)
    dfn.print_log(f"--> Pressure drop for {m_dot:0.2f} kg/s: {pressure_drop:0.3e} Pa")
    G = dfn.run_graph_flow(inflow, outflow, pressure_out + pressure_drop, pressure_out)
    return G, pressure_drop


def write_sample(filename, inputs, result):
    """ Write one dataset sample (inputs, curves, particle paths) to HDF5

    Parameters
    ----------
        filename : str
            output .h5 file

        inputs : dict
            scalar, string or small-array inputs (spacing, seed, p32, ...);
            stored as attributes of the "inputs" group

        result : dict
            output of pydfnworks.thermal.heat_extraction()

    Returns
    -------
        None
    """
    with h5py.File(filename, "w") as f:
        g = f.create_group("inputs")
        for k, v in inputs.items():
            g.attrs[k] = json.dumps(v) if isinstance(v, dict) else v
        c = f.create_group("curves")
        for k in ("t_eval", "F", "T_out", "Q"):
            c.create_dataset(k, data=result[k])
        c.attrs["units"] = "t_eval s, F -, T_out degC, Q W"
        p = f.create_group("particles")
        p.create_dataset("advect_time", data=result["advect_time"])
        p.create_dataset("beta", data=result["beta"])
        p.attrs["units"] = "advect_time s, beta s/m (sum 2 tau_i / b_i)"
        for k in ("m_dot", "flow_rate", "n_exit", "n_stuck"):
            f.attrs[k] = result[k]


def append_manifest(filename, row):
    """ Append one row (dict) to a CSV manifest, writing the header if new """
    new = not os.path.isfile(filename)
    with open(filename, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        if new:
            w.writeheader()
        w.writerow(row)

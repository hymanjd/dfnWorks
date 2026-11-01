#"""
#   :synopsis: Driver run file for heat extraction on a graph representation
#              of a DFN: graph flow, finite-slab TDRW with thermal
#              properties, and the Gringarten spacing surrogate
#   :version: 1.0
#   :maintainer: Jeffrey Hyman
#.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
#
# Built on examples/graph_transport. The DFN, graph flow, and graph
# transport steps are the same, except that the fracture permeability is
# raised from 2e-12 to 2e-9 m^2 (aperture about 0.15 mm, a stimulated
# network) so that the 1 MPa drop drives a flow of geothermal scale
# through this 400 m block; what is added is
#   1. the surrogate geometry read from the flow graph: the fractures on
#      the inflow face (count and mean diameter) and the DFN's
#      characteristic fracture spacing from compute_dQ (1/dQ and 1/p32),
#   2. graph particle transport with TDRW driven by thermal rather than
#      solute properties, using the Phase 2 parameter mapping in
#      pydfnworks.thermal.dfn_heat (matrix_porosity ->
#      (rho c)_rock/(rho c)_fluid, matrix_diffusivity -> rock thermal
#      diffusivity), which is the ingredient the DFN-based heat backend
#      will build on. The infinite-matrix model is used: the finite-slab
#      models sample matrix excursions event by event at a rate
#      proportional to matrix_porosity * matrix_diffusivity, which is
#      about 1e7 times larger for heat than for a solute, and a single
#      particle takes minutes (see the notes in dfn_heat.py),
#   3. the Gringarten surrogate and spacing optimizer on the
#      graph-derived geometry.
#"""

import os
import functools
import numpy as np
import h5py

from pydfnworks import *
from pydfnworks.thermal import (heat_extraction, optimize_spacing, time_to_drawdown,
                                plot_optimization, thermal_to_tdrw_parameters,
                                dimensionless_parameters, GRANITE_PROPS, WATER_PROPS,
                                SECONDS_PER_YEAR)

src_path = os.getcwd()
jobname = f"{src_path}/output"

DFN = DFNWORKS(jobname, ncpu=8)

DFN.params['domainSize']['value'] = [400, 50, 50]
DFN.params['h']['value'] = 1
DFN.params['domainSizeIncrease']['value'] = [5, 5, 5]
DFN.params['ignoreBoundaryFaces']['value'] = True
DFN.params['boundaryFaces']['value'] = [1, 1, 0, 0, 0, 0]
DFN.params['keepOnlyLargestCluster']['value'] = True
DFN.params['disableFram']['value'] = True

DFN.add_fracture_family(shape="ell",
                        distribution="tpl",
                        kappa=10,
                        p32=0.25,
                        aspect=1,
                        theta=0.0,
                        phi=0.0,
                        alpha=1.8,
                        min_radius=10.0,
                        max_radius=20.0,
                        hy_variable='permeability',
                        hy_function='constant',
                        hy_params={"mu": 2e-9})

DFN.add_fracture_family(shape="ell",
                        distribution="tpl",
                        kappa=10,
                        p32=0.25,
                        aspect=1,
                        theta=0.0,
                        phi=270.0,
                        alpha=1.8,
                        min_radius=10.0,
                        max_radius=20.0,
                        hy_variable='permeability',
                        hy_function='constant',
                        hy_params={"mu": 2e-9})

DFN.add_fracture_family(shape="ell",
                        distribution="tpl",
                        kappa=10,
                        p32=0.25,
                        aspect=1,
                        theta=45.0,
                        phi=0.0,
                        alpha=1.8,
                        min_radius=10.0,
                        max_radius=20.0,
                        hy_variable='permeability',
                        hy_function='constant',
                        hy_params={"mu": 3e-9})

DFN.make_working_directory(delete=True)
DFN.check_input()
DFN.create_network()

pressure_in = 2 * 10**6
pressure_out = 10**6
G = DFN.run_graph_flow("left", "right", pressure_in, pressure_out)

# ---------------------------------------------------------------------
# Surrogate geometry from the flow graph. Inlet vertices are intersections
# of a fracture with the inflow face; their 'frac' attribute is (id, 's').
# ---------------------------------------------------------------------
inlet_nodes = [v for v in G.nodes if G.nodes[v]['inletflag']]
inlet_fractures = np.unique([G.nodes[v]['frac'][0] for v in inlet_nodes]).astype(int)
n_inlet = len(inlet_fractures)
fracture_height = float(np.mean(2 * DFN.radii[inlet_fractures - 1, 2]))
graph_flow_rate = sum(G.edges[u, v]['vol_flow_rate']
                      for u in inlet_nodes for v in G.successors(u))
p32, dQ, _ = DFN.compute_dQ(G)
spacing_dfn = 1.0 / dQ
DFN.print_log(f"--> {n_inlet} fractures on the inflow face, mean diameter {fracture_height:0.1f} m; "
              f"flow {graph_flow_rate:0.2e} m^3/s; spacing 1/dQ = {spacing_dfn:0.2f} m")

# ---------------------------------------------------------------------
# Graph transport with thermal properties. The infinite-matrix TDRW model
# draws one closed-form retention time per edge; do not pass
# fracture_spacing, which switches to the event-by-event slab sampler.
# ---------------------------------------------------------------------
tdrw = thermal_to_tdrw_parameters(GRANITE_PROPS, WATER_PROPS)
DFN.run_graph_transport(G, 10**4, "thermal_partime", "thermal_frac_sequence",
                        tdrw_flag=True, tdrw_model="infinite",
                        matrix_porosity=tdrw["matrix_porosity"],
                        matrix_diffusivity=tdrw["matrix_diffusivity"])
with h5py.File("thermal_partime.hdf5", "r") as fp:
    advective = np.asarray(fp["Advective time [s]"])
    total = np.asarray(fp["Total travel time [s]"])
# the retention distribution is heavy-tailed; early quantiles control breakthrough
for label, times in (("advective", advective), ("thermally retarded", total)):
    DFN.print_log(f"--> {label} time 1, 10, 50 % quantiles: "
                  + ", ".join(f"{v:0.2e}" for v in np.percentile(times, [1, 10, 50])) + " s")

# ---------------------------------------------------------------------
# Surrogate and spacing optimization on the graph-derived geometry
# ---------------------------------------------------------------------
m_dot = WATER_PROPS["rho"] * graph_flow_rate  # the rate the network carries
T_res = 226.0  # deg C
T_inj = 40.0  # deg C
t_eval = np.linspace(0.0, 30.0 * SECONDS_PER_YEAR, 361)

Q_dfn = heat_extraction(spacing_dfn, m_dot, T_res, T_inj, GRANITE_PROPS,
                        WATER_PROPS, n_inlet, fracture_height, t_eval)
p = dimensionless_parameters(spacing_dfn, m_dot, GRANITE_PROPS, WATER_PROPS,
                             n_inlet, fracture_height)
DFN.print_log(f"--> Surrogate at 1/dQ spacing: beta = {p['beta']:0.3f}, piston time "
              f"{p['t_piston']/SECONDS_PER_YEAR:0.2f} yr, Q(0) = {Q_dfn[0]/1e6:0.2f} MW")

# stages fill the domain length: n = L / spacing; objective: time to 10 % drawdown
L = DFN.params['domainSize']['value'][0]
result = optimize_spacing(L / np.arange(2, int(L / 2.0) + 1),
                          m_dot, T_res, T_inj, GRANITE_PROPS, WATER_PROPS,
                          lambda spacing: int(round(L / spacing)), fracture_height, t_eval,
                          objective=functools.partial(time_to_drawdown, fraction=0.9))
DFN.print_log(f"--> Best spacing {result['spacing']:0.1f} m with {result['n_fractures']} stages: "
              f"time to 10 % drawdown {result['objective']/SECONDS_PER_YEAR:0.2f} yr")

np.savetxt("objective_vs_spacing.csv", np.c_[result["spacings"], result["objectives"]],
           delimiter=",", header="spacing [m], time to 10 % drawdown [s]", comments="# ")
np.savetxt("heat_rate_best.csv", np.c_[t_eval / SECONDS_PER_YEAR, result["Q"]],
           delimiter=",", header=f"time [yr], heat rate [W]; spacing {result['spacing']:0.2f} m",
           comments="# ")
plot_optimization(result, "spacing_optimization.png",
                  reference=(f"1/dQ spacing: {spacing_dfn:0.1f} m, n = {n_inlet}", Q_dfn),
                  seconds_per_year=SECONDS_PER_YEAR)
DFN.print_log("--> Wrote spacing_optimization.png, objective_vs_spacing.csv, heat_rate_best.csv")

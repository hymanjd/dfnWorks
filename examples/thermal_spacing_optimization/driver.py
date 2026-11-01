#"""
#   :synopsis: Driver run file for hydraulic fracture spacing optimization
#              with the pydfnworks.thermal heat-extraction surrogate on a
#              FORGE-representative DFN and injection/production well pair
#   :version: 1.0
#   :maintainer: Jeffrey Hyman
#.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
#
# Workflow
#   1. Generate a FORGE-representative DFN (Utah FORGE granitoid, ~2.3 km
#      depth): one set parallel to the hydraulic fractures (strike along
#      S_Hmax, near vertical; the model x axis is the lateral direction and
#      the normal of this set) and one natural set (strike ~300, dip ~70).
#   2. Place an injection lateral (16A(78)-32 analogue) and a production
#      lateral 100 m above it (16B(78)-32 analogue) with the well package,
#      and find the fractures each well intersects.
#   3. Take the surrogate geometry from the DFN: the number of fractures
#      cut by the injection lateral and their mean diameter (fracture
#      height), plus the lateral length inside the domain.
#   4. Run the Gringarten multiple-parallel-fracture surrogate at the
#      DFN-derived spacing, then optimize the spacing of stages filling
#      the lateral for 30-year cumulative heat, and write Q(t) for the
#      best spacing and the objective-vs-spacing curve.
#
# Reservoir numbers are representative of FORGE (T ~ 226 C, granitoid
# rho = 2750 kg/m^3, c = 790 J/kg/K, k = 3.05 W/m/K); the 600 m block is
# roughly a quarter of the 16A lateral, so the circulation rate is scaled
# down to 20 kg/s.
#"""

import os
import numpy as np

from pydfnworks import *
from pydfnworks.dfnGen.well_package.wells import get_segments
from pydfnworks.thermal import (heat_extraction, optimize_spacing, cumulative_heat,
                                plot_optimization, GRANITE_PROPS, WATER_PROPS,
                                SECONDS_PER_YEAR)

src_path = os.getcwd()
jobname = f"{src_path}/output"

DFN = DFNWORKS(jobname, ncpu=8)

DFN.params['domainSize']['value'] = [600, 600, 600]
DFN.params['h']['value'] = 5.0
DFN.params['domainSizeIncrease']['value'] = [20, 20, 20]
DFN.params['stopCondition']['value'] = 1
DFN.params['ignoreBoundaryFaces']['value'] = True
DFN.params['boundaryFaces']['value'] = [1, 1, 0, 0, 0, 0]
DFN.params['keepOnlyLargestCluster']['value'] = False
# reduced mesh only: enough for the well intersection search
DFN.params['visualizationMode']['value'] = True
DFN.params['seed']['value'] = 1

# Set 1: parallel to the hydraulic fractures (normal along the lateral, x)
DFN.add_fracture_family(shape="ell",
                        distribution="tpl",
                        kappa=25,
                        p32=0.03,
                        aspect=1,
                        theta=90.0,
                        phi=0.0,
                        alpha=2.4,
                        min_radius=60.0,
                        max_radius=200.0,
                        hy_variable="permeability",
                        hy_function="constant",
                        hy_params={"mu": 1e-11})

# Set 2: FORGE natural fracture set (strike ~300, dip ~70)
DFN.add_fracture_family(shape="ell",
                        distribution="tpl",
                        kappa=12,
                        p32=0.01,
                        aspect=1,
                        theta=70.0,
                        phi=60.0,
                        alpha=2.6,
                        min_radius=40.0,
                        max_radius=120.0,
                        hy_variable="permeability",
                        hy_function="constant",
                        hy_params={"mu": 1e-12})

DFN.make_working_directory(delete=True)
DFN.check_input()
DFN.create_network()

# Wells: injection lateral at z = -50 m, production lateral at z = +50 m,
# both along x (8.5 inch hole, r = 0.11 m)
inject_well = {"name": "inject", "filename": "well_inject.dat", "r": 0.11}
extract_well = {"name": "extract", "filename": "well_extract.dat", "r": 0.11}
wells = [inject_well, extract_well]
for well in wells:
    os.symlink(f"{src_path}/{well['filename']}", well['filename'])

# builds reduced_mesh.inp (visualization mode) and well_points.dat
DFN.find_well_intersection_points(wells)

# ---------------------------------------------------------------------
# Surrogate geometry from the DFN: fractures cut by the injection lateral
# ---------------------------------------------------------------------
_, _, inject_fractures = get_segments(f"well_{inject_well['name']}_intersect.inp")
inject_ids = np.unique(np.asarray(inject_fractures, dtype=int))
n_dfn = len(inject_ids)
if n_dfn < 2:
    DFN.print_log("Error. The injection lateral cuts fewer than two fractures; "
                  "increase p32 or the lateral length.", 'error')
fracture_height = float(np.mean(2 * DFN.radii[inject_ids - 1, 2]))
half_domain = 0.5 * np.asarray(DFN.params['domainSize']['value'], dtype=float)
well_xyz = np.clip(np.loadtxt(f"{src_path}/{inject_well['filename']}"), -half_domain, half_domain)
lateral_length = float(np.linalg.norm(well_xyz[-1] - well_xyz[0]))
spacing_dfn = lateral_length / (n_dfn - 1)
DFN.print_log(f"--> Injection lateral cuts {n_dfn} fractures, mean diameter "
              f"{fracture_height:0.1f} m, lateral {lateral_length:0.0f} m, "
              f"spacing {spacing_dfn:0.1f} m")

# ---------------------------------------------------------------------
# Heat extraction surrogate and spacing optimization
# ---------------------------------------------------------------------
m_dot = 20.0  # kg/s, scaled to the 600 m block
T_res = 226.0  # deg C
T_inj = 40.0  # deg C
t_eval = np.linspace(0.0, 30.0 * SECONDS_PER_YEAR, 361)

Q_dfn = heat_extraction(spacing_dfn, m_dot, T_res, T_inj, GRANITE_PROPS,
                        WATER_PROPS, n_dfn, fracture_height, t_eval)
DFN.print_log(f"--> DFN spacing: Q(0) = {Q_dfn[0]/1e6:0.2f} MW, "
              f"Q(30 yr) = {Q_dfn[-1]/1e6:0.2f} MW, "
              f"{cumulative_heat(t_eval, Q_dfn)/1e15:0.2f} PJ over 30 yr")

# stages fill the lateral: n = lateral_length / spacing, integer n
result = optimize_spacing(lateral_length / np.arange(2, int(lateral_length / 5.0) + 1),
                          m_dot, T_res, T_inj, GRANITE_PROPS, WATER_PROPS,
                          lambda spacing: int(round(lateral_length / spacing)),
                          fracture_height, t_eval, objective=cumulative_heat)
DFN.print_log(f"--> Best spacing {result['spacing']:0.1f} m with {result['n_fractures']} "
              f"stages: {result['objective']/1e15:0.2f} PJ over 30 yr "
              f"({result['spacings'].size} spacings evaluated)")

np.savetxt("objective_vs_spacing.csv", np.c_[result["spacings"], result["objectives"]],
           delimiter=",", header="spacing [m], cumulative heat 30 yr [J]", comments="# ")
np.savetxt("heat_rate_best.csv", np.c_[t_eval / SECONDS_PER_YEAR, result["Q"]],
           delimiter=",", header=f"time [yr], heat rate [W]; spacing {result['spacing']:0.2f} m",
           comments="# ")
plot_optimization(result, "spacing_optimization.png",
                  reference=(f"DFN-derived: {spacing_dfn:0.0f} m, n = {n_dfn}", Q_dfn),
                  seconds_per_year=SECONDS_PER_YEAR)
DFN.print_log("--> Wrote spacing_optimization.png, objective_vs_spacing.csv, heat_rate_best.csv")

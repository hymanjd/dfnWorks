#"""
#   :synopsis: Driver run file for a DFN heat-extraction dataset: parallel
#              hydraulic fractures at a range of spacings in a natural
#              fracture network, graph flow and transport, produced
#              temperature and heat rate per sample for downstream (ML) use
#   :version: 1.0
#   :maintainer: Jeffrey Hyman
#.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
#
# For each (spacing, seed):
#   1. dfnGen: parallel rectangular hydraulic fractures (normal y, spacing
#      D, spanning the domain in x and z) plus two natural fracture sets
#   2. graph flow from the left (-x) to the right (+x) face at the
#      pressure drop that delivers the target mass flow rate: water enters
#      every hydraulic fracture and flows along it, with the natural
#      network cross-connecting them
#   3. pydfnworks.thermal.heat_extraction: graph particle tracking gives
#      each path's advective time and exchange integral; the slab kernel
#      with half-width D/2 gives T_out(t) and Q(t)
#   4. write_sample: inputs, curves and particle paths to HDF5, one row in
#      manifest.csv
#
# Reservoir numbers are representative of Utah FORGE (226 C granitoid).
# A target mass flow rate is imposed on every sample (the pressure drop is
# solved for), so spacing and network are the only things that change.
#"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pydfnworks import *
from pydfnworks.thermal import (heat_extraction, add_hydraulic_fractures, graph_flow_for_rate,
                                write_sample, append_manifest, GRANITE_PROPS, WATER_PROPS,
                                SECONDS_PER_YEAR)

src_path = os.getcwd()
dataset_path = f"{src_path}/output/dataset"
os.makedirs(dataset_path, exist_ok=True)
manifest = f"{dataset_path}/manifest.csv"
if os.path.isfile(manifest):
    os.remove(manifest)

domain = [400.0, 200.0, 200.0]     # x: flow direction, y: hydraulic fracture normal
spacings = [10.0, 20.0, 40.0, 80.0]
seeds = [1, 2]
k_hydraulic = 1e-8                 # m^2, aperture ~0.35 mm by the cubic law
k_natural = 1e-12
p32_natural = 0.02                 # per set
m_dot_target = 10.0                # kg/s through the block; the pressure drop is solved for
pressure_out = 1.0e6
T_res, T_inj = 226.0, 40.0
t_eval = np.linspace(0.0, 30.0 * SECONDS_PER_YEAR, 361)
curves = {}

for spacing in spacings:
    for seed in seeds:
        name = f"spacing_{spacing:g}_seed_{seed}"
        os.chdir(src_path)
        DFN = DFNWORKS(f"{src_path}/output/{name}", ncpu=8)
        DFN.params['domainSize']['value'] = domain
        DFN.params['h']['value'] = 2.0
        DFN.params['domainSizeIncrease']['value'] = [5, 5, 5]
        DFN.params['stopCondition']['value'] = 1
        DFN.params['ignoreBoundaryFaces']['value'] = True
        DFN.params['boundaryFaces']['value'] = [1, 1, 0, 0, 0, 0]
        DFN.params['keepOnlyLargestCluster']['value'] = True
        DFN.params['disableFram']['value'] = True
        DFN.params['seed']['value'] = seed

        # natural sets: normals along x (cross the hydraulic planes, span y)
        # and along z (sub-horizontal), truncated power-law sizes
        for theta, phi in ((90.0, 0.0), (0.0, 0.0)):
            DFN.add_fracture_family(shape="ell", distribution="tpl", kappa=10, p32=p32_natural,
                                    aspect=1, theta=theta, phi=phi, alpha=2.2,
                                    min_radius=20.0, max_radius=60.0,
                                    hy_variable="permeability", hy_function="constant",
                                    hy_params={"mu": k_natural})
        positions = add_hydraulic_fractures(DFN, spacing, length=domain[1], radius=230.0,
                                            permeability=k_hydraulic, normal=(0, 1, 0))

        DFN.make_working_directory(delete=True)
        DFN.check_input()
        DFN.create_network()
        G, pressure_drop = graph_flow_for_rate(DFN, "left", "right", m_dot_target, WATER_PROPS,
                                               pressure_out)
        p32, dQ, _ = DFN.compute_dQ(G)

        result = heat_extraction(DFN, G, T_res, T_inj, t_eval, GRANITE_PROPS, WATER_PROPS,
                                 slab_half_width=spacing / 2, nparticles=10**4)
        inputs = {"spacing": spacing, "n_hydraulic": int(positions.size), "seed": seed,
                  "domain": domain, "k_hydraulic": k_hydraulic, "k_natural": k_natural,
                  "p32_natural_per_set": p32_natural, "m_dot_target": m_dot_target,
                  "pressure_drop": pressure_drop, "pressure_out": pressure_out,
                  "T_res": T_res, "T_inj": T_inj,
                  "rock_props": GRANITE_PROPS, "fluid_props": WATER_PROPS,
                  "num_fractures": int(DFN.num_frac), "p32_total": float(p32), "dQ": float(dQ)}
        write_sample(f"{dataset_path}/{name}.h5", inputs, result)
        append_manifest(manifest, {"sample": f"{name}.h5", "spacing": spacing, "seed": seed,
                                   "n_hydraulic": int(positions.size), "num_fractures": int(DFN.num_frac),
                                   "m_dot_kg_s": result["m_dot"], "pressure_drop_Pa": pressure_drop,
                                   "n_exit": result["n_exit"],
                                   "n_stuck": result["n_stuck"], "T_out_30yr": result["T_out"][-1],
                                   "Q0_MW": result["Q"][0] / 1e6, "Q_30yr_MW": result["Q"][-1] / 1e6})
        curves[(spacing, seed)] = result
        DFN.print_log(f"--> {name}: {DFN.num_frac} fractures, m_dot {result['m_dot']:0.2f} kg/s, "
                      f"Q(0) {result['Q'][0]/1e6:0.2f} MW, Q(30 yr) {result['Q'][-1]/1e6:0.2f} MW")

# figure: produced temperature and heat rate for every sample
os.chdir(src_path)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
colors = dict(zip(spacings, ("#7a1f2b", "#b23a48", "#d9775c", "#e8a87c")))
for (spacing, seed), r in curves.items():
    t = r["t_eval"] / SECONDS_PER_YEAR
    style = dict(color=colors[spacing], lw=2, ls="-" if seed == seeds[0] else "--")
    ax1.plot(t, r["T_out"], label=f"D = {spacing:g} m, seed {seed}", **style)
    ax2.plot(t, r["Q"] / 1e6, **style)
ax1.set_xlabel("time [yr]"); ax1.set_ylabel("produced temperature [C]")
ax2.set_xlabel("time [yr]"); ax2.set_ylabel("heat extraction rate [MW]")
ax1.legend(fontsize=8, frameon=False)
for ax in (ax1, ax2):
    ax.grid(alpha=0.25); ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{dataset_path}/dataset_curves.png", dpi=150)
print(f"--> Dataset written to {dataset_path} ({len(curves)} samples)")

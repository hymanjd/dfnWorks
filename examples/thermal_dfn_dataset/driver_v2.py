#"""
#   :synopsis: Graph-based EGS prototype: segmented lateral wells and 3 hydraulic
#              fractures, graph flow and particle tracking with native pydfnworks calls,
#              then produced temperature and heat rate from the particle paths
#   :version: 0.1
#   :maintainer: Jeffrey Hyman
#.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
#
#   usage: python driver_v2.py [runs/<run>.yaml]     (no file: make_yaml_file.DEFAULT_CONFIG)
#
#   Steps:
#     1. read the run settings (YAML run file on top of make_yaml_file.DEFAULT_CONFIG)
#     2. build the DFN: two segmented lateral wells + evenly spaced hydraulic fractures (HFs)
#     3. graph flow from the injection well's -x end segment to the production well's
#     4. advective graph particle tracking: travel time tau and exchange integral beta per particle
#     5. thermal model: cooled fraction F(t), produced temperature T_out(t), heat rate Q(t)
#     6. write particles.npz, thermal.npz, run_info.json; add a row to the master CSV
#     7. gather the outputs in results/<jobname>/ and delete the dfnWorks job folder
#"""

from pydfnworks import *
import os
import sys
import json
import numpy as np

from helper_functions import *

# ## Run settings: the YAML file given on the command line, on top of the defaults
run_file = sys.argv[1] if len(sys.argv) > 1 else None
cfg = load_run_config(run_file)
# short names for the setting groups used most below
hf, wl, fl, th = cfg["hydraulic_fractures"], cfg["wells"], cfg["flow"], cfg["thermal"]
num_hf = hf["number"]
p_in, p_out = fl["p_in"], fl["p_out"]
T_res, T_inj = th["T_res"], th["T_inj"]
t_end_years = th["t_end_years"]
nparticles = cfg["transport"]["nparticles"]

print(f"--> Run file: {run_file or '(defaults)'}")
print(f"--> {num_hf} HFs, radius {hf['radius']} m, aspect {hf['aspect']}")

# the run's dfnWorks job folder sits in the folder the driver is run from
src_path = os.getcwd()
jobname = f"{src_path}/{cfg['jobname']}"


# ## DFN: domain and dfnGen options
DFN = DFNWORKS(jobname, ncpu=cfg["ncpu"])

DFN.params['domainSize']['value'] = cfg["domain"]["size"]
DFN.params['domainCenter']['value'] = cfg["domain"]["center"]
DFN.params['h']['value'] = cfg["domain"]["h"]                 # minimum feature size; below the well segment overlap
DFN.params['orientationOption']['value'] = 1
DFN.params['stopCondition']['value'] = 1
DFN.params['seed']['value'] = cfg["domain"]["seed"]
DFN.params['ignoreBoundaryFaces']['value'] = True             # flow is well to well, not face to face
DFN.params['disableFram']['value'] = True
DFN.params['keepOnlyLargestCluster']['value'] = False


# ## Hydraulic fracture centers (same layout rule as ideal_HF)
# num_hf fractures spaced evenly over x in [-length/2, length/2], half a spacing from each end
domain_length = hf["length"]
hf_spacing = domain_length / num_hf
centers = [-0.5*domain_length + -0.5 * hf_spacing + hf_spacing * (i + 1) for i in range(int(num_hf))]

# ## Lateral wells (helper_functions.add_lateral_wells)
wells = add_lateral_wells(DFN, centers, src_path, lateral_y=wl["lateral_y"], lateral_x=wl["lateral_x"],
                          seg_length=wl["seg_length"], seg_overlap=wl["seg_overlap"],
                          seg_width=wl["seg_width"], aperture=wl["aperture"])
# fracture numbering: well segments first (injection well, then production well), then the HFs
n_seg = wells["n_seg"]
inj_frac, prod_frac = wells["inj_frac"], wells["prod_frac"]
hf_fracs = [wells["n_fracs"] + i + 1 for i in range(num_hf)]

# ## Hydraulic fractures: ellipses normal to x, centered between the two wells (y = 0)
for fracture_center in centers:
    DFN.add_user_fract(shape='ell',
                       radii=hf["radius"],
                       aspect_ratio=hf["aspect"],
                       translation=[fracture_center, 0, 0],
                       normal_vector=[1, 0, 0],
                       number_of_vertices=16,
                       aperture=hf["aperture"])



# ## Generate the network (dfnWorks changes into the job folder here)
DFN.make_working_directory(delete=True)
write_run_config(cfg, f"{jobname}/run_config.yaml")         # the merged settings this run used

DFN.print_domain_parameters()
DFN.check_input()
DFN.create_network()

# the mesh is only needed to paint the graph pressure onto the DFN (below)
if cfg["mesh"]:
    DFN.mesh_network()


# ## Graph flow: the -x end segments of the two laterals are the source and target
# intersection graph: nodes are fracture intersections, edges run along fractures.
# "left"/"right" are placeholders (nothing touches the faces); the source and target
# are then moved onto the well end segments.
G = DFN.create_graph("intersection", "left", "right")
G = DFN.add_fracture_source(G, [inj_frac])
G = DFN.add_fracture_target(G, [prod_frac])
G = DFN.run_graph_flow("injection_lateral", "production_lateral", p_in, p_out, G=G,
                       fluid_viscosity=fl["fluid_viscosity"],
                       conductance_model=fl["conductance_model"], mp_correction=fl["mp_correction"])

# every interior node should pass on what it receives
check_mass_balance(G, fl["fluid_viscosity"])

# paint the per-fracture mean of the intersection pressures onto reduced_mesh.inp
if cfg["mesh"]:
    map_graph_flow_onto_dfn(DFN, G, variable='pressure')

flow_rate = inlet_flow_rate(G)
print(f"--> Inlet flow rate: {flow_rate:0.4e} m^3/s at dP = {p_in - p_out:0.2e} Pa")

# effective fracture spacing for the thermal slab: geometric (1/P32) and flow-weighted
# (1/dQ), without the cased well segments; or the HF spacing
spacing = effective_spacing(DFN, G, exclude_fracs=() if wl["exchange_heat"] else wells["well_fracs"])
spacing["spacing_hf"] = hf_spacing
print(f"--> Effective spacing: 1/P32 = {spacing['spacing_p32']:0.1f} m, 1/dQ = {spacing['spacing_dQ']:0.1f} m, "
      f"HF spacing = {hf_spacing:0.1f} m (using thermal.spacing = {th['spacing']})")

# paint the flow channeling onto reduced_mesh.inp: each fracture's flow as a fraction of
# the doublet's flow (the per-fracture flow dQ is built from), linear and log10
if cfg["mesh"]:
    flow_fraction = spacing["Qf"] / flow_rate
    paint_fracture_values(DFN, "flow_fraction", flow_fraction)
    paint_fracture_values(DFN, "log10_flow_fraction", np.log10(np.maximum(flow_fraction, 1e-12)))

# ## Graph particle tracking (advective)
# edge area is flow-equivalent: aperture x the width implied by the Doolaeghe conductance
report_hf_widths(G, hf_fracs)

# Cased wells (wells.exchange_heat: false): particles move through the segments
# but exchange no heat there
no_exchange = None if wl["exchange_heat"] else wells["well_fracs"]
# flux-weighted release; purely advective (the heat exchange is applied analytically below)
particles = DFN.run_graph_transport(G, nparticles, "partime", None,
                                    initial_positions="flux", tdrw_flag=False,
                                    no_exchange_fractures=no_exchange)

# only particles that reach the production well count
exited = [p for p in particles if p.exit_flag]
tau = np.array([p.advect_time for p in exited])
beta = np.array([p.beta for p in exited])
# with one aperture everywhere, beta / (2 tau / b) is the fraction of travel time
# spent in fractures that exchange heat
beta_pp = 2.0 * tau / hf["aperture"]

print(f"\n--> Particles: {len(exited)} exited, {len(particles) - len(exited)} stuck")
for name, x in (("advect_time [s]", tau), ("beta [s/m]", beta), ("exchanging tau frac", beta / beta_pp)):
    print(f"--> {name:20s} min {x.min():0.3e}  median {np.median(x):0.3e}  max {x.max():0.3e}")

paths = particle_paths(exited, n_seg, hf_fracs)
shares = hf_shares(exited, hf_fracs)
report_flow_paths(paths, shares, hf_fracs, centers)

np.savez(f"{jobname}/particles.npz", advect_time=tau, beta=beta, paths=np.array(paths))
print(f"\n--> Particle data written to {jobname}/particles.npz")

# ## Thermal: produced temperature and heat rate from the particle (tau, beta)
# Rock between neighbouring HFs is shared, so by default each fracture draws on
# a slab of half-width hf_spacing / 2.
rock, fluid = th["rock"], th["fluid"]
if th["slab_half_width"] is not None:
    slab_half_width = th["slab_half_width"]
elif th["spacing"] in ("dQ", "p32", "hf"):
    slab_half_width = 0.5 * spacing[f"spacing_{th['spacing']}"]
else:
    raise ValueError(f"thermal.spacing must be dQ, p32 or hf, got {th['spacing']!r}")
t_eval = np.linspace(0.0, t_end_years * SECONDS_PER_YEAR, 12 * t_end_years + 1)   # monthly
F = thermal_breakthrough(tau, beta, t_eval, rock, fluid, slab_half_width)          # cooled fraction
m_dot = fluid["rho"] * flow_rate                                                    # mass flow rate [kg/s]
T_out = T_res - (T_res - T_inj) * F                                                 # produced temperature [C]
Q = m_dot * fluid["c"] * (T_res - T_inj) * (1.0 - F)                                # heat extraction rate [W]

print(f"\n--> Thermal: m_dot = {m_dot:0.3f} kg/s, slab half-width = {slab_half_width:0.1f} m")
print(f"    {'year':>5s} {'F':>7s} {'T_out [C]':>10s} {'Q [MW]':>9s}")
for year in (0.25, 0.5, 1, 2, 5, 10, 20, 30):
    if year > t_end_years:
        continue
    i = np.argmin(np.abs(t_eval - year * SECONDS_PER_YEAR))
    print(f"    {year:5.2f} {F[i]:7.3f} {T_out[i]:10.2f} {Q[i] / 1e6:9.3f}")

np.savez(f"{jobname}/thermal.npz", t_eval=t_eval, F=F, T_out=T_out, Q=Q, m_dot=m_dot,
         slab_half_width=slab_half_width, T_res=T_res, T_inj=T_inj)
print(f"--> Thermal curves written to {jobname}/thermal.npz")

# run inputs for plot_results.py
run_info = {"jobname": os.path.basename(jobname), "num_hf": num_hf, "hf_spacing": hf_spacing,
            "hf_radii": hf["radius"], "hf_aspect": hf["aspect"], "hf_centers": centers,
            "lateral_y": wells["lateral_y"], "n_seg": n_seg, "seg_length": wells["seg_length"],
            "seg_overlap": wells["seg_overlap"], "seg_width": wells["seg_width"], "aperture": hf["aperture"],
            "num_frac": int(DFN.num_frac), "inj_frac": inj_frac, "prod_frac": prod_frac,
            "p_in": p_in, "p_out": p_out, "flow_rate": float(flow_rate),
            "conductance_model": fl["conductance_model"], "mp_correction": fl["mp_correction"],
            "no_exchange": "none" if wl["exchange_heat"] else "well segments",
            "nparticles": nparticles, "n_exit": len(exited),
            "rock_props": rock, "fluid_props": fluid, "run_file": run_file}
with open(f"{jobname}/run_info.json", "w") as fp:
    json.dump(run_info, fp, indent=2)
print(f"--> Run inputs written to {jobname}/run_info.json")

# ## Total production: one row in the master CSV with every input setting (dotted
# names, so the swept ones are columns) and the run's production, for statistics
E = cumulative_heat(t_eval, Q)                                                      # [GWh]
E_no_cooling = Q[0] * t_eval[-1] / 3.6e12                                           # [GWh] outflow at T_res throughout


def T_out_at(year):
    """ Produced temperature at the output time nearest `year`; NaN past the end of the run """
    return float(T_out[np.argmin(np.abs(t_eval - year * SECONDS_PER_YEAR))]) if year <= t_end_years else np.nan


# inputs (every setting except the bookkeeping ones), then derived quantities and outputs
row = {"jobname": cfg["jobname"], "run_file": run_file or "",
       **flatten_config({k: v for k, v in cfg.items() if k not in ("jobname", "master_csv")}),
       "num_frac": int(DFN.num_frac), "hf_spacing": hf_spacing, "slab_half_width": slab_half_width,
       "p32": spacing["p32"], "dQ": spacing["dQ"], "spacing_p32_m": spacing["spacing_p32"],
       "spacing_dQ_m": spacing["spacing_dQ"],
       "flow_rate_m3_s": float(flow_rate), "m_dot_kg_s": float(m_dot),
       "E_total_GWh": float(E[-1]), "E_no_cooling_GWh": float(E_no_cooling),
       "E_fraction_of_no_cooling": float(E[-1] / E_no_cooling),
       "Q_initial_MW": float(Q[0]) / 1e6, "Q_end_MW": float(Q[-1]) / 1e6,
       "T_out_1yr_C": T_out_at(1), "T_out_10yr_C": T_out_at(10), "T_out_end_C": float(T_out[-1]),
       "tau_median_s": float(np.median(tau)), "beta_median_s_m": float(np.median(beta)),
       "hf_share_max": float(shares.max()),
       "hf_effective_number": float(1.0 / np.sum(shares**2)),   # 1: one HF carries everything; num_hf: even split
       "n_exit": len(exited), "n_stuck": len(particles) - len(exited)}
master_csv = os.path.join(src_path, cfg["master_csv"])
n_rows = append_master_csv(master_csv, row)
print(f"--> Total production {row['E_total_GWh']:0.2f} GWh over {t_end_years} years; "
      f"row added to {master_csv} ({n_rows} runs)")

# ## Gather the output data in results/<jobname>/ and remove the dfnWorks job folder
os.chdir(src_path)            # dfnWorks ran inside the job folder
gather_and_clean(jobname, os.path.join(src_path, cfg["results_dir"]), cleanup=cfg["cleanup"],
                 extra_files=[wells["polygon_file"]])

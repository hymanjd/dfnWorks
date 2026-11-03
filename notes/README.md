# Pickup notes: graph-based EGS heat extraction

Work moved from the `pydfnworks.thermal` dataset pipeline to a graph-based EGS
prototype that uses native dfnWorks calls. Two repos, both pushed 2 November 2026:

- `hymanjd/dfnWorks`, branch `feature/thermal-mpf-surrogate` (dfnWorks fixes, slides)
- `hymanjd/EGS_scripts`, branch `main`, folder `graph_based/` (the prototype)

Nothing is on `lanl/dfnWorks`.

## Where things are

| What | Where |
|---|---|
| Current driver (YAML-controlled) | `EGS_scripts/graph_based/driver_v2.py`; `driver_v3.py` adds background fractures |
| Defaults and sweep generator | `graph_based/make_yaml_file.py` (`DEFAULT_CONFIG`, `SWEEP`, `FIXED`, Latin hypercube) |
| Helpers (wells, thermal kernel, diagnostics, results) | `graph_based/helper_functions.py` |
| Plots | `graph_based/plot_results.py` (`<run>`, or `--compare 'lhs_hf_*'`) |
| Statistics | `graph_based/analysis.ipynb` (total thermal production vs swept inputs) |
| Gather finished runs into `results/` | `graph_based/gather_outputs.py` |
| Last sweep | `graph_based/runs/lhs_hf/` (30 runs: HF number 2-20, radius, aspect); summary `results/master_results.csv` |
| Old drivers | `graph_based/old/` |
| Slides | `dfnWorks/Documentation/graph_egs_slides/` (Beamer; `python make_figures.py`, then `pdflatex`) |

## Run

From `EGS_scripts/graph_based/`, with `python3.13` (pydfnworks installed from this
dfnWorks checkout):

    python3.13 driver_v2.py runs/base.yaml          # one run
    python3.13 make_yaml_file.py                    # write a sweep (edit SWEEP/FIXED first)
    bash runs/lhs_hf/run_all.sh                     # run it; progress in runs/lhs_hf/progress.log
    python3.13 plot_results.py --compare 'lhs_hf_*'

Each run gathers its outputs into `results/<jobname>/`, deletes the dfnWorks job
folder, and adds a row to `results/master_results.csv`. `run_all.sh` skips finished
runs, so rerunning resumes a sweep.

After changing pydfnworks, reinstall it into python3.13 (the install is a copy, not
a link; clear the stale `build/` first):

    cd pydfnworks && rm -rf build && python3.13 -m pip install -e . --break-system-packages --no-deps

## dfnWorks fixes on this branch (2 November)

- `add_fracture_target` accepts a bare int.
- Graph flow uses the same conductance for pressure, edge flow and area: node mass
  balance holds for the Doolaeghe model; edge area is flow-equivalent (`area_geometric`
  keeps the old value).
- `run_graph_transport(..., no_exchange_fractures=[...])`: no heat/matrix exchange on
  those fractures (cased wells).
- **User-fracture apertures follow dfnGen's insertion order (polygons first).** Before
  this, with polygon wells and ellipse HFs, the HF apertures landed on the first well
  segments. Any sweep that varied HF aperture before this fix is invalid (the old
  100-run sweep was deleted).
- Intersection-graph edges shorter than 1 mm are set to 1 mm (coincident intersection
  centers broke the flow solve with `KeyError: 'weight'`).

## State at the end of 2 November

- Thermal slab: `thermal.spacing` sets B = spacing/2 from `dQ` (default, 1/dQ from
  `DFN.compute_dQ`, cased wells excluded), `p32` (1/P32) or `hf` (HF spacing).
  `compute_dQ` normalizes by the whole domain volume, so 1/P32 is large for HF-only runs.
- Painting (`mesh: true`): per-fracture `pressure`, `flow_fraction` and
  `log10_flow_fraction` on `reduced_mesh.inp` (pydfnworks meshes in visual mode by
  default, so ~4000 fractures mesh in about a minute). Example: `runs/bg_pressure.yaml`
  -> `results/bg_pressure/reduced_mesh.inp`.
- Last background sweep (`driver_v3.py`, total_p32 = 0.288): 27 ok, 3 failed
  (`lhs_hf_11`, 20, 22) at the numbering check: `keepOnlyLargestCluster = True` dropped or
  renumbered an HF.
- Current sweeps barely cool (E / no-cooling >= 0.98): with HF aperture 1e-4 m and a
  1 MPa drop, m_dot is 0.07-0.5 kg/s, so production is flow-limited and the slab choice
  hardly matters yet.

## What to do next

1. Painting by node: give each mesh node the pressure of the nearest graph node on the
   same fracture (`add_variable_to_mesh(..., node_based=True)`); per-fracture means hide
   the pressure drop across each HF (every HF paints ~1.5 MPa).
2. `flow_fraction` is `compute_dQ`'s per-fracture Qf over the doublet flow; Qf sums flow
   across every intersection on a fracture, so values exceed 1 (wells up to ~4). Fine as a
   channeling ranking, not a true fraction: rename or normalize.
3. Cluster removal (`driver_v3.py`): keep all fractures (`keepOnlyLargestCluster = False`)
   and restrict the graph to the connected component holding source and target before
   `run_graph_flow`, so disconnected clusters never enter the solve (fixes the 3 failures).
4. Make cooling matter: target flow rate (one extra linear solve) or larger HF aperture /
   pressure drop; decouple HF number from spacing (sweep `hydraulic_fractures.length`).
5. Layout: producer at the +x end (`prod_frac = 2 n_seg`).
6. Validation: compare T_out(t) for two or three designs against the TOUGH UDFM runs in
   `EGS_scripts/ideal_HF/TOUGH_udfm/`.

## Carrying the Claude session over

The code is on the two forks; only the conversation is machine-local. From the
Claude desktop app: **Open in > Cloud** carries a summary and clones the branch
(uncommitted files do not travel), or use **Remote Control** while this Mac stays on.
A fresh session has what it needs in this file and the script headers.

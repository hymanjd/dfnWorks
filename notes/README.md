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

## What to do next

1. Slides: update `graph_egs_slides`, dropping the dfnWorks-development material
   (graph-flow fix slide, verification table rows, before/after fixes slide) and adding
   the sweep results.
2. Background fractures (`driver_v3.py`): keep `keepOnlyLargestCluster = False`, or the
   cluster removal can delete the wells and HFs. Run a sweep with it.
3. Layout: injection and production are both at the -x end, so the HF nearest the
   well ends takes most of the flow. Try the producer at the +x end
   (`prod_frac = 2 n_seg`).
4. Rock contact: compare the flow-equivalent width with a whole-fracture sweep
   (2 A_f / Q_f per fracture) as an upper bound; HF radius/aspect barely affect
   production under the current model.
5. Re-sweep HF aperture and pressure drop now that apertures are assigned correctly.

## Carrying the Claude session over

The code is on the two forks; only the conversation is machine-local. From the
Claude desktop app: **Open in > Cloud** carries a summary and clones the branch
(uncommitted files do not travel), or use **Remote Control** while this Mac stays on.
A fresh session has what it needs in this file and the script headers.

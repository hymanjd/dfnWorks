# EGS heat-extraction sweep (pydfnworks graph flow + thermal)

Runs a parameter sweep of an enhanced geothermal system (EGS) doublet and
asks which design parameters control total thermal production.

Each run builds a DFN with two segmented lateral wells (injection and
production, at y = -50 and y = +50 m by default), a set of evenly spaced
hydraulic fractures (HFs) between them, and four families of background
fractures. It then:

1. solves graph flow from the injection lateral to the production lateral,
2. tracks particles on the intersection graph (advective travel time `tau`
   and fracture-rock exchange integral `beta` for each particle),
3. turns the particle paths into the produced temperature `T_out(t)`, the heat
   rate `Q(t)` and the total production `E` over the run, using a
   finite-slab heat-exchange kernel.

No PFLOTRAN or FEHM is needed; meshing is optional and off in sweeps.
Reservoir numbers are representative of Utah FORGE (226 C granitoid, 40 C
injection, 30 years).

## Files

| File | Role |
|---|---|
| `make_yaml_file.py` | Defines the default run (`DEFAULT_CONFIG`) and the sweep; writes one YAML config per run |
| `driver.py` | Runs one simulation from a YAML config |
| `helper_functions.py` | Wells, config loading, thermal kernel, output gathering |
| `run_batch.sh` | Generates the configs, then runs `driver.py` on each in turn |
| `analysis.ipynb` | Statistics on the finished sweep (regression, random forest, best/worst designs) |

## Workflow

### 1. Set up the sweep

Edit the sweep settings in `make_yaml_file.py`:

- `SWEEP`: the parameters to vary, each as `dotted.key: (low, high, scale)`,
  where scale is `"log"`, `"linear"` or `"int"`. Currently: HF number (2-20),
  HF radius (60-100 m), HF aspect ratio (1-3) and HF aperture (1e-5 to 1e-3 m, log).
- `FIXED`: settings held for every run in this sweep (e.g. `mesh: False`).
- `N_samples`: number of Latin hypercube points; `lhs_seed`: the sampling seed.
- `seeds`: `None` gives each run its own DFN seed (`seed_base + index`), or a
  list of seeds to repeat every point with.
- `sweep_name`: names the config folder and each run's jobname (`<sweep_name>_<N>`).

Any key in `DEFAULT_CONFIG` can be swept or fixed. A typo or a bad value
stops the script before it writes any files.

Note on apertures: the sampled `hydraulic_fractures.aperture` is shared by all
HFs within a run. The well segments use `wells.aperture` (1e-3 m), and the
background fractures get a radius-correlated aperture, b = 1.7e-6 r^0.5
(about 8e-6 to 2e-5 m), set in `driver.py`.

### 2. Generate the configs

    python3.13 make_yaml_file.py

This writes to `runs/<sweep_name>/`:

- `config_1.yaml ... config_N.yaml`: each lists only the settings that differ
  from `DEFAULT_CONFIG` (fixed settings, sampled values, jobname, seed)
- `samples.csv`: the sampled values of every run
- `run_all.sh`: an alternative runner (see below)

Re-running the script overwrites the configs. Adding or removing a swept
parameter, or changing `N_samples` or `lhs_seed`, changes every sample, not
just the new one.

### 3. Run the sweep

From this folder:

    bash run_batch.sh

This regenerates the configs and then runs `driver.py` on each config in turn.
The loop count in `run_batch.sh` must match `N_samples`.

To run a single case:

    python3.13 driver.py runs/lhs_hf/config_7.yaml

With no argument, `driver.py` runs `DEFAULT_CONFIG` alone.

`runs/<sweep_name>/run_all.sh` is a more careful runner. It skips runs already
in `results/` (so an interrupted sweep resumes), carries on past failures,
logs progress to `runs/<sweep_name>/progress.log` and prints a time estimate.

### 4. Outputs

Each run writes its dfnWorks job folder `<jobname>/` and `<jobname>.log` in
this folder. When the run finishes, its outputs are gathered into
`results/<jobname>/`:

- `run_config.yaml`: the full merged settings the run used
- `run_info.json`: geometry, flow and particle summary
- `thermal.npz`: `t_eval`, cooled fraction `F`, `T_out`, `Q`, `m_dot`, slab half-width
- `particles.npz`, `graph_flow.hdf5`, `polygons.dat`, `<jobname>.log`
  (and `reduced_mesh.inp` if meshing was on)

The job folder is then deleted (`cleanup: True`), but only if every output
copied intact. A `<jobname>/` folder left behind in this folder means that run
failed or did not finish; its `.log` says why.

Each finished run also adds one row to `results/master_results.csv`: every
input setting (dotted names such as `hydraulic_fractures.aperture`) plus
derived and output quantities. These include `E_total_GWh`,
`E_fraction_of_no_cooling`, `Q_initial_MW`, `Q_end_MW`, `T_out` at 1 and 10 years
and at the end, median `tau` and `beta`, the HF flow shares and the
exited/stuck particle counts. Re-running a jobname replaces its row.

Start each new sweep with an empty `results/` (or a new `sweep_name`).
Otherwise its rows mix with the old sweep's in the master CSV.

### 5. Analyse

Open `analysis.ipynb` and run it from this folder. It reads
`results/master_results.csv`, finds the inputs that vary across runs and
relates them to total production `E_total_GWh`. It needs no edits for a
different sweep.

## Requirements

- pydfnworks from this branch (`cd pydfnworks && pip install -e .`)
- dfnGen, with its path set in `~/.dfnworksrc`
- Python packages: numpy, scipy, pandas, pyyaml; for the notebook also
  matplotlib and scikit-learn

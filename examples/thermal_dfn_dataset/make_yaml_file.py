"""
Sweep driver_v2.py runs with parameters drawn by Latin hypercube:
python3.13 make_yaml_file.py

N_samples points (scipy.stats.qmc): each of the N equal-width bins of every swept
parameter (in log10 for "log", linear for "linear" and "int") holds exactly one point.
Each point is crossed with seeds (if a list). Writes runs/<sweep_name>/config_N.yaml,
samples.csv and run_all.sh; then run: bash runs/<sweep_name>/run_all.sh

A config lists only what differs from DEFAULT_CONFIG (below): the fixed
settings below, the sampled values and the run's jobname and seed.

To set up a sweep, edit the settings below: DEFAULT_CONFIG (every run's baseline),
sweep_name, SWEEP (what to vary and over what range), FIXED (what to hold for this
sweep), N_samples, lhs_seed, seeds and python_exe.
"""
import os
import sys
import itertools
import yaml
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import qmc

from helper_functions import load_run_config, _merge, GRANITE_PROPS, WATER_PROPS

# Run inputs for driver_v2.py. A run's YAML file overrides any of these; keys it
# leaves out keep these values (helper_functions.load_run_config).
DEFAULT_CONFIG = {
    "jobname": "graph_sim_v2",
    "ncpu": 8,
    "mesh": True,                      # mesh_network, needed to paint the graph pressure on the DFN
    "master_csv": "results/master_results.csv",  # one row per run (inputs + total production), relative to where the driver runs
    "results_dir": "results",          # each run's output data is gathered in results/<jobname>/
    "cleanup": True,                   # delete the dfnWorks job folder once its outputs are gathered
    "domain": {"size": [500, 500, 250], "center": [0, 0, -150],
               "h": 0.5,               # must be below wells.seg_overlap
               "seed": 10},
    "hydraulic_fractures": {"number": 3,
                            "length": 300,     # the HFs are spaced evenly over x in [-length/2, length/2]
                            "radius": 70, "aspect": 2, "aperture": 1.0e-4},
    "wells": {"lateral_y": 50, "lateral_x": [-200, 200], "seg_length": 20,
              "seg_overlap": 2, "seg_width": 4, "aperture": 1.0e-3,
              "exchange_heat": False},  # False: particles cross the segments with no heat exchange
    "flow": {"p_in": 2.0e6, "p_out": 1.0e6, "fluid_viscosity": 8.9e-4,
             "conductance_model": "doolaeghe", "mp_correction": "divide"},
    "transport": {"nparticles": 10000},
    "thermal": {"T_res": 226.0, "T_inj": 40.0, "t_end_years": 30,
                "spacing": "dQ",          # slab half-width = spacing/2 from: dQ (1/dQ, flow channeling), p32 (1/P32), hf (HF spacing)
                "slab_half_width": None,  # [m]; a number overrides thermal.spacing
                "rock": GRANITE_PROPS, "fluid": WATER_PROPS},
}

# ## Sweep settings
sweep_name = "lhs_hf"       # names the output folder runs/<sweep_name>/ and each run's jobname

# Swept settings: dotted key in DEFAULT_CONFIG -> (low, high, scale)
#   "log": uniform in log10, "linear": uniform, "int": uniform then rounded
SWEEP = {
    "hydraulic_fractures.number":   (2, 20, "int"),
    "hydraulic_fractures.radius":   (60.0, 100.0, "linear"),   # y half-extent = radius * aspect must exceed wells.lateral_y
    "hydraulic_fractures.aspect":   (1.0, 3.0, "linear"),
    # "hydraulic_fractures.aperture": (1.0e-4, 1.0e-3, "log"),
    # "flow.p_in":                    (1.1e6, 1.0e7, "log"),      # p_out = 1e6 Pa
}

# Settings fixed for every run in the sweep (dotted key -> value)
FIXED = {
    "transport.nparticles": 10000,
    "mesh": False,                 # skip mesh_network: no pressure painting in a sweep
}

N_samples = 30              # number of LHS points
lhs_seed  = 7               # sampling seed (not the DFN seed)
seeds     = None            # None: own DFN seed per run, seed = seed_base + index; or a list to cross with
seed_base = 1000
python_exe = "python3.13"   # interpreter run_all.sh calls driver_v2.py with


# ## Sampling

def lhs_points(N=N_samples, sweep=SWEEP, seed=lhs_seed):
    """ N points from a Latin hypercube over the swept settings; dict key -> array """
    keys = list(sweep)
    # u: N x d points in the unit cube, one per row of N equal bins in every dimension
    u = qmc.LatinHypercube(d=len(keys), seed=seed).random(N)
    points = {}
    # map each unit-interval column onto its setting's range and scale
    for j, key in enumerate(keys):
        low, high, scale = sweep[key]
        if scale == "log":
            points[key] = 10.0 ** (np.log10(low) + u[:, j] * (np.log10(high) - np.log10(low)))
        elif scale == "linear":
            points[key] = low + u[:, j] * (high - low)
        elif scale == "int":
            # equal-width bins over [low - 0.5, high + 0.5] so every integer is equally likely
            points[key] = np.floor(low - 0.5 + u[:, j] * (high - low + 1) + 0.5).astype(int)
        else:
            sys.exit(f"Error: unknown scale '{scale}' for {key}; use log, linear or int")
    return points


# ## Building one run's settings

def set_key(config, dotted, value):
    """ Set config['a']['b'] for dotted = 'a.b', creating nested mappings """
    keys = dotted.split(".")
    node = config
    for k in keys[:-1]:
        node = node.setdefault(k, {})
    node[keys[-1]] = value


def make_config(idx, sample, seed, fixed=FIXED):
    """ One run's overrides: the fixed settings, the sampled values, jobname and seed.
    Derived settings (e.g. one parameter set from others) go here. """
    config = {"jobname": f"{sweep_name}_{idx}"}
    for key, value in fixed.items():
        set_key(config, key, value)
    for key, value in sample.items():
        # numpy scalars -> plain Python numbers, so the YAML stays readable
        set_key(config, key, value.item() if hasattr(value, "item") else value)
    set_key(config, "domain.seed", int(seed))
    return config


# ## Writing the sweep

def generate_configs(N=N_samples, sweep=SWEEP, lhs_seed=lhs_seed, seeds=seeds, seed_base=seed_base,
                     output_dir=None, manifest_path=None):
    """ One YAML config per (LHS point, seed); returns the files and the manifest """
    output_dir = Path(output_dir or Path("runs") / sweep_name)
    manifest_path = Path(manifest_path or output_dir / "samples.csv")

    points = lhs_points(N, sweep, lhs_seed)
    samples = [{k: points[k][i] for k in sweep} for i in range(N)]   # one dict per LHS point

    configs, records = [], []
    # LHS gives N points: cross each with the seeds, don't cross the parameters
    for idx, (sample, seed) in enumerate(itertools.product(samples, seeds or [None]), start=1):
        if seed is None:
            seed = seed_base + idx                 # own DFN realization per run
        config = make_config(idx, sample, seed)
        try:
            _merge(DEFAULT_CONFIG, config)         # unknown keys or bad numbers
        except (KeyError, TypeError) as e:
            sys.exit(f"Error: {e.args[0]}")
        configs.append(config)
        # manifest row: run index, jobname, seed and the sampled values
        records.append({"index": idx, "jobname": config["jobname"], "seed": int(seed),
                        **{k: config_value(config, k) for k in sweep}})

    # write files only after every config passed the check above
    output_dir.mkdir(parents=True, exist_ok=True)
    config_list = []
    for idx, config in enumerate(configs, start=1):
        filepath = output_dir / f"config_{idx}.yaml"
        with open(filepath, "w") as f:
            f.write(f"# {sweep_name} run {idx}; settings not listed take make_yaml_file.DEFAULT_CONFIG values\n")
            yaml.safe_dump(config, f, sort_keys=False, default_flow_style=None)
        load_run_config(filepath)                  # round trip: the file loads as written
        config_list.append(filepath)

    manifest = pd.DataFrame.from_records(records)
    manifest.to_csv(manifest_path, index=False)

    results_dir = fixed_results_dir()
    write_run_all(output_dir, [(fp, c["jobname"]) for fp, c in zip(config_list, configs)], results_dir)

    return config_list, manifest


def fixed_results_dir(fixed=FIXED):
    """ The results folder the sweep's runs gather into (FIXED can override the default) """
    return fixed.get("results_dir", DEFAULT_CONFIG["results_dir"])


# bash for run_all.sh; @...@ are filled in by write_run_all
RUN_ALL_HEAD = """#!/bin/bash
# @NAME@: @TOTAL@ runs of driver_v2.py. Run from the folder holding driver_v2.py:
#   bash @DIR@/run_all.sh
# Prints each run's progress and logs it to @DIR@/progress.log. A run whose
# @RESULTS@/<jobname>/thermal.npz exists is skipped, so rerunning resumes an
# interrupted sweep. A failed run is reported and the sweep goes on; its console
# output stays in <jobname>.out. A finished run's .out moves into @RESULTS@/<jobname>/.

PY=@PY@
LOG=@DIR@/progress.log
RESULTS=@RESULTS@
TOTAL=@TOTAL@
n_ok=0; n_fail=0; n_skip=0; t_runs=0; failed=""
START=$(date +%s)

hms() { printf "%dh%02dm%02ds" $(($1 / 3600)) $(($1 % 3600 / 60)) $(($1 % 60)); }

run() {
    # run <index> <jobname> <config file>
    local i=$1 job=$2 cfg=$3
    if [ -f "$RESULTS/$job/thermal.npz" ]; then
        echo "[$i/$TOTAL] $job: already in $RESULTS/, skipping" | tee -a "$LOG"
        n_skip=$((n_skip + 1)); return
    fi
    echo "[$i/$TOTAL] $job: started $(date '+%Y-%m-%d %H:%M:%S')  ($cfg)" | tee -a "$LOG"
    local t0=$(date +%s) status
    if $PY driver_v3.py "$cfg" > "$job.out" 2>&1; then
        status="ok"; n_ok=$((n_ok + 1))
        [ -d "$RESULTS/$job" ] && mv "$job.out" "$RESULTS/$job/"
    else
        status="FAILED (see $job.out)"; n_fail=$((n_fail + 1)); failed="$failed $job"
    fi
    local dt=$(( $(date +%s) - t0 ))
    t_runs=$((t_runs + dt))
    local n_run=$((n_ok + n_fail)) left=$((TOTAL - i))
    local eta=$(( left * t_runs / n_run ))
    echo "[$i/$TOTAL] $job: $status in $(hms $dt) | ok $n_ok, failed $n_fail, skipped $n_skip | elapsed $(hms $(( $(date +%s) - START ))), ~$(hms $eta) left" | tee -a "$LOG"
    if [ "$status" = "ok" ]; then
        grep -h "Total production" "$RESULTS/$job/$job.out" 2>/dev/null | sed 's/^/    /'
    fi
}

echo "=== @NAME@: $TOTAL runs, started $(date '+%Y-%m-%d %H:%M:%S') ===" | tee -a "$LOG"
"""

RUN_ALL_TAIL = """
echo "=== @NAME@ done $(date '+%Y-%m-%d %H:%M:%S'): ok $n_ok, failed $n_fail, skipped $n_skip, total $(hms $(( $(date +%s) - START ))) ===" | tee -a "$LOG"
[ -n "$failed" ] && echo "Failed:$failed" | tee -a "$LOG"
[ "$n_fail" -eq 0 ]
"""


def write_run_all(output_dir, runs, results_dir="results", name=None, py=None):
    """ Write output_dir/run_all.sh: runs the sweep's configs in turn with progress
    printing, a progress.log, skipping of finished runs and a summary at the end

    runs is a list of (config file, jobname), in run order.
    """
    output_dir = Path(output_dir)
    fill = {"@NAME@": name or output_dir.name, "@TOTAL@": str(len(runs)), "@DIR@": str(output_dir),
            "@RESULTS@": str(results_dir), "@PY@": py or python_exe}
    text = RUN_ALL_HEAD
    for token, value in fill.items():
        text = text.replace(token, value)
    for i, (filepath, jobname) in enumerate(runs, start=1):
        text += f"run {i} {jobname} {filepath}\n"
    tail = RUN_ALL_TAIL
    for token, value in fill.items():
        tail = tail.replace(token, value)
    script = output_dir / "run_all.sh"
    script.write_text(text + tail)
    os.chmod(script, 0o755)
    return script


def config_value(config, dotted):
    """ Read config['a']['b'] for dotted = 'a.b' """
    node = config
    for k in dotted.split("."):
        node = node[k]
    return node


if __name__ == "__main__":

    configs, manifest = generate_configs()
    print(f"Generated {len(configs)} config files in {configs[0].parent} and a manifest with {len(manifest)} rows:")
    print(manifest.to_string(index=False))
    print(f"Run them all with: bash {configs[0].parent / 'run_all.sh'}")

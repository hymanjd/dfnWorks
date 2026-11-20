#!/usr/bin/env python3
"""Does the wiggle option change how well accepted radii match the prescribed distribution?

Generates the same stochastic family many times (one run per seed) for every
combination of wiggle on/off and rFram on/off, then compares the radii of all
accepted fractures with the analytic distribution the family was given.

Only DFNGen runs (no meshing), so a full sweep takes a few minutes.

Usage
-----
    python3 run_sweep.py                      # 20 seeds, all distributions
    python3 run_sweep.py --seeds 50 --jobs 8
    python3 run_sweep.py --dists tpl exp --out /some/dir
    python3 run_sweep.py --npoly 120          # fixed fracture count, no P32 target

Outputs (in --out, default ./_sweep)
------------------------------------
    runs.csv      one row per DFNGen run
    summary.csv   one row per (distribution, rFram, wiggle)
    summary.md    the same table, readable

Columns
-------
    ks            Kolmogorov-Smirnov distance between accepted radii and the
                  prescribed distribution (0 = perfect match). "ks_run" is the
                  mean over seeds of the per-run distance; "ks_pooled" pools
                  the radii of every seed first, which removes most of the
                  sampling noise of a single small network.
    mean_ratio    mean accepted radius / mean of the prescribed distribution.
    r_smallest    smallest accepted radius, averaged over seeds. DFNGen sorts its
                  sampled radii largest first and stops at the P32 target, so
                  when fewer fractures are rejected the smallest sampled radii
                  are never inserted and this rises above the minimum radius.
    ks_on_vs_off  two-sample KS distance between pooled wiggle-on and
                  wiggle-off radii at the same rFram setting.
"""

import argparse
import csv
import math
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

MIN_RADIUS = 1.0
MAX_RADIUS = 10.0

# Family definitions. Each has an analytic CDF on [MIN_RADIUS, MAX_RADIUS].
FAMILIES = {
    "tpl": dict(distribution="tpl", alpha=1.8),
    "log_normal": dict(distribution="log_normal", log_mean=math.log(2.0), log_std=0.5),
    "exp": dict(distribution="exp", exp_mean=2.5),
}


def _norm_cdf(x):
    return 0.5 * (1.0 + np.vectorize(math.erf)(x / math.sqrt(2.0)))


def analytic_cdf(dist, r):
    """CDF of the prescribed radius distribution, truncated to [MIN_RADIUS, MAX_RADIUS]."""
    r = np.asarray(r, dtype=float)
    fam = FAMILIES[dist]

    if dist == "tpl":
        raw = lambda x: 1.0 - (MIN_RADIUS / x)**fam["alpha"]
    elif dist == "log_normal":
        raw = lambda x: _norm_cdf((np.log(x) - fam["log_mean"]) / fam["log_std"])
    elif dist == "exp":
        raw = lambda x: 1.0 - np.exp(-x / fam["exp_mean"])
    else:
        raise ValueError(dist)

    lo, hi = raw(np.array(MIN_RADIUS)), raw(np.array(MAX_RADIUS))
    return np.clip((raw(r) - lo) / (hi - lo), 0.0, 1.0)


def analytic_mean(dist):
    """Mean of the prescribed distribution, by numerical integration of 1 - CDF."""
    x = np.linspace(MIN_RADIUS, MAX_RADIUS, 20001)
    y = 1.0 - analytic_cdf(dist, x)
    return MIN_RADIUS + float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def ks_one_sample(radii, dist):
    r = np.sort(np.asarray(radii))
    n = len(r)
    cdf = analytic_cdf(dist, r)
    return float(max(np.max(np.arange(1, n + 1) / n - cdf), np.max(cdf - np.arange(0, n) / n)))


def ks_two_sample(a, b):
    a, b = np.sort(a), np.sort(b)
    grid = np.concatenate([a, b])
    return float(np.max(np.abs(np.searchsorted(a, grid, side="right") / len(a)
                               - np.searchsorted(b, grid, side="right") / len(b))))


# ---------------------------------------------------------------- one DFNGen run

def worker(dist, rfram, wiggle, seed, run_dir, p32, domain, h, npoly, list_increase):
    """Runs inside its own process: generate one network, no meshing."""
    from pydfnworks import DFNWORKS

    os.makedirs(run_dir, exist_ok=True)
    os.chdir(run_dir)
    DFN = DFNWORKS(os.path.join(run_dir, "output"), ncpu=1)
    DFN.params['domainSize']['value'] = [domain, domain, domain]
    DFN.params['h']['value'] = h
    DFN.params['domainSizeIncrease']['value'] = [0.5, 0.5, 0.5]
    # Keep every accepted fracture so cluster removal does not enter the comparison
    DFN.params['keepIsolatedFractures']['value'] = True
    DFN.params['ignoreBoundaryFaces']['value'] = True
    DFN.params['seed']['value'] = seed
    DFN.params['rFram']['value'] = bool(rfram)
    DFN.params['wiggle']['value'] = bool(wiggle)

    if list_increase >= 0:
        DFN.params['radiiListIncrease']['value'] = list_increase

    if npoly > 0:
        # Fixed number of fractures instead of a P32 target
        DFN.params['stopCondition']['value'] = 0
        DFN.params['nPoly']['value'] = npoly
        stop = dict(probability=1.0)
    else:
        stop = dict(p32=p32)

    DFN.add_fracture_family(shape="ell",
                            orientation_distribution="fisher",
                            kappa=1.0,
                            theta=0.0,
                            phi=0.0,
                            min_radius=MIN_RADIUS,
                            max_radius=MAX_RADIUS,
                            number_of_points=8,
                            hy_variable='permeability',
                            hy_function='constant',
                            hy_params={"mu": 1e-12},
                            **stop,
                            **FAMILIES[dist])
    DFN.make_working_directory(delete=True)
    DFN.check_input()
    DFN.create_network()


def read_run(run_dir):
    """Accepted radii and generation counters of one finished run, or None if it failed."""
    gen = os.path.join(run_dir, "output", "dfnGen_output")
    radii_file = os.path.join(gen, "radii.dat")
    report = os.path.join(run_dir, "output", "DFN_output.txt")

    if not (os.path.isfile(radii_file) and os.path.isfile(report)):
        return None

    # radii.dat lists every accepted fracture: xRadius yRadius family [R]
    radii = [float(line.split()[0]) for line in open(radii_file).readlines()[1:] if line.strip()]
    text = open(report).read()

    def grab(pattern, cast=int, default=0):
        m = re.search(pattern, text)
        return cast(m.group(1)) if m else default

    return dict(radii=np.array(radii),
                rejected=grab(r"Total Fractures Rejected: (\d+)"),
                retranslated=grab(r"Total Fractures Re-translated: (\d+)"),
                wiggled=grab(r"Total Fractures Wiggled: (\d+)"),
                p32=grab(r"Total Fracture Intensity \(P32\): ([0-9.eE+-]+)", float, float("nan")))


def launch(task, args):
    dist, rfram, wiggle, seed = task
    run_dir = os.path.join(args.out, f"{dist}_rfram{rfram}_wiggle{wiggle}", f"seed{seed}")
    os.makedirs(run_dir, exist_ok=True)
    cmd = [sys.executable, os.path.abspath(__file__), "--worker", dist, str(rfram), str(wiggle),
           str(seed), run_dir, str(args.p32), str(args.domain), str(args.h), str(args.npoly),
           str(args.radii_list_increase)]

    with open(os.path.join(run_dir, "run.log"), "w") as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                            env=dict(os.environ, MPLBACKEND="Agg")).returncode

    return task, run_dir, rc


# ---------------------------------------------------------------- sweep

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        dist, rfram, wiggle, seed, run_dir, p32, domain, h, npoly, list_increase = sys.argv[2:12]
        worker(dist, int(rfram), int(wiggle), int(seed), run_dir, float(p32), float(domain), float(h),
               int(npoly), float(list_increase))
        return

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", type=int, default=20, help="number of seeds per configuration")
    parser.add_argument("--dists", nargs="+", default=list(FAMILIES), choices=list(FAMILIES))
    parser.add_argument("--p32", type=float, default=1.0, help="target P32 of the family")
    parser.add_argument("--npoly", type=int, default=0,
                        help="insert exactly this many fractures instead of stopping at --p32")
    parser.add_argument("--radii-list-increase", type=float, default=-1.0,
                        help="override radiiListIncrease (default: leave the dfnWorks default of 0.1)")
    parser.add_argument("--domain", type=float, default=15.0, help="cubic domain edge length")
    parser.add_argument("--h", type=float, default=0.1)
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    parser.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "_sweep"))
    args = parser.parse_args()
    args.out = os.path.abspath(args.out)

    tasks = [(d, rf, w, s) for d in args.dists for rf in (0, 1) for w in (0, 1)
             for s in range(1, args.seeds + 1)]
    print(f"{len(tasks)} DFNGen runs, {args.jobs} at a time -> {args.out}")

    with ThreadPoolExecutor(args.jobs) as pool:
        finished = list(pool.map(lambda t: launch(t, args), tasks))

    rows, pooled, failed = [], {}, 0

    for (dist, rfram, wiggle, seed), run_dir, rc in finished:
        data = read_run(run_dir)

        if data is None or len(data["radii"]) == 0:
            failed += 1
            print(f"FAILED rc={rc}: {run_dir}")
            continue

        radii = data["radii"]
        rows.append(dict(dist=dist, rfram=rfram, wiggle=wiggle, seed=seed, n=len(radii),
                         ks=ks_one_sample(radii, dist),
                         mean_ratio=float(radii.mean()) / analytic_mean(dist),
                         r_smallest=float(radii.min()),
                         rejected=data["rejected"], retranslated=data["retranslated"],
                         wiggled=data["wiggled"], p32=data["p32"]))
        pooled.setdefault((dist, rfram, wiggle), []).append(radii)

    with open(os.path.join(args.out, "runs.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = []

    for dist in args.dists:
        for rfram in (0, 1):
            for wiggle in (0, 1):
                sel = [r for r in rows if (r["dist"], r["rfram"], r["wiggle"]) == (dist, rfram, wiggle)]

                if not sel:
                    continue

                allr = np.concatenate(pooled[(dist, rfram, wiggle)])
                other = pooled.get((dist, rfram, 0))
                col = lambda k: np.array([r[k] for r in sel], dtype=float)
                summary.append(dict(
                    dist=dist, rfram=rfram, wiggle=wiggle, runs=len(sel),
                    n_per_run=round(col("n").mean(), 1),
                    ks_run=round(col("ks").mean(), 4), ks_run_sd=round(col("ks").std(ddof=1) if len(sel) > 1 else 0.0, 4),
                    ks_pooled=round(ks_one_sample(allr, dist), 4),
                    mean_ratio=round(float(allr.mean()) / analytic_mean(dist), 4),
                    r_smallest=round(col("r_smallest").mean(), 3),
                    ks_on_vs_off=round(ks_two_sample(allr, np.concatenate(other)), 4) if wiggle and other else "",
                    retranslated=round(col("retranslated").mean(), 1),
                    rejected=round(col("rejected").mean(), 1),
                    wiggled=round(col("wiggled").mean(), 1),
                    p32=round(col("p32").mean(), 4)))

    with open(os.path.join(args.out, "summary.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    keys = list(summary[0])
    lines = ["| " + " | ".join(keys) + " |", "|" + "---|" * len(keys)]
    lines += ["| " + " | ".join(str(s[k]) for k in keys) + " |" for s in summary]
    table = "\n".join(lines)

    with open(os.path.join(args.out, "summary.md"), "w") as f:
        stop = f"nPoly = {args.npoly}" if args.npoly > 0 else f"p32 = {args.p32}"
        extra = f", radiiListIncrease = {args.radii_list_increase}" if args.radii_list_increase >= 0 else ""
        f.write(f"{stop}{extra}, domain = {args.domain}, h = {args.h}, radii in [{MIN_RADIUS}, {MAX_RADIUS}], "
                f"{args.seeds} seeds per row. Per-run values are means over seeds.\n\n{table}\n")

    print(table)
    print(f"\n{failed} failed run(s). Results in {args.out}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()

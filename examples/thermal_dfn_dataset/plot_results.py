#"""
#   :synopsis: Plot the graph-based EGS results written by driver_v1.py / driver_v2.py:
#                - plot_run: one run's six-panel figure (summary, produced temperature,
#                  heat rate, cumulative heat, travel times, share of flow per HF)
#                - compare_production: total heat production across many runs
#              A run folder is a dfnWorks job folder or its gathered copy in results/;
#              it needs thermal.npz (both plots) and particles.npz (plot_run).
#   :version: 0.1
#   :maintainer: Jeffrey Hyman
#.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
#
#   usage: python plot_results.py <run>                 -> <run folder>/results.png
#          python plot_results.py --compare [<run> ...] [--out production_comparison.png]
#
#   A <run> is a folder path, a run name looked up in results/ (lhs_hf_1), or a
#   pattern matched in results/ (quote it so the shell leaves it alone):
#          python plot_results.py lhs_hf_1
#          python plot_results.py --compare 'lhs_hf_*'      # one sweep
#          python plot_results.py --compare                 # every run in results/
#"""

import os
import re
import sys
import glob
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from helper_functions import thermal_breakthrough, cumulative_heat, GRANITE_PROPS, WATER_PROPS, SECONDS_PER_YEAR

RESULTS_DIR = "results"     # where driver_v2.py gathers each run (DEFAULT_CONFIG results_dir)

# ## Figure style
# Categorical slots in fixed order (one per HF), plus text and grid inks
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
REFERENCE = "#8a8984"

plt.rcParams.update({
    "text.usetex": False,       # importing pydfnworks turns LaTeX on, which eats '%'
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": TEXT_2, "axes.titlecolor": TEXT,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.labelsize": 9, "xtick.color": TEXT_2, "ytick.color": TEXT_2,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 2, "legend.frameon": False, "legend.fontsize": 8,
})

# ## One run

def plot_run(jobdir):
    """ Six-panel figure for one run: summary, temperature, heat rate,
    cumulative heat, travel times and HF shares; writes <jobdir>/results.png """
    th = np.load(os.path.join(jobdir, "thermal.npz"))
    pa = np.load(os.path.join(jobdir, "particles.npz"))
    # run inputs, written by driver_v1.py / driver_v2.py; older runs do not have them
    info_file = os.path.join(jobdir, "run_info.json")
    info = json.load(open(info_file)) if os.path.isfile(info_file) else None

    years = th["t_eval"] / SECONDS_PER_YEAR
    T_res, T_inj = float(th["T_res"]), float(th["T_inj"])
    tau, beta, paths = pa["advect_time"], pa["beta"], pa["paths"]

    # HF each particle went through (paths look like 'inj-HF1-prod')
    hf_of = np.array([int(m.group(1)) if (m := re.search(r"HF(\d+)", str(p))) else 0 for p in paths])
    hfs = sorted(set(hf_of) - {0})
    share = {k: np.mean(hf_of == k) for k in hfs}

    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    ax_info, ax_T, ax_Q, ax_E, ax_tau, ax_hf = axes.ravel()

    # Produced temperature, with each HF's share of the mixed outflow
    ax_T.plot(years, th["T_out"], color=TEXT, label="mixed outflow")
    for k in hfs:
        sel = hf_of == k
        F_k = thermal_breakthrough(tau[sel], beta[sel], th["t_eval"], GRANITE_PROPS, WATER_PROPS,
                                   float(th["slab_half_width"]))
        T_k = T_res - (T_res - T_inj) * F_k
        c = SERIES[(k - 1) % len(SERIES)]
        ax_T.plot(years, T_k, color=c, linewidth=1.5, label=f"HF{k} ({100 * share[k]:.1f}% of flow)")
        ax_T.annotate(f"HF{k}", (years[-1], T_k[-1]), xytext=(4, 0), textcoords="offset points",
                      va="center", fontsize=8, color=TEXT_2)
    for T_ref, name in ((T_res, "T_res"), (T_inj, "T_inj")):
        ax_T.axhline(T_ref, color=REFERENCE, linewidth=1, linestyle=(0, (4, 3)))
        ax_T.annotate(f"{name} {T_ref:.0f} °C", (0, T_ref), xytext=(4, 3), textcoords="offset points",
                      fontsize=8, color=TEXT_2)
    ax_T.set_title("Produced temperature")
    ax_T.set_xlabel("time [years]")
    ax_T.set_ylabel("T_out [°C]")
    ax_T.set_xlim(0, years[-1] * 1.06)
    ax_T.legend(loc="center right", bbox_to_anchor=(0.94, 0.40))

    # Heat extraction rate
    ax_Q.plot(years, th["Q"] / 1e6, color=SERIES[0])
    ax_Q.set_title("Heat extraction rate")
    ax_Q.set_xlabel("time [years]")
    ax_Q.set_ylabel("Q [MW]")
    ax_Q.set_xlim(0, years[-1])
    ax_Q.set_ylim(bottom=0)
    ax_Q.annotate(f"m_dot = {float(th['m_dot']):.2f} kg/s\nslab half-width = {float(th['slab_half_width']):.0f} m",
                  (0.98, 0.95), xycoords="axes fraction", ha="right", va="top", fontsize=8, color=TEXT_2)

    # Cumulative heat extracted, against the no-cooling limit m_dot c (T_res - T_inj) t
    t = th["t_eval"]
    Q = th["Q"]
    E = cumulative_heat(t, Q)                # GWh
    E_max = Q[0] * t / 3.6e12
    ax_E.plot(years, E_max, color=REFERENCE, linewidth=1, linestyle=(0, (4, 3)))
    ax_E.annotate("no cooling", (years[-1], E_max[-1]), xytext=(-4, -4), textcoords="offset points",
                  ha="right", va="top", fontsize=8, color=TEXT_2)
    ax_E.plot(years, E, color=SERIES[0])
    ax_E.annotate(f"{E[-1]:.1f} GWh ({100 * E[-1] / E_max[-1]:.0f}% of no cooling)", (years[-1], E[-1]),
                  xytext=(-4, 6), textcoords="offset points", ha="right", fontsize=8, color=TEXT)
    ax_E.set_title("Cumulative heat extracted")
    ax_E.set_xlabel("time [years]")
    ax_E.set_ylabel("E [GWh thermal]")
    ax_E.set_xlim(0, years[-1])
    ax_E.set_ylim(0, E_max[-1] * 1.05)


    def at_year(x, year):
        """ Value of x at the output time nearest `year` """
        return x[np.argmin(np.abs(years - year))]


    # Run summary as text; the full inputs need run_info.json, otherwise only what
    # thermal.npz and particles.npz hold
    lines = []
    if info is not None:
        lines += [f"run  {info['jobname']}", "",
                  "Geometry",
                  f"  {info['num_hf']} HFs, spacing {info['hf_spacing']:.0f} m, radius {info['hf_radii']} m, "
                  f"aspect {info['hf_aspect']}",
                  f"  HF x = " + ", ".join(f"{c:.0f}" for c in info['hf_centers']) + " m",
                  f"  laterals at y = ±{info['lateral_y']} m, {info['n_seg']} segments each",
                  f"  segments {info['seg_length']} m x {info['seg_width']} m, overlap {info['seg_overlap']} m",
                  f"  aperture {info['aperture']:.1e} m, {info['num_frac']} fractures",
                  f"  source fracture {info['inj_frac']}, target fracture {info['prod_frac']}", "",
                  "Flow",
                  f"  p_in - p_out = {(info['p_in'] - info['p_out']) / 1e6:.2f} MPa",
                  f"  Q_vol = {info['flow_rate']:.3e} m^3/s"]
        if "conductance_model" in info:
            lines += [f"  conductance {info['conductance_model']}, mp_correction {info['mp_correction']}"]
    else:
        lines += [f"run  {os.path.basename(os.path.abspath(jobdir))}",
                  "  (no run_info.json; rerun driver_v2.py for inputs)", "",
                  "Flow"]
    lines += [f"  m_dot = {float(th['m_dot']):.2f} kg/s", "",
              "Particles",
              f"  {tau.size} exited" + (f" of {info['nparticles']}" if info is not None else ""),
              f"  median travel time {np.median(tau) / 3600:.2f} h",
              f"  median beta {np.median(beta):.2e} s/m", "",
              "Thermal",
              f"  T_res {T_res:.0f} °C, T_inj {T_inj:.0f} °C, slab half-width {float(th['slab_half_width']):.0f} m",
              f"  T_out at 1 / 10 / {years[-1]:.0f} yr: " + " / ".join(
                  f"{at_year(th['T_out'], y):.0f}" for y in (1, 10, years[-1])) + " °C",
              f"  E at {years[-1]:.0f} yr: {E[-1]:.1f} GWh"]
    ax_info.axis("off")
    ax_info.set_title("Summary")
    ax_info.text(0.0, 1.0, "\n".join(lines), transform=ax_info.transAxes, va="top", ha="left",
                 fontsize=8, family="monospace", color=TEXT, linespacing=1.4)

    # Travel time CDF per HF (log time axis: the HFs differ by orders of magnitude)
    for k in hfs:
        t_k = np.sort(tau[hf_of == k]) / 3600.0
        c = SERIES[(k - 1) % len(SERIES)]
        ax_tau.step(t_k, np.arange(1, t_k.size + 1) / t_k.size, where="post", color=c, label=f"HF{k} (n = {t_k.size})")
        ax_tau.annotate(f"HF{k}", (t_k[-1], 1.0), xytext=(4, -6), textcoords="offset points",
                        fontsize=8, color=TEXT_2)
    ax_tau.set_xscale("log")
    ax_tau.set_title("Particle travel time")
    ax_tau.set_xlabel("advective travel time [hours]")
    ax_tau.set_ylabel("cumulative fraction of particles")
    ax_tau.set_ylim(0, 1.02)
    ax_tau.legend(loc="lower left", bbox_to_anchor=(0.1, 0.03))

    # Share of particles through each HF
    x = np.arange(len(hfs))
    bars = ax_hf.bar(x, [share[k] for k in hfs], width=0.6,
                     color=[SERIES[(k - 1) % len(SERIES)] for k in hfs], edgecolor=SURFACE, linewidth=2)
    for b, k in zip(bars, hfs):
        ax_hf.annotate(f"{100 * share[k]:.1f}%", (b.get_x() + b.get_width() / 2, b.get_height()),
                       xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8, color=TEXT)
    ax_hf.set_xticks(x, [f"HF{k}" for k in hfs])
    ax_hf.set_title("Share of particles through each HF")
    ax_hf.set_ylabel("fraction of particles")
    ax_hf.set_ylim(0, 1.1)
    ax_hf.grid(axis="x", visible=False)

    out = os.path.join(jobdir, "results.png")
    fig.savefig(out, dpi=150)
    print(f"--> Figure written to {out}")


# ## Many runs

def compare_production(jobdirs, out="production_comparison.png"):
    """ Compare total heat production across runs

    Left: cumulative heat E(t) for every run. Middle: total heat at the end of
    each run, ranked. Right: total heat against the mass flow rate. Up to eight
    runs get their own colours; with more, every run is drawn in one colour and
    only the highest and lowest totals are highlighted. Writes the figure and a
    CSV of the totals next to it (same name, .csv).
    """
    runs = []
    for d in jobdirs:
        if not os.path.isdir(d):
            continue                       # e.g. lhs_hf_1.log matched by lhs_hf_*
        if not os.path.isfile(os.path.join(d, "thermal.npz")):
            print(f"--> Skipping {d}: no thermal.npz (run unfinished or failed)")
            continue
        th = np.load(os.path.join(d, "thermal.npz"))
        E = cumulative_heat(th["t_eval"], th["Q"])
        runs.append({"name": os.path.basename(os.path.normpath(d)), "years": th["t_eval"] / SECONDS_PER_YEAR,
                     "E": E, "E_end": E[-1], "m_dot": float(th["m_dot"]),
                     "T_out_end": float(th["T_out"][-1]), "Q_end": float(th["Q"][-1]) / 1e6})
    if not runs:
        sys.exit("Error: no runs to compare")
    order = sorted(range(len(runs)), key=lambda i: runs[i]["E_end"], reverse=True)   # largest total first
    few = len(runs) <= len(SERIES)          # few runs: one colour each; many: highlight the extremes
    hi, lo = order[0], order[-1]

    def color(i):
        """ Run i's colour: its own slot when few runs, else blue/orange for the extremes, grey otherwise """
        if few:
            return SERIES[i]
        return SERIES[0] if i == hi else SERIES[1] if i == lo else GRID

    fig, (ax_E, ax_bar, ax_m) = plt.subplots(1, 3, figsize=(16, 5.5), constrained_layout=True,
                                             gridspec_kw={"width_ratios": [1.3, 1, 1]})
    t_end = max(r["years"][-1] for r in runs)

    # cumulative heat; highlighted runs drawn last so they sit on top
    for i in sorted(range(len(runs)), key=lambda i: (not few) and i in (hi, lo)):
        r = runs[i]
        emphasized = few or i in (hi, lo)
        ax_E.plot(r["years"], r["E"], color=color(i), linewidth=2 if emphasized else 1,
                  label=r["name"] if emphasized else None, zorder=3 if emphasized else 2)
        if emphasized and (not few or len(runs) <= 4):   # direct labels up to 4 runs; legend beyond
            ax_E.annotate(f"{r['name']}  {r['E_end']:.1f} GWh", (r["years"][-1], r["E"][-1]),
                          xytext=(4, 0), textcoords="offset points", va="center", fontsize=8, color=TEXT_2)
    ax_E.set_title(f"Cumulative heat extracted ({len(runs)} runs)")
    ax_E.set_xlabel("time [years]")
    ax_E.set_ylabel("E [GWh thermal]")
    ax_E.set_xlim(0, t_end * (1.25 if (not few or len(runs) <= 4) else 1.02))   # room for end labels
    ax_E.set_ylim(bottom=0)
    if few:
        ax_E.legend(loc="upper left")

    # total heat at the end, ranked: one labelled bar per run when they fit,
    # otherwise a rank plot with only the top and bottom three labelled
    if len(runs) <= 30:
        y = np.arange(len(runs))
        ax_bar.barh(y, [runs[i]["E_end"] for i in order], height=0.7,
                    color=[color(i) if (few or i in (hi, lo)) else REFERENCE for i in order],
                    edgecolor=SURFACE, linewidth=1)
        ax_bar.set_yticks(y, [runs[i]["name"] for i in order], fontsize=8 if len(runs) <= 20 else 6)
        ax_bar.invert_yaxis()
        if len(runs) <= 20:
            for k, i in enumerate(order):
                ax_bar.annotate(f"{runs[i]['E_end']:.1f}", (runs[i]["E_end"], k), xytext=(3, 0),
                                textcoords="offset points", va="center", fontsize=8, color=TEXT)
        ax_bar.set_xlabel("E [GWh thermal]")
        ax_bar.grid(axis="y", visible=False)
    else:
        rank = np.arange(1, len(runs) + 1)
        E_sorted = np.array([runs[i]["E_end"] for i in order])
        ax_bar.plot(rank, E_sorted, color=REFERENCE, linewidth=1.5, marker="o", markersize=3)
        for k in list(range(3)) + list(range(len(runs) - 3, len(runs))):
            i = order[k]
            emphasized = i in (hi, lo)
            ax_bar.plot(rank[k], E_sorted[k], marker="o", markersize=7 if emphasized else 5, linestyle="none",
                        color=color(i) if emphasized else TEXT_2, markeredgecolor=SURFACE, markeredgewidth=1.5)
            # bottom three sit close together: stack their labels upward
            dy = 12 * (len(runs) - 1 - k) if k >= len(runs) - 3 else 0
            ax_bar.annotate(f"{runs[i]['name']}  {E_sorted[k]:.1f}", (rank[k], E_sorted[k]),
                            xytext=(6, dy), textcoords="offset points", va="center", fontsize=8, color=TEXT_2)
        ax_bar.set_xlabel("rank (1 = most heat)")
        ax_bar.set_ylabel("E [GWh thermal]")
        ax_bar.set_xlim(0, len(runs) * 1.3)
        ax_bar.set_ylim(bottom=0)
    ax_bar.set_title(f"Total heat at {t_end:.0f} years")

    # total heat against mass flow rate
    for i, r in enumerate(runs):
        emphasized = few or i in (hi, lo)
        ax_m.plot(r["m_dot"], r["E_end"], marker="o", markersize=8 if emphasized else 6, linestyle="none",
                  color=color(i) if emphasized else REFERENCE, markeredgecolor=SURFACE, markeredgewidth=1.5,
                  zorder=3 if emphasized else 2)
        if emphasized and (not few or len(runs) <= 4):
            ax_m.annotate(r["name"], (r["m_dot"], r["E_end"]), xytext=(6, 0), textcoords="offset points",
                          va="center", fontsize=8, color=TEXT_2)
    ax_m.set_title("Total heat against mass flow rate")
    m_dots = np.array([r["m_dot"] for r in runs])
    if m_dots.min() > 0 and m_dots.max() / m_dots.min() > 30:
        ax_m.set_xscale("log")              # flow rates spanning decades
    ax_m.set_xlabel("m_dot [kg/s]")
    ax_m.set_ylabel("E [GWh thermal]")
    ax_m.set_ylim(bottom=0)

    fig.savefig(out, dpi=150)
    plt.close(fig)
    csv = os.path.splitext(out)[0] + ".csv"
    with open(csv, "w") as fp:
        fp.write("run,E_end_GWh,m_dot_kg_s,T_out_end_C,Q_end_MW\n")
        for i in order:
            r = runs[i]
            fp.write(f"{r['name']},{r['E_end']:.6g},{r['m_dot']:.6g},{r['T_out_end']:.6g},{r['Q_end']:.6g}\n")
    print(f"--> Figure written to {out}; totals in {csv}")
    print(f"    highest: {runs[hi]['name']} {runs[hi]['E_end']:.1f} GWh, lowest: {runs[lo]['name']} {runs[lo]['E_end']:.1f} GWh")


# ## Command line

def resolve_runs(names, results_dir=RESULTS_DIR):
    """ Run folders for the given names: folder paths, run names in results_dir,
    or glob patterns tried as given and inside results_dir. No names: every run
    folder in results_dir. Folders are returned in natural order (run_2 before run_10). """
    if not names:
        names = ["*"]
    found = []
    for name in names:
        if os.path.isdir(name):
            found.append(name)
            continue
        matches = [d for d in glob.glob(name) + glob.glob(os.path.join(results_dir, name)) if os.path.isdir(d)]
        if not matches:
            print(f"--> No run folder matches '{name}' (looked here and in {results_dir}/)")
        found += matches
    natural = lambda d: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", os.path.normpath(d))]
    return sorted(dict.fromkeys(os.path.normpath(d) for d in found), key=natural)


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "--compare":
        out = "production_comparison.png"
        if "--out" in args:
            k = args.index("--out")
            out = args[k + 1]
            args = args[:k] + args[k + 2:]
        compare_production(resolve_runs(args[1:]), out)
    elif len(args) == 1:
        runs = resolve_runs(args)
        if len(runs) != 1:
            sys.exit(f"Error: '{args[0]}' should name one run; it matches {len(runs)}")
        plot_run(runs[0])
    else:
        sys.exit("usage: python plot_results.py <run>\n"
                 "       python plot_results.py --compare [<run> ...] [--out production_comparison.png]\n"
                 "  <run>: a folder, a run name in results/ (lhs_hf_1), or a quoted pattern ('lhs_hf_*')")

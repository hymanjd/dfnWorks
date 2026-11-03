"""Figures for graph_egs_prototype.tex, from the run outputs in data/

    data/current  one run: 3 HFs (radius 70 m, aspect 2), the geometry and results slides
    data/sweep    30-run Latin hypercube sweep (HF number 2-20, radius, aspect): master_results.csv

usage: python make_figures.py      (writes figs/)
"""
import os
import re
import sys
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SECONDS_PER_YEAR = 365.25 * 24 * 3600.0
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
TEXT, TEXT_2, GRID, REFERENCE = "#0b0b0b", "#52514e", "#e4e3df", "#8a8984"
plt.rcParams.update({
    "text.usetex": False, "font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": TEXT_2,
    "axes.titlecolor": TEXT, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "xtick.color": TEXT_2, "ytick.color": TEXT_2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
    "legend.frameon": False, "legend.fontsize": 8,
})
os.makedirs("figs", exist_ok=True)


def read_polygons(filename):
    polys = []
    for line in open(filename).read().splitlines()[1:]:
        pts = re.findall(r"\{([^}]*)\}", line)
        if pts:
            polys.append(np.array([[float(v) for v in p.split(",")] for p in pts]))
    return polys


# ---------------------------------------------------------------- geometry
info = json.load(open("data/current/run_info.json"))
polys = read_polygons("data/current/dfn_polygons.dat")
n_seg = info["n_seg"]
wells, hfs = polys[:2 * n_seg], polys[2 * n_seg:]

fig, (ax, ax_x) = plt.subplots(1, 2, figsize=(9.2, 3.6), gridspec_kw={"width_ratios": [3.2, 1]},
                               constrained_layout=True)
for k, P in enumerate(hfs):
    ax.plot([P[0, 0]] * 2, [P[:, 1].min(), P[:, 1].max()], color=SERIES[0], linewidth=3, solid_capstyle="butt")
    ax.annotate(f"HF{k + 1}", (P[0, 0], P[:, 1].max()), xytext=(0, 4), textcoords="offset points",
                ha="center", fontsize=8, color=TEXT)
for i, P in enumerate(wells):
    well_color = SERIES[1] if i < n_seg else SERIES[2]
    if np.ptp(P[:, 1]) > 0.1:      # horizontal segment: a strip in plan view
        ax.fill(P[:, 0], P[:, 1], color=well_color, alpha=0.55, linewidth=0)
    else:                          # vertical segment: a line in plan view
        ax.plot(P[:, 0], P[:, 1], color=well_color, linewidth=1.2)
for frac, name, color in ((info["inj_frac"], "injection", SERIES[1]), (info["prod_frac"], "production", SERIES[2])):
    P = polys[frac - 1]
    ax.plot(P[:, 0].min(), P[:, 1].mean(), marker="o", markersize=8, color=color,
            markeredgecolor="white", markeredgewidth=1.5, zorder=5)
    ax.annotate(f"{name} (fracture {frac})", (P[:, 0].min(), P[:, 1].mean()), xytext=(-4, 7),
                textcoords="offset points", ha="left", va="bottom", fontsize=8, color=TEXT)
ax.set_aspect("equal")
ax.set_xlim(-240, 230)
ax.set_ylim(-160, 165)
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_title("Plan view: laterals, HFs, source and target")

# X cross-section of one well at a segment overlap
s = 0.5 * info["seg_width"]
ax_x.plot([-s, s], [0, 0], color=SERIES[1], linewidth=3)
ax_x.plot([0, 0], [-s, s], color=SERIES[1], linewidth=3, alpha=0.55)
ax_x.annotate("normal (0,0,1)", (s, 0), xytext=(-2, -10), textcoords="offset points", ha="right", fontsize=8, color=TEXT_2)
ax_x.annotate("normal (0,1,0)", (0, s), xytext=(4, -4), textcoords="offset points", ha="left", fontsize=8, color=TEXT_2)
ax_x.set_aspect("equal")
ax_x.set_xlim(-3.2, 3.2)
ax_x.set_ylim(-3.2, 3.2)
ax_x.set_xlabel("y - y_well [m]")
ax_x.set_ylabel("z [m]")
ax_x.set_title("Well cross-section")
fig.savefig("figs/geometry.pdf")
plt.close(fig)


# ---------------------------------------------------------------- sweep
sw = pd.read_csv("data/sweep/master_results.csv")
inputs = [("hydraulic_fractures.number", "HF number"), ("hydraulic_fractures.radius", "HF radius [m]"),
          ("hydraulic_fractures.aspect", "HF aspect")]
fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.0), constrained_layout=True, sharey=True)
for ax, (c, name) in zip(axes, inputs):
    ax.scatter(sw[c], sw["E_total_GWh"], s=16, color=SERIES[0], edgecolor="white", linewidth=0.6, zorder=2)
    if sw[c].nunique() <= 25 and np.allclose(sw[c], sw[c].round()):
        means = sw.groupby(c)["E_total_GWh"].mean()
    else:
        q = pd.qcut(sw[c], q=6, duplicates="drop")
        means = sw.groupby(q, observed=True).agg(x=(c, "mean"), e=("E_total_GWh", "mean")).set_index("x")["e"]
    ax.plot(means.index, means.values, color=SERIES[1], linewidth=2, marker="o", markersize=3, zorder=3)
    rho = stats.spearmanr(sw[c], sw["E_total_GWh"]).statistic
    ax.set_title(f"{name}  (ρ = {rho:+.2f})")
    ax.set_xlabel(name)
axes[0].set_ylabel("total production E [GWh]")
fig.savefig("figs/sweep.pdf")
plt.close(fig)

# regression numbers quoted on the sweep slide
X = np.column_stack([np.ones(len(sw))] + [sw[c].values for c, _ in inputs])
y = sw["E_total_GWh"].values
coef, *_ = np.linalg.lstsq(X, y, rcond=None)
resid = y - X @ coef
dof = len(y) - X.shape[1]
se = np.sqrt(np.diag(resid @ resid / dof * np.linalg.inv(X.T @ X)))
p = 2 * stats.t.sf(np.abs(coef / se), dof)
r2 = 1 - resid @ resid / np.sum((y - y.mean()) ** 2)
print(f"--> sweep: {len(sw)} runs, E {y.min():.1f}-{y.max():.1f} GWh, linear R^2 = {r2:.3f}")
for (c, name), b, pp in zip(inputs, coef[1:], p[1:]):
    print(f"    {name:14s} coef {b:+.3f} GWh per unit, p = {pp:.2g}")

print("--> figs/geometry.pdf, figs/sweep.pdf")
print("--> results: python <EGS_scripts>/graph_based/plot_results.py data/current; cp data/current/results.png figs/")

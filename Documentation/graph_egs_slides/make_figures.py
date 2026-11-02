"""Figures for graph_egs_prototype.tex, from the run outputs in data/

    data/current  driver_v1.py as committed (hf_aspect = 2), after the dfnWorks fixes
    data/before   hf_aspect = 1, before the fixes (Karra-form edge flows, intersection-wide
                  edges, heat exchange in the well segments)
    data/after    hf_aspect = 1, after the fixes

usage: python make_figures.py      (writes figs/)
"""
import os
import re
import sys
import json
import numpy as np
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


# ---------------------------------------------------------------- before / after
def curves(run):
    th = np.load(f"data/{run}/thermal.npz")
    t = th["t_eval"]
    Q = th["Q"]
    E = np.concatenate([[0.0], np.cumsum(0.5 * (Q[1:] + Q[:-1]) * np.diff(t))]) / 3.6e12
    return t / SECONDS_PER_YEAR, th["T_out"], E, float(th["m_dot"])


fig, (ax_T, ax_E) = plt.subplots(1, 2, figsize=(9.2, 3.4), constrained_layout=True)
for run, label, c in (("before", "before fixes", SERIES[1]), ("after", "after fixes", SERIES[0])):
    yr, T, E, m = curves(run)
    ax_T.plot(yr, T, color=c, label=f"{label} (m_dot {m:.2f} kg/s)")
    ax_E.plot(yr, E, color=c, label=label)
    ax_E.annotate(f"{E[-1]:.0f} GWh", (yr[-1], E[-1]), xytext=(-4, 5), textcoords="offset points",
                  ha="right", fontsize=8, color=TEXT)
ax_T.set_title("Produced temperature")
ax_T.set_xlabel("time [years]")
ax_T.set_ylabel("T_out [°C]")
ax_T.set_xlim(0, 30)
ax_T.legend(loc="upper right")
ax_E.set_title("Cumulative heat extracted")
ax_E.set_xlabel("time [years]")
ax_E.set_ylabel("E [GWh thermal]")
ax_E.set_xlim(0, 30)
ax_E.set_ylim(bottom=0)
ax_E.legend(loc="upper left")
fig.savefig("figs/before_after.pdf")
plt.close(fig)

print("--> figs/geometry.pdf, figs/before_after.pdf")
print("--> results: python <EGS_scripts>/graph_based/plot_results.py data/current; cp data/current/results.png figs/")

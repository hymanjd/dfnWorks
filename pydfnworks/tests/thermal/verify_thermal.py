"""
verify_thermal.py
-----------------
Verification of pydfnworks.thermal (DFN-based heat extraction).

Tests
-----
  1. Slab kernel limits: beta_G -> inf is the semi-infinite erfc solution
     (Lauwerier 1955); beta_G -> 0 is a step at t_D = beta_G with unit
     heat balance
  2. Slab kernel at finite beta_G against an independent finite-difference
     solution of the slab-conduction / fracture-advection equations
  3. Single path: thermal_breakthrough() with infinite rock equals the
     Lauwerier closed form erfc(phi sqrt(D) beta / (2 sqrt(t - tau))), and
     with finite slabs equals the slab kernel
  4. Consistency with the graph-transport TDRW sampler: the Monte Carlo
     CDF of the infinite-matrix return time (tdrw.infinite) over one
     edge matches thermal_breakthrough() for that edge
  5. Superposition: two groups of paths give the weighted average, and
     stuck/empty inputs are rejected
  6. Dataset helpers: hydraulic fracture placement, HDF5 round trip,
     manifest

Run with:  python verify_thermal.py   (exit code 0 if all pass)
"""

import os
import sys
import tempfile
import numpy as np
from scipy.special import erfc, erfcinv
from scipy.linalg import solve_banded

from pydfnworks.thermal import (slab_drawdown, thermal_breakthrough, thermal_to_tdrw_parameters,
                                hydraulic_fracture_positions, write_sample, append_manifest,
                                GRANITE_PROPS, WATER_PROPS, SECONDS_PER_YEAR)
from pydfnworks.dfnGraph.transport.tdrw.infinite import t_diff_unlimited

_trapezoid = getattr(np, "trapezoid", None) or np.trapz
P = thermal_to_tdrw_parameters(GRANITE_PROPS, WATER_PROPS)
PHI, D = P["matrix_porosity"], P["matrix_diffusivity"]


def test_slab_kernel_limits():
    t_D = np.logspace(-3, 3, 60)
    err_inf = np.max(np.abs(slab_drawdown(t_D, np.inf) - erfc(0.5 / np.sqrt(t_D))))
    err_big = np.max(np.abs(slab_drawdown(t_D, 1e5) - erfc(0.5 / np.sqrt(t_D))))
    # small beta_G: step at t_D = beta_G; integral of (1 - T_wD) over t_D equals beta_G
    bg = 0.3
    t = np.linspace(0, 6 * bg, 3001)
    heat = _trapezoid(1 - slab_drawdown(t, bg), t) / bg
    print(f"  |T_wD - erfc|: beta_G = inf {err_inf:.1e}, 1e5 {err_big:.1e}; "
          f"beta_G = {bg}: heat / piston = {heat:.4f}")
    assert err_inf < 1e-12 and err_big < 1e-4 and abs(heat - 1) < 1e-2


def fd_outlet_drawdown(beta_G, t_end, nz, nx, nt):
    """ Finite differences for: slab d theta/dt = d2 theta/dx2 on 0 < x < beta_G,
    insulated at x = beta_G; fracture d theta_w/dz = d theta/dx|0+, theta_w(0) = 1 """
    dx, dz = beta_G / nx, 1.0 / nz
    theta = np.zeros((nz + 1, nx + 1))
    c = dz / (2 * dx)
    T, times, out = 0.0, [0.0], [0.0]
    for dt in list(np.geomspace(t_end / nt * 1e-3, t_end / nt, 100)) + [t_end / nt] * nt:
        if T >= t_end - 1e-12:
            break
        r = dt / (2 * dx**2)
        ab = np.zeros((4, nx + 1))
        ab[2, :], ab[1, 1:], ab[3, :-1] = 1 + 2 * r, -r, -r
        ab[3, nx - 1] = -2 * r
        ab[2, 0], ab[1, 1], ab[0, 2] = 1 + 3 * c, -4 * c, c
        ab_in = ab.copy()
        ab_in[2, 0], ab_in[1, 1], ab_in[0, 2] = 1.0, 0.0, 0.0
        theta_w = 1.0
        for j in range(nz + 1):
            old = theta[j]
            rhs = old.copy()
            rhs[1:-1] = old[1:-1] + r * (old[:-2] - 2 * old[1:-1] + old[2:])
            rhs[-1] = old[-1] + 2 * r * (old[-2] - old[-1])
            rhs[0] = 1.0 if j == 0 else theta_w
            theta[j] = solve_banded((1, 2), ab_in if j == 0 else ab, rhs)
            theta_w = theta[j, 0]
        T += dt
        times.append(T)
        out.append(theta_w)
    return np.array(times), np.array(out)


def test_slab_kernel_finite_difference():
    errs = []
    for nz, nx, nt in ((50, 100, 1500), (100, 200, 3000)):
        t, fd = fd_outlet_drawdown(1.0, 6.0, nz, nx, nt)
        m = t > 0.05
        errs.append(np.max(np.abs(fd[m] - slab_drawdown(t[m], 1.0))))
    print(f"  beta_G = 1: max |FD - Stehfest| coarse {errs[0]:.1e}, fine {errs[1]:.1e}")
    assert errs[1] < errs[0] and errs[1] < 1e-2


def test_single_path():
    tau, b = 2.0e4, 3.0e-4            # advective time [s], aperture [m]
    beta = 2 * tau / b                # exchange integral of one edge
    t = tau + np.logspace(5, 10, 80)
    F = thermal_breakthrough([tau], [beta], t)
    F_ref = erfc(PHI * np.sqrt(D) * beta / (2 * np.sqrt(t - tau)))
    err_inf = np.max(np.abs(F - F_ref))
    B = 5.0
    t_s = (PHI * np.sqrt(D) * beta)**2
    F_B = thermal_breakthrough([tau], [beta], t, slab_half_width=B)
    err_B = np.max(np.abs(F_B - slab_drawdown((t - tau) / t_s, B / (PHI * D * beta))))
    print(f"  infinite rock: max |F - Lauwerier| = {err_inf:.1e}; slab B = {B} m: "
          f"beta_G = {B/(PHI*D*beta):.2f}, max |F - kernel| = {err_B:.1e}; F = 0 before tau: {np.all(thermal_breakthrough([tau],[beta],[0, tau]) == 0)}")
    assert err_inf < 1e-12 and err_B < 1e-12
    assert np.all(F_B >= F - 0.02) and F_B[-1] > F[-1]   # finite slabs cool sooner (ringing tolerance near the step)


def test_matches_tdrw_sampler():
    tau, b = 2.0e4, 3.0e-4
    a = PHI * np.sqrt(D) / b           # tdrw.infinite parameterization
    rng = np.random.default_rng(0)
    samples = np.array([t_diff_unlimited(a, tau, xi) for xi in rng.uniform(size=20000)])
    t = tau + np.logspace(5, 10, 40)
    F_mc = np.array([(samples <= ti - tau).mean() for ti in t])
    F = thermal_breakthrough([tau], [2 * tau / b], t)
    err = np.max(np.abs(F - F_mc))
    print(f"  20000 TDRW infinite-matrix samples: max |F - F_MC| = {err:.3f}")
    assert err < 0.02


def test_superposition_and_errors():
    t = np.logspace(5, 10, 30)
    tau = np.array([1e4, 5e4])
    beta = np.array([1e8, 3e8])
    F_each = [thermal_breakthrough([tau[i]], [beta[i]], t) for i in range(2)]
    F_w = thermal_breakthrough(tau, beta, t, weights=[1.0, 3.0])
    assert np.max(np.abs(F_w - (F_each[0] + 3 * F_each[1]) / 4)) < 1e-12
    for bad in ([], ):
        try:
            thermal_breakthrough(bad, bad, t)
            raise AssertionError("empty particle set accepted")
        except ValueError:
            pass
    print("  weighted average of two paths reproduced; empty particle set rejected")


def test_dataset_helpers():
    pos = hydraulic_fracture_positions(100.0, 10.0)
    assert pos.size == 10 and abs(pos[0] + 45.0) < 1e-12 and abs(np.diff(pos).mean() - 10) < 1e-12
    try:
        hydraulic_fracture_positions(10.0, 20.0)
        raise AssertionError("spacing larger than the domain accepted")
    except ValueError:
        pass
    t = np.linspace(0, 30 * SECONDS_PER_YEAR, 50)
    result = {"t_eval": t, "F": np.linspace(0, 1, 50), "T_out": np.linspace(226, 40, 50),
              "Q": np.linspace(1e7, 0, 50), "m_dot": 20.0, "flow_rate": 0.02,
              "advect_time": np.array([1e4, 2e4]), "beta": np.array([1e8, 2e8]), "n_exit": 2, "n_stuck": 0}
    inputs = {"spacing": 10.0, "seed": 1, "rock_props": GRANITE_PROPS, "label": "test"}
    with tempfile.TemporaryDirectory() as d:
        write_sample(os.path.join(d, "s.h5"), inputs, result)
        import h5py, json
        with h5py.File(os.path.join(d, "s.h5")) as f:
            assert f["inputs"].attrs["spacing"] == 10.0
            assert json.loads(f["inputs"].attrs["rock_props"])["k"] == GRANITE_PROPS["k"]
            assert np.allclose(f["curves/Q"][...], result["Q"]) and f.attrs["n_exit"] == 2
            assert f["particles/beta"].shape == (2,)
        m = os.path.join(d, "manifest.csv")
        append_manifest(m, {"sample": "s.h5", "spacing": 10.0})
        append_manifest(m, {"sample": "t.h5", "spacing": 20.0})
        assert open(m).read().count("\n") == 3
    print("  10 fractures at 10 m in 100 m, first at -45 m; HDF5 round trip and manifest OK")


TESTS = [
    ("Slab kernel limits (erfc, piston)", test_slab_kernel_limits),
    ("Slab kernel vs finite differences", test_slab_kernel_finite_difference),
    ("Single path vs Lauwerier and slab kernel", test_single_path),
    ("Consistency with the TDRW infinite-matrix sampler", test_matches_tdrw_sampler),
    ("Superposition and input errors", test_superposition_and_errors),
    ("Dataset helpers", test_dataset_helpers),
]

if __name__ == "__main__":
    n_fail = 0
    for i, (label, test) in enumerate(TESTS, 1):
        print(f"\nTest {i}: {label}")
        try:
            test()
            print("  PASS")
        except Exception as err:
            n_fail += 1
            print(f"  FAIL: {type(err).__name__}: {err}")
    print(f"\n{len(TESTS) - n_fail} of {len(TESTS)} tests passed")
    sys.exit(1 if n_fail else 0)

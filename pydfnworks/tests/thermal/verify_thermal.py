"""
verify_thermal.py
-----------------
Verification of pydfnworks.thermal (Gringarten MPF surrogate, DFN backend
interface, spacing optimizer).

Tests
-----
  1. Wide-spacing limit equals the semi-infinite single-fracture erfc
     solution (Lauwerier; Bodvarsson and Tsang 1982)
  2. Close spacing: never more heat than wide spacing, and much earlier
     deep drawdown
  3. Energy balance at close spacing: cumulative heat = heat in the slabs
  4. Stehfest inversion matches the GEOPHIRES-X mpmath route
  5. Finite spacing: a finite-difference solution of the governing
     equations converges to the inverted solution at beta = 1
  6. Physical bounds and runtime
  7. Interface: both backends reject the same bad inputs; the DFN stub
     raises NotImplementedError; thermal-to-TDRW mapping
  8. Optimizer: objectives, grid, stage-count coupling, backend swap

Run with:  python verify_thermal.py   (exit code 0 if all pass)
The test_* functions also run under pytest.
"""

import sys
import time
import functools
import numpy as np
from scipy.special import erfc
from scipy.linalg import solve_banded

from pydfnworks.thermal import (mpf_gringarten as mpf, dfn_heat, optimize_spacing,
                                cumulative_heat, time_to_drawdown,
                                GRANITE_PROPS, WATER_PROPS, SECONDS_PER_YEAR)

ROCK, FLUID = GRANITE_PROPS, WATER_PROPS
M_DOT, T_RES, T_INJ, N, H = 60.0, 226.0, 40.0, 20, 60.0
T_EVAL = np.linspace(0.0, 30 * SECONDS_PER_YEAR, 361)
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


def Q_of(spacing, t_eval=T_EVAL, n=N, backend=mpf.heat_extraction):
    return backend(spacing, M_DOT, T_RES, T_INJ, ROCK, FLUID, n, H, t_eval)


def test_wide_spacing_limit():
    p = mpf.dimensionless_parameters(1e6, M_DOT, ROCK, FLUID, N, H)
    alpha = ROCK["k"] / (ROCK["rho"] * ROCK["c"])
    q_w = M_DOT / (FLUID["rho"] * N * H)
    T_D = np.zeros_like(T_EVAL)
    T_D[1:] = erfc(ROCK["k"] * H / (FLUID["rho"] * FLUID["c"] * q_w
                                    * np.sqrt(alpha * T_EVAL[1:])))
    Q_ref = M_DOT * FLUID["c"] * (T_RES - T_INJ) * (1 - T_D)
    err = np.max(np.abs(Q_of(1e6) - Q_ref)) / Q_ref[0]
    t_D = np.logspace(-3, 3, 50)
    err_D = np.max(np.abs(mpf.dimensionless_drawdown(t_D, p["beta"])
                          - mpf.single_fracture_drawdown(t_D)))
    print(f"  beta = {p['beta']:.1e}: |Q - Q_erfc|/Q0 {err:.1e}, |T_wD - erfc| {err_D:.1e}")
    assert err < 1e-5 and err_D < 1e-5


def test_close_vs_wide():
    Q_close, Q_wide = Q_of(5.0), Q_of(100.0)
    t = np.linspace(0, 2 * SECONDS_PER_YEAR, 4001)
    t75 = [time_to_drawdown(t, Q_of(s, t), 0.25) for s in (5.0, 100.0)]
    print(f"  max(Q_close - Q_wide)/Q0 {np.max(Q_close - Q_wide)/Q_close[0]:.1e}; "
          f"E ratio {_trapezoid(Q_close, T_EVAL)/_trapezoid(Q_wide, T_EVAL):.3f}; "
          f"t_75% {t75[0]/86400:.0f} d vs {t75[1]/86400:.0f} d")
    assert np.max(Q_close - Q_wide) < 1e-3 * Q_close[0]
    assert _trapezoid(Q_close, T_EVAL) < 0.5 * _trapezoid(Q_wide, T_EVAL)
    assert t75[0] < 0.8 * t75[1]


def test_energy_balance():
    for spacing in (2.0, 5.0):
        p = mpf.dimensionless_parameters(spacing, M_DOT, ROCK, FLUID, N, H)
        t = np.linspace(0, 6 * p["t_piston"], 2001)
        ratio = _trapezoid(Q_of(spacing, t), t) / (p["heat_in_place"] * (T_RES - T_INJ))
        print(f"  spacing {spacing} m, beta {p['beta']:.2f}: E / E_rock = {ratio:.4f}")
        assert abs(ratio - 1) < 1e-2


def test_matches_geophires_mpmath():
    from mpmath import exp, invertlaplace, sqrt, tanh, workdps
    t_D = np.logspace(-1, 3, 30)
    for beta in (0.5, 1.5, 6.0):
        fp = lambda s: (1 / s) * exp(-sqrt(s) * tanh(beta * sqrt(s)))
        with workdps(15):  # GEOPHIRES-X default precision
            ref = np.array([float(invertlaplace(fp, t, method="stehfest")) for t in t_D])
        err = np.max(np.abs(mpf.dimensionless_drawdown(t_D, beta) - ref))
        print(f"  beta {beta}: max |stehfest - mpmath| = {err:.1e}")
        assert err < 1e-3


def fd_outlet_drawdown(beta, t_end, nz, nx, nt):
    """ Finite-difference solution of the dimensionless MPF equations:
    slab  d theta/dt = d^2 theta/dx^2 on 0 < x < beta, insulated at x = beta;
    fracture d theta_w/dz = d theta/dx |_0+, theta_w(0) = 1.
    Crank-Nicolson in t per slab column, implicit wall coupling along z. """
    dx, dz = beta / nx, 1.0 / nz
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


def test_finite_difference():
    errs = []
    for nz, nx, nt in ((50, 100, 1500), (100, 200, 3000)):
        t, fd = fd_outlet_drawdown(1.0, 6.0, nz, nx, nt)
        m = t > 0.05
        errs.append(np.max(np.abs(fd[m] - mpf.dimensionless_drawdown(t[m], 1.0))))
    print(f"  beta = 1: max |FD - Stehfest| coarse {errs[0]:.1e}, fine {errs[1]:.1e}")
    assert errs[1] < errs[0] < 1e-1 and errs[1] < 1e-2


def test_bounds_and_runtime():
    Q0 = M_DOT * FLUID["c"] * (T_RES - T_INJ)
    for spacing in (1.0, 20.0, 1e4):
        Q = Q_of(spacing)
        assert abs(Q[0] - Q0) < 1e-9 * Q0 and np.all(Q >= 0) and np.all(Q <= Q0 * (1 + 1e-12))
    tic = time.perf_counter()
    for _ in range(200):
        Q_of(20.0)
    per_call = (time.perf_counter() - tic) / 200
    print(f"  Q(0) = {Q0/1e6:.2f} MW; bounds hold; {per_call*1e6:.0f} us per call")
    assert per_call < 0.1


def test_interface():
    base = [20.0, M_DOT, T_RES, T_INJ, ROCK, FLUID, N, H, T_EVAL]
    bad = [(0, -1.0), (1, 0.0), (3, T_RES), (6, 0), (4, {"rho": 1, "c": 1}),
           (8, np.array([-1.0, 1.0]))]
    for backend in (mpf.heat_extraction, dfn_heat.heat_extraction):
        for i, value in bad:
            args = list(base)
            args[i] = value
            try:
                backend(*args)
                raise AssertionError(f"{backend.__module__} accepted arg {i} = {value}")
            except ValueError:
                pass
    try:
        dfn_heat.heat_extraction(*base)
        raise AssertionError("stub did not raise")
    except NotImplementedError as err:
        assert "run_graph_transport" in str(err)
    tdrw = dfn_heat.thermal_to_tdrw_parameters(ROCK, FLUID)
    assert abs(tdrw["matrix_diffusivity"] - ROCK["k"] / (ROCK["rho"] * ROCK["c"])) < 1e-15
    print(f"  6 bad inputs rejected by both backends; stub raises; "
          f"TDRW porosity-like {tdrw['matrix_porosity']:.3f}, D {tdrw['matrix_diffusivity']:.2e}")


def test_optimizer():
    t, Q = np.array([0.0, 1.0, 2.0, 3.0]), np.array([10.0, 10.0, 5.0, 0.0])
    assert abs(cumulative_heat(t, Q) - 20.0) < 1e-12
    assert abs(time_to_drawdown(t, Q, 0.9) - 1.2) < 1e-12
    assert time_to_drawdown(t, np.full(4, 10.0)) == 3.0
    common = (M_DOT, T_RES, T_INJ, ROCK, FLUID)
    # fixed count: wider spacing extracts more over 30 yr
    r = optimize_spacing(np.geomspace(2, 200, 15), *common, N, H, T_EVAL)
    assert r["spacings"].size == 15 and r["spacing"] == 200.0
    assert np.all(np.diff(r["objectives"]) >= -1e-6 * r["objective"])
    # stages filling a lateral: closest spacing delays drawdown most
    L = 300.0
    r = optimize_spacing(L / np.arange(2, 31), *common, lambda s: int(round(L / s)), H, T_EVAL,
                         objective=functools.partial(time_to_drawdown, fraction=0.9))
    assert r["spacing"] == L / 30 and r["n_fractures"] == 30
    try:
        optimize_spacing([5.0, 50.0], *common, N, H, T_EVAL, heat_extraction=dfn_heat.heat_extraction)
        raise AssertionError("backend swap did not raise")
    except NotImplementedError:
        pass
    print("  objectives OK; fixed-count best at upper bound; fixed-lateral best at closest; "
          "backend swap raises NotImplementedError")


TESTS = [
    ("Wide-spacing erfc limit (Bodvarsson-Tsang)", test_wide_spacing_limit),
    ("Close vs wide spacing", test_close_vs_wide),
    ("Energy balance at close spacing", test_energy_balance),
    ("Stehfest vs GEOPHIRES-X mpmath", test_matches_geophires_mpmath),
    ("Finite-difference cross-check at beta = 1", test_finite_difference),
    ("Physical bounds and runtime", test_bounds_and_runtime),
    ("Interface, DFN stub, TDRW mapping", test_interface),
    ("Spacing optimizer", test_optimizer),
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

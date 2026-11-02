# pydfnworks.thermal: status and pickup notes

*Branch `feature/thermal-mpf-surrogate` on hymanjd/dfnWorks, 2 November 2026*

## What is on the branch

`pydfnworks/pydfnworks/thermal/` is the DFN-based heat-extraction package
that produces per-sample `T_out(t)` and `Q(t)` for downstream (ML) use. The
earlier analytic Gringarten surrogate, optimizer, and slides were removed
on 2 November 2026 at the project lead's direction; the DFN model is the
product.

- `dfn_heat.py`: `heat_extraction(dfn, G, T_res, T_inj, t_eval, ...,
  slab_half_width=D/2)` runs purely advective graph particle tracking and
  applies the finite-slab fracture-rock exchange kernel analytically to
  each particle from its advective time `tau` and exchange integral
  `beta = sum 2 tau_i / b_i` (which `Particle` accumulates without TDRW).
  Laplace kernel `exp[-phi sqrt(D s) beta tanh(B sqrt(s/D))]`, inverted
  with the shared Stehfest routine (18 terms); `thermal_breakthrough()`
  sums the particle curves with flux weights; `slab_drawdown()` is the
  kernel; `thermal_to_tdrw_parameters()` is the heat-as-solute mapping.
- `dataset.py`: `add_hydraulic_fractures()` (parallel user rectangles at
  spacing D), `graph_flow_for_rate()` (pressure drop that delivers a
  target mass flow rate), `write_sample()` (HDF5: `inputs/`, `curves/`,
  `particles/`), `append_manifest()`.
- `pydfnworks/tests/thermal/verify_thermal.py`: 6 tests, all passing:
  erfc and piston limits of the kernel, finite-difference cross-check,
  single path vs the Lauwerier closed form and vs the kernel, Monte Carlo
  consistency with the graph-transport infinite-matrix TDRW sampler,
  superposition, HDF5 round trip.
- `examples/thermal_dfn_dataset/driver.py`: 400 x 200 x 200 m block,
  hydraulic fractures normal to y, two natural TPL sets, flow left to
  right at a target 10 kg/s, spacings 10 to 80 m x 2 seeds, one HDF5 per
  sample plus `manifest.csv` and a figure.
- Slides: `Documentation/thermal_dfn_dataset_method.pptx`, built by
  `Documentation/thermal_slides/build_deck.py` (equations typeset by
  `render_eqs.py` with MacTeX `latex` + `dvipng`; needs python-pptx).

## Open issues found on the last run (2 November 2026)

The example runs end to end for 6 of 8 samples but the produced
temperature stays at `T_res` for 30 years in every sample, and one sample
fails in dfnGen. Both need fixing before the dataset is meaningful.

1. **Exchange integral is ~50x too large on the intersection graph.** For
   10 m spacing, 20 hydraulic fractures and 0.5 kg/s each, the parallel-
   plate estimate is `beta = 2 L H / Q_f = 3.2e8 s/m` (L = 400 m,
   H = 200 m), giving a slab piston time of about 26 years. The stored
   particles have median `beta = 1.6e10 s/m` (1st percentile 9.4e9),
   median `t_s` of 3e6 years and median `beta_G` of 4e-4, so no cooling
   within 30 years. Cause: the pipe-network edges carry the flow through
   an area `b x w_edge` with `w_edge` ~ 40 m (median edge area 1.4e-2 m^2
   at b = 3.46e-4 m) rather than the full 200 m fracture height, so edge
   velocities are low and `tau_i / b_i` is high; detours through natural
   fractures with `b = 3.5e-6 m` add more. Options: compute `beta` at the
   fracture level from the per-fracture flow partition
   (`beta_f = 2 A_f / Q_f` along the fracture sequence the particle
   records in `frac_seq`), or correct the edge widths to the fracture
   extent; then re-verify against the single-fracture limit using a DFN
   of one user rectangle, which should reproduce the Lauwerier curve.
2. **Sparse connectivity at wide spacing.** `spacing_80_seed_1` (two
   planes 80 m apart) kept only 3 of 133 fractures after dfnGen's
   largest-cluster removal and graph flow then found no inlet nodes
   ("Error. There are no nodes in the inlet", which exits). The driver
   should either guarantee connectivity (higher natural p32 or longer
   natural fractures), request fewer samples at wide spacing, or catch
   the failure and skip the sample; `keepOnlyLargestCluster=False` is
   not safe for graph flow without checking that every node reaches the
   boundaries.
3. **Stehfest at small beta_G** is adequate (step at `t_D = beta_G`
   within 0.1 of a 30-digit Talbot reference for beta_G = 0.02 to 0.1;
   the reference itself diverges below 0.02). Not the cause of item 1.
4. Minor: `numpy` floats passed to `add_user_fract` are written verbatim
   into dfnGen's input (`np.float64(...)`), which dfnGen parses as zeros;
   `add_hydraulic_fractures` now converts to plain floats.

## Picking up on another machine

    git clone git@github.com:hymanjd/dfnWorks.git && cd dfnWorks
    git checkout feature/thermal-mpf-surrogate
    cd pydfnworks && pip install -e . && python tests/thermal/verify_thermal.py
    cd ../examples/thermal_dfn_dataset && python driver.py

dfnGen must be built and named in `~/.dfnworksrc` (`dfnworks_PATH`);
LaGriT, PFLOTRAN and FEHM are not needed for this example. The DFN
outputs of the last run are not committed (they regenerate in seconds).

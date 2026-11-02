# DFN heat-extraction dataset (pydfnworks.thermal)

Builds a set of DFN heat-extraction samples for downstream modelling.
Each sample places parallel hydraulic fractures at spacing D (normal y,
spanning the 400 x 200 x 200 m block in x and z) in a natural fracture
network of two truncated power-law sets, solves graph flow from the left
to the right face, tracks particles on the flow graph, and converts each
particle's advective time and fracture-rock exchange integral into the
produced temperature T_out(t) and heat rate Q(t) with the finite-slab
kernel of half-width D/2 (`pydfnworks.thermal.heat_extraction`). No
meshing, PFLOTRAN or FEHM.

`driver.py` loops over spacings 10, 20, 40, 80 m and two seeds and writes
to `output/dataset/`:

- `spacing_<D>_seed_<s>.h5`: group `inputs` (spacing, seed, permeabilities,
  pressures, temperatures, properties, fracture count, p32, dQ as
  attributes), group `curves` (`t_eval`, `F`, `T_out`, `Q`), group
  `particles` (`advect_time`, `beta`), and file attributes `m_dot`,
  `flow_rate`, `n_exit`, `n_stuck`
- `manifest.csv`: one row per sample with the inputs and headline outputs
- `dataset_curves.png`: T_out(t) and Q(t) for every sample

Reservoir numbers are representative of Utah FORGE (226 C granitoid at
about 2.3 km). Every sample is run at the same target mass flow rate
(10 kg/s); `graph_flow_for_rate` solves for the pressure drop that delivers
it, which the sample records, so spacing and the network are the only
things that differ between samples.

## Running on another machine

    git fetch origin && git checkout feature/thermal-mpf-surrogate
    cd pydfnworks && pip install -e . && python tests/thermal/verify_thermal.py
    cd ../examples/thermal_dfn_dataset && python driver.py

Needs dfnGen only (set in `~/.dfnworksrc`). Eight samples take a few
minutes; add spacings and seeds to the lists at the top of `driver.py`.

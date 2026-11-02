# Pickup notes: pydfnworks.thermal (DFN heat-extraction dataset)

Branch `feature/thermal-mpf-surrogate` on the fork `hymanjd/dfnWorks`.
Last push 2 November 2026 (`f8dfbc67`). Nothing is on `lanl/dfnWorks`.

## Where things are

| What | Where |
|---|---|
| Package | `pydfnworks/pydfnworks/thermal/` (`dfn_heat.py`, `dataset.py`) |
| Tests | `pydfnworks/tests/thermal/verify_thermal.py` (6 tests, all pass) |
| Example | `examples/thermal_dfn_dataset/driver.py` + README |
| Status, open defects, numbers | `Documentation/thermal_status.md` |
| Slides | `Documentation/thermal_dfn_dataset_method.pptx`; sources in `Documentation/thermal_slides/` |
| Sphinx page | `Documentation/sphinx-docs/source/pydfnThermal.rst` |

## Set up on a new machine

    git clone git@github.com:hymanjd/dfnWorks.git && cd dfnWorks
    git checkout feature/thermal-mpf-surrogate
    cd pydfnworks && pip install -e . && python tests/thermal/verify_thermal.py

dfnGen must be built (`DFNGen/DFNGen`) and `~/.dfnworksrc` must name the
checkout in `dfnworks_PATH`. LaGriT, PFLOTRAN and FEHM are not needed for
the thermal example. `pip install -e .` must point at this checkout, not
another dfnWorks clone, or `import pydfnworks.thermal` fails.

Run the example (not required to pick up the work; it is the thing to fix):

    cd examples/thermal_dfn_dataset && python driver.py

## What to do next

Read `Documentation/thermal_status.md` first. Two defects from the last
run, in priority order:

1. The particle exchange integral `beta = sum 2 tau_i / b_i` on the
   intersection graph is about 50x the parallel-plate value (median
   1.6e10 s/m vs 3.2e8 expected at 10 m spacing), so no sample cools
   within 30 years. Pipe-network edges carry flow through `b x ~40 m`,
   not the fracture height. Candidate fix: fracture-level `beta` from the
   per-fracture flow partition along the particle's `frac_seq`, or
   corrected edge widths; then re-verify on a one-user-rectangle DFN
   against the Lauwerier curve (the test file has the closed form).
2. At 80 m spacing dfnGen's largest-cluster removal left 3 of 133
   fractures and graph flow exited with "no nodes in the inlet". Make the
   driver guarantee connectivity or skip the sample.

## Carrying the Claude session over

The code is fully on the fork; only the conversation is machine-local.
Options, from the Claude desktop app:

- **Open in > Cloud** on the session: carries a summary of the
  conversation and clones this branch. Uncommitted files do not travel.
- **Remote Control** (`/remote-control`): drive the session from another
  device while the original Mac stays on.
- **Copy the transcript**: on the original Mac the session is
  `~/.claude/projects/-Users-hymanjd-src-dfnWorks-jdhdev/d06dd10f-6bd4-44ea-a2d5-720a8dd77b85.jsonl`,
  and the memory notes are in the `memory/` folder of the project
  directory under `~/.claude/projects/` named after the original scratch
  workspace. Copy both to the same paths on the other Mac and run
  `claude --resume d06dd10f-6bd4-44ea-a2d5-720a8dd77b85`.

A fresh session with none of that still has what it needs: this file,
`Documentation/thermal_status.md`, and the module docstrings describe the
model, the state, and the next step.

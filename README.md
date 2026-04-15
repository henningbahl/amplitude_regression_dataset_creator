# amplitude_regression_dataset_creator

Generate parton-level events for an arbitrary MadGraph5_aMC@NLO process and
compute the tree-level matrix element |M|² for every event using the
`allmatrix2py` Python interface (see
[FAQ-General-4](https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/FAQ-General-4)).

## Usage

```
./generate_amp_dataset.py --process "<MG process>" [options]
```

On the first invocation, point the script at your MadGraph binary with
`--mg-path`. The path is cached in `.mg_path` next to the script and reused
on subsequent calls.

```
./generate_amp_dataset.py --mg-path /path/to/MG/bin/mg5_aMC --process "q q > z g"
```

Arguments:

- `--process STR` *(required)* — the MadGraph generate string, e.g.
  `"q q > z g"`, `"q q > z g g"`, `"g g > t t~"`. The leading
  `generate ` keyword is added automatically. The shortcut
  `q = u d s c u~ d~ s~ c~` is always defined.
- `--n-events N` *(optional, default 10000)* — number of unweighted events
  to generate.
- `--mg-path PATH` *(optional)* — absolute path to `mg5_aMC`. Saved to
  `.mg_path` and reused on subsequent runs. Required on the first call
  if no `MG/bin/mg5_aMC` exists next to the script.
- `--m-inv lo,hi` — cut on the invariant mass of all final-state
  particles (in GeV). Either bound may be omitted, e.g. `--m-inv 400,`.
- `--costheta lo,hi` — cut on a Lorentz-invariant scattering-angle
  proxy `(u - t)/(u + t)`. Only valid for 2 → 2 processes.
  **Warning:** this expression equals the true cos θ* only if at least
  one of the two final-state particles is massless. For two massive
  final states it is only a proxy.
- `--workdir PATH` — output directory (default `./output/proc_<slug>`,
  where `<slug>` is derived from the final-state of the process string).
- `--plot` — produce a hexbin plot of `(m_inv, cos θ)` colored by mean
  |M|² (2 → 2 only). Same massless-final-state caveat as `--costheta`.

The `m_inv` and `costheta` cuts are applied **at generation time** via
a patched `dummy_cuts` function in MadGraph (see
[FAQ-General-1](https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/FAQ-General-1)),
so the unweighted events are sampled directly inside the requested
phase-space window. Default jet cuts are: `etaj = -1` (no eta cut),
`drjj = 0.3`.

## Output

In the working directory (`output/proc_<slug>/` by default):

- `proc_<slug>.npy` — array of shape `(N_events, 4*n_part + 2)`. Each
  row is `(E0,px0,py0,pz0, E1,px1,py1,pz1, ..., alpha_s, me2)` with
  particles in MadGraph LHE order (initial state first).
- `proc_<slug>_minv_costheta.png` — hexbin plot (if `--plot`).
- `madgraph.log` — full MadGraph stdout/stderr (event generation,
  `make allmatrix2py.so`, etc.).
- `madevent/`, `standalone/` — the two MG output directories.

## Requirements

- MadGraph5_aMC@NLO; pass its `mg5_aMC` binary via `--mg-path` once
  (cached in `.mg_path`), or place a symlink at `./MG/bin/mg5_aMC`.
- Python with `numpy`, `matplotlib`
- `meson` and `ninja` (needed by `numpy.f2py`'s meson backend on
  Python ≥ 3.12) — install via `pip install meson ninja`
- `gfortran`

## Example

```
./generate_amp_dataset.py \
    --mg-path /opt/MG5_aMC_v3_5_4/bin/mg5_aMC \
    --process "q q > z g" \
    --n-events 100000 \
    --m-inv 400,500 \
    --costheta=-0.5,0.5 \
    --plot
```

Generates 100k `q q > Z g` events with `m(Zg) ∈ [400, 500] GeV` and
`|cos θ| < 0.5` (cuts enforced inside MadGraph), computes |M|² per
event, and writes `output/proc_z_g/proc_z_g.npy` plus the hexbin plot.

A subsequent run can drop `--mg-path` and try a different process:

```
./generate_amp_dataset.py --process "g g > t t~" --n-events 50000
```

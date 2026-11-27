# Periodic slab inputs and candidate structures

The surface is a 120-atom ZnO m-plane slab. QE calculations use PBE, 60/480 Ry
wavefunction/density cutoffs, a 2 × 2 × 1 k-point mesh and 0.003 Ry cold smearing.
The archived PBE labels contain no D3 correction. `slab_layers.json` identifies
the fixed layers. The cell height is 27.3209 Å; this is not the empty vacuum gap.

## O85 benchmark

The early matched-geometry test follows O85 and Zn20. Three water-branch
structures have 158 atoms; three release-branch structures have 170. Energies
are compared within each branch, never between the two compositions.

```bash
python slab/o85_benchmark/compare_forces.py --water-jobs 7055417 7055417 7055962 --release-jobs 7055418 7056013 7055961 --mace-job 7055419
```

This parses QE energies and forces, checks the QE/MACE coordinate identity,
and writes `outputs/benchmark_comparison.csv`. The release observable in the
MACE pilot is the sum of forces on the whole departing complex, not the force
on Zn alone. The water observable is H158 projected along S133 → O85.

`build_structures.py` constructs the archived steric starting models from the
relaxed bare slab and published complex coordinates. It does not optimize a
reaction path:

```bash
python slab/o85_benchmark/build_structures.py slab/o85_benchmark/outputs/qe-relax-7055128.out --out /tmp/zno-starting-structures
```

The `source/` directory holds the starting POSCAR and the attributed published
complex geometry. Existing `.in` files are the primary record if numerical
optimizer differences affect newly generated steric placements.

## Zn20–O80 candidates

The proposed same-site mechanism uses **O80**, not O85. `o80_candidates/inputs/`
contains the starting models; `bound.extxyz` is the later water-formed,
Zn-bound candidate used in the proposal. `last-evaluated-frames.extxyz` retains
an earlier A/B/C snapshot; its B frame is superseded by `bound.extxyz` for the
contact analysis. C is the separately relaxed detached candidate.

```bash
python slab/inspect_candidates.py
```

All atom labels are one-based: Zn20 is the tracked Zn; O80 is the selected
water oxygen; O65/O75 are neighboring lattice O atoms; N145/N159 are the two
MeIm donor atoms. The same 170 atoms are retained in the O80 models.
The script applies periodic minimum-image distances and writes
`candidate_contacts.json`. Independent relaxations do not prove a connected
path. These structures are not validated minima or transition states.

## Running QE

Pseudopotentials are not bundled. `pseudopotentials.json` lists exact filenames
and SHA-256 values for SSSP 1.3.0 PBE efficiency. Obtain the matching files from
[SSSP](https://www.materialscloud.org/discover/sssp/table/efficiency), place them
in a `pseudo/` directory under the calculation working directory, and compare
the hashes before running. Inputs use relative `./pseudo` and `./scratch` paths.
For example, from `slab/o85_benchmark/` after providing `pseudo/`:

```bash
mkdir -p scratch
pw.x -in structures/water1.in > water1.out
```

The archived `outputs/*.in` record minor electronic-convergence adjustments
used for the completed calculations; use these for an exact settings audit.

## Cell-height check

`vacuum_check/` records a fixed-geometry comparison of release1 and release3
with Lz increased from 27.3209 to 37.3209 Å. Their relative energy changes by
−0.000614 eV. The two complex-normal forces change by +0.015107 and −0.021014
eV Å⁻¹. This checks those two geometries only; it does not establish convergence
of a transition state, the entire path or the MACE neighbor graph.

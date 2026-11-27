# ZnO dissolution: preliminary calculations

**Harry Linfeng Hu · USC Chemistry · November 2026**

This repository contains the calculations and analysis supporting my qualifying
proposal on ligand-assisted ZnO dissolution. The chemical question is whether
zinc remains difficult to release after water has formed at the surface.
The proposed comparison follows one Zn atom through water formation and zinc
release on a periodic ZnO m-plane slab.

The files here document the preliminary work: cage thermochemistry, published
kinetics, candidate slab structures, and a small MACE pilot. The two surface
transition states have not yet been validated. Endpoint energy differences and
partially relaxed structures must not be read as activation barriers.

## Contents

| Directory | What it contains | Start here |
| --- | --- | --- |
| `cage/` | Gas and implicit-acetonitrile thermochemistry for 2 + PhSH → 3 + H₂O; dispersion checks | [Methods and energy references](cage/README.md) |
| `geometry/` | State-2 alignment, ligand orientations and D3 pair decomposition | [Geometry comparison](geometry/README.md) |
| `kinetics/` | Published SI marker coordinates, fits and Table S1 consistency check | [Data provenance](kinetics/README.md) |
| `slab/` | QE inputs, O85 benchmark, Zn20–O80 candidates and a cell-height check | [Structures and atom numbering](slab/README.md) |
| `mace/` | 117-frame dataset, training/evaluation scripts, archived predictions and path diagnostics | [Pilot design and limitations](mace/README.md) |

## Reproduce the analysis

Python 3.11 or later is required (tested with 3.12.14; package versions are in
`requirements-tested.txt`). From this directory:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python verify.py
```

This checks the archive, rebuilds the gas thermochemistry from the saved mobile
Hessians, aggregates the PCM species energies, refits the published markers,
recomputes geometric contacts and D3 terms, and verifies the MACE split and error
tables. Results are written beside each analysis. It does **not** submit jobs,
run Q-Chem or Quantum ESPRESSO, or retrain MACE.

To repeat an individual calculation, use the commands in its module README.
Q-Chem inputs require Q-Chem; QE inputs require matching pseudopotentials.
Model training requires a separate CUDA/MACE environment and the upstream
MACE-MH-1 checkpoint. The lightweight analysis environment does not install those
programs or model weights.

## Main observations

At 303.15 K, the local cage calculations give estimated reaction free energies
of **+113.89 kJ mol⁻¹ in gas** and **+82.80 kJ mol⁻¹ in implicit acetonitrile**
after adding D3(BJ) single-point corrections to the respective B3LYP geometries
and retaining their B3LYP thermal corrections. Water formation remains uphill
in both environments. The corresponding values without D3 are +133.04 and
+127.90 kJ mol⁻¹. The published gas value, +102.8 kJ mol⁻¹, is a separate
literature reference, not the baseline used to calculate our solvent shift.

The slab optimizations provide candidate structures with water formed and Zn
still attached, and with Zn detached. Their stability and connection through
transition states remain open questions. The early MACE test improves some
force errors but worsens the water-branch relative-energy error; it does not
establish accurate barriers or a net computational speedup.

## Sources and reuse

The starting mechanism, cage coordinates and experimental plot markers come
from Koskela et al., *Chemical Science* (2022),
[DOI: 10.1039/D1SC06667F](https://doi.org/10.1039/D1SC06667F).
The experimental points are extracted from the published figures, not original
instrument records. See [THIRD_PARTY.md](THIRD_PARTY.md) for attribution and
[provenance.json](provenance.json) for source paths and checksums.
Paths beginning `archive/` identify original records; they are not dependencies
needed to run this repository.

Inputs and results were selected from a larger working archive. Scheduler
logs, failed environment setup, duplicate snapshots, manuscript drafts and
personal exam records are omitted. Program outputs retain their scientific
content; machine-specific paths are normalized where needed. Changes are
recorded in the provenance manifest. No blanket license is assigned to this
collection; third-party materials retain their original terms.

Please cite this archive using `CITATION.cff` and cite the original paper when
using its data or structures.

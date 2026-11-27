# MACE pilot

The pilot starts from MACE-MH-1, retaining the native multihead architecture and
selecting `oc20_usemppbe`. The training script updates shared parameters using
PBE energies and forces from QE; it is not a D3-trained potential.

## Data and archived results

There are 117 unique frames after removal of three duplicates: 84 training,
25 validation and 8 test frames. The split follows parent trajectories where
available. The eight test frames mainly change MeIm orientation and provide
limited coverage of bond-breaking and proton-sharing regions. Reference
structures used to define relative energies are not extra independent test
observations. The manifest records their origin and duplicate aliases.

```bash
python mace/check_dataset.py
```

This checks frame identities, split counts and trajectory separation, then
recomputes errors from `reference/test_metrics.csv`.

| Observable | Pretrained MAE | Fine-tuned MAE |
| --- | ---: | ---: |
| Water relative energy / eV | 0.1683 | 0.2914 |
| Water directional force / eV Å⁻¹ | 0.4090 | 0.1066 |
| Release relative energy / eV | 0.3060 | 0.1692 |
| Release directional force / eV Å⁻¹ | 0.4743 | 0.4706 |

Energy improvement is not uniform. No barrier accuracy or net speedup has been
established. The water training branch concerns O85; the later proposed
reaction site is Zn20–O80.

## Repeat training and inference

The archived training run used Python 3.11 and CUDA. The baseline inference
record identifies mace-torch 0.3.15 and the upstream checkpoint SHA-256
`ec00a2705854622fbbd898ccfb7701072fcd674709102d009fb919c1b8cc5dde`. Install PyTorch, ASE, e3nn and a
compatible `mace-torch` version in a separate environment, and obtain MACE-MH-1
from its official release. Model weights are not redistributed here. The full
original package lockfile was not archived, so cross-version numerical identity
of a fresh training run is not claimed.

```bash
python mace/train_native_head.py mace/data /path/to/mace-mh-1.model /tmp/zno-trained
python mace/evaluate_native.py mace/data /path/to/mace-mh-1.model /tmp/zno-trained/zno_mh1_pilot.model /tmp/zno-evaluation
```

Use a trusted upstream checkpoint; these scripts load the serialized model
object. Both output directories must be new. The seed is 42, learning rate
2 × 10⁻⁵, batch size 2 and maximum training length 28 epochs. The archived best
checkpoint was epoch 27. Composition-specific constant offsets are fitted on
the training set. The loss weights are 0.5 for total energy, 3.0 for force
components and 0.15 for the reaction-direction force. The validation score is
a numerical sum of energy and directional-force MAEs with different units;
it is a pilot selection rule, not a dimensionless physical error measure.

The original structures use periodic boundaries in all three directions.
Some extended release structures can acquire neighbors across the nominal
vacuum. This remains a limitation of the pilot and must be checked before
reusing it for an isolated surface. Changing PBC flags changes the model input;
it must not be described as reproducing the original run.

## O80 path diagnostic

`path_pilot/` contains the archived path coverage and residual-force tables.
At eight interior-image QE evaluations, the MACE-initialized path had a maximum
residual of 0.254 eV Å⁻¹ versus 3.260 for the interpolated path. Neither reached
the 0.05 eV Å⁻¹ target. The MACE arm has one additional recorded cycle at ten
evaluations; it is not a matched-cost comparison. The tables are diagnostics,
not a converged NEB result or evidence of overall computational savings.

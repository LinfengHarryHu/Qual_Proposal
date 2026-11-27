# Cage water-formation thermochemistry

Reaction: **2 + PhSH → 3 + H₂O**. All terms in the main comparison refer to
303.15 K. Each reaction difference is the sum of product energies minus the sum
of reactant energies. These are endpoint thermodynamic estimates, not barriers.

## Recalculate

```bash
python cage/si_single_points/ledger.py
python cage/summarize.py
python cage/check_dispersion.py
```

The ledger reconstructs mobile-region harmonic frequencies from the saved gas
Hessians, checks atom ordering against the inputs, and calculates the gas
thermal terms at 30 and 75 °C. `summarize.py` reads the 30 °C row, parses the
native PCM/D3 outputs, and writes `comparison.json`.

| Treatment | ΔG at 30 °C / kJ mol⁻¹ |
| --- | ---: |
| Published gas reference, Table 1 | +102.80 |
| Local gas, B3LYP | +133.04 |
| C-PCM on local gas structures; gas thermal terms retained | +110.67 |
| Local structures reoptimized in C-PCM, B3LYP | +127.90 |
| Local gas structures + D3(BJ); B3LYP thermal terms retained | +113.89 |
| Local C-PCM structures + D3(BJ); B3LYP thermal terms retained | +82.80 |

The fixed-structure estimate is 133.04390 − 22.37787 = 110.66603 kJ mol⁻¹.
The earlier archive used 75 °C and reported 110.98; `results_75C.json` retains
that record. The present script uses only its electronic solvent shift and
applies the 30 °C gas thermal correction. Rounding the input gas value to
133.04 first would give 110.66; the script retains full precision.

## Electronic structure and thermal corrections

The local calculations use all-electron B3LYP/def2-SVP in Q-Chem 6.3, an SG-1
integration grid and a constrained cage backbone. Twenty cage atoms are fixed:
5–24 in state 2 and 4–23 in state 3 (one-based numbering). Gas thermal terms for
these cages use the mobile Cartesian Hessian, without assigning whole-cage
translation or rotation. Free PhSH and water include molecular translation,
rotation and vibration; water's rotational symmetry number is corrected to 2.

The solvent calculations use C-PCM acetonitrile (dielectric constant 36.60),
SwiG surface discretization and UFF radii scaled by 1.20. The native outputs
call the polarization-inclusive electronic energy “Total Free Energy”; this
quantity alone is not the thermochemical Gibbs free energy. `summarize.py`
adds the archived vibrational/molecular thermal correction.

`pcm/thermochemistry.json` contains the accepted frequencies and the individual
ZPE, thermal enthalpy, entropy, symmetry and standard-state terms. Cage terms
use the mobile region. The free molecules receive the same 1 atm → 1 M
correction (+8.10071 kJ mol⁻¹ at 303.15 K), which cancels between PhSH and H₂O.
The PCM reaction thermal correction is +7.86989 kJ mol⁻¹. Full finite-difference
force logs for the PCM Hessians are omitted; this package reproduces aggregation
from the archived frequencies/thermal ledger, not their original force jobs.

D3(BJ) uses B3LYP parameters with the three-body ATM term off. Gas corrections
were evaluated with simple-dftd3 1.3.1; PCM corrections were obtained from native
Q-Chem single points. Coordinates and frequencies were **not** reoptimized with
D3. The gas correction is −19.15378 kJ mol⁻¹ and the PCM correction is about
−45.09 kJ mol⁻¹. D3 therefore changes both the endpoint difference and the
apparent solvent shift for these particular structures.

State 3 in PCM met the historical maximum-gradient criterion at a saved frame
(2.73 × 10⁻⁴ Eh bohr⁻¹), but the optimizer did not terminate normally. Its
mobile-region frequencies were positive. This is why its ΔG is labeled an
estimate. These constraints and a single local geometry limit comparison with
an unconstrained molecular free energy or a conformational ensemble.

## Literature comparison

`si_baseline/` and `si_single_points/` hold single points on the published SI
coordinates. Their electronic energies can be audited separately from the
locally optimized geometries. Transferring the local thermal correction to
those SI/D3 energies gives +111.49 kJ mol⁻¹; it is **not** the +113.89 local-gas
value. Table 1 does not explicitly resolve whether the reported energies
include the D3(BJ) check described in the methods. No paper-to-local difference
is treated as a pure solvent effect.

The `.in` files record the actual method settings. For an individual archived
single point, for example:

```bash
qchem -nt 8 cage/pcm/state2.in state2.out
```

Rebuilding the full optimization/frequency calculation requires the constrained
optimization and mobile-Hessian procedure, not just this single-point input.

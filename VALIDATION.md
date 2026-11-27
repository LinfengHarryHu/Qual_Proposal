# Validation record

2026-11-26 · Python 3.12.14 · macOS

The analysis also passed from an isolated copy outside the original workspace,
launched with `/tmp` as the working directory. The lightweight analysis completed
successfully using the versions in
`requirements-tested.txt`:

- Gas mobile Hessians reproduced the archived 30 and 75 °C thermal ledger.
- Four gas and four PCM D3 calculations matched the archived dispersion terms.
- Product/reactant atom counts balanced in both cage environments.
- PCM species energies and thermal terms reproduced the archived reaction ΔG.
- Geometry alignment and pairwise dispersion reproduced the recorded metrics.
- Re-extraction from the source SI reproduced all 38 marker records exactly.
- QE/MACE comparison checked matching coordinates before comparing energies/forces.
- O80 candidate contacts, 84/25/8 split membership and test errors were checked.
- PNG author metadata reads Harry Linfeng Hu; pixel content is unchanged.

Quantum-chemistry jobs, MACE training and OVITO rendering were not rerun.
The archived PCM thermal ledger was used; its original finite-difference
force calculations were not reconstructed. A fresh training environment has
not been validated by these analysis checks.

Run `python verify.py` to repeat the included checks. It writes a short
machine-readable `verification.json` after success. Source-PDF re-extraction
is a separate command in `kinetics/README.md` because the PDF is not bundled.

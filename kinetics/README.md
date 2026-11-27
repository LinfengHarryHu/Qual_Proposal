# Published kinetics: extraction and numerical checks

`marker_points.csv` contains vector-marker centers from Figures S9–S12 and
S14–S15 of the Koskela supporting information. Each record retains its PDF page,
drawing index and page coordinates. Two labeled ticks calibrate each axis;
ligand-dependence panels are fitted in the printed log₁₀ coordinates.
These are extracted published markers, not raw experimental observations or
independent replicate measurements.

```bash
python kinetics/analyze.py
```

The four concentration-versus-time slopes are approximately 3.09, 4.80, 7.65
and 14.84 μmol L⁻¹ s⁻¹ at 30, 40, 55 and 75 °C. Table S1 prints 1.48 × 10⁻⁶
at 75 °C, whereas Figure S12 and the table's transformed ordinate support
approximately 1.48 × 10⁻⁵ mol L⁻¹ s⁻¹. Both records are retained in the reference
JSON; the fit follows the plotted points.

The script also reproduces the paper's numerical convention R ln(rh/kBT).
Because r is a volumetric product flux rather than a first-order rate constant,
the fitted intercept is not a molecular activation entropy without a site and
standard-state model. The numerical apparent enthalpy is 27.445 kJ mol⁻¹.
OLS uncertainty from plot markers does not include experimental systematic
errors, digitization uncertainty or sample-to-sample variation.

There is an unresolved absolute-concentration issue: the nominal 0.446 mmol
ZnO in 4.5 mL corresponds to 99.1 mM, below the largest extracted S10 value
(about 135 mM). The archive preserves the published axes rather than silently
rescaling them. Normalized temperature comparisons additionally assume a common
reacting-site concentration across temperatures.

To repeat the extraction, obtain the SI PDF associated with
[the paper](https://doi.org/10.1039/D1SC06667F), then run:

```bash
python kinetics/extract_markers.py /path/to/SC-013-D1SC06667F-s012.pdf --output extracted.csv
```

The script checks SHA-256
`ec14ab5c79462b2d9211a0f3ebf3d54489cf2493de77209b0f06b48eb9db5f51`.
A different PDF must be inspected before reusing drawing indices.

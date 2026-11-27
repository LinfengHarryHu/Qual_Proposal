"""Fit the published plot markers and check the Table S1 rate convention."""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.constants import R, h, k
from scipy.stats import linregress

ROOT = Path(__file__).resolve().parent


def main():
    with (ROOT/'marker_points.csv').open() as f:
        rows=list(csv.DictReader(f))
    fits={}
    for name in ('S9','S10','S11','S12','S14','S15'):
        group=[r for r in rows if r['figure']==name]
        fit=linregress([float(r['x']) for r in group],[float(r['y']) for r in group])
        fits[name]={'n':len(group),'slope':fit.slope,'intercept':fit.intercept,
                    'r_squared':fit.rvalue**2,'temperature_C':int(group[0]['temperature_c'])}
    names=('S9','S10','S11','S12')
    temperatures=np.array([fits[n]['temperature_C']+273.15 for n in names])
    rates=np.array([fits[n]['slope'] for n in names])
    # Reproduce the paper's numerical flux convention, not a molecular rate constant.
    fit=linregress(1/temperatures,R*np.log(rates*h/(k*temperatures)))
    result={'fits':fits,'apparent_enthalpy_kJ_mol':-fit.slope/1000,
            'flux_convention_entropy_J_mol_K':fit.intercept,
            'rate_75C_from_printed_ordinate':float(np.exp(-338.64/R)*k*348.15/h)}
    reference=json.loads((ROOT/'reference.json').read_text())
    for name in fits:
        if not np.isclose(fits[name]['slope'],reference['fits'][name]['slope'],rtol=1e-10):
            raise ValueError(f'Fit changed: {name}')
    (ROOT/'fits.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()

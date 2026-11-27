"""Reconstruct the 303.15 K endpoint comparison from archived species data."""
import csv
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
HARTREE_KJ = 2625.49963948
SIGNS = {'state2': -1, 'phsh': -1, 'state3': 1, 'water': 1}


def property_value(text, label):
    values = re.findall(re.escape(label) + r'\s*=?\s*([-+0-9.Ee]+)', text)
    if not values:
        raise ValueError(f'Missing {label}')
    return float(values[-1])


def main():
    thermal = json.loads((ROOT/'pcm/thermochemistry.json').read_text())['completed']
    manifest = json.loads((ROOT/'reference/pcm_manifest.json').read_text())['species']
    energies, dispersion = {}, {}
    for name in SIGNS:
        text = (ROOT/'pcm'/f'{name}.out').read_text()
        if 'Thank you very much for using Q-Chem' not in text or 'Doing -D3(BJ)' not in text:
            raise ValueError(f'Incomplete D3 calculation: {name}')
        # The final PCM energy includes polarization and the D3 correction.
        energies[name] = property_value(text, 'Total Free Energy (H0 + V/2 + non-elec)')
        dispersion[name] = property_value(text, '-D3 energy without 3body term')
        if abs(energies[name] - dispersion[name] - manifest[name]['baseline_e_eh']) > 5e-7:
            raise ValueError(f'Non-D3 baseline changed: {name}')
    delta = lambda values: sum(SIGNS[s]*values[s] for s in SIGNS)
    electronic = delta(energies)*HARTREE_KJ
    correction = delta({s:thermal[s]['g_correction_kj_mol'] for s in SIGNS})
    with (ROOT/'si_single_points/reaction_ledger.csv').open() as f:
        gas = next(r for r in csv.DictReader(f) if float(r['T_C']) == 30)
    gas0 = float(gas['optimized_no_d3_delta_G_kJ_mol'])
    gas_d3 = json.loads((ROOT/'gas/dispersion.json').read_text())
    fixed = json.loads((ROOT/'fixed_structure_pcm/results_75C.json').read_text())
    shift = delta({s:float(fixed['species'][s]['pcm_Eh'])-float(fixed['species'][s]['gas_Eh']) for s in SIGNS})*HARTREE_KJ
    values = {
        'published_gas': float(gas['published_delta_G_kJ_mol']),
        'own_gas_no_d3': gas0,
        'fixed_structure_pcm_no_d3_estimate': gas0+shift,
        'reoptimized_pcm_no_d3_estimate': electronic-delta(dispersion)*HARTREE_KJ+correction,
        'own_gas_d3_estimate': gas0+gas_d3['reaction_d3bj_shift_kJ_per_mol'],
        'reoptimized_pcm_d3_estimate': electronic+correction,
    }
    result = {'temperature_K':303.15, 'reaction':'2 + PhSH -> 3 + H2O',
              'units':'kJ mol^-1', 'delta_G':values,
              'pcm_d3_delta_E':electronic, 'pcm_thermal_correction':correction,
              'fixed_structure_pcm_shift':shift}
    (ROOT/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    for name, value in values.items():print(f'{name:42s} {value:9.3f} kJ mol^-1')


if __name__ == '__main__':
    main()

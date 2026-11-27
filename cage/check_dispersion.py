"""Check gas and PCM dispersion terms at the archived coordinates."""
import csv
import json
from pathlib import Path
import numpy as np
from ase.io import read
from dftd3.interface import DispersionModel, RationalDampingParam

ROOT=Path(__file__).resolve().parent


def main():
    gas=json.loads((ROOT/'gas/dispersion.json').read_text())['species']
    with (ROOT/'reference/pcm_species_results.csv').open() as f:
        pcm={r['species']:float(r['d3bj_e_eh']) for r in csv.DictReader(f)}
    parameters=RationalDampingParam(method='b3lyp',atm=False)
    for environment in ('gas','pcm'):
        for species in gas:
            suffix='_opt.xyz' if environment=='gas' else '.xyz'
            atoms=read(ROOT/environment/(species+suffix))
            model=DispersionModel(atoms.numbers,atoms.positions/0.529177210903)
            energy=float(model.get_dispersion(parameters,grad=False)['energy'])
            expected=gas[species]['d3bj_energy_hartree'] if environment=='gas' else pcm[species]
            assert abs(energy-expected)<1e-8,(environment,species,energy,expected)
    print('Eight gas/PCM D3 single points agree with the archive.')


if __name__=='__main__':
    main()

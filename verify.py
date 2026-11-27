"""Run the archived analyses and check their scientific and provenance invariants."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from collections import Counter
import numpy as np
from ase.io import read
from PIL import Image

ROOT=Path(__file__).resolve().parent


def run(*args):
    subprocess.run([sys.executable,*args],cwd=ROOT,check=True)


def main():
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        digest,name=line.split('  ',1)
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    run('cage/si_single_points/ledger.py')
    run('cage/summarize.py')
    run('cage/check_dispersion.py')
    run('geometry/inspect_state2.py')
    run('kinetics/analyze.py')
    run('slab/o85_benchmark/compare_forces.py','--water-jobs','7055417','7055417','7055962',
        '--release-jobs','7055418','7056013','7055961','--mace-job','7055419')
    run('slab/inspect_candidates.py')
    run('mace/check_dataset.py')
    comparison=json.loads((ROOT/'cage/comparison.json').read_text())['delta_G']
    reference=json.loads((ROOT/'cage/reference/pcm_reaction_summary.json').read_text())
    assert abs(comparison['reoptimized_pcm_d3_estimate']-reference['pcm_d3bj_delta_G_estimate_kj_mol'])<1e-5
    assert abs(comparison['own_gas_d3_estimate']-reference['gas_d3bj_same_optimized_geometry_delta_G_estimate_kj_mol'])<1e-5
    for environment in ('gas','pcm'):
        balance=Counter()
        for species,sign in [('state2',-1),('phsh',-1),('state3',1),('water',1)]:
            suffix='_opt.xyz' if environment=='gas' else '.xyz'
            atoms=read(ROOT/'cage'/environment/(species+suffix))
            for element,n in Counter(atoms.get_chemical_symbols()).items():balance[element]+=sign*n
        assert all(n==0 for n in balance.values()),balance
    actual=json.loads((ROOT/'geometry/metrics.json').read_text())
    expected=json.loads((ROOT/'geometry/reference.json').read_text())
    for key in ('phenyl_plane_angle_degrees','MeIm_plane_angle_degrees','D3_pcm_minus_gas_kJ_mol'):
        assert np.isclose(actual[key],expected[key],atol=1e-7),key
    assert Image.open(ROOT/'geometry/state2_overlay.png').info['Author']=='Harry Linfeng Hu'
    report={'status':'passed','modules':5,'checks':['file checksums','gas Hessian reconstruction',
            'PCM energy accounting','reaction atom balance','geometry and D3 decomposition',
            'published-marker fits','matched QE/MACE coordinates','O80 contacts','dataset splits and error tables','PNG author'],
            'not_run':['new quantum calculations','MACE training','OVITO rendering']}
    (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('All archive checks passed.')


if __name__=='__main__':
    main()

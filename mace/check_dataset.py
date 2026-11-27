"""Check split membership and recalculate held-out errors from archived predictions."""
import csv
import json
from pathlib import Path
import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parent


def main():
    manifest=json.loads((ROOT/'data/manifest.json').read_text())
    splits={s:read(ROOT/'data'/f'{s}.extxyz',index=':') for s in ('train','valid','test')}
    names={s:{a.info['frame'] for a in atoms} for s,atoms in splits.items()}
    for s,atoms in splits.items():
        assert len(atoms)==len(names[s])==manifest['counts'][s]
        assert all(a.info['partition']==s for a in atoms)
        assert all(np.isfinite(a.arrays['REF_forces']).all() for a in atoms)
    assert not any(names[a]&names[b] for a,b in [('train','valid'),('train','test'),('valid','test')])
    trajectories={s:{r.get('parent_trajectory') for r in manifest['frames']
                     if r['partition']==s and r.get('parent_trajectory')} for s in names}
    assert not trajectories['train']&trajectories['valid']
    with (ROOT/'reference/test_metrics.csv').open() as f:rows=list(csv.DictReader(f))
    reference=json.loads((ROOT/'reference/summary.json').read_text())['groups']
    metrics={}
    for group in ('water','release'):
        subset=[r for r in rows if r['group']==group];assert len(subset)==4
        metrics[group]={}
        for quantity in ('relative_eV','directional_eV_A'):
            for model in ('original','tuned'):
                mae=float(np.mean([abs(float(r[f'{model}_{quantity}'])-float(r[f'qe_{quantity}'])) for r in subset]))
                assert abs(mae-reference[group][quantity][model+'_mae'])<1e-10
                metrics[group][model+'_'+quantity+'_MAE']=mae
    (ROOT/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    print('Split counts:',{s:len(a) for s,a in splits.items()})
    print(json.dumps(metrics,indent=2))


if __name__=='__main__':
    main()

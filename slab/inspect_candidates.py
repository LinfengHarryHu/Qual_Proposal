"""Measure local contacts in the archived Zn20-O80 candidates."""
import json
from pathlib import Path
from ase.io import read

ROOT=Path(__file__).resolve().parent
PAIRS=((20,65),(20,75),(20,80),(80,157),(80,158),(20,145),(20,159))


def main():
    prior=read(ROOT/'o80_candidates/last-evaluated-frames.extxyz',index=':')
    frames={a.info['case']:a for a in prior}
    frames['B']=read(ROOT/'o80_candidates/bound.extxyz')
    rows={}
    for name,atoms in frames.items():
        assert len(atoms)==170 and atoms[19].symbol=='Zn' and atoms[79].symbol=='O'
        rows[name]={f'{i}-{j}':atoms.get_distance(i-1,j-1,mic=True) for i,j in PAIRS}
    (ROOT/'candidate_contacts.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':
    main()

"""Evaluate the pretrained and validation-selected MACE models on held-out ASE inputs."""

import csv
import json
from pathlib import Path
import sys

from ase.geometry import find_mic
from ase.io import read
from mace.calculators import MACECalculator
import numpy as np
import torch


def direction_force(atoms, forces):
    if atoms.info["reaction"] == "release":
        return float(forces[[19, *range(120,156), *range(158,170)],2].sum())
    vector = find_mic(atoms.positions[84]-atoms.positions[132], atoms.cell, atoms.pbc)[0]
    return float(forces[157] @ (vector/np.linalg.norm(vector)))


def predict(atoms, calculators):
    energies = {"qe": float(atoms.info["REF_energy"])}
    forces = {"qe": direction_force(atoms, atoms.arrays["REF_forces"])}
    for key, calculator in calculators.items():
        atoms.calc = calculator
        energies[key] = float(atoms.get_potential_energy())
        forces[key] = direction_force(atoms, atoms.get_forces())
    if not (np.isfinite(list(energies.values())).all() and np.isfinite(list(forces.values())).all()):
        raise ValueError(f"Nonfinite prediction: {atoms.info['frame']}")
    return energies, forces


def main(data, foundation, model, out):
    out.mkdir()
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA unavailable for the held-out evaluation")
    calculators = {
        "original": MACECalculator(model_paths=foundation, device="cuda", default_dtype="float64",
                                   head="oc20_usemppbe"),
        "tuned": MACECalculator(model_paths=model, device="cuda", default_dtype="float64",
                                head="oc20_usemppbe"),
    }
    refs = {a.info["reaction"]: predict(a, calculators)[0] for a in read(data/"references.extxyz", index=":")}
    samples = read(data/"test.extxyz", index=":")
    if len(samples) != 8 or any(a.info["partition"] != "test" for a in samples):
        raise ValueError("Evaluation split differs from the 8-frame hold-out")
    rows = []
    for atoms in samples:
        group = atoms.info["reaction"]
        energy, force = predict(atoms, calculators)
        relative = {key: value-refs[group][key] for key,value in energy.items()}
        rows.append(dict(group=group, frame=atoms.info["frame"], qe_job=atoms.info["qe_source"],
                         **{f"{k}_relative_eV":v for k,v in relative.items()},
                         **{f"{k}_directional_eV_A":v for k,v in force.items()}))
    report = {"reference_structures":{"water":"water1", "release":"release1"},
              "force_observables":{"water":"H158 along S133 to O85", "release":"sum of complex forces along +z"},
              "checkpoint":str(model), "groups":{}}
    for group in ("water","release"):
        subset = [r for r in rows if r["group"] == group]
        assert len(subset)==4
        metrics={"n": len(subset)}
        for observable in ("relative_eV","directional_eV_A"):
            errors = {model: float(np.mean([abs(r[f"{model}_{observable}"]-r[f"qe_{observable}"]) for r in subset]))
                      for model in calculators}
            metrics[observable] = {
                "original_mae": errors["original"], "tuned_mae":errors["tuned"],
                "absolute_improvement":errors["original"]-errors["tuned"],
                "improvement_percent":100*(1-errors["tuned"]/errors["original"])
            }
        metrics["direction_sign_matches"]={model:sum(r[f"{model}_directional_eV_A"]*r["qe_directional_eV_A"]>0 for r in subset)
                                            for model in calculators}
        report["groups"][group]=metrics
    with (out/"test_metrics.csv").open("w", newline="") as stream:
        writer=csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    (out/"summary.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(report,indent=2),flush=True)


if __name__ == "__main__":
    make_data, base, tuned, output=sys.argv[1:]
    main(Path(make_data),base,tuned,Path(output))

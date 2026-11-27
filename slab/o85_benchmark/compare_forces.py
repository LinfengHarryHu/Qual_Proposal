"""Compare PBE and MACE energies/forces on identical periodic frames."""

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import re

import numpy as np
from scipy.constants import physical_constants

from build_structures import NORMAL, unit, wrap_displacement

RY_EV = physical_constants["Rydberg constant times hc in eV"][0]
BOHR_A = physical_constants["Bohr radius"][0] * 1e10
NUMBER = r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?"
BASE = Path(__file__).resolve().parent


def check_qe_geometry(name, path, symbols, xyz, cell):
    text = path.read_text()
    lines = text.split("ATOMIC_POSITIONS angstrom\n", 1)[1].split("K_POINTS", 1)[0].strip().splitlines()
    rows = [line.split() for line in lines]
    lattice = text.split("CELL_PARAMETERS angstrom\n", 1)[1].splitlines()[:3]
    if ([row[0] for row in rows] != symbols
            or not np.array_equal(np.array([row[1:] for row in rows], dtype=float), xyz)
            or not np.array_equal(np.array([row.split() for row in lattice], dtype=float), np.diag(cell))):
        raise ValueError(f"{name}: QE geometry and MACE input differ")


def qe_result(path, symbols):
    text = path.read_text()
    for marker in ("convergence has been achieved", "JOB DONE.", "Forces acting on atoms"):
        if marker not in text:
            raise ValueError(f"{path}: missing {marker}")
    e = float(re.search(rf"^!\s*total energy\s*=\s*({NUMBER}) Ry", text, re.M)[1])
    block = text.split("Forces acting on atoms (cartesian axes, Ry/au):", 1)[1].split("Total force", 1)[0]
    rows = re.findall(rf"atom\s+(\d+) type\s+(\d+)\s+force =\s*({NUMBER})\s*({NUMBER})\s*({NUMBER})", block)
    expected_type = {el: i for i, el in enumerate(("Zn", "O", "S", "C", "H", "N"), 1)}
    if [int(row[0]) for row in rows] != list(range(1, len(symbols)+1)):
        raise ValueError(f"{path}: force atom indices differ")
    if [int(row[1]) for row in rows] != [expected_type[el] for el in symbols]:
        raise ValueError(f"{path}: force atom types differ")
    forces = np.asarray([row[2:] for row in rows], dtype=float) * (RY_EV / BOHR_A)
    smearing = re.search(rf"smearing contrib\.\s*\(-TS\)\s*=\s*({NUMBER}) Ry", text)
    return {"energy": e * RY_EV, "forces": forces,
            "smearing_eV": None if smearing is None else float(smearing[1])*RY_EV}


def read_xyz(path):
    lines = path.read_text().splitlines()
    if 'pbc="T T T"' not in lines[1]:
        raise ValueError(f"{path}: missing periodic cell")
    lattice = np.array(re.search(r'Lattice="([^"]+)"', lines[1])[1].split(), dtype=float).reshape(3, 3)
    rows = [line.split() for line in lines[2:]]
    symbols = [row[0] for row in rows]
    xyz = np.array([row[1:] for row in rows], dtype=float)
    return symbols, xyz, lattice.diagonal()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--water-jobs", nargs=3)
    parser.add_argument("--release-jobs", nargs=3)
    parser.add_argument("--mace-job", required=True)
    args = parser.parse_args()
    groups = [(name, jobs) for name, jobs in (("water", args.water_jobs), ("release", args.release_jobs)) if jobs]
    if not groups:
        parser.error("provide at least one QE series")
    metadata = json.loads((BASE / "structures/metadata.json").read_text())
    if metadata["preview_on_unrelaxed_POSCAR"]:
        raise ValueError("Only relaxed slab structures are accepted")
    results = BASE / "outputs"
    mace = json.loads((results / f"mace-series-{args.mace_job}.json").read_text())
    if mace["head"] != "oc20_usemppbe":
        raise ValueError("Unexpected MACE head")
    labels = {name: value-1 for name, value in metadata["labels"].items()}
    movers = [i-1 for i in metadata["release_movers"]]
    all_rows = []
    for group, jobs in groups:
        reference = None
        for point, job in enumerate(jobs, 1):
            name = f"{group}{point}"
            path = BASE / "structures" / (name + ".xyz")
            symbols, xyz, cell = read_xyz(path)
            if sha256(path.read_bytes()).hexdigest() != mace["structures"][name]["xyz_sha256"]:
                raise ValueError(f"{name}: MACE structure hash differs")
            check_qe_geometry(name, results/f"{name}-{job}.in", symbols, xyz, cell)
            qe = qe_result(results/f"{name}-{job}.out", symbols)
            classical = mace["structures"][name]
            forces = [qe["forces"], np.array(classical["forces_eV_per_A"])]
            energies = [qe["energy"], classical["energy_eV"]]
            if reference is None:
                reference = energies
            row = {"frame": name, "qe_job": job, "mace_job": args.mace_job,
                   "qe_relative_eV": energies[0]-reference[0],
                   "mace_relative_eV": energies[1]-reference[1],
                   "force_component_rmse_eV_per_A": float(np.sqrt(np.mean((forces[0]-forces[1])**2))),
                   "qe_smearing_eV": qe["smearing_eV"]}
            extra = metadata["frames"][name]
            if group == "water":
                s, o, h = labels["incoming_S"], labels["water_O"], labels["transferring_H"]
                direction = unit(wrap_displacement(xyz[o]-xyz[s], cell))
                row.update({"coordinate_A": extra["S_H_A"]-extra["O_H_A"],
                            "key_force": "H along S to O",
                            "companion_force": "water O outward"})
                first = [force[h] @ direction for force in forces]
                second = [force[o] @ NORMAL for force in forces]
            else:
                zn = labels["anchor_Zn"]
                row.update({"coordinate_A": extra["anchor_Zn_height_A"],
                            "key_force": "Zn outward",
                            "companion_force": "whole Zn complex outward"})
                first = [force[zn] @ NORMAL for force in forces]
                second = [force[movers].sum(axis=0) @ NORMAL for force in forces]
            for method, a, b in zip(("qe", "mace"), first, second):
                row[f"{method}_key_force_eV_per_A"] = float(a)
                row[f"{method}_companion_force_eV_per_A"] = float(b)
            all_rows.append(row)

    name = "benchmark_comparison" if len(groups) == 2 else f"{groups[0][0]}_comparison"
    with (results / f"{name}.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    for row in all_rows:
        print(f"{row['frame']:9s} ΔE(QE/MACE) {row['qe_relative_eV']:8.3f} {row['mace_relative_eV']:8.3f} eV;"
              f" key F {row['qe_key_force_eV_per_A']:8.3f} {row['mace_key_force_eV_per_A']:8.3f} eV/Å")


if __name__ == "__main__":
    main()

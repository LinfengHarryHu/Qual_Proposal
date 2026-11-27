"""Generate matched frozen geometries for QE and MACE from a relaxed slab."""

import argparse
import json
from pathlib import Path
import re

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parent
BASELINE = ROOT / "source"
SI = ROOT / "source"
ZN, WATER_O = 19, 84  # 0-based positions in the original slab
NORMAL = np.array([0.0, 0.0, 1.0])
H_BOND = .98
# Molecular product (s011) fragment ordering from its published coordinates.
PHS_FIRST = [2, *range(4, 15)]
PHS_SECOND = [3, *range(15, 26)]
MEIM_FIRST = [32, 26, 27, 28, 29, 30, 31, 33, 46, 47, 48, 49]
MEIM_SECOND = [40, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45]
COVALENT = {"Zn": 1.22, "O": .66, "S": 1.05, "C": .76, "H": .31, "N": .71}
RADII = {"Zn": 1.39, "O": 1.52, "S": 1.80, "C": 1.70, "H": 1.20, "N": 1.55}
PSEUDO = {
    "Zn": (65.38, "Zn_pbe_v1.uspp.F.UPF"),
    "O": (15.999, "O.pbe-n-kjpaw_psl.0.1.UPF"),
    "S": (32.06, "s_pbe_v1.4.uspp.F.UPF"),
    "C": (12.011, "C.pbe-n-kjpaw_psl.1.0.0.UPF"),
    "H": (1.008, "H.pbe-rrkjus_psl.1.0.0.UPF"),
    "N": (14.007, "N.pbe-n-radius_5.UPF"),
}


def unit(vector):
    return vector / np.linalg.norm(vector)


def wrap_displacement(displacement, cell):
    frac = displacement / cell
    return (frac - np.rint(frac)) * cell


def distance(first, second, cell):
    return np.linalg.norm(wrap_displacement(first - second, cell))


def xyz_source(number):
    path = SI / f"SC-013-D1SC06667F-s{number:03d}.xyz"
    lines = path.read_text().splitlines()
    rows = [row.split() for row in lines[2:] if row.strip()]
    if int(lines[0]) != len(rows):
        raise ValueError(f"Invalid XYZ {path}")
    return [row[0] for row in rows], np.array([row[1:] for row in rows], dtype=float)


def read_slab(path, trial):
    reference = (BASELINE / "POSCAR").read_text().splitlines()
    matrix = np.array([row.split() for row in reference[2:5]], dtype=float)
    cell = matrix.diagonal()
    if not np.allclose(matrix, np.diag(cell)):
        raise ValueError("The seed generator expects an orthorhombic cell")
    reference_pos = np.array([row.split()[:3] for row in reference[8:128]], dtype=float) * cell
    if trial:
        return cell, reference_pos
    text = path.read_text()
    if "bfgs converged" not in text.lower() or "JOB DONE." not in text:
        raise ValueError("The QE slab optimization lacks a completion/convergence marker")
    block = text.rsplit("Begin final coordinates", 1)[1].split("End final coordinates", 1)[0]
    section = re.search(r"^ATOMIC_POSITIONS\s*\((crystal|angstrom)\)\s*\n([\s\S]+)", block, re.M)
    if not section:
        raise ValueError("Missing final slab coordinates")
    rows = [line.split() for line in section.group(2).splitlines() if line.strip()]
    if [row[0] for row in rows] != ["Zn"] * 60 + ["O"] * 60:
        raise ValueError("Final slab atom order differs")
    final = np.array([row[1:4] for row in rows], dtype=float)
    if section.group(1) == "crystal":
        final *= cell
    fixed = np.array(json.loads((ROOT / "slab_layers.json").read_text())["fixed_atoms"]) - 1
    if abs(wrap_displacement(final[fixed] - reference_pos[fixed], cell)).max() > 1e-5:
        raise ValueError("Frozen slab layers moved")
    return cell, final


def rotate_vector(source, target):
    source, target = unit(source), unit(target)
    cross = np.cross(source, target)
    norm = np.linalg.norm(cross)
    if norm < 1e-10:
        if np.dot(source, target) > 0:
            return Rotation.identity()
        basis = np.eye(3)[int(np.argmin(abs(source)))]
        return Rotation.from_rotvec(np.pi * unit(np.cross(source, basis)))
    return Rotation.from_rotvec(np.arctan2(norm, np.dot(source, target)) * cross / norm)


def overlap(left, right, limits, cell):
    distances = np.linalg.norm(wrap_displacement(left[:, None] - right[None, :], cell), axis=-1)
    return float(np.sum(np.maximum(0.0, limits - distances) ** 2))


def make_limits(left, right):
    return .70 * (np.array([RADII[a] for a in left])[:, None] +
                  np.array([RADII[a] for a in right])[None, :])


def water_group(slab, cell):
    elements, product = xyz_source(11)
    selected = [i - 1 for i in PHS_FIRST + PHS_SECOND + MEIM_FIRST]
    ligand_elements = [elements[i] for i in selected]
    centered = product - product[0]
    template = centered[selected]
    # Ligand-bound surface Zn is initialized at a small outward displacement.
    lattice_start = slab.copy()
    lattice_start[ZN] += .70 * NORMAL
    site = lattice_start[ZN]
    oxygen = site + wrap_displacement(slab[WATER_O] - site, cell)
    fixed = np.delete(lattice_start, WATER_O, axis=0)
    fixed_symbols = ["Zn"] * 60 + ["O"] * 59
    cross_limits = make_limits(ligand_elements + ["H", "H"], fixed_symbols)
    detached_o_limits = make_limits(["O"], ligand_elements)
    branch = slice(12, 24)
    atom_groups = (slice(0, 12), branch, slice(24, 36))
    group_limits = {(a, b): make_limits(ligand_elements[atom_groups[a]], ligand_elements[atom_groups[b]])
                    for a, b in ((0, 1), (0, 2), (1, 2))}
    hydrogen_limits = make_limits(["H", "H"], ligand_elements)
    # The donor points above the lattice; scanning the internal roll lets each ligand
    # face away from the solid while preserving the published ligand conformation.
    def oriented_fragment(group, azimuth, height):
        direction = np.array([np.cos(azimuth), np.sin(azimuth), height])
        direction = unit(direction)
        initial = rotate_vector(template[group.start], direction).apply(template[group])
        rolls = [Rotation.from_rotvec(t * direction).apply(initial)
                 for t in np.linspace(0, 2*np.pi, 24, endpoint=False)]
        outward = unit(direction - direction[2]*NORMAL)
        return max(rolls, key=lambda x: .85*x[:, 2].min() + .95*(x[1:] @ outward).mean()) + site

    all_symbols = ["Zn"] * 60 + ["O"] * 60 + ligand_elements + ["H", "H"]
    size = len(all_symbols)
    pairs = np.triu_indices(size, 1)
    radii = np.array([COVALENT[element] for element in all_symbols])
    hard = np.maximum(.74, .80*(radii[pairs[0]] + radii[pairs[1]]))
    groups = np.zeros(size, dtype=int)
    groups[120:132], groups[132:144], groups[144:156] = 1, 2, 3
    groups[[WATER_O, 156, 157]] = 4
    pair_vdw = np.array([RADII[e] for e in all_symbols])
    separated = groups[pairs[0]] != groups[pairs[1]]
    intentional = (pairs[0] == 132) & (pairs[1] == 157)  # S-H being transferred
    for_resolved_contact = separated & ~intentional
    hard[for_resolved_contact] = np.maximum(hard[for_resolved_contact], .64*(
        pair_vdw[pairs[0][for_resolved_contact]] + pair_vdw[pairs[1][for_resolved_contact]]))

    def score_geometry(ligands, lift):
        geometries = []
        total_score = 0.0
        expansion = unit(ligands[12] - site)*.52
        for fraction in (0., .5, 1.):
            local = ligands.copy()
            local[branch] += (1.0 - fraction)*expansion
            side = oxygen - local[12]
            side[2] = 0.0
            acceptor = oxygen + fraction*lift*NORMAL + .42*fraction**2*unit(side)
            vector = local[12] - acceptor
            length = np.linalg.norm(vector)
            axis = unit(vector)
            lower, upper = ((2.90, 4.05), (2.48, 3.25), (2.76, 3.95))[int(2*fraction)]
            total_score += 9000*max(lower-length, length-upper, 0)**2
            proton_length = (length - 1.355, length/2, H_BOND)[int(2*fraction)]
            if fraction == 1.0:
                total_score += 200*max(2.76 - length, 0)**2
            proton = acceptor + proton_length*axis
            tangent = unit(NORMAL - np.dot(axis, NORMAL)*axis)
            angle = np.deg2rad(104.5)
            previous_h = acceptor + H_BOND*(np.cos(angle)*axis + np.sin(angle)*tangent)
            ads = np.vstack([local, previous_h, proton])
            total_score += overlap(ads, fixed, cross_limits, cell)
            total_score += overlap(acceptor[None, :], local, detached_o_limits, cell)
            contacts = hydrogen_limits.copy()
            contacts[-1, 12] = 1.30 if fraction < 1.0 else 1.55
            total_score += overlap(ads[36:], local, contacts, cell)
            for a, b in group_limits:
                total_score += overlap(local[atom_groups[a]], local[atom_groups[b]], group_limits[a, b], cell)
            total_score += 2*(max(2.68-length, 0)**2 + max(length-3.30, 0)**2)
            lattice = lattice_start.copy()
            lattice[WATER_O] = acceptor
            geometry = np.vstack([lattice, ads])
            separations = np.linalg.norm(wrap_displacement(geometry[pairs[0]]-geometry[pairs[1]], cell), axis=1)
            total_score += 620*np.sum(np.maximum(0.0, hard-separations)**2)
            geometries.append(geometry)
        return total_score, geometries

    acceptor_azimuth = np.arctan2(*(oxygen - site)[[1, 0]])
    candidates = []
    for offset in range(-75, 76, 15):
        azimuth = acceptor_azimuth + np.deg2rad(offset)
        for height in (.65, .95, 1.30, 1.60):
            for handedness in (-1, 1):
                arm = handedness*2*np.pi/3
                ligands = np.vstack([oriented_fragment(atom_groups[0], azimuth+arm, 1.1),
                                     oriented_fragment(atom_groups[1], azimuth, height),
                                     oriented_fragment(atom_groups[2], azimuth-arm, 1.1)])
                for lift in (2.45, 2.8, 3.15):
                    score, _ = score_geometry(ligands, lift)
                    candidates.append((score, ligands, lift))
    polished = []
    # Jointly resolve ligand-ligand clashes without using either model's energies.
    for _, initial, lift in sorted(candidates, key=lambda row:row[0])[:10]:
        ligands = initial.copy()
        for _ in range(2):
            for group in atom_groups:
                axis = unit(ligands[group.start] - site)
                best_roll = None
                for roll in np.linspace(0, 2*np.pi, 16, endpoint=False):
                    frame = ligands.copy()
                    frame[group] = Rotation.from_rotvec(roll*axis).apply(ligands[group]-site) + site
                    score, geometries = score_geometry(frame, lift)
                    if best_roll is None or score < best_roll[0]:
                        best_roll = score, frame, geometries
                score, ligands, geometries = best_roll
        polished.append((score, ligands, lift))

    best = None
    for _, initial, lift in sorted(polished, key=lambda row: row[0])[:4]:
        def rotated(x):
            blocks = [Rotation.from_rotvec(x[3*i:3*i+3]).apply(initial[g]-site)+site
                      for i, g in enumerate(atom_groups)]
            return np.vstack(blocks)

        def objective(x):
            return score_geometry(rotated(x), lift)[0] + .20*(x @ x)

        result = minimize(objective, np.zeros(9), method="L-BFGS-B",
                          bounds=[(-.65, .65)]*9, options={"maxiter": 75, "ftol": 1e-6})
        score, geometries = score_geometry(rotated(result.x), lift)
        if best is None or score < best[0]:
            best = score, geometries
    return all_symbols, best[1], best[0]


def release_group(water_symbols, product, cell):
    elements, solution = xyz_source(11)
    selection = [i-1 for i in PHS_FIRST + PHS_SECOND + MEIM_FIRST + MEIM_SECOND]
    ligand_symbols = [elements[i] for i in selection]
    assert ligand_symbols[:36] == water_symbols[120:156]
    template = (solution - solution[0])[selection]
    site = product[ZN]
    # Precoordinated product geometry near the surface, then rigid lifting.
    offsets = (1.85, 2.55, 3.60)
    solid_idx = [i for i in range(120) if i not in (ZN, WATER_O)]
    water_idx = [WATER_O, 156, 157]
    solid = product[solid_idx]
    liquid = product[water_idx]
    solid_symbols = [water_symbols[i] for i in solid_idx]
    slab_limits = make_limits(ligand_symbols, solid_symbols)
    water_limits = make_limits(["O", "H", "H"], solid_symbols)
    solvent_limits = make_limits(["O", "H", "H"], ["Zn"] + ligand_symbols)
    ligand_groups = np.arange(48)//12
    ii, jj = np.triu_indices(48, 1)
    inter = ligand_groups[ii] != ligand_groups[jj]
    radii = np.array([RADII[el] for el in ligand_symbols])
    limit_internal = .66*(radii[ii] + radii[jj])

    def ligand_score(ligands):
        penalty = 0.0
        for shift in offsets:
            shifted = ligands+site+shift*NORMAL
            penalty += overlap(shifted, solid, slab_limits, cell)
            penalty += 350*overlap(shifted, solid, .90*slab_limits, cell)
        distances = np.linalg.norm(wrap_displacement(ligands[ii]-ligands[jj], cell), axis=1)
        penalty += 220*np.sum(np.maximum(0, limit_internal[inter]-distances[inter])**2)
        return penalty

    rng = np.random.default_rng(2011)
    initial = [rot.apply(template) for rot in Rotation.random(260, random_state=rng)]
    trial = sorted(initial, key=ligand_score)[:14]
    finalists = []
    for ligands in trial:
        ligands = ligands.copy()
        for _ in range(2):
            for start in (0, 12, 24, 36):
                group = slice(start, start+12)
                axis = unit(ligands[start])
                rolls = []
                for roll in np.linspace(0, 2*np.pi, 17, endpoint=False):
                    candidate = ligands.copy()
                    candidate[group] = Rotation.from_rotvec(roll*axis).apply(ligands[group])
                    rolls.append((ligand_score(candidate), candidate))
                _, ligands = min(rolls, key=lambda pair: pair[0])
        finalists.append((ligand_score(ligands), ligands))
    best = None
    for ligand_penalty, ligands in sorted(finalists, key=lambda pair: pair[0])[:5]:
        complex_points = np.vstack([site, site+ligands])
        for radius in (0., .9, 1.8, 2.8, 3.6):
            for phi in np.linspace(0, 2*np.pi, 24, endpoint=False):
                delta = np.array([radius*np.cos(phi), radius*np.sin(phi), .2])
                moved = liquid + delta
                score = ligand_penalty + .025*radius**2
                score += overlap(moved, solid, water_limits, cell)
                for lift in offsets:
                    complex_shifted = complex_points + lift*NORMAL
                    score += overlap(moved, complex_shifted, solvent_limits, cell)
                    score += 350*overlap(moved, complex_shifted, .91*solvent_limits, cell)
                if best is None or score < best[0]:
                    best = (score, ligands, moved, delta)

    geometries = []
    movers = [ZN, *range(120, 156), *range(158, 170)]
    for shift in offsets:
        frame = np.vstack([product, site+best[1][36:]])
        frame[120:156] = site+best[1][:36]
        frame[water_idx] = best[2]
        frame[movers] += shift*NORMAL
        geometries.append(frame)
    return water_symbols + ligand_symbols[36:], geometries, movers, best[0], best[3]


def verify_contacts(symbols, frame, cell):
    pairs = np.triu_indices(len(frame), 1)
    lengths = np.linalg.norm(wrap_displacement(frame[pairs[0]] - frame[pairs[1]], cell), axis=1)
    elements = np.array(symbols)
    involves_hydrogen = (elements[pairs[0]] == "H") | (elements[pairs[1]] == "H")
    radii = np.array([COVALENT[symbol] for symbol in symbols])
    minimum = np.maximum(.74, .80*(radii[pairs[0]] + radii[pairs[1]]))
    bad = np.flatnonzero(lengths < minimum)
    if bad.size:
        i = bad[0]
        pair = [int(pairs[j][i]+1) for j in (0, 1)]
        raise ValueError(f"Short contact {pair}, {lengths[i]:.3f} Å")
    return float(lengths.min()), float(lengths[~involves_hydrogen].min())


def write_geometry(output, key, symbols, geometry, cell):
    rows = "\n".join(f"{symbol:2s} {x: .10f} {y: .10f} {z: .10f}"
                     for symbol, (x, y, z) in zip(symbols, geometry))
    lattice = " ".join(f"{value:.10f}" for value in np.diag(cell).ravel())
    (output / f"{key}.xyz").write_text(f'{len(symbols)}\nLattice="{lattice}" Properties=species:S:1:pos:R:3 pbc="T T T"\n{rows}\n')
    species = "\n".join(f"{symbol} {PSEUDO[symbol][0]} {PSEUDO[symbol][1]}" for symbol in PSEUDO)
    cell_card = "\n".join(" ".join(f"{v:.10f}" for v in row) for row in np.diag(cell))
    (output / f"{key}.in").write_text(f"""&CONTROL
 calculation='scf', prefix='{key}', pseudo_dir='./pseudo',
 outdir='./scratch', tprnfor=.true., disk_io='low',
/
&SYSTEM
 ibrav=0, nat={len(symbols)}, ntyp=6, ecutwfc=60, ecutrho=480, input_dft='PBE',
 occupations='smearing', smearing='mv', degauss=0.003, nosym=.true., noinv=.true.,
/
&ELECTRONS
 conv_thr=1.0d-8, electron_maxstep=140, mixing_beta=0.35,
/
ATOMIC_SPECIES
{species}
CELL_PARAMETERS angstrom
{cell_card}
ATOMIC_POSITIONS angstrom
{rows}
K_POINTS automatic
2 2 1 0 0 0
""")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("relax_output", type=Path, nargs="?")
    parser.add_argument("--trial", action="store_true")
    parser.add_argument("--out", type=Path, default=ROOT / "structures")
    args = parser.parse_args()
    if not args.trial and not args.relax_output:
        parser.error("a converged QE relaxation log is required")
    cell, slab = read_slab(args.relax_output, args.trial)
    water_symbols, water, score_water = water_group(slab, cell)
    release_symbols, release, movers, score_release, water_shift = release_group(water_symbols, water[-1], cell)
    args.out.mkdir(parents=True, exist_ok=True)
    metadata = {"preview_on_unrelaxed_POSCAR": args.trial, "geometry_selection": "steric-only deterministic",
                "water_overlap_score": score_water, "release_overlap_score": score_release,
                "release_spectator_water_displacement_A": water_shift.tolist(),
                "labels": {"anchor_Zn": ZN+1, "water_O": WATER_O+1, "first_S": 121,
                           "incoming_S": 133, "first_OH_H": 157, "transferring_H": 158},
                "release_movers": [i+1 for i in movers], "frames": {}}
    for group, symbols, geometries, expected in (("water", water_symbols, water, 158),
                                                  ("release", release_symbols, release, 170)):
        for i, frame in enumerate(geometries, 1):
            assert len(symbols) == len(frame) == expected
            shortest, shortest_heavy = verify_contacts(symbols, frame, cell)
            key = f"{group}{i}"
            write_geometry(args.out, key, symbols, frame, cell)
            details = {"n_atoms": len(frame), "min_pair_A": shortest, "min_heavy_pair_A": shortest_heavy,
                       "S_H_A": distance(frame[132], frame[157], cell),
                       "O_H_A": distance(frame[WATER_O], frame[157], cell),
                       "water_O_nearest_Zn_A": min(distance(frame[WATER_O], pos, cell) for pos in frame[:60]),
                       "anchor_Zn_height_A": float(frame[ZN, 2] - slab[ZN, 2]),
                       "Zn_closest_lattice_O_A": min(distance(frame[ZN], frame[j], cell) for j in range(60, 120) if j != WATER_O)}
            metadata["frames"][key] = details
            print(key, {k: round(v, 3) for k, v in details.items()})
    (args.out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()

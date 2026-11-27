"""Rebuild harmonic thermochemistry from the mobile Cartesian Hessian."""

import argparse
import hashlib
import json
import math
import re

import numpy as np

from analyze_thermo import (
    GAS_CONSTANT,
    HARTREE_TO_KJ,
    KCAL_TO_KJ,
    ROOT,
    TEMPERATURE,
    final_coordinates,
    property_value,
    read_output,
)


MASSES = {"H": 1.00783, "C": 12.0, "N": 14.00307,
          "O": 15.99491, "S": 31.97207, "Zn": 63.92910}
ATOMIC_NUMBERS = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "Zn": 30}
WAVENUMBER_FACTOR = 5140.48714  # sqrt(Eh / (bohr^2 amu)) / (2 pi c), in cm^-1
WAVENUMBER_TO_KJ = 0.0119626565639
PARTIAL_JOBS = {"state2": ("7049753", 558), "state3": ("7049754", 606)}


def frequencies_from_weighted_hessian(matrix):
    eigenvalues = np.linalg.eigvalsh(matrix)
    return np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * WAVENUMBER_FACTOR


def printed_modes(text):
    return np.asarray([
        float(value)
        for row in re.findall(r"^\s*Frequency:\s*(.+)$", text, re.M)
        for value in row.split()
    ])


def active_hessian(species):
    job, electrons = PARTIAL_JOBS[species]
    folder = ROOT / f"job-{job}"
    text = (folder / f"{species}-{job}.out").read_text()
    if "Thank you very much for using Q-Chem" not in text:
        raise ValueError(f"{species}: partial Hessian did not finish")
    alpha, beta = re.search(r"There are\s+(\d+) alpha and\s+(\d+) beta electrons", text).groups()
    if int(alpha) != electrons // 2 or int(beta) != electrons // 2:
        raise ValueError(f"{species}: unexpected electronic state")

    alist = (ROOT / "partial_hessian" / f"{species}.in").read_text().split("$alist\n", 1)[1].split("$end", 1)[0]
    mobile = np.asarray(alist.split(), dtype=int) - 1
    coords = final_coordinates(text)
    reference_rows = (ROOT / f"{species}_opt.xyz").read_text().splitlines()[2:]
    reference = [row.split() for row in reference_rows]
    frozen = sorted(set(range(len(reference))) - set(mobile))
    order = [*mobile, *frozen]
    # PHESS places the selected atoms first in the actual Hessian coordinate order.
    if [a for a, _ in coords] != [reference[index][0] for index in order]:
        raise ValueError(f"{species}: Q-Chem atom order does not match the active-first input mapping")
    positions = np.asarray([xyz for _, xyz in coords])
    reference_xyz = np.asarray([reference[index][1:4] for index in order], dtype=float)
    actual_distances = np.linalg.norm(positions[:, None] - positions[None, :], axis=2)
    reference_distances = np.linalg.norm(reference_xyz[:, None] - reference_xyz[None, :], axis=2)
    if not np.allclose(actual_distances, reference_distances, rtol=0.0, atol=1e-7):
        raise ValueError(f"{species}: geometry reported in the Hessian calculation differs")

    part = text.split("Hessian Limited to Following Atoms", 1)[1]
    part = part.strip().split("\n\n", 1)[0].splitlines()[1:]
    rows = [line.split() for line in part if line.strip()]
    listed_z = [int(row[1]) for row in rows]
    listed_xyz = np.asarray([row[2:5] for row in rows], dtype=float)
    if (listed_z != [ATOMIC_NUMBERS[atom] for atom, _ in coords[:len(mobile)]]
            or not np.allclose(listed_xyz, positions[:len(mobile)], rtol=0.0, atol=6e-6)):
        raise ValueError(f"{species}: Q-Chem Hessian active atom order differs from $alist")

    cartesian = np.fromfile(folder / f"{species}-{job}.hessian.bin", dtype="<f8")
    if cartesian.size == (3 * len(coords)) ** 2:
        dim = 3 * len(coords)
        mobile_dimensions = 3 * len(mobile)
        matrix = cartesian.reshape(dim, dim)[:mobile_dimensions, :mobile_dimensions]
    elif cartesian.size == (3 * len(mobile)) ** 2:
        dim = 3 * len(mobile)
        matrix = cartesian.reshape(dim, dim)
    else:
        raise ValueError(f"{species}: unexpected Hessian length {cartesian.size}")
    if not np.isfinite(matrix).all() or np.max(np.abs(matrix - matrix.T)) > 1e-6:
        raise ValueError(f"{species}: invalid or non-symmetric Cartesian Hessian")

    mass = np.repeat([MASSES[atom] for atom, _ in coords[:len(mobile)]], 3)
    weighted = (matrix + matrix.T) / (2 * np.sqrt(mass[:, None] * mass[None, :]))
    frequencies = frequencies_from_weighted_hessian(weighted)
    vib, thermo = re.split(
        r"STANDARD THERMODYNAMIC QUANTITIES AT\s+348\.15\s+K\s+AND\s+1\.00\s+ATM",
        text, maxsplit=1,
    )
    printed = printed_modes(vib.rsplit("VIBRATIONAL ANALYSIS", 1)[1])
    if printed.size and not np.allclose(frequencies[-printed.size:], printed, atol=0.06, rtol=0.0):
        raise ValueError(f"{species}: independent frequencies disagree with Q-Chem's printed spectrum")
    print(f"{species}: {len(mobile)} mobile atoms; {dim}² saved values; active frequencies "
          f"{frequencies.min():.2f} to {frequencies.max():.2f} cm⁻¹; "
          f"nonpositive: {frequencies[frequencies <= 0].round(2)}; printed: {printed.size}")
    return {
        "energy": property_value(text, "Total energy ="),
        "frequencies": frequencies,
        "printed": printed,
        "vendor_h": KCAL_TO_KJ * property_value(thermo, "Vibrational Enthalpy"),
        "vendor_s": KCAL_TO_KJ / 1000 * property_value(thermo, "Vibrational Entropy"),
    }


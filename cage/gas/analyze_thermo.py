"""Validate archived Q-Chem output and assemble the 348.15 K reaction free energy."""

import argparse
import itertools
import math
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
TEMPERATURE = 348.15
HARTREE_TO_KJ = 2625.499639
KCAL_TO_KJ = 4.184
GAS_CONSTANT = 0.00831446261815324  # kJ mol^-1 K^-1
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
SPECIES = {
    "state2": ("7048719", 49, 558, range(5, 25), -1),
    "phsh": ("7048721", 13, 58, (), -1),
    "state3": ("7048720", 59, 606, range(4, 24), 1),
    "water": ("7048721", 3, 10, (), 1),
}


def property_value(text, label):
    match = re.search(r"^\s*" + re.escape(label) + r"\s*:?\s*(" + NUMBER + ")", text, re.M)
    if match is None:
        raise ValueError(f"Missing Q-Chem output property {label}")
    return float(match.group(1))


def final_coordinates(freq_job):
    table = freq_job.split("Standard Nuclear Orientation (Angstroms)", 1)[1]
    coordinates = []
    for line in table.splitlines():
        tokens = line.split()
        if (len(tokens) == 5 and tokens[0].isdigit()
                and re.fullmatch(r"[A-Z][a-z]?", tokens[1])):
            coordinates.append((tokens[1], tuple(map(float, tokens[2:]))))
        elif coordinates:
            break
    return coordinates


def initial_coordinates(species):
    raw = (ROOT / f"{species}.in").read_text()
    start = raw.split("$molecule\n", 1)[1].split("$end", 1)[0]
    rows = [row.split() for row in start.strip().splitlines()[1:]]
    return [(row[0], tuple(map(float, row[1:]))) for row in rows]


def read_output(species, job_id):
    _, n_atoms, electrons, frozen, sign = SPECIES[species]
    output_path = ROOT / f"job-{job_id}" / f"{species}-{job_id}.out"
    sections = re.split(r"^Running Job \d+ of \d+[^\n]*", output_path.read_text(), flags=re.M)
    if len(sections) not in (3, 4):
        raise ValueError(f"{output_path.name}: unrecognized OPT/FREQ job sequence")
    opt, freq_job = sections[-2:]
    if "OPTIMIZATION CONVERGED" not in opt or "Thank you very much for using Q-Chem" not in freq_job:
        raise ValueError(f"{output_path.name}: incomplete OPT/FREQ")

    electron_pairs = re.findall(r"There are\s+(\d+) alpha and\s+(\d+) beta electrons", freq_job)
    if electron_pairs != [(str(electrons // 2), str(electrons // 2))]:
        raise ValueError(f"{output_path.name}: unexpected electron count {electron_pairs}")

    start = initial_coordinates(species)
    end = final_coordinates(freq_job)
    if len(end) != n_atoms or [row[0] for row in end] != [row[0] for row in start]:
        raise ValueError(f"{output_path.name}: atom order or count changed")
    frozen_delta = max((
        abs(math.dist(start[i - 1][1], start[j - 1][1])
            - math.dist(end[i - 1][1], end[j - 1][1]))
        for i, j in itertools.combinations(frozen, 2)
    ), default=0.0)
    if frozen_delta > 1.0e-6:
        raise ValueError(f"{output_path.name}: frozen framework changed by {frozen_delta:.3g} Å")

    vib, thermo = re.split(
        r"STANDARD THERMODYNAMIC QUANTITIES AT\s+348\.15\s+K\s+AND\s+1\.00\s+ATM",
        freq_job, maxsplit=1,
    )
    latest_analysis = vib.rsplit("VIBRATIONAL ANALYSIS", 1)[1]
    modes = [
        float(value)
        for row in re.findall(r"^\s*Frequency:\s*(.+)$", latest_analysis, re.M)
        for value in row.split()
    ]
    if not modes:
        raise ValueError(f"{output_path.name}: missing frequencies")
    imaginary = int(re.search(r"has\s+(\d+)\s+Imaginary Frequencies", thermo).group(1))
    frequencies = [v for v in modes if v != 0.0]
    if not frozen and len(frequencies) != 3 * n_atoms - 6:
        raise ValueError(f"{output_path.name}: unexpected gas-phase mode count")

    suffix = "Vib." if frozen else "Total"
    h_anchor = property_value(thermo, "Vibrational Enthalpy" if frozen else "Total Enthalpy")
    s_anchor = property_value(thermo, "Vibrational Entropy" if frozen else "Total Entropy")
    gradients = gradient_components(opt, frozen)
    return {
        "species": species,
        "job_id": job_id,
        "sign": sign,
        "e": property_value(freq_job, "Total energy ="),
        "e_opt": property_value(opt, "Final energy is"),
        "initial_gradient": gradients[0],
        "final_gradient": gradients[-1],
        "modes": frequencies,
        "imaginary": imaginary,
        "frozen_delta": frozen_delta,
        "h_anchor": h_anchor,
        "s_anchor": s_anchor,
        "h_qrrho": property_value(thermo, f"QRRHO-{suffix} Enthalpy"),
        "s_qrrho": property_value(thermo, f"QRRHO-{suffix} Entropy"),
        "h_gas": property_value(thermo, "Total Enthalpy"),
        "s_gas": property_value(thermo, "Total Entropy"),
        "sigma": property_value(thermo, "Rotational Symmetry Number is"),
        "coords": end,
    }


def gradient_components(opt, frozen):
    if frozen:
        return [float(value) for value in re.findall(
            r"^\s*Gradient\s+(" + NUMBER + r")\s+" + NUMBER + r"\s+(?:YES|NO)\b",
            opt, re.M,
        )]
    return [
        float(value)
        for value in re.findall(r"Max(?:imum)? gradient component\s*=\s*(" + NUMBER + ")", opt, re.I)
    ]


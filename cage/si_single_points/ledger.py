"""Compare fixed-SI D3(BJ) single points and transfer archived thermal corrections."""

import csv
import math
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent
BASELINE = ROOT.parent / "si_baseline"
sys.path.insert(0, str(ROOT.parent / "gas"))

from analyze_thermo import (  # noqa: E402
    GAS_CONSTANT,
    HARTREE_TO_KJ,
    KCAL_TO_KJ,
    TEMPERATURE,
    property_value,
    read_output,
)
from validate_active_hessian import WAVENUMBER_TO_KJ, active_hessian  # noqa: E402

SPECIES = {
    "state2": (-1, 558, "s003"),
    "phsh": (-1, 58, "s005"),
    "state3": (1, 606, "s010"),
    "water": (1, 10, "s009"),
}
PUBLISHED = {30: 102.8, 75: 102.0}
NATIVE_JOB = "7079211"


def parse_single_point(path, electrons):
    output = path.read_text()
    required = ("Convergence criterion met", "Thank you very much for using Q-Chem")
    if not all(marker in output for marker in required):
        raise ValueError(f"{path}: incomplete calculation")
    pairs = re.findall(r"There are\s+(\d+) alpha and\s+(\d+) beta electrons", output)
    if pairs != [(str(electrons // 2), str(electrons // 2))]:
        raise ValueError(f"{path}: unexpected electron count {pairs}")
    basis = int(re.search(r"There are\s+\d+ shells and\s+(\d+) basis functions", output)[1])
    return output, property_value(output, "Total energy ="), basis


def electronic_rows():
    rows = []
    for name, (sign, electrons, si_id) in SPECIES.items():
        original = (BASELINE / f"{name}.in").read_bytes()
        native_input = (ROOT / "returned" / f"{name}.in").read_bytes()
        setting = b"DFT_D            D3_BJ\n"
        if native_input.count(setting) != 1 or native_input.replace(setting, b"") != original:
            raise ValueError(f"{name}: D3(BJ) is not the sole input change")
        raw, uncorrected, n_basis = parse_single_point(
            BASELINE / "job-7048629" / f"{name}-7048629.out", electrons
        )
        output, corrected, native_basis = parse_single_point(
            ROOT / "returned" / f"{name}-{NATIVE_JOB}.out", electrons
        )
        if "Doing -D3(BJ) with" not in output or "Adding empirical dispersion" in raw:
            raise ValueError(f"{name}: unexpected dispersion setup")
        dispersion = property_value(output, "-D3 energy without 3body term =")
        if native_basis != n_basis or abs(corrected - uncorrected - dispersion) > 3e-8:
            raise ValueError(f"{name}: total-energy or basis-count comparison failed")
        rows.append({
            "species": name,
            "si_xyz": f"{si_id}.xyz",
            "reaction_sign": sign,
            "electrons": electrons,
            "basis_functions": n_basis,
            "no_d3_Eh": uncorrected,
            "native_d3bj_total_Eh": corrected,
            "native_d3bj_term_Eh": dispersion,
            "weighted_d3bj_term_kJ_mol": sign * dispersion * HARTREE_TO_KJ,
        })
    return rows


def vibration(frequencies, temperature):
    quantum = WAVENUMBER_TO_KJ * np.asarray(frequencies)
    if not np.all(np.isfinite(quantum) & (quantum > 0)):
        raise ValueError("Harmonic thermochemistry requires positive finite frequencies")
    x = quantum / (GAS_CONSTANT * temperature)
    occupation = 1 / np.expm1(x)
    zero_point = float(0.5 * np.sum(quantum))
    enthalpy = zero_point + float(np.sum(quantum * occupation))
    entropy = float(GAS_CONSTANT * np.sum(x * occupation - np.log(-np.expm1(-x))))
    return zero_point, enthalpy, entropy


def thermal_rows(temperature, clusters, gas):
    rows = []
    for name, (sign, _, _) in SPECIES.items():
        if name in clusters:
            record = clusters[name]
            electronic = record["energy"]
            zpe, h, s_vib = vibration(record["frequencies"], temperature)
            s = s_vib
        else:
            record = gas[name]
            electronic = record["e"]
            if record["imaginary"] != 0:
                raise ValueError(f"{name}: archived gas-phase frequency has imaginary modes")
            zpe, h_now, s_vib = vibration(record["modes"], temperature)
            _, h_ref, s_ref = vibration(record["modes"], TEMPERATURE)
            # Both gas molecules are nonlinear: translation, pV and rotation yield 4R.
            h = record["h_anchor"] * KCAL_TO_KJ + h_now - h_ref + 4 * GAS_CONSTANT * (temperature - TEMPERATURE)
            s = record["s_anchor"] * KCAL_TO_KJ / 1000 + s_vib - s_ref
            s += 4 * GAS_CONSTANT * math.log(temperature / TEMPERATURE)
            if name == "water":
                s -= GAS_CONSTANT * math.log(2 / record["sigma"])
        rows.append({
            "T_K": temperature,
            "species": name,
            "reaction_sign": sign,
            "optimized_Eh": electronic,
            "ZPE_kJ_mol": zpe,
            "H_thermal_excluding_ZPE_kJ_mol": h - zpe,
            "H_thermal_kJ_mol": h,
            "S_vibrational_J_mol_K": s_vib * 1000,
            "S_translation_rotation_J_mol_K": (s - s_vib) * 1000,
            "S_J_mol_K": s * 1000,
            "minus_TS_kJ_mol": -temperature * s,
        })
    return rows


def reaction_sum(rows, field):
    return math.fsum(row["reaction_sign"] * row[field] for row in rows)


def write_csv(filename, rows):
    with (ROOT / filename).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    electronic = electronic_rows()
    delta_si = HARTREE_TO_KJ * reaction_sum(electronic, "no_d3_Eh")
    delta_si_d3 = HARTREE_TO_KJ * reaction_sum(electronic, "native_d3bj_total_Eh")
    d3_shift = HARTREE_TO_KJ * reaction_sum(electronic, "native_d3bj_term_Eh")
    if abs(delta_si_d3 - delta_si - d3_shift) > 0.001:
        raise ValueError("Native D3(BJ) reaction-energy accounting failed")
    clusters = {name: active_hessian(name) for name in ("state2", "state3")}
    gas = {name: read_output(name, "7048721") for name in ("phsh", "water")}
    reaction, thermal = [], []
    for degrees, paper_g in PUBLISHED.items():
        rows = thermal_rows(degrees + 273.15, clusters, gas)
        thermal.extend(rows)
        delta_opt = HARTREE_TO_KJ * reaction_sum(rows, "optimized_Eh")
        h_corr = reaction_sum(rows, "H_thermal_kJ_mol")
        minus_ts = reaction_sum(rows, "minus_TS_kJ_mol")
        g_corr = h_corr + minus_ts
        reaction.append({
            "T_C": degrees,
            "SI_no_d3_delta_E_kJ_mol": delta_si,
            "SI_native_d3bj_shift_kJ_mol": d3_shift,
            "SI_native_d3bj_delta_E_kJ_mol": delta_si_d3,
            "optimized_no_d3_delta_E_kJ_mol": delta_opt,
            "optimized_delta_ZPE_kJ_mol": reaction_sum(rows, "ZPE_kJ_mol"),
            "optimized_delta_H_thermal_excluding_ZPE_kJ_mol": reaction_sum(rows, "H_thermal_excluding_ZPE_kJ_mol"),
            "optimized_delta_H_minus_delta_E_kJ_mol": h_corr,
            "optimized_delta_S_J_mol_K": reaction_sum(rows, "S_J_mol_K"),
            "optimized_minus_T_delta_S_kJ_mol": minus_ts,
            "optimized_delta_G_thermal_kJ_mol": g_corr,
            "optimized_no_d3_delta_G_kJ_mol": delta_opt + g_corr,
            "SI_no_d3_with_optimized_thermo_delta_G_kJ_mol": delta_si + g_corr,
            "SI_d3bj_with_optimized_thermo_delta_G_kJ_mol": delta_si_d3 + g_corr,
            "published_delta_G_kJ_mol": paper_g,
            "SI_d3bj_with_optimized_thermo_minus_published_kJ_mol": delta_si_d3 + g_corr - paper_g,
            "thermal_term_needed_to_match_paper_at_SI_d3bj_kJ_mol": paper_g - delta_si_d3,
        })
    if abs(reaction[-1]["optimized_no_d3_delta_G_kJ_mol"] - 133.3566) > 0.001:
        raise ValueError("Recomputed 75 °C undispersed ΔG differs from the archive")
    write_csv("species_electronic.csv", electronic)
    write_csv("species_thermal.csv", thermal)
    write_csv("reaction_ledger.csv", reaction)
    for row in reaction:
        print(" | ".join(f"{key}={value:.5f}" for key, value in row.items()))


if __name__ == "__main__":
    main()

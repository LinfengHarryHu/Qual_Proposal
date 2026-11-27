"""Align the two state 2 endpoints and save traceable native OVITO views."""

from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
from zoneinfo import ZoneInfo

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import numpy as np
from dftd3.interface import DispersionModel, RationalDampingParam

ROOT = Path(__file__).resolve().parent
ANALYSIS = ROOT.parent / "cage"
SOURCES = {
    "gas": ANALYSIS / "gas/state2_opt.xyz",
    "pcm": ANALYSIS / "pcm/state2.xyz",
}
NUMBERS = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "Zn": 30}
FROZEN = np.arange(4, 24)
PHENYL = np.array([25, 26, 28, 30, 32, 34])
MEIM_RING = np.array([37, 38, 39, 40, 41])
GROUPS = {
    "oxide_and_hydroxyl_H": list(range(24)) + [36],
    "thiophenolate": list(range(24, 36)),
    "MeIm": list(range(37, 49)),
}
TO_KJ = 2625.49963948
ELEMENT_COLORS = {
    "Zn": (0.53, 0.60, 0.65), "O": (0.79, 0.25, 0.21),
    "S": (0.85, 0.68, 0.08), "C": (0.25, 0.28, 0.30),
    "N": (0.16, 0.39, 0.76), "H": (0.87, 0.88, 0.90),
}
RADII = {"Zn": 0.45, "O": 0.29, "S": 0.36, "C": 0.23, "N": 0.25, "H": 0.14}
GAS_COLOR = (0.87, 0.40, 0.10)
PCM_COLOR = (0.08, 0.46, 0.79)


def read_xyz(path):
    lines = path.read_text().splitlines()
    rows = [line.split() for line in lines[2:] if line.strip()]
    assert len(rows) == int(lines[0]) == 49
    return [row[0] for row in rows], np.array([row[1:] for row in rows], float)


def rmsd(first, second, ids):
    return float(np.sqrt(np.mean(np.sum((first[ids] - second[ids]) ** 2, axis=1))))


def plane_angle(first, second, ids):
    normals = [np.linalg.svd(x[ids] - x[ids].mean(axis=0))[2][-1] for x in (first, second)]
    return float(np.degrees(np.arccos(np.clip(abs(normals[0] @ normals[1]), 0, 1))))


def dispersion(elements, positions):
    model = DispersionModel(np.array([NUMBERS[e] for e in elements]), positions / 0.529177210903)
    parameters = RationalDampingParam(method="b3lyp", atm=False)
    energy = float(model.get_dispersion(parameters, grad=False)["energy"])
    pairs = model.get_pairwise_dispersion(parameters)["additive pairwise energy"]
    assert abs(pairs.sum() - energy) < 1e-10
    pieces = {}
    names = list(GROUPS)
    for i, left in enumerate(names):
        for right in names[i:]:
            value = pairs[np.ix_(GROUPS[left], GROUPS[right])].sum()
            if left != right:
                value += pairs[np.ix_(GROUPS[right], GROUPS[left])].sum()
            pieces[f"{left} / {right}"] = float(value * TO_KJ)
    assert abs(sum(pieces.values()) - energy * TO_KJ) < 1e-7
    return {"total_Eh": energy, "total_kJ_mol": energy * TO_KJ, "pair_groups_kJ_mol": pieces}


def bonds(elements, positions):
    cutoffs = {
        frozenset(("Zn", "O")): 2.30, frozenset(("Zn", "S")): 2.65,
        frozenset(("Zn", "N")): 2.50, frozenset(("O", "H")): 1.20,
        frozenset(("S", "C")): 1.95, frozenset(("C", "C")): 1.65,
        frozenset(("C", "N")): 1.65, frozenset(("C", "H")): 1.25,
    }
    return [(j, i) for i in range(49) for j in range(i)
            if np.linalg.norm(positions[i] - positions[j]) < cutoffs.get(frozenset((elements[i], elements[j])), 0)]


def make_scene(mode, elements, frames, camera):
    import ovito
    from ovito.data import DataCollection, ParticleType
    from ovito.pipeline import Pipeline, StaticSource
    from ovito.vis import BondsVis, TachyonRenderer, TextLabelOverlay, Viewport
    from PySide6.QtCore import Qt
    data = DataCollection()
    points, atom_ids, colors, radii, topology = [], [], [], [], set()
    included = ("gas", "pcm") if mode == "overlay" else (mode,)
    shared = {}
    for environment in included:
        mapping = {}
        for index, element in enumerate(elements):
            if mode == "overlay" and environment == "pcm" and index in FROZEN:
                mapping[index] = shared[index]
                continue
            mapping[index] = len(points)
            if environment == "gas" and index in FROZEN:
                shared[index] = len(points)
            points.append(frames[environment][index])
            atom_ids.append(index + 1 + (1000 if environment == "pcm" and mode == "overlay" else 0))
            color = ELEMENT_COLORS[element]
            if mode == "overlay":
                color = (0.66, 0.69, 0.71) if index in FROZEN else (GAS_COLOR if environment == "gas" else PCM_COLOR)
            colors.append(color)
            radii.append(RADII[element])
        for a, b in bonds(elements, frames[environment]):
            topology.add(tuple(sorted((mapping[a], mapping[b]))))
    particles = data.create_particles(count=len(points))
    particles.create_property("Position", data=points)
    particles.create_property("Particle Identifier", data=atom_ids)
    particles.create_property("Particle Type", data=np.ones(len(points), int)).types.append(ParticleType(id=1, name="State 2"))
    particles.create_property("Color", data=colors)
    particles.create_property("Radius", data=radii)
    bs = particles.create_bonds()
    bs.create_property("Topology", data=sorted(topology))
    bs.vis.width = 0.09
    bs.vis.coloring_mode = BondsVis.ColoringMode.ByParticle
    pipeline = Pipeline(source=StaticSource(data=data))
    pipeline.add_to_scene()
    viewport = ovito.scene.viewports.active_vp
    viewport.type = Viewport.Type.Ortho
    viewport.camera_dir = -np.array(camera["eye"])
    viewport.camera_up = camera["up"]
    viewport.camera_pos = camera["position"]
    viewport.fov = camera["fov"]
    viewport.overlays.clear()
    title = "Gas (orange) + acetonitrile (blue)" if mode == "overlay" else ("Gas phase" if mode == "gas" else "Implicit acetonitrile")
    viewport.overlays.append(TextLabelOverlay(text=title, alignment=Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop,
                                            font_size=0.035, text_color=(0.08, 0.08, 0.08)))
    ovito.scene.viewports.maximized_vp = viewport
    ovito.scene.render_settings.size = (1400, 1100)
    ovito.scene.render_settings.background_color = (1, 1, 1)
    viewport.preview_mode = True
    ovito.scene.save(str(ROOT / f"state2_{mode}.ovito"))
    viewport.render_image(filename=str(ROOT / f"state2_{mode}.png"), size=(1400, 1100),
                          background=(1, 1, 1), renderer=TachyonRenderer(ambient_occlusion=True))
    from PIL import Image, PngImagePlugin
    image_path = ROOT / f"state2_{mode}.png"
    raster = Image.open(image_path)
    metadata = PngImagePlugin.PngInfo()
    metadata.add_text("Author", "Harry Linfeng Hu")
    raster.save(image_path, pnginfo=metadata)
    pipeline.remove_from_scene()


def main():
    elements, gas = read_xyz(SOURCES["gas"])
    other_elements, pcm = read_xyz(SOURCES["pcm"])
    assert elements == other_elements
    gas_center, pcm_center = gas[FROZEN].mean(axis=0), pcm[FROZEN].mean(axis=0)
    u, _, vt = np.linalg.svd((pcm[FROZEN] - pcm_center).T @ (gas[FROZEN] - gas_center))
    rotation = u @ np.diag((1, 1, np.linalg.det(u @ vt))) @ vt
    pcm_aligned = (pcm - pcm_center) @ rotation + gas_center
    assert rmsd(gas, pcm_aligned, FROZEN) < 1e-7
    with (ROOT / "state2_pcm_aligned.xyz").open("w") as handle:
        handle.write("49\nRigid alignment to gas using frozen atoms 5-24\n")
        for element, point in zip(elements, pcm_aligned):
            handle.write(f"{element} {point[0]:.10f} {point[1]:.10f} {point[2]:.10f}\n")
    metrics = {
        "timestamp_local": datetime.now(ZoneInfo("America/Los_Angeles")).isoformat(timespec="seconds"),
        "source_files": {key: {"path": str(path.relative_to(ROOT.parent)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                         for key, path in SOURCES.items()},
        "frozen_atom_ids": list(range(5, 25)),
        "frozen_rmsd_A": rmsd(gas, pcm_aligned, FROZEN),
        "phenyl_heavy_rmsd_A": rmsd(gas, pcm_aligned, PHENYL),
        "phenyl_plane_angle_degrees": plane_angle(gas, pcm_aligned, PHENYL),
        "MeIm_ring_rmsd_A": rmsd(gas, pcm_aligned, MEIM_RING),
        "MeIm_plane_angle_degrees": plane_angle(gas, pcm_aligned, MEIM_RING),
        "geometry": {},
        "dispersion": {},
    }
    for label, positions in (("gas", gas), ("pcm", pcm_aligned)):
        ds = np.linalg.norm(positions[PHENYL, None, :] - positions[None, FROZEN, :], axis=2)
        metrics["geometry"][label] = {
            "Zn1_S25_A": float(np.linalg.norm(positions[0] - positions[24])),
            "Zn1_N38_A": float(np.linalg.norm(positions[0] - positions[37])),
            "phenyl_center_to_cage_center_A": float(np.linalg.norm(positions[PHENYL].mean(axis=0) - positions[:24].mean(axis=0))),
            "minimum_phenyl_C_to_frozen_cage_atom_A": float(ds.min()),
        }
        metrics["dispersion"][label] = dispersion(elements, positions)
    metrics["D3_pcm_minus_gas_kJ_mol"] = metrics["dispersion"]["pcm"]["total_kJ_mol"] - metrics["dispersion"]["gas"]["total_kJ_mol"]
    metrics["D3_pair_group_differences_kJ_mol"] = {
        name: metrics["dispersion"]["pcm"]["pair_groups_kJ_mol"][name] - value
        for name, value in metrics["dispersion"]["gas"]["pair_groups_kJ_mol"].items()
    }
    gas_archive = json.loads((ANALYSIS / "gas/dispersion.json").read_text())
    assert metrics["source_files"]["gas"]["sha256"] == gas_archive["species"]["state2"]["xyz_sha256"]
    assert abs(metrics["dispersion"]["gas"]["total_Eh"] - gas_archive["species"]["state2"]["d3bj_energy_hartree"]) < 1e-10
    import csv
    with (ANALYSIS / "reference/pcm_species_results.csv").open() as handle:
        native = next(row for row in csv.DictReader(handle) if row["species"] == "state2")
    assert abs(metrics["dispersion"]["pcm"]["total_Eh"] - float(native["d3bj_e_eh"])) < 1e-8
    metrics["D3_native_crosscheck"] = "Pairwise sums reproduce archived gas D3 and native Q-Chem PCM D3 within output precision."
    eye = np.array((0.25, -0.1, 1.0)); eye /= np.linalg.norm(eye)
    right = np.cross((0, 1, 0), eye); right /= np.linalg.norm(right)
    up = np.cross(eye, right)
    basis = np.column_stack((right, up))
    projection = np.vstack((gas, pcm_aligned)) @ basis
    center = (projection.min(axis=0) + projection.max(axis=0)) / 2
    focus = basis @ center
    camera = {"eye": eye.tolist(), "up": up.tolist(), "position": (focus + eye * 50).tolist(),
              "fov": float(max(np.ptp(projection[:, 1]), np.ptp(projection[:, 0]) / (1400 / 1100)) * 0.60)}
    metrics["camera"] = camera
    if "--render" in __import__("sys").argv:
        for mode in ("gas", "pcm", "overlay"):
            make_scene(mode, elements, {"gas": gas, "pcm": pcm_aligned}, camera)
    (ROOT / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps({key: metrics[key] for key in ("frozen_rmsd_A", "phenyl_plane_angle_degrees", "D3_pcm_minus_gas_kJ_mol", "D3_pair_group_differences_kJ_mol")}, indent=2))


if __name__ == "__main__":
    main()

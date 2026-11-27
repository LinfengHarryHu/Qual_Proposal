"""Reanalyze public SI plot markers and demonstrate a conditional rate ambiguity."""

import argparse
import csv
import hashlib
import json
import os
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "qual_zno_mpl"))

import pymupdf as fitz
import numpy as np
from scipy.constants import R, h, k
from scipy.stats import linregress, t



ROOT = Path(__file__).resolve().parent
SI = ROOT / "source_si.pdf"
SI_SHA256 = "ec14ab5c79462b2d9211a0f3ebf3d54489cf2493de77209b0f06b48eb9db5f51"


@dataclass(frozen=True)
class Panel:
    figure: str
    page: int
    x_ticks: tuple
    y_ticks: tuple
    marker_range: tuple
    published_slope: float
    x_quantity: str
    y_quantity: str
    temperature_c: int


# Each pair maps (PDF drawing sequence number, printed tick value).
TIME = "time_s"
PRODUCT = "released_complex_mol_per_L"
LOG_RATE = "log10_rate_mol_per_L_per_s"
PANELS = (
    Panel("S9", 6, ((54, 0), (58, 20000)), ((65, 0), (69, .08)), (76, 81),
          3.09e-6, TIME, PRODUCT, 30),
    Panel("S10", 6, ((99, 0), (102, 30000)), ((108, 0), (111, .15)), (117, 124),
          4.80e-6, TIME, PRODUCT, 40),
    Panel("S11", 7, ((52, 0), (54, 20000)), ((59, 0), (62, .15)), (68, 75),
          7.65e-6, TIME, PRODUCT, 55),
    Panel("S12", 7, ((93, 0), (95, 10000)), ((100, 0), (103, .15)), (109, 116),
          1.48e-5, TIME, PRODUCT, 75),
    Panel("S14", 8, ((108, -.7), (112, -.3)), ((119, -5.6), (122, -5.0)), (128, 131),
          -.03, "log10_MeIm_mol_per_L", LOG_RATE, 40),
    Panel("S15", 9, ((58, -.7), (62, -.3)), ((69, -5.8), (73, -5.0)), (80, 83),
          -.12, "log10_PhSH_mol_per_L", LOG_RATE, 40),
)


def marker_positions(page, panel):
    paths = {path["seqno"]: path for path in page.get_drawings()}
    x_pdf, x_data, y_pdf, y_data = [], [], [], []
    for axis, pairs, coords, labels in [
        ("x", panel.x_ticks, x_pdf, x_data), ("y", panel.y_ticks, y_pdf, y_data),
    ]:
        for seq, value in pairs:
            path = paths[seq]
            rect = path["rect"]
            if path["type"] != "s" or len(path["items"]) != 1:
                raise ValueError(f"{panel.figure}: expected a single-stroke tick")
            if (rect.width if axis == "x" else rect.height) > 0.001:
                raise ValueError(f"{panel.figure}: unexpected tick orientation")
            coords.append(rect.x0 if axis == "x" else rect.y0)
            labels.append(value)

    rows = []
    start, end = panel.marker_range
    for seq in range(start, end + 1):
        path = paths[seq]
        rect = path["rect"]
        if (path["type"] != "f" or len(path["items"]) != 4
                or not 9 < rect.width < 10 or abs(rect.width - rect.height) > .01):
            raise ValueError(f"{panel.figure}: expected a filled round data marker at {seq}")
        px, py = (rect.x0 + rect.x1) / 2, (rect.y0 + rect.y1) / 2
        x = x_data[0] + (px - x_pdf[0]) * (x_data[1] - x_data[0]) / (x_pdf[1] - x_pdf[0])
        y = y_data[0] + (py - y_pdf[0]) * (y_data[1] - y_data[0]) / (y_pdf[1] - y_pdf[0])
        rows.append({"figure": panel.figure, "pdf_page": panel.page, "pdf_path_seq": seq,
                     "pdf_x_pt": px, "pdf_y_pt": py, "temperature_c": panel.temperature_c,
                     "x_quantity": panel.x_quantity, "x": x, "y_quantity": panel.y_quantity, "y": y})
    return sorted(rows, key=lambda row: row["x"])


def summarize(panel, rows):
    observed_x, observed_y = ([row[axis] for row in rows] for axis in ("x", "y"))
    fit = linregress(observed_x, observed_y)
    half_width = t.ppf(.975, len(rows) - 2) * fit.stderr
    return {"panel": asdict(panel), "n": len(rows), "slope": float(fit.slope),
            "intercept": float(fit.intercept), "r_squared": float(fit.rvalue**2),
            "slope_standard_error": float(fit.stderr),
            "ols_only_95pct_slope_interval": [float(fit.slope-half_width), float(fit.slope+half_width)]}


def main():
    parser = argparse.ArgumentParser(description='Extract vector markers from the published SI PDF.')
    parser.add_argument('pdf', type=Path)
    parser.add_argument('--output', type=Path, default=ROOT/'marker_points.csv')
    args = parser.parse_args()
    if hashlib.sha256(args.pdf.read_bytes()).hexdigest() != SI_SHA256:
        raise ValueError('PDF checksum differs; drawing indices need to be checked again.')
    rows = []
    with fitz.open(args.pdf) as doc:
        for panel in PANELS:
            rows.extend(marker_positions(doc[panel.page-1], panel))
    with args.output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    print(f'Wrote {len(rows)} markers to {args.output}')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""testing_breakup_matplotlib.py

Matplotlib version of the breakup gamma vs d_child/d_parent plot.

This script mirrors the MATLAB breakup test as closely as practical:
- scatter of gamma vs dRatio, colored by case
- per-case binned mean trend line
- fixed x-window and log x-axis by default
- legend below the axes

Data source order:
1. breakup_analysis_by_case.mat from the latest plot data directory
2. breakup_events.csv from Testing/test_outputs as a fallback

The fallback exists so the plot can still be generated on a Python install
without SciPy, as long as the CSV export exists.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import zipfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    try:
        import seaborn as sns  # type: ignore
    except Exception:
        sns = None
except ImportError as exc:  # pragma: no cover - environment-dependent.
    raise SystemExit(
        "This script needs numpy and matplotlib.\n"
        "Install them in the Python you are using, for example:\n"
        "  py -3.13 -m pip install numpy matplotlib\n"
        "Optional extras:\n"
        "  py -3.13 -m pip install scipy seaborn\n"
        "SciPy is only needed if you want to load the MAT file instead of the\n"
        "CSV fallback; seaborn is only used for a nicer style when available."
    ) from exc

try:  # Optional dependency used when the MAT file is available.
    from scipy.io import loadmat as scipy_loadmat  # type: ignore
    from scipy.io.matlab import mat_struct  # type: ignore
except Exception:  # pragma: no cover - SciPy is optional.
    scipy_loadmat = None
    mat_struct = None


MATLAB_LINES = [
    "#0072BD",  # blue
    "#D95319",  # orange
    "#EDB120",  # yellow
    "#7E2F8E",  # purple
    "#77AC30",  # green
    "#4DBEEE",  # cyan
]

DEFAULT_MAT_DIR = r"C:\Users\kbsanjayvasanth\Downloads\plot_data_mat"
DEFAULT_CSV = Path(__file__).resolve().parent / "test_outputs" / "breakup_events.csv"
DEFAULT_OUT_DIR = Path(__file__).resolve().parent / "test_outputs" / "BreakupAnalysisMatplotlib"


@dataclass
class BreakupEvent:
    gamma: float
    dRatio: float


@dataclass
class BreakupCase:
    case_name: str
    Re: float
    kD: float
    events: List[BreakupEvent]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a matplotlib version of the breakup gamma vs dRatio plot."
    )
    parser.add_argument(
        "--mat-dir",
        default=DEFAULT_MAT_DIR,
        help="Directory containing breakup_analysis_by_case.mat",
    )
    parser.add_argument(
        "--csv",
        default=str(DEFAULT_CSV),
        help="Fallback breakup_events.csv file",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUT_DIR),
        help="Output directory for the matplotlib figure",
    )
    parser.add_argument(
        "--xscale",
        choices=("log", "linear"),
        default="log",
        help="X-axis scale for d_child/d_parent",
    )
    parser.add_argument(
        "--xmin",
        type=float,
        default=0.11,
        help="Lower x-axis limit",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=1.4,
        help="Upper x-axis limit",
    )
    parser.add_argument(
        "--marker-size",
        type=float,
        default=24.0,
        help="Scatter marker area in points^2",
    )
    parser.add_argument(
        "--marker-alpha",
        type=float,
        default=0.25,
        help="Scatter transparency",
    )
    parser.add_argument(
        "--max-trend-bins",
        type=int,
        default=12,
        help="Maximum number of binned trend segments per case",
    )
    parser.add_argument(
        "--min-bin-count",
        type=int,
        default=5,
        help="Minimum points required per bin to draw trend points",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Figure DPI",
    )
    args = parser.parse_args()

    configure_matplotlib()

    mat_dir = resolve_existing_path(args.mat_dir)
    csv_path = resolve_existing_path(args.csv)
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    cases, source_label = load_breakup_cases(mat_dir, csv_path)
    if not cases:
        raise RuntimeError("No breakup cases could be loaded.")

    total_events = sum(len(case.events) for case in cases)
    print(f"Loaded {len(cases)} cases from {source_label}")
    print(f"Total breakup events: {total_events}")
    print()
    print(f"{'Case':<12} {'k/d':>8} {'N events':>10}")
    print("-" * 32)
    for case in cases:
        print(f"{case.case_name:<12} {case.kD:>8.2f} {len(case.events):>10d}")

    fig_path = plot_breakup_gamma_vs_dratio(
        cases=cases,
        out_dir=out_dir,
        xscale=args.xscale,
        xlim=(args.xmin, args.xmax),
        marker_size=args.marker_size,
        marker_alpha=args.marker_alpha,
        max_trend_bins=args.max_trend_bins,
        min_bin_count=args.min_bin_count,
        dpi=args.dpi,
    )
    print()
    print(f"Saved matplotlib figure to: {fig_path}")
    return 0


def configure_matplotlib() -> None:
    # Prefer a cleaner seaborn-like look, but stay fully functional if the
    # seaborn package is not installed.
    if sns is not None:
        sns.set_theme(style="whitegrid", context="talk")
    else:
        try:
            plt.style.use("seaborn-v0_8-whitegrid")
        except Exception:
            plt.style.use("default")

    plt.rcParams.update(
        {
            "figure.facecolor": "#f8fafc",
            "axes.facecolor": "white",
            "axes.edgecolor": "0.25",
            "axes.linewidth": 1.4,
            "axes.labelsize": 24,
            "xtick.labelsize": 18,
            "ytick.labelsize": 18,
            "legend.fontsize": 12,
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.width": 1.4,
            "ytick.major.width": 1.4,
            "xtick.minor.width": 1.0,
            "ytick.minor.width": 1.0,
            "grid.color": "0.90",
            "grid.linewidth": 0.9,
            "grid.alpha": 1.0,
        }
    )


def build_case_palette(n_cases: int) -> List[Any]:
    if n_cases <= 0:
        return []
    if sns is not None:
        return list(sns.color_palette("colorblind", n_colors=n_cases))
    cmap = plt.get_cmap("tab10")
    if n_cases <= cmap.N:
        return [cmap(i) for i in range(n_cases)]
    return [cmap(i % cmap.N) for i in range(n_cases)]


def resolve_existing_path(path_like: str) -> Path:
    """Resolve Windows-style paths on Windows or WSL-style paths on Linux."""

    path = Path(path_like).expanduser()
    if path.exists():
        return path

    text = str(path_like)
    match = re.match(r"^([A-Za-z]):[\\/](.*)$", text)
    if match:
        drive = match.group(1).lower()
        rest = match.group(2).replace("\\", "/")
        candidate = Path(f"/mnt/{drive}/{rest}")
        if candidate.exists():
            return candidate

    return path


def load_breakup_cases(mat_dir: Path, csv_path: Path) -> Tuple[List[BreakupCase], str]:
    """Load breakup cases from the MAT file, falling back to the CSV export."""

    mat_file = mat_dir / "breakup_analysis_by_case.mat"
    if mat_file.exists():
        cases = load_breakup_cases_from_mat(mat_file)
        if cases is not None:
            return cases, str(mat_file)

    if csv_path.exists():
        cases = load_breakup_cases_from_csv(csv_path)
        return cases, str(csv_path)

    raise FileNotFoundError(
        "Could not find breakup_analysis_by_case.mat or breakup_events.csv."
    )


def load_breakup_cases_from_csv(csv_path: Path) -> List[BreakupCase]:
    """Parse the combined CSV export produced by write_breakup_analysis_csv."""

    grouped: Dict[Tuple[str, float, float], List[BreakupEvent]] = {}
    with csv_path.open("r", newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            case_name = str(row.get("Case", "")).strip() or "Case"
            re_val = _safe_float(row.get("Re"))
            k_val = _safe_float(row.get("kD"))
            gamma = _safe_float(row.get("gamma"))
            d_ratio = _safe_float(row.get("dRatio"))
            if not (math.isfinite(gamma) and math.isfinite(d_ratio)):
                continue
            key = (case_name, re_val, k_val)
            grouped.setdefault(key, []).append(BreakupEvent(gamma=gamma, dRatio=d_ratio))

    cases: List[BreakupCase] = []
    for (case_name, re_val, k_val), events in grouped.items():
        cases.append(BreakupCase(case_name=case_name, Re=re_val, kD=k_val, events=events))
    return cases


def load_breakup_cases_from_mat(mat_file: Path) -> Optional[List[BreakupCase]]:
    if scipy_loadmat is None or mat_struct is None:
        return None

    mat = scipy_loadmat(str(mat_file), squeeze_me=True, struct_as_record=False)
    if "allBreakup" not in mat:
        raise KeyError(f"'allBreakup' not found in {mat_file}")

    raw_cases = mat["allBreakup"]
    cases_py = _mat_to_python(raw_cases)
    case_list = _ensure_list(cases_py)

    cases: List[BreakupCase] = []
    for case_obj in case_list:
        if not isinstance(case_obj, dict):
            continue
        case_name = str(case_obj.get("caseName", "")).strip() or "Case"
        re_val = _safe_float(case_obj.get("Re", float("nan")))
        k_val = _safe_float(case_obj.get("kD", float("nan")))
        events_raw = _ensure_list(case_obj.get("events", []))
        events: List[BreakupEvent] = []
        for ev in events_raw:
            if not isinstance(ev, dict):
                continue
            gamma = _safe_float(ev.get("gamma"))
            d_ratio = _safe_float(ev.get("dRatio"))
            if math.isfinite(gamma) and math.isfinite(d_ratio):
                events.append(BreakupEvent(gamma=gamma, dRatio=d_ratio))
        cases.append(BreakupCase(case_name=case_name, Re=re_val, kD=k_val, events=events))

    return cases


def _mat_to_python(obj: Any) -> Any:
    """Convert scipy.io.loadmat output into plain Python containers."""

    if mat_struct is not None and isinstance(obj, mat_struct):
        return {field: _mat_to_python(getattr(obj, field)) for field in obj._fieldnames}

    if isinstance(obj, np.ndarray):
        if obj.dtype == object:
            items = [_mat_to_python(item) for item in obj.flat]
            if obj.size == 1:
                return items[0]
            return items

        if obj.dtype.kind in {"U", "S"}:
            if obj.size == 1:
                return str(obj.reshape(-1)[0])
            return "".join(str(x) for x in obj.reshape(-1))

        if obj.size == 1:
            return obj.reshape(-1)[0].item()

        return obj

    if isinstance(obj, np.generic):
        return obj.item()

    return obj


def load_breakup_cases_from_xlsx(xlsx_path: Path) -> List[BreakupCase]:
    """Parse the Excel export produced by write_breakup_analysis_xlsx."""

    with zipfile.ZipFile(str(xlsx_path), "r") as zf:
        shared_strings = _load_shared_strings(zf)
        sheet_targets = _load_sheet_targets(zf)

        cases: List[BreakupCase] = []
        for sheet_name, target in sheet_targets:
            rows = _read_sheet_rows(zf, target, shared_strings)
            if not rows:
                continue
            headers = [str(cell).strip() for cell in rows[0]]
            header_index = {name: idx for idx, name in enumerate(headers) if name}
            case_name, re_val, k_val, events = _rows_to_breakup_case(
                sheet_name=sheet_name,
                rows=rows[1:],
                header_index=header_index,
            )
            cases.append(
                BreakupCase(
                    case_name=case_name or sheet_name,
                    Re=re_val,
                    kD=k_val,
                    events=events,
                )
            )

    return cases


def _load_shared_strings(zf: zipfile.ZipFile) -> List[str]:
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []

    ns = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    strings: List[str] = []
    for si in root.findall("a:si", ns):
        text_parts = []
        for node in si.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t"):
            if node.text:
                text_parts.append(node.text)
        strings.append("".join(text_parts))
    return strings


def _load_sheet_targets(zf: zipfile.ZipFile) -> List[Tuple[str, str]]:
    ns = {
        "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    }
    workbook_root = ET.fromstring(zf.read("xl/workbook.xml"))
    rels_root = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))

    rel_map: Dict[str, str] = {}
    for rel in rels_root:
        rel_id = rel.attrib.get("Id")
        target = rel.attrib.get("Target")
        if rel_id and target:
            rel_map[rel_id] = target

    targets: List[Tuple[str, str]] = []
    sheets_node = workbook_root.find("a:sheets", ns)
    if sheets_node is None:
        return targets

    for sheet in sheets_node:
        sheet_name = sheet.attrib.get("name", "Sheet")
        rel_id = sheet.attrib.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
        target = rel_map.get(rel_id or "", "")
        if not target:
            continue
        if not target.startswith("xl/"):
            target = "xl/" + target.lstrip("/")
        targets.append((sheet_name, target))
    return targets


def _read_sheet_rows(
    zf: zipfile.ZipFile,
    target: str,
    shared_strings: Sequence[str],
) -> List[List[str]]:
    ns = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    root = ET.fromstring(zf.read(target))
    sheet_data = root.find("a:sheetData", ns)
    if sheet_data is None:
        return []

    rows: List[List[str]] = []
    for row_el in sheet_data.findall("a:row", ns):
        cells: Dict[int, str] = {}
        max_idx = -1
        for cell in row_el.findall("a:c", ns):
            ref = cell.attrib.get("r", "")
            col_idx = _excel_col_to_index(ref)
            cells[col_idx] = _read_xlsx_cell(cell, shared_strings)
            max_idx = max(max_idx, col_idx)
        if max_idx < 0:
            continue
        row = [""] * (max_idx + 1)
        for col_idx, value in cells.items():
            row[col_idx] = value
        rows.append(row)
    return rows


def _read_xlsx_cell(cell: ET.Element, shared_strings: Sequence[str]) -> str:
    cell_type = cell.attrib.get("t", "")
    if cell_type == "inlineStr":
        is_node = cell.find("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}is")
        if is_node is None:
            return ""
        text_parts = []
        for node in is_node.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t"):
            if node.text:
                text_parts.append(node.text)
        return "".join(text_parts)

    value_node = cell.find("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}v")
    if value_node is None or value_node.text is None:
        return ""

    raw = value_node.text
    if cell_type == "s":
        try:
            idx = int(raw)
            if 0 <= idx < len(shared_strings):
                return shared_strings[idx]
        except Exception:
            return raw
        return raw

    if cell_type == "b":
        return "1" if raw.strip() == "1" else "0"

    return raw


def _excel_col_to_index(cell_ref: str) -> int:
    letters = []
    for ch in cell_ref:
        if ch.isalpha():
            letters.append(ch.upper())
        else:
            break
    if not letters:
        return 0
    idx = 0
    for ch in letters:
        idx = idx * 26 + (ord(ch) - ord("A") + 1)
    return idx - 1


def _rows_to_breakup_case(
    sheet_name: str,
    rows: Sequence[Sequence[str]],
    header_index: Dict[str, int],
) -> Tuple[str, float, float, List[BreakupEvent]]:
    case_name = sheet_name
    re_val = float("nan")
    k_val = float("nan")
    events: List[BreakupEvent] = []

    idx_case = header_index.get("Case")
    idx_re = header_index.get("Re")
    idx_kd = header_index.get("kD")
    idx_gamma = header_index.get("gamma")
    idx_dratio = header_index.get("dRatio")

    for row in rows:
        if not any(str(cell).strip() for cell in row):
            continue
        if idx_case is not None and idx_case < len(row):
            maybe_case = str(row[idx_case]).strip()
            if maybe_case:
                case_name = maybe_case
        if math.isnan(re_val) and idx_re is not None and idx_re < len(row):
            re_val = _safe_float(row[idx_re])
        if math.isnan(k_val) and idx_kd is not None and idx_kd < len(row):
            k_val = _safe_float(row[idx_kd])
        if idx_gamma is None or idx_dratio is None:
            continue
        if idx_gamma >= len(row) or idx_dratio >= len(row):
            continue
        gamma = _safe_float(row[idx_gamma])
        d_ratio = _safe_float(row[idx_dratio])
        if math.isfinite(gamma) and math.isfinite(d_ratio):
            events.append(BreakupEvent(gamma=gamma, dRatio=d_ratio))

    return case_name, re_val, k_val, events


def _safe_float(value: Any) -> float:
    try:
        if value is None:
            return float("nan")
        if isinstance(value, str):
            text = value.strip()
            if not text:
                return float("nan")
            return float(text)
        return float(value)
    except Exception:
        return float("nan")


def _ensure_list(value: Any) -> List[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return [value.item()]
        return list(value.reshape(-1))
    return [value]


def plot_breakup_gamma_vs_dratio(
    cases: Sequence[BreakupCase],
    out_dir: Path,
    xscale: str,
    xlim: Tuple[float, float],
    marker_size: float,
    marker_alpha: float,
    max_trend_bins: int,
    min_bin_count: int,
    dpi: int,
) -> Path:
    x_lower, x_upper = sorted((float(xlim[0]), float(xlim[1])))
    if xscale not in ("log", "linear"):
        xscale = "log"
    if xscale == "log" and x_lower <= 0:
        raise ValueError("Log x-axis requires xlim lower bound to be positive.")

    all_x: List[float] = []
    all_y: List[float] = []
    for case in cases:
        for ev in case.events:
            if math.isfinite(ev.gamma) and math.isfinite(ev.dRatio):
                if xscale != "log" or ev.dRatio > 0:
                    all_x.append(ev.dRatio)
                    all_y.append(ev.gamma)

    if not all_x:
        raise RuntimeError("No finite breakup events found for plotting.")

    all_x_arr = np.asarray(all_x, dtype=float)
    all_y_arr = np.asarray(all_y, dtype=float)

    visible_all = (all_x_arr >= x_lower) & (all_x_arr <= x_upper)
    if not np.any(visible_all):
        visible_all = np.ones_like(all_x_arr, dtype=bool)

    y_visible = all_y_arr[visible_all]
    y_min, y_max = _robust_y_limits(y_visible)

    fig, ax = plt.subplots(figsize=(9.2, 8.8), dpi=dpi)
    fig.patch.set_facecolor("#f8fafc")
    ax.set_facecolor("white")
    ax.set_axisbelow(True)

    if xscale == "log":
        ax.set_xscale("log")

    # A light density backdrop helps with the crowded lower half of the plot
    # without overwhelming the case-specific colors.
    try:
        ax.hexbin(
            all_x_arr,
            all_y_arr,
            gridsize=(52, 36),
            extent=(x_lower, x_upper, y_min, y_max),
            bins="log",
            cmap="Greys",
            mincnt=1,
            linewidths=0.0,
            alpha=0.22,
            xscale="log" if xscale == "log" else "linear",
            zorder=0,
        )
    except Exception:
        # Hexbin is a nice-to-have; if anything goes wrong, keep the plot
        # working with the scatter + trend lines alone.
        pass

    legend_handles: List[Line2D] = []
    legend_labels: List[str] = []
    palette = build_case_palette(len(cases))

    for idx, case in enumerate(cases):
        color = palette[idx % len(palette)] if palette else MATLAB_LINES[idx % len(MATLAB_LINES)]
        x_vals: List[float] = []
        y_vals: List[float] = []
        for ev in case.events:
            if not (math.isfinite(ev.gamma) and math.isfinite(ev.dRatio)):
                continue
            if xscale == "log" and ev.dRatio <= 0:
                continue
            if ev.dRatio < x_lower or ev.dRatio > x_upper:
                continue
            x_vals.append(ev.dRatio)
            y_vals.append(ev.gamma)

        if x_vals:
            ax.scatter(
                x_vals,
                y_vals,
                s=marker_size,
                c=[color],
                alpha=marker_alpha,
                edgecolors="none",
                linewidths=0.0,
                zorder=2,
            )

            x_trend, y_trend = build_binned_mean_trend(
                np.asarray(x_vals, dtype=float),
                np.asarray(y_vals, dtype=float),
                max_trend_bins=max_trend_bins,
                min_bin_count=min_bin_count,
                xscale=xscale,
                xlim=(x_lower, x_upper),
            )
            if len(x_trend) >= 2:
                ax.plot(
                    x_trend,
                    y_trend,
                    color=color,
                    lw=2.8,
                    marker="o",
                    markersize=3.5,
                    zorder=3,
                )

            legend_handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="none",
                    markerfacecolor=color,
                    markeredgecolor="none",
                    markersize=16,
                )
            )
            legend_labels.append(f"k/d = {case.kD:.2f}")

    ax.axhline(0, color="0.55", ls="--", lw=1.8, zorder=0)
    ax.set_xlim(x_lower, x_upper)
    ax.set_ylim(y_min, y_max)

    ax.set_xlabel(r"$d_{\mathrm{child}}/d_{\mathrm{parent}}$", labelpad=18)
    ax.set_ylabel(r"$\gamma = (x_{\mathrm{child}} - x_{\mathrm{parent}})\,/\,d$", labelpad=18)
    ax.tick_params(which="major", direction="in", length=12, width=2.0, top=True, right=True)
    ax.tick_params(which="minor", direction="in", length=6, width=1.2, top=True, right=True)
    ax.minorticks_on()
    ax.grid(True, which="major", color="0.90", linewidth=0.9)
    if xscale == "log":
        ax.grid(True, which="minor", axis="x", color="0.94", linewidth=0.5)
    for spine in ax.spines.values():
        spine.set_linewidth(1.4)
        spine.set_color("0.25")

    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="lower center",
            bbox_to_anchor=(0.0, 0.0, 1.0, 0.12),
            mode="expand",
            ncol=max(1, len(legend_handles)),
            frameon=True,
            fancybox=True,
            framealpha=0.92,
            facecolor="white",
            edgecolor="0.85",
            handlelength=1.0,
            handletextpad=0.45,
            columnspacing=1.2,
            borderaxespad=0.0,
        )

    # Reserve space for the single-row legend while keeping a square canvas.
    fig.subplots_adjust(left=0.13, right=0.98, top=0.98, bottom=0.22)

    if sns is not None:
        sns.despine(ax=ax, trim=True)

    out_dir.mkdir(parents=True, exist_ok=True)
    base_name = f"Breakup_gamma_vs_dRatio_python_nice_{xscale}"
    out_path = out_dir / f"{base_name}.png"
    fig.savefig(out_path, dpi=dpi, facecolor="white")
    plt.close(fig)
    return out_path


def build_binned_mean_trend(
    x_vals: np.ndarray,
    y_vals: np.ndarray,
    max_trend_bins: int,
    min_bin_count: int,
    xscale: str,
    xlim: Tuple[float, float],
) -> Tuple[np.ndarray, np.ndarray]:
    valid = np.isfinite(x_vals) & np.isfinite(y_vals)
    x_vals = x_vals[valid]
    y_vals = y_vals[valid]

    if xscale == "log":
        positive = x_vals > 0
        x_vals = x_vals[positive]
        y_vals = y_vals[positive]

    if x_vals.size < 2 * min_bin_count:
        return np.array([]), np.array([])

    n_bins = min(max_trend_bins, max(1, int(math.floor(x_vals.size / float(min_bin_count)))))
    if n_bins < 1:
        return np.array([]), np.array([])

    x0 = max(float(xlim[0]), float(np.min(x_vals)))
    x1 = min(float(xlim[1]), float(np.max(x_vals)))
    if not math.isfinite(x0) or not math.isfinite(x1) or x0 >= x1:
        return np.array([]), np.array([])

    if xscale == "log":
        edges = np.logspace(math.log10(x0), math.log10(x1), n_bins + 1)
    else:
        edges = np.linspace(x0, x1, n_bins + 1)

    edges = np.unique(edges)
    if edges.size < 2:
        return np.array([]), np.array([])

    x_trend: List[float] = []
    y_trend: List[float] = []
    for bin_idx in range(edges.size - 1):
        left = edges[bin_idx]
        right = edges[bin_idx + 1]
        if bin_idx < edges.size - 2:
            mask = (x_vals >= left) & (x_vals < right)
        else:
            mask = (x_vals >= left) & (x_vals <= right)
        if int(np.sum(mask)) < min_bin_count:
            continue
        x_trend.append(float(np.mean(x_vals[mask])))
        y_trend.append(float(np.mean(y_vals[mask])))

    return np.asarray(x_trend), np.asarray(y_trend)


def _robust_y_limits(y_vals: np.ndarray) -> Tuple[float, float]:
    finite = y_vals[np.isfinite(y_vals)]
    if finite.size == 0:
        return -0.1, 0.1

    lo, hi = np.percentile(finite, [1, 99])
    lo = min(float(lo), 0.0)
    hi = max(float(hi), 0.0)
    if not math.isfinite(lo) or not math.isfinite(hi) or lo >= hi:
        span = float(np.max(np.abs(finite)))
        if not math.isfinite(span) or span <= 0:
            span = 1.0
        lo, hi = -0.1 * span, 0.1 * span

    pad = 0.08 * (hi - lo)
    return lo - pad, hi + pad


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarize stage-1 neutral calibration time-series outputs."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional


NUMERIC_FIELDS = [
    "N",
    "Nc_over_K",
    "pi_neutral",
    "Ne_pi",
    "Ne_pi_over_K",
    "Ne_pi_over_Nc",
    "density_factor",
    "births_total",
    "OMU_count",
    "AMU_males",
    "dominant_males",
    "mean_age",
]


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return [dict(row) for row in csv.DictReader(f, delimiter="\t")]


def to_float(value: object) -> Optional[float]:
    try:
        if value is None:
            return None
        text = str(value).strip()
        if text == "" or text.upper() in {"NA", "NAN", "NULL"}:
            return None
        x = float(text)
        if math.isnan(x):
            return None
        return x
    except Exception:
        return None


def mean(values: Iterable[float]) -> str:
    vals = list(values)
    if not vals:
        return "NA"
    return f"{statistics.fmean(vals):.12g}"


def sd(values: Iterable[float]) -> str:
    vals = list(values)
    if len(vals) <= 1:
        return "0"
    return f"{statistics.stdev(vals):.12g}"


def read_manifest(path: Path) -> Dict[Path, Dict[str, str]]:
    if not path.exists():
        return {}
    rows = read_tsv(path)
    out: Dict[Path, Dict[str, str]] = {}
    for row in rows:
        p = Path(row["neutral_calibration_out"])
        out[p] = row
    return out


def summarize_one(path: Path, meta: Dict[str, str], stable_points: int) -> Dict[str, str]:
    rows = read_tsv(path)
    if not rows:
        raise ValueError(f"空文件: {path}")
    stable = rows[-stable_points:] if len(rows) > stable_points else rows

    first = rows[0]
    k = meta.get("K") or first.get("K", "NA")
    replicate = meta.get("replicate") or "NA"

    out: Dict[str, str] = {
        "K": k,
        "replicate": replicate,
        "file": str(path),
        "records_total": str(len(rows)),
        "records_stable": str(len(stable)),
        "cycle_min": stable[0].get("cycle", "NA"),
        "cycle_max": stable[-1].get("cycle", "NA"),
    }

    for field in NUMERIC_FIELDS:
        vals = [x for x in (to_float(row.get(field)) for row in stable) if x is not None]
        out[field + "_mean"] = mean(vals)
        out[field + "_sd"] = sd(vals)
    return out


def write_rows(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="汇总阶段1 neutral calibration 输出。")
    parser.add_argument("--outdir", default="outputs/neutral_calibration")
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--stable-points", type=int, default=20, help="每个replicate取最后多少个输出点作为稳定窗口。")
    parser.add_argument("--by-rep-out", default=None)
    parser.add_argument("--summary-out", default=None)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    manifest_path = Path(args.manifest) if args.manifest else outdir / "neutral_calibration_manifest.tsv"
    manifest = read_manifest(manifest_path)

    files = sorted(manifest.keys()) if manifest else sorted(outdir.glob("*.neutral_calibration.tsv"))
    if not files:
        raise SystemExit(f"没有找到 neutral calibration 输出: {outdir}")

    by_rep: List[Dict[str, str]] = []
    for path in files:
        if not path.exists():
            continue
        by_rep.append(summarize_one(path, manifest.get(path, {}), args.stable_points))

    if not by_rep:
        raise SystemExit("manifest存在，但没有任何可读取的neutral calibration输出。")

    by_rep_fields = [
        "K", "replicate", "file", "records_total", "records_stable", "cycle_min", "cycle_max",
    ]
    for field in NUMERIC_FIELDS:
        by_rep_fields.extend([field + "_mean", field + "_sd"])

    by_rep_out = Path(args.by_rep_out) if args.by_rep_out else outdir / "neutral_calibration_by_rep.tsv"
    write_rows(by_rep_out, by_rep, by_rep_fields)

    grouped: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in by_rep:
        grouped[row["K"]].append(row)

    summary_rows: List[Dict[str, str]] = []
    for k, rows in sorted(grouped.items(), key=lambda item: float(item[0])):
        srow: Dict[str, str] = {"K": k, "n_replicates": str(len(rows))}
        for field in [
            "N_mean",
            "Nc_over_K_mean",
            "pi_neutral_mean",
            "Ne_pi_mean",
            "Ne_pi_over_K_mean",
            "Ne_pi_over_Nc_mean",
        ]:
            vals = [x for x in (to_float(row.get(field)) for row in rows) if x is not None]
            srow[field.replace("_mean", "") + "_mean"] = mean(vals)
            srow[field.replace("_mean", "") + "_sd"] = sd(vals)
        summary_rows.append(srow)

    summary_fields = [
        "K", "n_replicates",
        "N_mean", "N_sd",
        "Nc_over_K_mean", "Nc_over_K_sd",
        "pi_neutral_mean", "pi_neutral_sd",
        "Ne_pi_mean", "Ne_pi_sd",
        "Ne_pi_over_K_mean", "Ne_pi_over_K_sd",
        "Ne_pi_over_Nc_mean", "Ne_pi_over_Nc_sd",
    ]
    summary_out = Path(args.summary_out) if args.summary_out else outdir / "neutral_calibration_summary.tsv"
    write_rows(summary_out, summary_rows, summary_fields)

    print(f"已写出 replicate 汇总: {by_rep_out}")
    print(f"已写出 K 汇总: {summary_out}")


if __name__ == "__main__":
    main()

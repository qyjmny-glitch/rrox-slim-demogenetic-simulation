#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Estimate K5200 from observed ancient Ne using stage-1 neutral calibration."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return [dict(row) for row in csv.DictReader(f, delimiter="\t")]


def nearest_tested_k(test_ks: List[float], estimate_mid: float) -> float:
    return min(test_ks, key=lambda k: abs(k - estimate_mid))


def main() -> None:
    parser = argparse.ArgumentParser(description="根据Ne_pi/K反推K5200候选范围。")
    parser.add_argument("--summary", default="outputs/neutral_calibration/neutral_calibration_summary.tsv")
    parser.add_argument("--observed-ne", type=float, default=None, help="古代观测Ne点估计。")
    parser.add_argument("--observed-ne-min", type=float, default=None, help="古代观测Ne区间下限。")
    parser.add_argument("--observed-ne-max", type=float, default=None, help="古代观测Ne区间上限。")
    parser.add_argument("--out", default="outputs/neutral_calibration/K5200_estimate_from_observed_Ne.tsv")
    args = parser.parse_args()

    if args.observed_ne is not None:
        ne_min = args.observed_ne
        ne_max = args.observed_ne
    else:
        if args.observed_ne_min is None or args.observed_ne_max is None:
            raise SystemExit("请提供 --observed-ne 或同时提供 --observed-ne-min 和 --observed-ne-max")
        ne_min = args.observed_ne_min
        ne_max = args.observed_ne_max

    rows = read_tsv(Path(args.summary))
    if not rows:
        raise SystemExit(f"summary为空: {args.summary}")

    tested_ks = [float(row["K"]) for row in rows]
    out_rows: List[Dict[str, str]] = []
    for row in rows:
        ratio = float(row["Ne_pi_over_K_mean"])
        if not math.isfinite(ratio) or ratio <= 0:
            continue
        k_min = ne_min / ratio
        k_max = ne_max / ratio
        k_mid = (k_min + k_max) / 2.0
        out_rows.append({
            "ratio_source_K": row["K"],
            "Ne_obs_min": f"{ne_min:.12g}",
            "Ne_obs_max": f"{ne_max:.12g}",
            "Ne_pi_over_K_mean": f"{ratio:.12g}",
            "K5200_min": f"{k_min:.12g}",
            "K5200_max": f"{k_max:.12g}",
            "nearest_tested_K": f"{nearest_tested_k(tested_ks, k_mid):.12g}",
        })

    if not out_rows:
        raise SystemExit("没有可用的Ne_pi_over_K_mean正值。")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "ratio_source_K",
        "Ne_obs_min",
        "Ne_obs_max",
        "Ne_pi_over_K_mean",
        "K5200_min",
        "K5200_max",
        "nearest_tested_K",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)
    print(f"已写出K5200估计结果: {out_path}")


if __name__ == "__main__":
    main()

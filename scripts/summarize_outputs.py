#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarize_outputs.py

汇总SLiM forward输出。

输入：
  outputs/forward/*.forward_demography.tsv
  outputs/forward/*.forward_genetic.tsv
  outputs/forward/*.forward_events.tsv

输出：
  summary/replicate_extinction_summary.tsv
  summary/extinction_probability_by_scenario.tsv
  summary/extinction_time_distribution.tsv
  summary/final_genetic_summary.tsv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


def infer_replicate_from_name(path: Path) -> Optional[int]:
    m = re.search(r"rep(\d+)", path.name)
    if not m:
        return None
    return int(m.group(1))


def read_tsv_if_exists(path: Path) -> Optional[pd.DataFrame]:
    if not path.exists() or path.stat().st_size == 0:
        return None
    return pd.read_csv(path, sep="\t")


def summarize_demography_file(path: Path) -> Dict[str, object]:
    df = pd.read_csv(path, sep="\t")

    if df.empty:
        raise ValueError(f"demography文件为空: {path}")

    df = df.sort_values(["cycle", "bp"], ascending=[True, False]).reset_index(drop=True)
    first = df.iloc[0]
    last = df.iloc[-1]

    replicate = infer_replicate_from_name(path)
    scenario = str(first.get("scenario", "NA"))
    genetic_mode = str(first.get("genetic_mode", "NA"))

    demographic_rows = df[df["N"] == 0]
    demographic_extinct = len(demographic_rows) > 0
    demographic_extinction_bp = int(demographic_rows.iloc[0]["bp"]) if demographic_extinct else None

    functional_extinct = False
    functional_extinction_bp = None
    if "functional_extinct" in df.columns:
        rows = df[df["functional_extinct"].astype(int) == 1]
        functional_extinct = len(rows) > 0
        functional_extinction_bp = int(rows.iloc[0]["bp"]) if functional_extinct else None

    quasi_extinct = False
    quasi_extinction_bp = None
    if "quasi_extinct" in df.columns:
        rows = df[df["quasi_extinct"].astype(int) == 1]
        quasi_extinct = len(rows) > 0
        quasi_extinction_bp = int(rows.iloc[0]["bp"]) if quasi_extinct else None

    return {
        "file": str(path),
        "replicate": replicate,
        "scenario": scenario,
        "genetic_mode": genetic_mode,
        "demographic_extinct": int(demographic_extinct),
        "demographic_extinction_bp": demographic_extinction_bp,
        "functional_extinct": int(functional_extinct),
        "functional_extinction_bp": functional_extinction_bp,
        "quasi_extinct": int(quasi_extinct),
        "quasi_extinction_bp": quasi_extinction_bp,
        "start_N": int(first["N"]),
        "final_bp": int(last["bp"]),
        "final_N": int(last["N"]),
        "min_N": int(df["N"].min()),
        "max_N": int(df["N"].max()),
        "mean_N": float(df["N"].mean()),
    }


def summarize_genetic_file(path: Path) -> Dict[str, object]:
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        raise ValueError(f"genetic文件为空: {path}")

    df = df.sort_values(["cycle", "bp"], ascending=[True, False]).reset_index(drop=True)
    last_nonzero = df[df["N"] > 0]
    row = last_nonzero.iloc[-1] if len(last_nonzero) else df.iloc[-1]

    replicate = infer_replicate_from_name(path)
    scenario = str(row.get("scenario", "NA"))
    genetic_mode = str(row.get("genetic_mode", "NA"))

    keep = {
        "file": str(path),
        "replicate": replicate,
        "scenario": scenario,
        "genetic_mode": genetic_mode,
        "last_genetic_bp": int(row["bp"]),
        "last_genetic_N": int(row["N"]),
    }

    for col in [
        "sample_n",
        "mean_genetic_fitness",
        "mean_2B",
        "mean_het_del",
        "mean_hom_del",
        "mean_weak_het",
        "mean_moderate_het",
        "mean_strong_het",
        "mean_vstrong_het",
        "segregating_del_muts",
        "fixed_del_muts",
    ]:
        if col in row.index:
            keep[f"last_{col}"] = row[col]

    return keep


def add_time_bin(bp: Optional[float]) -> str:
    if pd.isna(bp):
        return "not_extinct"
    bp = float(bp)
    if bp >= 4000:
        return "5200-4000_BP"
    if bp >= 2500:
        return "3999-2500_BP"
    if bp >= 1000:
        return "2499-1000_BP"
    if bp >= 300:
        return "999-300_BP"
    return "299-0_BP"


def main() -> None:
    parser = argparse.ArgumentParser(description="汇总SLiM forward输出。")
    parser.add_argument("--forward-dir", default="outputs/forward")
    parser.add_argument("--outdir", default="outputs/summary")
    args = parser.parse_args()

    forward_dir = Path(args.forward_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    demo_files = sorted(forward_dir.glob("*.forward_demography.tsv"))
    gen_files = sorted(forward_dir.glob("*.forward_genetic.tsv"))

    if not demo_files:
        raise FileNotFoundError(f"在{forward_dir}中未找到*.forward_demography.tsv文件")

    demo_summary = pd.DataFrame([summarize_demography_file(p) for p in demo_files])
    demo_summary["demographic_extinction_bin"] = demo_summary["demographic_extinction_bp"].apply(add_time_bin)
    demo_summary["functional_extinction_bin"] = demo_summary["functional_extinction_bp"].apply(add_time_bin)
    demo_summary["quasi_extinction_bin"] = demo_summary["quasi_extinction_bp"].apply(add_time_bin)

    demo_summary_path = outdir / "replicate_extinction_summary.tsv"
    demo_summary.to_csv(demo_summary_path, sep="\t", index=False)

    group_cols = ["scenario", "genetic_mode"]

    probability = (
        demo_summary
        .groupby(group_cols, dropna=False)
        .agg(
            n_replicates=("file", "count"),
            demographic_extinction_probability=("demographic_extinct", "mean"),
            functional_extinction_probability=("functional_extinct", "mean"),
            quasi_extinction_probability=("quasi_extinct", "mean"),
            median_demographic_extinction_bp=("demographic_extinction_bp", "median"),
            median_functional_extinction_bp=("functional_extinction_bp", "median"),
            median_quasi_extinction_bp=("quasi_extinction_bp", "median"),
            mean_final_N=("final_N", "mean"),
            min_final_N=("final_N", "min"),
        )
        .reset_index()
    )

    probability_path = outdir / "extinction_probability_by_scenario.tsv"
    probability.to_csv(probability_path, sep="\t", index=False, float_format="%.6g")

    time_dist = (
        demo_summary
        .groupby(group_cols + ["demographic_extinction_bin"], dropna=False)
        .size()
        .reset_index(name="n")
    )
    time_dist_path = outdir / "extinction_time_distribution.tsv"
    time_dist.to_csv(time_dist_path, sep="\t", index=False)

    if gen_files:
        gen_summary = pd.DataFrame([summarize_genetic_file(p) for p in gen_files])
        gen_summary_path = outdir / "final_genetic_summary.tsv"
        gen_summary.to_csv(gen_summary_path, sep="\t", index=False)
    else:
        gen_summary_path = None

    print("已写出:")
    print(demo_summary_path)
    print(probability_path)
    print(time_dist_path)
    if gen_summary_path:
        print(gen_summary_path)


if __name__ == "__main__":
    main()

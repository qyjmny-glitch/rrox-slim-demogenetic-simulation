#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarize_outputs.py

Summarize SLiM forward outputs for the Rhinopithecus roxellana model.

This version is compatible with the updated forward_demography.tsv format:
  - demography files may contain a `stage` column;
  - when `stage` exists, only `stage == post_mortality` rows are used for
    extinction probability, final N, mean N, and time-bin summaries;
  - paternity QA columns are summarized when present:
        births_total, births_own_OMU, births_other_OMU, births_AMU.

Inputs:
  outputs/forward/*.forward_demography.tsv
  outputs/forward/*.forward_genetic.tsv
  outputs/forward/*.forward_events.tsv

Outputs:
  summary/replicate_extinction_summary.tsv
  summary/extinction_probability_by_scenario.tsv
  summary/extinction_time_distribution.tsv
  summary/final_genetic_summary.tsv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Optional

import pandas as pd


PATERNITY_COLS = [
    "births_total",
    "births_own_OMU",
    "births_other_OMU",
    "births_AMU",
]


def infer_replicate_from_name(path: Path) -> Optional[int]:
    m = re.search(r"rep(\d+)", path.name)
    if not m:
        return None
    return int(m.group(1))


def read_tsv_if_exists(path: Path) -> Optional[pd.DataFrame]:
    if not path.exists() or path.stat().st_size == 0:
        return None
    return pd.read_csv(path, sep="\t")


def select_demography_stage(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """
    Use post_mortality rows for final demographic summaries.

    The updated forward script can output both pre_breeding and post_mortality rows
    for the same cycle/bp. Extinction summaries must use the annual final state,
    i.e. post_mortality. Older files without stage are kept for backward compatibility.
    """
    if "stage" not in df.columns:
        out = df.copy()
        out["stage_used_for_summary"] = "no_stage_column"
        return out

    post = df[df["stage"].astype(str) == "post_mortality"].copy()
    if post.empty:
        raise ValueError(
            f"{path}: has a stage column, but no rows with stage == 'post_mortality'."
        )

    post["stage_used_for_summary"] = "post_mortality"
    return post


def safe_int(value):
    if pd.isna(value):
        return None
    return int(value)


def safe_sum(df: pd.DataFrame, col: str):
    if col not in df.columns:
        return None
    return float(pd.to_numeric(df[col], errors="coerce").fillna(0).sum())


def safe_fraction(num, den):
    if num is None or den is None or den == 0:
        return None
    return num / den


def summarize_demography_file(path: Path) -> Dict[str, object]:
    raw = pd.read_csv(path, sep="\t")
    if raw.empty:
        raise ValueError(f"demography file is empty: {path}")

    raw = raw.sort_values(["cycle", "bp"], ascending=[True, False]).reset_index(drop=True)
    df = select_demography_stage(raw, path)
    df = df.sort_values(["cycle", "bp"], ascending=[True, False]).reset_index(drop=True)

    first = df.iloc[0]
    last = df.iloc[-1]

    replicate = infer_replicate_from_name(path)
    scenario = str(first.get("scenario", "NA"))
    genetic_mode = str(first.get("genetic_mode", "NA"))

    demographic_rows = df[df["N"] == 0]
    demographic_extinct = len(demographic_rows) > 0
    demographic_extinction_bp = safe_int(demographic_rows.iloc[0]["bp"]) if demographic_extinct else None

    functional_extinct = False
    functional_extinction_bp = None
    if "functional_extinct" in df.columns:
        rows = df[pd.to_numeric(df["functional_extinct"], errors="coerce").fillna(0).astype(int) == 1]
        functional_extinct = len(rows) > 0
        functional_extinction_bp = safe_int(rows.iloc[0]["bp"]) if functional_extinct else None

    quasi_extinct = False
    quasi_extinction_bp = None
    if "quasi_extinct" in df.columns:
        rows = df[pd.to_numeric(df["quasi_extinct"], errors="coerce").fillna(0).astype(int) == 1]
        quasi_extinct = len(rows) > 0
        quasi_extinction_bp = safe_int(rows.iloc[0]["bp"]) if quasi_extinct else None

    births_total = safe_sum(df, "births_total")
    births_own = safe_sum(df, "births_own_OMU")
    births_other = safe_sum(df, "births_other_OMU")
    births_amu = safe_sum(df, "births_AMU")

    return {
        "file": str(path),
        "replicate": replicate,
        "scenario": scenario,
        "genetic_mode": genetic_mode,
        "stage_used_for_summary": str(first.get("stage_used_for_summary", "NA")),
        "n_raw_rows": len(raw),
        "n_summary_rows": len(df),
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
        "births_total_sum": births_total,
        "births_own_OMU_sum": births_own,
        "births_other_OMU_sum": births_other,
        "births_AMU_sum": births_amu,
        "births_own_OMU_fraction": safe_fraction(births_own, births_total),
        "births_other_OMU_fraction": safe_fraction(births_other, births_total),
        "births_AMU_fraction": safe_fraction(births_amu, births_total),
    }


def summarize_genetic_file(path: Path) -> Dict[str, object]:
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        raise ValueError(f"genetic file is empty: {path}")

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
    parser = argparse.ArgumentParser(description="Summarize SLiM forward outputs.")
    parser.add_argument("--forward-dir", default="outputs/forward")
    parser.add_argument("--outdir", default="outputs/summary")
    args = parser.parse_args()

    forward_dir = Path(args.forward_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    demo_files = sorted(forward_dir.glob("*.forward_demography.tsv"))
    gen_files = sorted(forward_dir.glob("*.forward_genetic.tsv"))

    if not demo_files:
        raise FileNotFoundError(f"No *.forward_demography.tsv files found in {forward_dir}")

    demo_summary = pd.DataFrame([summarize_demography_file(p) for p in demo_files])
    demo_summary["demographic_extinction_bin"] = demo_summary["demographic_extinction_bp"].apply(add_time_bin)
    demo_summary["functional_extinction_bin"] = demo_summary["functional_extinction_bp"].apply(add_time_bin)
    demo_summary["quasi_extinction_bin"] = demo_summary["quasi_extinction_bp"].apply(add_time_bin)

    demo_summary_path = outdir / "replicate_extinction_summary.tsv"
    demo_summary.to_csv(demo_summary_path, sep="\t", index=False)

    group_cols = ["scenario", "genetic_mode"]

    agg_dict = dict(
        n_replicates=("file", "count"),
        demographic_extinction_probability=("demographic_extinct", "mean"),
        functional_extinction_probability=("functional_extinct", "mean"),
        quasi_extinction_probability=("quasi_extinct", "mean"),
        median_demographic_extinction_bp=("demographic_extinction_bp", "median"),
        median_functional_extinction_bp=("functional_extinction_bp", "median"),
        median_quasi_extinction_bp=("quasi_extinction_bp", "median"),
        mean_final_N=("final_N", "mean"),
        min_final_N=("final_N", "min"),
        total_births=("births_total_sum", "sum"),
        total_births_own_OMU=("births_own_OMU_sum", "sum"),
        total_births_other_OMU=("births_other_OMU_sum", "sum"),
        total_births_AMU=("births_AMU_sum", "sum"),
    )

    probability = demo_summary.groupby(group_cols, dropna=False).agg(**agg_dict).reset_index()

    probability["births_own_OMU_fraction"] = probability["total_births_own_OMU"] / probability["total_births"].replace({0: pd.NA})
    probability["births_other_OMU_fraction"] = probability["total_births_other_OMU"] / probability["total_births"].replace({0: pd.NA})
    probability["births_AMU_fraction"] = probability["total_births_AMU"] / probability["total_births"].replace({0: pd.NA})

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

    print("Wrote:")
    print(demo_summary_path)
    print(probability_path)
    print(time_dist_path)
    if gen_summary_path:
        print(gen_summary_path)


if __name__ == "__main__":
    main()

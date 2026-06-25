#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate or run stage-1 neutral calibration commands.

This runner keeps rrox_burnin.slim and rrox_forward.slim untouched. It drives
slim/rrox_neutral_calibration.slim with pure-neutral settings to estimate
pi_neutral, Ne_pi, Ne_pi/K and Nc/K for candidate K5200 values.
"""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from run_burnin_grid import (
    build_slim_command,
    load_all_params,
    load_param_descriptions,
    parse_scalar_value,
    read_param_descriptions_psv,
    slim_literal,
    write_command_script,
)


STAGE1_DESCRIPTIONS = {
    "NEUTRAL_CALIBRATION_OUT": "阶段1中性校准输出表，记录K、N、pi_neutral、Ne_pi、Ne_pi/K和Nc/K。",
    "RUN_MODE": "阶段1runner固定生成：neutral_calibration。",
    "GENETIC_MODE": "阶段1固定为neutral；所有新突变均为中性突变。",
    "K5200": "候选5200 BP环境承载量，用于估计Ne_pi/K。",
    "INITIAL_N": "阶段1初始普查数量；默认等于候选K。",
    "NEUTRAL_FRACTION": "阶段1固定为1.0，只输入中性突变。",
    "DELETERIOUS_FRACTION": "阶段1固定为0.0，不输入有害突变。",
}


def read_stage1_params(path: Path) -> Dict[str, object]:
    params: Dict[str, object] = {}
    with path.open("r", encoding="utf-8-sig") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = [x.strip() for x in line.split("|")]
            if fields[0] == "parameter":
                continue
            if len(fields) < 5:
                raise ValueError(f"{path}: 参数行格式错误: {line}")
            name, value, typ, required, _desc = fields[:5]
            if required.lower() == "yes" and value == "":
                raise ValueError(f"{path}: 必填参数为空: {name}")
            if value != "":
                params[name] = parse_scalar_value(value, typ)
    return params


def parse_k_values(value: object) -> List[float]:
    out: List[float] = []
    for token in str(value).split(","):
        token = token.strip()
        if token:
            out.append(float(token))
    if not out:
        raise ValueError("K_LIST不能为空")
    return out


def format_k_label(k: float) -> str:
    if abs(k - round(k)) < 1e-9:
        return str(int(round(k)))
    return str(k).replace(".", "p")


def write_manifest(path: Path, rows: Iterable[Dict[str, object]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "K",
        "replicate",
        "seed",
        "out_prefix",
        "neutral_calibration_out",
        "burnin_trees_out",
        "burnin_social_state_out",
        "command",
        "status",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="生成或运行阶段1 neutral calibration 命令。")
    parser.add_argument("--config-dir", default="config")
    parser.add_argument("--stage1-config", default="config/stage1_neutral_calibration_params.psv")
    parser.add_argument("--slim-bin", default=None)
    parser.add_argument("--slim-script", default=None)
    parser.add_argument("--k-values", default=None, help="覆盖K_LIST，例如500,1000,2000")
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--replicate-start", type=int, default=None)
    parser.add_argument("--seed-base", type=int, default=None)
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--write-sh", default=None)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    config_dir = Path(args.config_dir)
    stage1_path = Path(args.stage1_config)
    stage1 = read_stage1_params(stage1_path)

    base_params = load_all_params(config_dir)
    descriptions = load_param_descriptions(config_dir)
    descriptions.update(read_param_descriptions_psv(stage1_path))
    descriptions.update(STAGE1_DESCRIPTIONS)

    slim_bin = args.slim_bin or str(stage1.get("SLIM_BIN", "slim"))
    slim_script = Path(args.slim_script or str(stage1.get("SLIM_SCRIPT", "slim/rrox_neutral_calibration.slim")))
    outdir = Path(args.outdir or str(stage1.get("OUTDIR", "outputs/neutral_calibration")))
    write_sh = Path(args.write_sh or str(stage1.get("WRITE_SH", "outputs/logs/run_stage1_neutral_calibration.sh")))
    k_values = parse_k_values(args.k_values if args.k_values is not None else stage1["K_LIST"])
    replicates = args.replicates if args.replicates is not None else int(stage1.get("REPLICATES", 1))
    replicate_start = args.replicate_start if args.replicate_start is not None else int(stage1.get("REPLICATE_START", 1))
    seed_base = args.seed_base if args.seed_base is not None else int(stage1.get("SEED_BASE", 700000))

    outdir.mkdir(parents=True, exist_ok=True)
    command_specs: List[Tuple[List[str], Dict[str, object], str]] = []
    manifest_rows: List[Dict[str, object]] = []

    for rep in range(replicate_start, replicate_start + replicates):
        for k in k_values:
            k_label = format_k_label(k)
            label = f"K{k_label}_rep{rep:03d}"
            params = dict(base_params)
            params["RUN_MODE"] = "neutral_calibration"
            params["GENETIC_MODE"] = "neutral"
            params["REPLICATE_ID"] = rep
            params["RANDOM_SEED"] = seed_base + rep + int(round(k)) * 1000
            params["K5200"] = k

            initial_n_mode = str(stage1.get("INITIAL_N_MODE", "equal_to_K"))
            if initial_n_mode == "equal_to_K":
                params["INITIAL_N"] = int(round(k))
            elif initial_n_mode == "from_config":
                params["INITIAL_N"] = int(params["INITIAL_N"])
            else:
                raise ValueError(f"未知INITIAL_N_MODE: {initial_n_mode}")

            params["BURNIN_YEARS"] = int(stage1["BURNIN_YEARS"])
            params["BURNIN_LOG_INTERVAL"] = int(stage1["BURNIN_LOG_INTERVAL"])
            params["SAMPLE_SIZE"] = int(stage1["SAMPLE_SIZE"])
            params["MU_RATE"] = float(stage1["MU_RATE"])
            params["SAVE_TREES"] = int(stage1["SAVE_TREES"])
            params["NEUTRAL_FRACTION"] = 1.0
            params["DELETERIOUS_FRACTION"] = 0.0

            params["OUT_PREFIX"] = str(outdir / label)
            params["NEUTRAL_CALIBRATION_OUT"] = str(outdir / f"{label}.neutral_calibration.tsv")
            params["BURNIN_TREES_OUT"] = str(outdir / f"{label}.trees")
            params["BURNIN_SOCIAL_STATE_OUT"] = str(outdir / f"{label}.social_state.tsv")

            cmd = build_slim_command(slim_bin, slim_script, params)
            command_specs.append((cmd, params, label))
            manifest_rows.append({
                "K": k,
                "replicate": rep,
                "seed": params["RANDOM_SEED"],
                "out_prefix": params["OUT_PREFIX"],
                "neutral_calibration_out": params["NEUTRAL_CALIBRATION_OUT"],
                "burnin_trees_out": params["BURNIN_TREES_OUT"],
                "burnin_social_state_out": params["BURNIN_SOCIAL_STATE_OUT"],
                "command": " ".join(shlex.quote(x) for x in cmd),
                "status": "written",
            })

    write_command_script(write_sh, command_specs, descriptions)
    manifest_path = outdir / "neutral_calibration_manifest.tsv"
    write_manifest(manifest_path, manifest_rows)

    print(f"已写出命令脚本: {write_sh}")
    print(f"已写出manifest: {manifest_path}")
    print(f"命令数量: {len(command_specs)}")

    if args.execute:
        for i, (cmd, _params, label) in enumerate(command_specs, start=1):
            print(f"[{i}/{len(command_specs)}] {label}")
            subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()

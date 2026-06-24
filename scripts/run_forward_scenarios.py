#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_forward_scenarios.py

生成并可选择执行SLiM forward命令。

流程设计：
  - 每个replicate先有一个selected burn-in和一个neutral burn-in。
  - M1-M4 selected共用同一replicate的selected burn-in .trees。
  - M1-M4 neutral共用同一replicate的neutral burn-in .trees。
  - M3/M4使用human_K_factor和hunting文件。
  - M1/M2虽然会传入文件路径用于脚本校验，但不施加人类压力或狩猎。
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

from run_burnin_grid import (
    build_slim_command,
    load_all_params,
    load_param_descriptions,
    write_command_script,
)


def parse_levels(spec: str) -> List[str]:
    return [x.strip() for x in spec.split(",") if x.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description="运行或写出SLiM forward情景命令。")
    parser.add_argument("--config-dir", default="config")
    parser.add_argument("--slim-script", default="slim/rrox_forward.slim")
    parser.add_argument("--slim-bin", default="slim")
    parser.add_argument("--input-dir", default="input")
    parser.add_argument("--burnin-dir", default="outputs/burnin")
    parser.add_argument("--outdir", default="outputs/forward")
    parser.add_argument("--logs-dir", default="outputs/logs")
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--replicate-start", type=int, default=1)
    parser.add_argument("--seed-base", type=int, default=900000)
    parser.add_argument("--genetic-modes", default="selected,neutral")
    parser.add_argument("--scenarios", default="M1,M2,M3,M4")
    parser.add_argument("--human-levels", default="medium", help="人类压力水平，如low,medium,high，多个值用英文逗号分隔。")
    parser.add_argument("--hunting-levels", default="high", help="狩猎压力水平，如low,medium,high，多个值用英文逗号分隔。")
    parser.add_argument("--climate-level", default=None, help="若设置，则使用climate_K_factor_<level>.tsv。")
    parser.add_argument("--k-label", default=None, help="可选K标签，用于匹配burn-in tree文件名。")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--write-sh", default="outputs/logs/run_forward_scenarios.sh")
    args = parser.parse_args()

    config_dir = Path(args.config_dir)
    slim_script = Path(args.slim_script)
    input_dir = Path(args.input_dir)
    burnin_dir = Path(args.burnin_dir)
    outdir = Path(args.outdir)
    logs_dir = Path(args.logs_dir)

    outdir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    base_params = load_all_params(config_dir)
    base_params["RUN_MODE"] = "forward"
    descriptions = load_param_descriptions(config_dir)

    genetic_modes = parse_levels(args.genetic_modes)
    scenarios = parse_levels(args.scenarios)
    human_levels = parse_levels(args.human_levels)
    hunting_levels = parse_levels(args.hunting_levels)

    if args.climate_level:
        climate_file = input_dir / f"climate_K_factor_{args.climate_level}.tsv"
    else:
        climate_file = input_dir / "climate_K_factor.tsv"

    if not climate_file.exists():
        print(f"WARNING: 气候文件目前不存在: {climate_file}")

    command_specs: List[Tuple[List[str], Dict[str, object], str]] = []

    for rep in range(args.replicate_start, args.replicate_start + args.replicates):
        for genetic_mode in genetic_modes:
            if args.k_label:
                burnin_tree = burnin_dir / f"burnin_rep{rep:03d}_{genetic_mode}_{args.k_label}.trees"
                burnin_social_state = burnin_dir / f"burnin_rep{rep:03d}_{genetic_mode}_{args.k_label}.social_state.tsv"
            else:
                # 匹配run_burnin_grid.py在包含K时的默认标签格式。
                # 如果burn-in文件名不同，请使用--k-label或调整这个模板。
                k_int = int(round(float(base_params["K5200"])))
                burnin_tree = burnin_dir / f"burnin_rep{rep:03d}_{genetic_mode}_K{k_int}.trees"
                burnin_social_state = burnin_dir / f"burnin_rep{rep:03d}_{genetic_mode}_K{k_int}.social_state.tsv"

            for scenario in scenarios:
                if scenario not in {"M1", "M2", "M3", "M4"}:
                    raise ValueError(f"未知情景: {scenario}")

                # M1/M2不使用human/hunting，但仍传入文件路径用于SLiM脚本校验。
                level_pairs = [("none", "none")]
                if scenario in {"M3", "M4"}:
                    level_pairs = [(h, q) for h in human_levels for q in hunting_levels]

                for human_level, hunting_level in level_pairs:
                    params = dict(base_params)
                    params["SCENARIO"] = scenario
                    params["GENETIC_MODE"] = genetic_mode
                    params["REPLICATE_ID"] = rep
                    params["RANDOM_SEED"] = (
                        args.seed_base
                        + rep
                        + {"M1": 10000, "M2": 20000, "M3": 30000, "M4": 40000}[scenario]
                        + (0 if genetic_mode == "selected" else 500000)
                    )
                    params["BURNIN_TREES_IN"] = str(burnin_tree)
                    params["BURNIN_SOCIAL_STATE_IN"] = str(burnin_social_state)
                    params["CLIMATE_FILE"] = str(climate_file)

                    if scenario in {"M3", "M4"}:
                        params["HUMAN_FILE"] = str(input_dir / f"human_K_factor_{human_level}.tsv")
                        params["HUNTING_FILE"] = str(input_dir / f"hunting_{hunting_level}.tsv")
                        pressure_label = f"human-{human_level}_hunt-{hunting_level}"
                    else:
                        params["HUMAN_FILE"] = str(input_dir / "human_K_factor_medium.tsv")
                        params["HUNTING_FILE"] = str(input_dir / "hunting_high.tsv")
                        pressure_label = "no-human"

                    label = f"rep{rep:03d}_{genetic_mode}_{scenario}_{pressure_label}"
                    params["OUT_PREFIX"] = str(outdir / label)

                    cmd = build_slim_command(args.slim_bin, slim_script, params)
                    command_specs.append((cmd, params, label))

    sh_path = Path(args.write_sh)
    write_command_script(sh_path, command_specs, descriptions)

    print(f"已写出命令脚本: {sh_path}")
    print(f"命令数量: {len(command_specs)}")

    if args.execute:
        for i, (cmd, _params, _label) in enumerate(command_specs, start=1):
            print(f"[{i}/{len(command_specs)}] {' '.join(shlex.quote(x) for x in cmd)}")
            subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()

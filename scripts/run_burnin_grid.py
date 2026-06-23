#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_burnin_grid.py

生成并可选择执行SLiM burn-in命令。

设计：
  - 每个replicate分别生成selected和neutral两类burn-in运行。
  - selected和neutral burn-in共用生活史、社会结构、基因组和人口学参数；
    GENETIC_MODE控制是否施加选择。
  - 输出的.trees文件供M1-M4 forward情景使用。
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


PARAM_CONFIG_FILES = [
    "core_params.psv",
    "life_history_params.psv",
    "social_params.psv",
    "genetic_params.psv",
]

NON_CONFIG_PARAM_DESCRIPTIONS = {
    "RUN_MODE": "runner生成：SLiM运行阶段，burnin表示burn-in，forward表示前向情景模拟。",
    "GENETIC_MODE": "runner生成或覆盖：遗传模式；selected表示包含选择，neutral表示中性对照。",
    "REPLICATE_ID": "runner生成：重复编号，用于区分独立模拟重复。",
    "RANDOM_SEED": "runner生成：随机数种子，用于保证每次重复可追踪。",
    "OUT_PREFIX": "runner生成：输出文件名前缀，各类结果文件会在此基础上添加后缀。",
    "BURNIN_TREES_OUT": "runner生成：burn-in结束后保存的tree-sequence文件路径。",
    "BURNIN_TREES_IN": "runner生成：forward阶段读取的burn-in tree-sequence文件路径。",
    "SCENARIO": "runner生成：forward情景编号，M1基线，M2气候，M3人类活动，M4联合情景。",
    "CLIMATE_FILE": "runner生成：气候适宜栖息地K因子逐年表路径。",
    "HUMAN_FILE": "runner生成：人类活动K因子逐年表路径。",
    "HUNTING_FILE": "runner生成：狩猎死亡率逐年表路径。",
    "CHROM_COUNT": "由genome_structure.psv派生：染色体或连锁群数量。",
    "GENE_COUNT": "由genome_structure.psv派生：压缩基因数量。",
    "GENE_LENGTH": "由genome_structure.psv派生：每个压缩基因的长度，单位bp。",
    "INTERGENIC_RECOMB_RATE": "由genome_structure.psv派生：基因间隔区重组率。",
}


def parse_scalar_value(value: str, typ: str):
    value = value.strip()
    typ = typ.strip().lower()

    if value == "NA":
        return value

    if typ == "int":
        return int(value)
    if typ == "float":
        return float(value)
    if typ == "bool":
        if value in {"1", "true", "True", "yes", "YES"}:
            return 1
        if value in {"0", "false", "False", "no", "NO"}:
            return 0
        raise ValueError(f"无法解析bool值: {value}")
    return value


def read_params_psv(path: Path) -> Dict[str, object]:
    """
    读取对齐PSV参数表，列为：
      parameter | value | type | required | description
    """
    params: Dict[str, object] = {}

    with path.open("r", encoding="utf-8") as f:
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


def read_param_descriptions_psv(path: Path) -> Dict[str, str]:
    """
    从对齐PSV参数表读取description列，列为：
      parameter | value | type | required | description
    """
    descriptions: Dict[str, str] = {}

    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue

            fields = [x.strip() for x in line.split("|")]
            if fields[0] == "parameter":
                continue

            if len(fields) < 5:
                raise ValueError(f"{path}: 参数行格式错误: {line}")

            name, _value, _typ, _required, desc = fields[:5]
            descriptions[name] = desc

    return descriptions


def read_genome_structure_psv(path: Path) -> Dict[str, object]:
    """
    读取genome_structure.psv并派生：
      CHROM_COUNT
      GENE_COUNT
      GENE_LENGTH
      INTERGENIC_RECOMB_RATE
    """
    rows = []
    with path.open("r", encoding="utf-8") as f:
        header = None
        for raw in f:
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = [x.strip() for x in line.split("|")]
            if header is None:
                header = fields
                continue
            if len(fields) < len(header):
                raise ValueError(f"{path}: 行格式错误: {line}")
            rows.append(dict(zip(header, fields)))

    if not rows:
        raise ValueError(f"{path}: 未找到基因组结构行")

    chrom_count = len(rows)
    gene_count = sum(int(r["gene_count"]) for r in rows)
    gene_lengths = sorted({int(r["gene_length_bp"]) for r in rows})
    intergenic_rates = sorted({r["intergenic_recomb_rate"] for r in rows})

    if len(gene_lengths) != 1:
        raise ValueError("当前runner要求所有行使用同一个gene_length_bp")
    if len(intergenic_rates) != 1:
        raise ValueError("当前runner要求所有行使用同一个intergenic_recomb_rate")

    return {
        "CHROM_COUNT": chrom_count,
        "GENE_COUNT": gene_count,
        "GENE_LENGTH": gene_lengths[0],
        "INTERGENIC_RECOMB_RATE": intergenic_rates[0],
    }


def slim_literal(value: object) -> str:
    """
    将Python值转换为SLiM -d字面量。
    字符串必须为SLiM加引号。
    """
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.12g}"

    text = str(value)
    # 已经是科学计数法形式的数值字符串？
    try:
        float(text)
        if any(c in text for c in ".eE"):
            return text
    except ValueError:
        pass

    escaped = text.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def build_slim_command(slim_bin: str, slim_script: Path, params: Dict[str, object]) -> List[str]:
    cmd = [slim_bin]
    for key in sorted(params.keys()):
        cmd.extend(["-d", f"{key}={slim_literal(params[key])}"])
    cmd.append(str(slim_script))
    return cmd


def load_all_params(config_dir: Path) -> Dict[str, object]:
    params: Dict[str, object] = {}

    for filename in PARAM_CONFIG_FILES:
        path = config_dir / filename
        if path.exists():
            params.update(read_params_psv(path))
        else:
            print(f"WARNING: 缺少可选config文件: {path}")

    genome_path = config_dir / "genome_structure.psv"
    if genome_path.exists():
        params.update(read_genome_structure_psv(genome_path))
    else:
        print(f"WARNING: 缺少可选genome_structure文件: {genome_path}")

    return params


def load_param_descriptions(config_dir: Path) -> Dict[str, str]:
    descriptions: Dict[str, str] = {}

    for filename in PARAM_CONFIG_FILES:
        path = config_dir / filename
        if path.exists():
            descriptions.update(read_param_descriptions_psv(path))

    descriptions.update(NON_CONFIG_PARAM_DESCRIPTIONS)
    return descriptions


def write_command_script(
    sh_path: Path,
    command_specs: Iterable[Tuple[List[str], Dict[str, object], str]],
    descriptions: Dict[str, str],
) -> None:
    command_specs = list(command_specs)
    sh_path.parent.mkdir(parents=True, exist_ok=True)

    with sh_path.open("w", encoding="utf-8") as f:
        f.write("#!/usr/bin/env bash\nset -euo pipefail\n\n")
        for cmd, _params, _label in command_specs:
            f.write(" ".join(shlex.quote(x) for x in cmd) + "\n")

        f.write("\n")
        f.write("# ============================================================\n")
        f.write("# 人工QA参数说明\n")
        f.write("# config/*.psv 的 description 是普通参数解释的单一来源。\n")
        f.write("# runner派生参数在生成脚本中补充解释。\n")
        f.write("# 下方列出每条命令实际传给 SLiM 的 -d 参数和值。\n")
        f.write("# ============================================================\n")

        for i, (_cmd, params, label) in enumerate(command_specs, start=1):
            f.write("\n")
            f.write(f"# ---- 命令 {i}: {label} ----\n")
            for key in sorted(params.keys()):
                desc = descriptions.get(
                    key,
                    "未在config/*.psv的description中找到；请补充description或确认该参数是否仍被SLiM使用。",
                )
                f.write(f"# {key} = {slim_literal(params[key])}\n")
                f.write(f"#   {desc}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="运行或写出SLiM burn-in网格命令。")
    parser.add_argument("--config-dir", default="config")
    parser.add_argument("--slim-script", default="slim/rrox_burnin.slim")
    parser.add_argument("--slim-bin", default="slim")
    parser.add_argument("--outdir", default="outputs/burnin")
    parser.add_argument("--logs-dir", default="outputs/logs")
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--replicate-start", type=int, default=1)
    parser.add_argument("--seed-base", type=int, default=100000)
    parser.add_argument("--genetic-modes", default="selected,neutral")
    parser.add_argument("--k-values", default=None, help="可选的K5200网格，多个值用英文逗号分隔。")
    parser.add_argument("--execute", action="store_true", help="实际执行生成的命令。")
    parser.add_argument("--write-sh", default="outputs/logs/run_burnin_grid.sh")
    args = parser.parse_args()

    config_dir = Path(args.config_dir)
    slim_script = Path(args.slim_script)
    outdir = Path(args.outdir)
    logs_dir = Path(args.logs_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    base_params = load_all_params(config_dir)
    base_params["RUN_MODE"] = "burnin"
    descriptions = load_param_descriptions(config_dir)

    genetic_modes = [x.strip() for x in args.genetic_modes.split(",") if x.strip()]
    if args.k_values:
        k_values = [float(x.strip()) for x in args.k_values.split(",") if x.strip()]
    else:
        k_values = [float(base_params.get("K5200", 0))]

    if not k_values or k_values[0] <= 0:
        raise ValueError("必须在config或--k-values中提供有效的K5200")

    command_specs: List[Tuple[List[str], Dict[str, object], str]] = []

    for rep in range(args.replicate_start, args.replicate_start + args.replicates):
        for genetic_mode in genetic_modes:
            for k in k_values:
                params = dict(base_params)
                params["GENETIC_MODE"] = genetic_mode
                params["REPLICATE_ID"] = rep
                params["RANDOM_SEED"] = args.seed_base + rep + (0 if genetic_mode == "selected" else 500000)
                params["K5200"] = k
                params["INITIAL_N"] = int(round(k))
                label = f"rep{rep:03d}_{genetic_mode}_K{int(round(k))}"
                params["OUT_PREFIX"] = str(outdir / label)
                params["BURNIN_TREES_OUT"] = str(outdir / f"burnin_{label}.trees")

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

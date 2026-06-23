#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepare_yearly_inputs.py

为川金丝猴（Rhinopithecus roxellana）SLiM模型生成逐年输入文件。

输出：
  climate_K_factor*.tsv
  human_K_factor_low/medium/high.tsv
  hunting_low/medium/high.tsv

流程设计：
  - SLiM直接读取逐年TSV文件。
  - bp必须覆盖START_BP..END_BP，通常为5200..0。
  - climate_K_factor基于相对于5200 BP的适宜栖息地面积。
  - human_K_factor基于耕地面积压力。
  - hunting死亡率基于聚落数量。

默认考古表使用总流程中提供的数据。
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


DEFAULT_ARCH_TABLE = [
    # start_bp, end_bp, 文化阶段, 聚落数量, 耕地面积（万hm2）
    (6900, 5500, "Yangshao", 9, 0.06),
    (5400, 4000, "Majiayao", 372, 1.99),
    (4000, 3600, "Qijia", 349, 2.01),
    (3400, 2500, "Xindian_Siwa", 68, 0.11),
    (2152, 1730, "Han", 56, 2.10),
    (1332, 1043, "Tang", 63, 5.18),
    (990, 674, "Song", 90, 7.87),
    (679, 582, "Yuan", 34, 5.67),
    (582, 306, "Ming", 240, 7.45),
    (306, 39, "Qing", 310, 12.09),
    (38, 1, "Republic", 528, 34.43),
]


def parse_named_floats(spec: str) -> Dict[str, float]:
    """
    解析如下格式的字符串：
      low:0.8,medium:0.5,high:0.2
    """
    out: Dict[str, float] = {}
    if not spec.strip():
        return out

    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"命名浮点参数格式错误: {item!r}")
        name, value = item.split(":", 1)
        out[name.strip()] = float(value.strip())

    return out


def expected_bp_values(start_bp: int, end_bp: int) -> List[int]:
    if start_bp < end_bp:
        raise ValueError("start_bp必须大于等于end_bp")
    return list(range(start_bp, end_bp - 1, -1))


def validate_yearly_bp(df: pd.DataFrame, start_bp: int, end_bp: int, file_label: str) -> None:
    expected = set(expected_bp_values(start_bp, end_bp))
    observed = set(df["bp"].astype(int))

    if len(df) != len(expected):
        raise ValueError(f"{file_label}: 预期{len(expected)}行，实际{len(df)}行")
    if df["bp"].duplicated().any():
        dup = df.loc[df["bp"].duplicated(), "bp"].head().tolist()
        raise ValueError(f"{file_label}: bp存在重复值，例如{dup}")
    missing = sorted(expected - observed, reverse=True)
    if missing:
        raise ValueError(f"{file_label}: 缺少bp值，例如{missing[:10]}")
    extra = sorted(observed - expected, reverse=True)
    if extra:
        raise ValueError(f"{file_label}: 存在范围外bp值，例如{extra[:10]}")


def build_yearly_arch_series(
    start_bp: int,
    end_bp: int,
    value_column: str,
    arch_table: List[Tuple[int, int, str, int, float]],
) -> pd.DataFrame:
    """
    根据考古表构建逐年序列。
    文化阶段内部保持常数，间隔年份使用线性插值。
    0 BP使用1 BP的数值填充。
    """
    years = np.array(expected_bp_values(start_bp, end_bp), dtype=int)
    df = pd.DataFrame({"bp": years, value_column: np.nan})

    for row in arch_table:
        p_start, p_end, _stage, settlements, cropland = row
        value = settlements if value_column == "settlements" else cropland

        s = min(p_start, start_bp)
        e = max(p_end, end_bp)
        if s < end_bp or e > start_bp:
            continue

        mask = (df["bp"] <= s) & (df["bp"] >= e)
        # 后面的行覆盖边界重叠年份；这与之前处理耕地面积的方式一致。
        df.loc[mask, value_column] = float(value)

    # 对间隔年份插值。由于bp递减，先按bp升序进行线性插值。
    asc = df.sort_values("bp").copy()
    asc[value_column] = asc[value_column].interpolate(method="linear", limit_direction="both")
    out = asc.sort_values("bp", ascending=False).reset_index(drop=True)

    if out[value_column].isna().any():
        raise ValueError(f"未能填满{value_column}的所有逐年数值")

    return out


def write_human_k_files(
    outdir: Path,
    start_bp: int,
    end_bp: int,
    fmin_values: Dict[str, float],
    arch_table: List[Tuple[int, int, str, int, float]],
) -> List[Path]:
    base = build_yearly_arch_series(start_bp, end_bp, "cropland_area", arch_table)

    area_5200 = float(base.loc[base["bp"] == start_bp, "cropland_area"].iloc[0])
    area_max = float(base["cropland_area"].max())
    if area_max <= area_5200:
        raise ValueError("cropland_area_max必须大于start_bp时的cropland_area")

    pressure = (base["cropland_area"] - area_5200) / (area_max - area_5200)
    base["cropland_pressure"] = pressure.clip(lower=0.0, upper=1.0)

    paths: List[Path] = []
    summary_rows = []

    for level, fmin in fmin_values.items():
        if not (0 < fmin <= 1):
            raise ValueError(f"fmin必须在(0,1]范围内，当前为{level}:{fmin}")

        beta_h = -math.log(fmin)

        out = base.copy()
        out["human_K_factor"] = np.exp(-beta_h * out["cropland_pressure"])

        if not ((out["human_K_factor"] > 0) & (out["human_K_factor"] <= 1)).all():
            raise ValueError(f"{level}的human_K_factor超出范围")

        out = out[["bp", "cropland_area", "cropland_pressure", "human_K_factor"]]
        path = outdir / f"human_K_factor_{level}.tsv"
        out.to_csv(path, sep="\t", index=False, float_format="%.6f")
        paths.append(path)

        summary_rows.append({
            "level": level,
            "fmin": fmin,
            "beta_h": beta_h,
            "cropland_area_start_bp": area_5200,
            "cropland_area_max": area_max,
            "file": path.name,
        })

    summary = pd.DataFrame(summary_rows)
    summary_path = outdir / "human_K_factor_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False, float_format="%.12g")
    paths.append(summary_path)

    return paths


def write_hunting_files(
    outdir: Path,
    start_bp: int,
    end_bp: int,
    q_anchor_values: Dict[str, float],
    arch_table: List[Tuple[int, int, str, int, float]],
    anchor_settlements: float,
    write_check: bool,
) -> List[Path]:
    base = build_yearly_arch_series(start_bp, end_bp, "settlements", arch_table)

    paths: List[Path] = []
    summary_rows = []

    for level, q_anchor in q_anchor_values.items():
        if not (0 <= q_anchor < 1):
            raise ValueError(f"q_anchor必须在[0,1)范围内，当前为{level}:{q_anchor}")

        out_full = base.copy()
        out_full["annual_hunting_mortality"] = q_anchor * out_full["settlements"] / anchor_settlements

        if not ((out_full["annual_hunting_mortality"] >= 0) & (out_full["annual_hunting_mortality"] < 1)).all():
            raise ValueError(f"{level}的annual_hunting_mortality超出范围")

        out = out_full[["bp", "annual_hunting_mortality"]]
        path = outdir / f"hunting_{level}.tsv"
        out.to_csv(path, sep="\t", index=False, float_format="%.6f")
        paths.append(path)

        if write_check:
            check_path = outdir / f"hunting_{level}_with_settlements_check.tsv"
            out_full[["bp", "settlements", "annual_hunting_mortality"]].to_csv(
                check_path, sep="\t", index=False, float_format="%.6f"
            )
            paths.append(check_path)

        summary_rows.append({
            "level": level,
            "q_anchor_Majiayao": q_anchor,
            "anchor_settlements": anchor_settlements,
            "q_Yangshao_9": q_anchor * 9 / anchor_settlements,
            "q_Qijia_349": q_anchor * 349 / anchor_settlements,
            "q_Republic_528": q_anchor * 528 / anchor_settlements,
            "file": path.name,
        })

    summary = pd.DataFrame(summary_rows)
    summary_path = outdir / "hunting_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False, float_format="%.12g")
    paths.append(summary_path)

    return paths


def write_climate_files(
    outdir: Path,
    climate_area_file: Path,
    start_bp: int,
    end_bp: int,
    beta_values: Dict[str, float],
) -> List[Path]:
    """
    输入文件必须包含：
      bp<TAB>suitable_area_km2

    100年间隔的数值按阶梯规则扩展为逐年数值：
      5200-5101使用5200时间片
      5100-5001使用5100时间片
      ...
    """
    area_df = pd.read_csv(climate_area_file, sep="\t")
    required = {"bp", "suitable_area_km2"}
    if not required.issubset(area_df.columns):
        raise ValueError(f"{climate_area_file}必须包含列: bp, suitable_area_km2")

    area_df["bp"] = area_df["bp"].astype(int)
    area_df = area_df.sort_values("bp", ascending=False).reset_index(drop=True)

    if start_bp not in set(area_df["bp"]):
        raise ValueError(f"气候面积文件必须包含start_bp={start_bp}")

    a_start = float(area_df.loc[area_df["bp"] == start_bp, "suitable_area_km2"].iloc[0])
    if a_start <= 0:
        raise ValueError("start_bp处的suitable_area_km2必须大于0")

    area_by_bp = dict(zip(area_df["bp"], area_df["suitable_area_km2"]))

    paths: List[Path] = []
    summary_rows = []

    for level, beta_c in beta_values.items():
        if beta_c <= 0:
            raise ValueError(f"beta_c必须大于0，当前为{level}:{beta_c}")

        rows = []
        for bp in expected_bp_values(start_bp, end_bp):
            step_bp = (bp // 100) * 100
            if step_bp > start_bp:
                step_bp = start_bp
            if step_bp < end_bp:
                step_bp = end_bp
            if step_bp not in area_by_bp:
                raise ValueError(f"缺少bp={step_bp}的气候时间片")

            ratio = float(area_by_bp[step_bp]) / a_start
            factor = min(1.0, ratio ** beta_c)
            if factor <= 0 or factor > 1:
                raise ValueError(f"bp={bp}处的climate_K_factor超出范围: {factor}")

            rows.append({"bp": bp, "climate_K_factor": factor})

        out = pd.DataFrame(rows)
        validate_yearly_bp(out, start_bp, end_bp, f"climate_K_factor_{level}")
        path = outdir / f"climate_K_factor_{level}.tsv"
        out.to_csv(path, sep="\t", index=False, float_format="%.6f")
        paths.append(path)

        summary_rows.append({
            "level": level,
            "beta_c": beta_c,
            "area_start_bp": a_start,
            "file": path.name,
        })

    summary = pd.DataFrame(summary_rows)
    summary_path = outdir / "climate_K_factor_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False, float_format="%.12g")
    paths.append(summary_path)

    return paths


def main() -> None:
    parser = argparse.ArgumentParser(
        description="为SLiM生成逐年气候、人类活动和狩猎输入文件。"
    )
    parser.add_argument("--outdir", default="input", help="输出目录。")
    parser.add_argument("--start-bp", type=int, default=5200)
    parser.add_argument("--end-bp", type=int, default=0)
    parser.add_argument(
        "--human-fmin",
        default="low:0.8,medium:0.5,high:0.2",
        help="命名fmin值，例如low:0.8,medium:0.5,high:0.2。",
    )
    parser.add_argument(
        "--hunting-q-anchor",
        default="low:0.010,medium:0.025,high:0.050",
        help="命名马家窑q锚点，例如low:0.010,medium:0.025,high:0.050。",
    )
    parser.add_argument("--anchor-settlements", type=float, default=372.0)
    parser.add_argument(
        "--climate-area-file",
        default=None,
        help="可选TSV文件，需包含bp和suitable_area_km2列。",
    )
    parser.add_argument(
        "--climate-beta",
        default="low:0.5,medium:1.0,high:1.5",
        help="命名beta_c值，例如low:0.5,medium:1.0,high:1.5。",
    )
    parser.add_argument(
        "--write-check-files",
        action="store_true",
        help="写出扩展检查文件，例如hunting_*_with_settlements_check.tsv。",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    paths: List[Path] = []

    paths.extend(write_human_k_files(
        outdir=outdir,
        start_bp=args.start_bp,
        end_bp=args.end_bp,
        fmin_values=parse_named_floats(args.human_fmin),
        arch_table=DEFAULT_ARCH_TABLE,
    ))

    paths.extend(write_hunting_files(
        outdir=outdir,
        start_bp=args.start_bp,
        end_bp=args.end_bp,
        q_anchor_values=parse_named_floats(args.hunting_q_anchor),
        arch_table=DEFAULT_ARCH_TABLE,
        anchor_settlements=args.anchor_settlements,
        write_check=args.write_check_files,
    ))

    if args.climate_area_file is not None:
        paths.extend(write_climate_files(
            outdir=outdir,
            climate_area_file=Path(args.climate_area_file),
            start_bp=args.start_bp,
            end_bp=args.end_bp,
            beta_values=parse_named_floats(args.climate_beta),
        ))

    print("已生成文件:")
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()

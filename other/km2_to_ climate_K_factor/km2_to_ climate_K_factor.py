import pandas as pd
import numpy as np

beta_c = 1.0
start_bp = 5200
end_bp = 0

area_df = pd.read_csv("climate_area_100yr.tsv", sep="\t")
area_df = area_df.sort_values("bp", ascending=False)

A_5200 = area_df.loc[area_df["bp"] == 5200, "suitable_area_km2"].iloc[0]

area_df["climate_K_factor"] = (area_df["suitable_area_km2"] / A_5200) ** beta_c
area_df["climate_K_factor"] = area_df["climate_K_factor"].clip(upper=1.0)

if (area_df["climate_K_factor"] <= 0).any():
    raise ValueError("ERROR: climate_K_factor must be > 0")

rows = []

for bp in range(start_bp, end_bp - 1, -1):
    # 100年阶梯：5200–5101用5200，5100–5001用5100
    step_bp = (bp // 100) * 100

    if step_bp > start_bp:
        step_bp = start_bp

    if step_bp < end_bp:
        step_bp = end_bp

    # 对于5200、5100、5000...正好匹配时间片
    matched = area_df.loc[area_df["bp"] == step_bp, "climate_K_factor"]

    if matched.empty:
        raise ValueError(f"ERROR: missing climate time slice for step_bp = {step_bp}")

    rows.append({
        "bp": bp,
        "climate_K_factor": round(float(matched.iloc[0]), 6)
    })

out_df = pd.DataFrame(rows)

# 硬性检查
if out_df["bp"].duplicated().any():
    raise ValueError("ERROR: duplicated bp in climate_K_factor.tsv")

expected_bp = set(range(5200, -1, -1))
observed_bp = set(out_df["bp"])
missing_bp = expected_bp - observed_bp

if missing_bp:
    raise ValueError(f"ERROR: missing bp values: {sorted(list missing_bp)[:10]}")

if not ((out_df["climate_K_factor"] > 0) & (out_df["climate_K_factor"] <= 1)).all():
    raise ValueError("ERROR: climate_K_factor must be in (0, 1]")

out_df.to_csv("climate_K_factor.tsv", sep="\t", index=False)
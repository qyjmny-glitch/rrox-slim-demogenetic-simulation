# 阶段1：neutral calibration

本阶段用于估计川金丝猴 nonWF + OMU/AMU 社会结构模型中 `K5200` 与中性核苷酸多样性推断的等效 `Ne_pi` 之间的换算关系。

## 文件

| 文件 | 作用 |
|---|---|
| `slim/rrox_neutral_calibration.slim` | 阶段1专用 SLiM 脚本；复制 burn-in 的生活史和社会结构，但遗传输入固定为纯中性 |
| `config/stage1_neutral_calibration_params.psv` | 阶段1专用参数表 |
| `scripts/run_stage1_neutral_calibration.py` | 生成或运行 K 网格 neutral calibration 命令 |
| `scripts/summarize_stage1_neutral_calibration.py` | 汇总每个 K 和 replicate 的稳定窗口统计 |
| `scripts/estimate_K5200_from_observed_Ne.py` | 根据古代观测 Ne 区间反推 K5200 |

`slim/rrox_burnin.slim` 和 `slim/rrox_forward.slim` 保留为已验证主脚本，不在阶段1中修改。

## 生成命令但不运行

```bash
python3 scripts/run_stage1_neutral_calibration.py \
  --config-dir config \
  --stage1-config config/stage1_neutral_calibration_params.psv \
  --write-sh outputs/logs/run_stage1_neutral_calibration.sh
```

检查 `outputs/logs/run_stage1_neutral_calibration.sh` 末尾的 `-d` 参数说明后再运行：

```bash
bash outputs/logs/run_stage1_neutral_calibration.sh
```

## 输出

每个 K × replicate 会输出：

```text
outputs/neutral_calibration/K500_rep001.neutral_calibration.tsv
outputs/neutral_calibration/K500_rep001.trees
outputs/neutral_calibration/K500_rep001.social_state.tsv
```

`*.neutral_calibration.tsv` 的核心字段：

| 字段 | 含义 |
|---|---|
| `K` | 当前候选 K5200 |
| `N` | 当前普查种群数量 Nc |
| `Nc_over_K` | Nc/K |
| `pi_neutral` | 根据中性突变频率计算的核苷酸多样性 |
| `Ne_pi` | `pi_neutral / (4 * MU_RATE)` |
| `Ne_pi_over_K` | Ne_pi/K，用于反推 K5200 |
| `Ne_pi_over_Nc` | Ne_pi/Nc |

## 汇总

```bash
python3 scripts/summarize_stage1_neutral_calibration.py \
  --outdir outputs/neutral_calibration \
  --stable-points 20
```

输出：

```text
outputs/neutral_calibration/neutral_calibration_by_rep.tsv
outputs/neutral_calibration/neutral_calibration_summary.tsv
```

## 根据古代 Ne 反推 K5200

例如古代 Ne 区间为 10000 到 40000：

```bash
python3 scripts/estimate_K5200_from_observed_Ne.py \
  --summary outputs/neutral_calibration/neutral_calibration_summary.tsv \
  --observed-ne-min 10000 \
  --observed-ne-max 40000 \
  --out outputs/neutral_calibration/K5200_estimate_from_observed_Ne.tsv
```

## 注意

阶段1不使用气候、人类活动、狩猎死亡、M1-M4 情景、ROH-like 或 selected DFE 输出。它只用于建立 `K -> pi_neutral -> Ne_pi` 的校准关系。

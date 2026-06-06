# 进度看板

> 更新规则：完成一个子任务就把 `[ ]` 改成 `[x]`，在后面加上完成日期和备注。
> DeepSeek 和 Claude Opus 都可以编辑此文件。

## 状态图例
- `[ ]` 待开始
- `[~]` 进行中
- `[x]` 已完成
- `[!]` 有问题/阻塞

---

## DeepSeek 任务 🟢

| # | 任务 | 状态 | 完成日 | 产出文件 | 备注 |
|---|---|---|---|---|---|
| D1 | 环境与仓库初始化 | [x] | 06-06 | `.gitignore`, `Makefile`, `README.md`, `src/physical_constants.py` | 目录结构已创建 |
| D2 | M1 Patched Conic | [x] | 06-06 | `src/patched_conic.py` | rp=0.2AU 验证通过, Δv=16.37 km/s 与钱学森一致 |
| D3 | M2 N体积分器 | [x] | 06-06 | `src/nbody.py` | 2体基准 PASS (rel_err=1.97e-5), 能量监控就绪 |
| D4 | M3 Horizons 验证 | [x] | 06-06 | `src/horizons_validate.py` | 三层回退(astroquery/proxy/analytic), THETA=295.3deg |
| D5 | M4 月球助推解析 | [x] | 06-06 | (Opus C1 已完成完整实现) | 不需要再单独做 |
| D6 | M6 扫描框架 | [x] | 06-06 | `src/optimize_launch.py` | coarse/fine/refine 三阶段, --quick/--fine 模式 |
| D7 | M8 可视化 | [x] | 06-06 | `src/visualize.py` | 5类图: 轨迹/能量/dv曲线/误差/热力图 |
| D8 | README 草稿 | [x] | 06-06 | `README.md` | 初版已完成，Opus 可润色 |

## Claude Opus 任务 🔴

| # | 任务 | 状态 | 完成日 | 产出文件 | 备注 |
|---|---|---|---|---|---|
| C1 | M4 数值+验证 | [x] | 06-06 | `src/lunar_flyby.py` | 矢量飞掠模型, N体验证误差<1% |
| C2 | M5 单点求解 | [x] | 06-06 | `src/single_trajectory.py` | 完整 pipeline, 约束全通, rp/日期扫描 |
| C3 | M6 结果分析调参 | [x] | 06-06 | `data/scan_2026_*.json` | 全年扫描755解, 最优day107 |
| C4 | M7 灵敏度 | [x] | 06-06 | `src/sensitivity.py` | 四参数偏导数, 月相±2°抵消节省 |
| C5 | LaTeX 报告 | [x] | 06-06 | `report.tex`, `report.pdf` | 13页, XeLaTeX编译通过 |
| C6 | AI-Agent.md | [x] | 06-06 | `AI-Agent.md` | 完整AI使用记录 |
| C7 | README 终稿 | [x] | 06-06 | `README.md` | 模块状态+关键结果+运行说明 |
| C8 | 选做项 O5 动画 | [x] | 06-06 | `figures/trajectory.mp4` | 45s 轨迹动画, ffmpeg编码 |

## 阻塞项

> 任何卡住的事情写在这里，包括需要对方先完成的依赖。

| 阻塞方 | 被什么阻塞 | 解决条件 |
|---|---|---|
| ~~C1~~ | ~~D3~~ | ✅ D3 已完成，阻塞解除 |
| ~~C2~~ | ~~D4+D5~~ | ✅ C1 已完成矢量飞掠模型，C2 直接整合 |
| ~~C3~~ | ~~D6~~ | ✅ D6 已完成，阻塞解除 |

---

*最后更新：2026-06-06 — 全部任务完成 (DeepSeek D1-D8 ✅ + Claude Opus C1-C8 ✅)*

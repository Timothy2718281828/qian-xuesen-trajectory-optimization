# AI-Agent.md — AI 工具使用记录

> 本文档记录本项目中所有 AI 辅助工具的使用方式、交互摘要和产出。

## 使用的 AI 工具

| 工具 | 模型 | 用途 |
|---|---|---|
| DeepSeek (Cursor) | DeepSeek v4-pro | 代码量大的基础模块开发 |
| Claude Opus (Cursor) | Claude Opus 4.6 | 跨模块联调、物理诊断、报告撰写 |

## 协作机制

两个 AI 代理通过三个共享文件协同工作：

- **`WORK_DIVISION.md`**：固定分工定义（只读参考）
- **`PROGRESS.md`**：进度看板（两方都更新，`[ ]` → `[x]`）
- **`HANDOFF.md`**：交接日志（每次切换时写一条，记录完成内容、产出文件、已知问题、建议）

工作流程：启动任一 AI → 读取 PROGRESS + HANDOFF → 执行任务 → 更新 PROGRESS → 写 HANDOFF → 切换到另一 AI。

---

## DeepSeek 使用记录

### 第 1 轮

**Prompt 摘要**：
> 先读 PROGRESS.md 和 HANDOFF.md，然后做 DeepSeek 的 D1-D3, D8 任务。

**完成任务**：

| 任务 | 产出 | 要点 |
|---|---|---|
| D1: 环境初始化 | `.gitignore`, `Makefile`, `README.md`, `physical_constants.py` | 统一单位制(km/s/rad)，目录结构 |
| D2: M1 Patched Conic | `patched_conic.py` | rp=0.2AU 验证 Δv=16.37 km/s，与钱学森一致 |
| D3: M2 N体积分器 | `nbody.py` | Velocity-Verlet, 二体基准 rel_err=1.97e-5 |
| D8: README 草稿 | `README.md` | 初版项目说明 |

**已知问题**：月球助推 `moon_flyby_turn()` 的 Δv 方向有误（标量模型局限），留给 Opus 修正。

### 第 2 轮

**Prompt 摘要**：
> 先读 PROGRESS.md 和 HANDOFF.md，继续做 D4, D6, D7。

**完成任务**：

| 任务 | 产出 | 要点 |
|---|---|---|
| D4: M3 Horizons 验证 | `horizons_validate.py` | 三层回退(astroquery/proxy/analytic)，THETA=295.3° |
| D6: M6 扫描框架 | `optimize_launch.py` | coarse/fine/refine 三阶段，--quick 模式验证通过 |
| D7: M8 可视化 | `visualize.py` | 5 类图表（轨迹/能量/dv曲线/误差/热力图） |

**重要发现**：
- Horizons 代理不可用（网络限制），采用解析星历回退
- 月球初始相位 θ₀ = 295.3° 由 DE430 平经度公式校准
- 28 天快速扫描确认 Moon assist 节省 ~0.05-0.08 km/s

---

## Claude Opus 使用记录

### 第 1 轮

**Prompt 摘要**：
> 认真阅读 PROJECT_SPEC.pdf、WORK_DIVISION.md、PROGRESS.md、HANDOFF.md，明确要完成的任务。然后开始 C1 和 C2。

**完成任务**：

| 任务 | 产出 | 要点 |
|---|---|---|
| C1: M4 数值+验证 | `lunar_flyby.py` | 完整矢量飞掠模型，修复了 DeepSeek 标量模型的方向错误 |
| C2: M5 单点求解 | `single_trajectory.py` | 完整 pipeline，约束 C1-C5 全通过，含 N体验证 |

**关键 Debug 过程**：

1. **DeepSeek 的 `moon_flyby_turn` 方向 bug**：标量模型无法正确处理 2D 矢量旋转，导致月球助推反而增加 Δv。Opus 重写为完整矢量模型（双曲线入射→旋转→出射→坐标变换），修复后节省量变为正值（~0.05 km/s）。

2. **N体验证 14% 误差**：初始将火箭放在月球精确位置，触发 `nbody.py` 的 `r < 1e-6` 距离保护，月球引力被跳过。修正：放在 SOI 边界。

3. **N体验证 400% 误差 + KeyError**：入射侧 SOI 初始化不准确（需要精确的撞击参数 b），改为出射侧初始化（用 v_inf_out），验证日心转移段精度。最终误差 < 1%。

### 第 2 轮

**Prompt 摘要**：
> DeepSeek 已完成，继续推进 C3-C8。

**完成任务**：

| 任务 | 产出 | 要点 |
|---|---|---|
| C3: M6 精细调参 | `data/scan_2026_*.json` | 全年扫描 755 可行解，最优 day 107, 节省 0.029 km/s |
| C4: M7 灵敏度 | `sensitivity.py`, `data/sensitivity.json` | 四参数偏导数，发现月相 ±2° 即抵消节省 |
| C5: LaTeX 报告 | `report.tex`, `report.pdf` | 13 页完整报告，XeLaTeX 编译通过 |
| C6: AI-Agent.md | 本文件 | AI 使用记录 |
| C7: README 终稿 | `README.md` | 更新模块状态和关键结果 |

---

## 人工干预

| 干预内容 | 影响 |
|---|---|
| 提供 PROJECT_SPEC.pdf 和分工文件 | 定义任务分配 |
| 指示"继续推进项目" | 触发 AI 读取最新 PROGRESS/HANDOFF 并继续 |
| 无代码层面的手动修改 | 所有代码由 AI 生成并自测 |

## AI 生成代码的验证

| 验证方式 | 结果 |
|---|---|
| M1 解析解 vs 钱学森教材 | Δv 偏差 < 0.1% ✅ |
| M2 二体基准 | 位置误差 1.97e-5 < 1e-4 ✅ |
| M4 N体 vs 解析近日点 | 误差 < 1% ✅ |
| 约束 C1-C5 | 全部通过 ✅ |
| LaTeX 编译 | XeLaTeX 无错误 ✅ |

## AI 的局限与教训

1. **DeepSeek 的标量模型缺陷**：速度快但在复杂矢量物理推理上出错（Δv 方向），需要 Opus 修正。这说明"快速生成"和"深度推理"需要搭配使用。

2. **N体验证的陷阱**：初始条件设置非常敏感——将火箭放在天体精确位置会触发距离保护，放在 SOI 边界的入射侧需要精确的撞击参数。最终选择出射侧初始化作为务实的验证策略。

3. **月球助推的物理现实**：两个 AI 一开始都预期月球助推有显著节省效果，但数值结果一致表明节省 < 0.5%。这是正确的物理结论（月球质量太小），而非代码 bug。AI 能客观接受反直觉的结果。

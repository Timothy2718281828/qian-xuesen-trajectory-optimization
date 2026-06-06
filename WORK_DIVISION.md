# Project Task Division —— DeepSeek + Claude Opus Collaboration

## 项目概览

| 项目 | 详情 |
|---|---|
| 主题 | 钱学森问题扩展求解 —— 含月球引力助推的绕日返回轨道优化 |

## 分工原则

- **DeepSeek v4-pro**：速度快、产能高 → 负责代码量大、逻辑清晰、公式标准化的基础模块
- **Claude Opus**：推理深、写作好 → 负责跨模块联调、物理诊断、LaTeX 报告、长文档

---

## 协同机制

### 三个共享文件

| 文件 | 用途 | 谁更新 |
|---|---|---|
| `WORK_DIVISION.md` | 分工定义（本文件）| 只读参考，基本不变 |
| `PROGRESS.md` | 进度看板，`[ ]` → `[x]` | **两人都更新**，完成任务时打勾 |
| `HANDOFF.md` | 交接日志 | **两人都更新**，每次切换时写一条 |

### 工作流程

```
你启动 DeepSeek
  ↓
1. 先读 PROGRESS.md → 看自己的任务哪些 [ ] 
2. 先读 HANDOFF.md → 看上次 Opus 有没有留言/坑
  ↓
3. 干一个或多个模块
  ↓
4. 更新 PROGRESS.md → [x] 已完成项
5. 写 HANDOFF.md → 记录本次做了什么、文件在哪、有什么坑
  ↓
你切换到 Claude Opus
  ↓
Opus 重复 1-5
```

### 启动 Prompt 模板

**启动 DeepSeek 时：**
> 先读 PROGRESS.md 和 HANDOFF.md，然后继续做 DeepSeek 的未完成任务。

**启动 Claude Opus 时：**
> 先读 PROGRESS.md 和 HANDOFF.md，然后继续做 Claude Opus 的未完成任务。DeepSeek 已完成的部分在 src/ 下，每个 .py 文件头部有接口文档。

---

## DeepSeek 负责板块 🟢

### D1. 环境与仓库初始化
- [ ] 目录结构（`src/`, `data/`）
- [ ] `.gitignore`（排除 build/、*.aux、*.log、大中间数据）
- [ ] `Makefile`（`make all` 一键生成 PDF+图+视频，`make clean` 清理）
- [ ] 读取 week9 参考代码（`JPL.py`, `eclipse_predict.py`, `Euler_vs_Verlet.py`, `jpl_forward.py`）
- [ ] 读取 `report.tex` §1–§5 了解已有框架

### D2. M1 — Patched Conic 解析解
- [ ] 将钱学森书 §3–§5 的逐步推导转为 Python 函数
- [ ] 输入：$r_p$（近日点距）、$r_1$（地球轨道半径）、$K_s$（太阳引力常数）
- [ ] 输出：椭圆轨道根数 + 各段 $\Delta v$
- [ ] 验证 $r_p = 0.2$ AU 算例，偏差 ≤ 0.1%
- [ ] 文件：`src/patched_conic.py`

### D3. M2 — N 体数值积分器
- [ ] Velocity-Verlet 积分器（主步长 3600s，精细段 60s）
- [ ] Sun-Earth-Moon-Rocket 四体加速度计算
- [ ] 二体基准验证（$\mu=4\pi^2$，无量纲化）：1 年位置相对误差 ≤ 10⁻⁴
- [ ] 能量监控（每 1000 步检查相对误差，记录日志）
- [ ] 文件：`src/nbody.py`

### D4. M3 — JPL Horizons 星历对比
- [ ] 通过代理获取 Sun-Earth-Moon 一年星历
- [ ] 用 N 体积分器相同初值传播，逐日比对
- [ ] 确保所有天体位置误差 ≤ 6000 km
- [ ] 注意 `CENTER='@10'` 而非 `CENTER='10'`（否则会得到地面站坐标）
- [ ] 文件：`src/horizons_validate.py`

### D5. M4 解析部分 — 月球引力助推公式
- [ ] 双曲线飞掠解析公式
- [ ] SOI 边界坐标变换（参考 Vallado §12.4 / Curtis §8.10）
- [ ] 文件：`src/lunar_flyby.py`（解析部分）

### D6. M6 扫描框架（骨架）
- [ ] 外层：遍历 2026 年 365 天
- [ ] 内层：$(r_m, r_p)$ 参数空间网格搜索
- [ ] 约束检查（C1-C5）
- [ ] 输出 $\Delta v_{total}(t_0)$ 数据
- [ ] 文件：`src/optimize_launch.py`（框架）

### D7. M8 — 可视化
- [ ] 最优轨迹 2D 黄道面图
- [ ] 能量守恒监控曲线
- [ ] $\Delta v_{total}(t_0)$ 全年扫描曲线
- [ ] 与 Horizons 误差对比图
- [ ] 文件：`src/visualize.py`

### D8. README.md 草稿
- [ ] 项目说明
- [ ] 环境要求
- [ ] 5 分钟内可复现的运行说明

---

## Claude Opus 负责板块 🔴

### C1. M4 数值部分 + 解析 vs 数值验证
- [ ] 在 N 体模拟中实际飞掠月球
- [ ] 与 D5 的解析结果对比
- [ ] 差异分析：定位是积分器问题、坐标变换问题还是初值问题
- [ ] 更新 `src/lunar_flyby.py`

### C2. M5 — 单点轨迹求解（最难的整合模块）
- [ ] 整合 M1 + M2 + M4 的所有模块
- [ ] 实现完整 pipeline：
  1. 地球表面 → 月球 SOI 边界
  2. 月球飞掠转弯
  3. 日心转移轨道 → 近日点
  4. 返回地球 + 再入速度匹配
- [ ] 约束全部满足（C1-C5）
- [ ] 输出该日期的 $\Delta v_{total}$ 及分项明细
- [ ] 与无月球助推（纯 patched conic）对比节能比例
- [ ] 文件：`src/single_trajectory.py`

### C3. M6 结果分析 + 精细调参
- [ ] 在 DeepSeek 的扫描框架基础上
- [ ] 粗扫 → 精化（二分 / 梯度）
- [ ] 分析是否有物理上异常的"最优解"
- [ ] 确定最终 $t_0^*$、$r_m^*$、$r_p^*$、$\Delta v_{total}^*$
- [ ] 输出优化后的数据和图表

### C4. M7 — 灵敏度分析
- [ ] $r_m$ 飞掠距离偏移影响
- [ ] 发射日期 ±N 天
- [ ] 积分步长 h 收敛性
- [ ] 文件：`src/sensitivity.py`

### C5. 报告 LaTeX（最长的单一文件）
- [ ] 填充 `report.tex` 全部章节
- [ ] 公式规范（正确使用数学模式、交叉引用）
- [ ] 图表嵌入与标注
- [ ] 逻辑连贯：问题 → 方法 → 实现 → 结果 → 分析
- [ ] 确保 ≤ 40 页
- [ ] 确保 XeLaTeX 编译通过
- [ ] 参考文献格式规范

### C6. AI-Agent.md
- [ ] 记录 DeepSeek 的使用方式、Prompt 摘要、产出
- [ ] 记录 Claude Opus 的使用方式、Prompt 摘要、产出
- [ ] 记录任何人工修正或干预

### C7. README.md 终稿
- [ ] 在 D8 草稿基础上润色完善
- [ ] 确保 `make all` 可复现

### C8. 选做项（选做）
- [ ] O1 — 3D 扩展（加入月球 5.1° 倾角）
- [ ] O2 — 相对论修正（近日点 ~3μ²/(c²r⁴)）
- [ ] O3 — 自实现 Newton-Raphson / Lambert 求解器
- [ ] O4 — 多次引力助推（Earth→Moon→Venus→Sun→Earth）
- [ ] O5 — 轨迹动画 MP4（30-60 秒）
- [ ] O6 — Tkinter / Web 实时交互
- [ ] O7 — 其他创新

---

## 执行顺序

```
阶段 0: DeepSeek D1          → 环境 + 仓库骨架
阶段 1: DeepSeek D2, D3      → M1 解析解 + M2 积分器（两个独立模块可并行）
阶段 2: DeepSeek D4          → M3 星历验证（依赖 M2）
阶段 3: DeepSeek D5 + D7     → M4 解析 + M8 可视化框架
阶段 4: Claude Opus C1       → M4 数值验证（依赖 D5 + D3）
阶段 5: Claude Opus C2       → M5 单点求解（依赖所有前序模块）
阶段 6: DeepSeek D6          → M6 扫描框架（与 C2 可部分重叠）
阶段 7: Claude Opus C3       → M6 结果分析（依赖 D6 跑完）
阶段 8: Claude Opus C4       → M7 灵敏度（依赖 C3）
阶段 9: Claude Opus C5, C6, C7 → 报告 + 文档
阶段 X: Claude Opus C8       → 选做项（穿插进行）
```

---

## 交接约定

### DeepSeek → Claude Opus 的接口

每个 DeepSeek 模块完成后，在对应 `.py` 文件头部写清楚：

```python
"""
模块：<模块名>
功能：<一句话>
输入：
  - param1: 类型, 含义
  - param2: 类型, 含义
输出：
  - result1: 类型, 含义
调用示例：
  >>> from xxx import yyy
  >>> result = yyy(param1, param2)
已知限制：
  - <当前版本未处理的情况>
"""
```

所有可调用函数集中在 `src/` 下，Claude Opus 直接 `import` 使用。

### 数据产物

- 扫描结果、星历缓存等中间数据统一放在 `data/`，JSON 格式
- 文件名清晰：`scan_2026_coarse.json`, `horizons_cache_2026.json` 等

### 物理量单位约定（全项目统一）

| 物理量 | 单位 |
|---|---|
| 长度 | km |
| 时间 | s |
| 速度 | km/s |
| 质量 | 无量纲（$\mu$ 制） |
| 角度 | rad |
| 参考系 | J2000 黄道面 |

---

## Token 预算参考

| 负责方 | 预估 token |
|---|---|
| DeepSeek 全部代码 | 150K - 250K |
| Claude Opus 代码 + 分析 | 120K - 200K |
| Claude Opus 报告 | 100K - 150K |
| **合计** | **370K - 600K** |

---

*最后更新：2026-06-06*

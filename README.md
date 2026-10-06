# 三轴机床 NURBS 轨迹规划与 NSGA-II 多目标优化

[![License: GPL v3 or later](https://img.shields.io/badge/License-GPLv3%2B-blue.svg)](LICENSE)

面向三轴机床 XYZ 路径的 MATLAB 研究项目，以加工时间、双向轮廓误差和振动指标为三个目标，使用 NSGA-II 搜索 Pareto 解。

几何层采用单位权重五次 B 样条，时间层采用全局七次 B 样条。控制点数量、速度节点数量和染色体维度随输入路径自动确定。

## 功能

- CSV 离散路径输入，五次 B 样条几何拟合和双向轮廓误差检查。
- 综合曲率、轴速度、几何跃度和工艺速度约束。
- LMSC 限速点检测、SPU 速度规划单元划分和七次柔性变速。
- XYZ 速度、加速度、跃度检查，以及可选 Snap 和共振约束。
- 全局时间样条光顺、NSGA-II 多目标优化与断点续算。
- Pareto 候选解筛选，轨迹、约束、振动指标和图形导出。

## 快速开始

需要 MATLAB。项目已有 MATLAB R2026a 的历史验证记录，详见 [算法说明](docs/ALGORITHM.md#验证状态)。公开版本包含核心 MATLAB 源码和示例路径；主流程使用项目自己的样条实现，无需安装原工作目录中的 NURBS、GeoPDEs 或 MATLAB MCP 工具。

```sh
git clone https://github.com/1resister/NSGA-II-NURBS.git
```

在 MATLAB 中将当前文件夹设为克隆得到的仓库根目录，然后执行：

```matlab
cd(fullfile('SR5NURBS', 'SR5NURBS'))
test_speed_planning(false)   % 合成用例及相关回归断言
main_function                % 完整 NSGA-II 优化
```

完整回归还包含实际输入路径和临时 Pareto 导出：

```matlab
test_speed_planning(true)
```

配置集中在 [`trajectory_config.m`](SR5NURBS/SR5NURBS/trajectory_config.m)。当前默认种群规模为 20、迭代次数为 20；可调整 `cfg.nsga.pop` 和 `cfg.nsga.gen`。计算成本取决于路径和约束，完整优化通常比单条轨迹评价耗时更长。

## 输入与结果

默认输入为 [`5stars.csv`](SR5NURBS/SR5NURBS/5stars.csv)，文件名保留原接口；另提供 [`square_20mm_spacing_1mm.csv`](SR5NURBS/SR5NURBS/square_20mm_spacing_1mm.csv)。CSV 至少包含三列 `X,Y,Z`，坐标单位为毫米，时间单位为秒。通过 `cfg.files.input_csv` 切换路径。

优化结束后，可查看并导出 Pareto 候选解：

```matlab
refresh_pareto_candidates
export_selected_solution(1)
```

结果写入代码目录中的 `output/solution_001/`，包括路径、位置、速度、加速度、跃度、Snap、Crackle、约束状态、振动指标及 PNG 图形。运行结果和 `chromosome_*.txt` 检查点由程序生成，不纳入版本控制。

使用相同配置和算法恢复中断的优化：

```matlab
resume_optimization
```

修改路径、约束或目标算法后，应重新运行 `main_function`，避免混用旧检查点。

## 项目结构

```text
SR5NURBS/SR5NURBS/   核心 MATLAB 源码、回归用例和示例 CSV
docs/ALGORITHM.md    算法细节、参数、导出字段和历史验证记录
LICENSE              GNU GPL v3 许可证全文
THIRD_PARTY_NOTICES.md  第三方代码来源及许可声明
```

这是从当前开发源码整理的独立公开版本。工作目录中的备份、答辩材料、参考论文、工具服务、第三方工具箱副本和历史运行结果未包含在本仓库中。

## 开源许可

项目新增代码和修改按 **GNU General Public License v3.0 or later**（`GPL-3.0-or-later`）发布。第三方文件保留原作者的版权和许可声明：NSGA-II 基础实现采用 BSD 2-Clause，兼容函数 `findspan.m` 采用 GPL v3 or later。具体文件和完整声明见 [第三方声明](THIRD_PARTY_NOTICES.md)。

欢迎通过 GitHub Issues 提供复现步骤、输入路径及配置，或通过 Pull Requests 提交改进。

# 三轴机床 NURBS 轨迹规划与 NSGA-II 多目标优化

## 项目目标

本项目面向三轴机床 XYZ 路径的时间、轮廓误差和振动三目标优化。当前采用“两层光顺”：几何层使用单位权重五次 B 样条，时间层使用全局七次 B 样条；控制点数、速度节点数及染色体维度会根据输入路径和误差目标自动确定。完整流程为：

```text
双向轮廓误差约束下的 C4 五次 B 样条几何拟合
  -> 综合速度约束
  -> LMSC（局部最小速度约束点）
  -> SPU（速度规划单元）
  -> 七次多项式反向/正向可达性扫描
  -> SPU 内加速/匀速/减速规划
  -> SPU 内部速度上界检查
  -> XYZ 单轴速度/加速度/跃度检查
  -> 超限邻域局部降速并重新规划
  -> 全局七次 B 样条 S(t) 二次光顺与约束缩放
  -> 固定 Ts 时间采样和结果导出
```

三个优化目标保持不变：

- `Time`：总加工时间 `time(end)`。
- `Error`：输入折线到光顺曲线、光顺曲线到输入折线两种距离的最大值。
- `Vibration`：未人为削峰的 X、Y、Z 三轴跃度 RMS、Snap、Crackle 与可选机床共振的组合指标。

振动目标先分别计算三轴分量，再合成为一个标量供 NSGA-II 使用：

```text
Vx = sqrt(mean(Jx^2))
Vy = sqrt(mean(Jy^2))
Vz = sqrt(mean(Jz^2))
Vrms = sqrt((wx*Vx)^2 + (wy*Vy)^2 + (wz*Vz)^2)
Vibration = Vrms + snap_weight*Ts*RMS(XYZ snap)
          + crackle_weight*Ts^2*RMS(XYZ crackle)
          + resonance_objective
```

其中 `snap=d(jerk)/dt`，`crackle=d(snap)/dt`。默认 `wx=wy=wz=1`、`peak_weight=0`，因此 XYZ 跃度峰值不参与振动目标，但仍会计算并导出供诊断。每个轴的 RMS、峰值和贡献率都会单独保存。若某一轴对机床振动更敏感，可增大对应权重，例如 `cfg.objectives.vibration_axis_weights = [1 1 1.5]` 会加强对 Z 轴振动的惩罚。

## 目录与入口

MATLAB 代码目录：

```text
SR5NURBS/SR5NURBS
```

主要入口：

- `main_function.m`：运行完整 NSGA-II 优化。
- `evaluate_trajectory.m`：评价一条按当前路径动态确定维度的染色体。
- `export_selected_solution(solution_id)`：从 `chromosome_final.txt` 导出指定 Pareto 解，调用方式保持不变。
- `refine_selected_lmsc_jerk_boundary.m`：计算各内部连接点的局部最大可行跃度，并只评价染色体给定的共享利用率。
- `test_speed_planning.m`：速度规划回归测试。

## 决策变量

每次启动时，`configure_dynamic_dimensions.m` 先分析当前输入路径，并在整次 NSGA-II 优化开始前固定本次染色体维度：

| 变量 | 数量 | 默认范围 | 实际作用 |
|---|---:|---|---|
| `x(1:Nw)` | 动态 `12～60` | 当前固定为 `[1, 1]` | 五次 B 样条兼容权重槽；固定单位权重可避免有理权重重新引入方向振荡 |
| `x(Nw+1:Nw+Nv)` | 动态 `7～25` | `[0.6, 1.8]` | 沿弧长均匀布置的速度倍率节点，只调节局部速度上界 |
| `x(Nw+Nv+1)` | 固定 1 | 当前固定为 `[0, 0]` | 旧版共享端点跃度利用率兼容槽；当前由全局时间样条取代 |

`Nv` 个速度倍率通过 PCHIP 插值得到 `speed_scale(s)`：

```text
v_scale(s)     = Vmax * speed_scale(s)
v_candidate(s) = min(v_limit(s), v_scale(s))
```

动态控制点数从 `min_control_points` 开始递增。对每个候选数量，`fit_adaptive_nurbs_path.m` 默认在显著尖角使用简单内部节点；五次 B 样条因此在内部节点保持 `C4` 连续。控制点通过固定首末点的最小二乘拟合得到，并使用三阶差分 `0.003`、四阶差分 `0.001` 的正则项抑制方向振荡。随后采用单位权重生成基准曲线，并同时计算“原始点到连续曲线”和“曲线采样点到原始折线”的最大距离；达到 `max_contour_error` 后停止增加，否则使用 `max_control_points`。动态速度节点数同时考虑最终控制点数和显著转角数量，并限制在 `min_speed_nodes～max_speed_nodes` 内。两种数量只在优化开始前确定一次，优化过程中保持不变。

速度倍率节点已经与 B 样条内部节点向量解耦，速度优化不会再改变拟合好的基准几何。新种群的第一个个体采用单位权重、单位速度倍率和 `rho=0`；其余个体只在尚未固定的变量范围内随机生成。

当前误差目标不再只检查输入采样点。`nurbs_contour_error.m` 一方面在参数域定位最近区间并对输入点执行连续曲线投影，另一方面检查光顺曲线到原始折线的反向距离，防止曲线在稀疏采样点之间鼓出。

弧长重参数化保留 `cfg.sampling.num_points` 个等弧长基础点，并自动加入所有 NURBS 节点和相邻节点区间中点对应的弧长位置。尖角附近的成簇节点因此会触发局部加密，曲率和 XYZ 动态约束不会因为固定 600 点恰好漏过狭窄过渡区而被低估。

速度倍率大于 1 时仍不能突破物理综合速度上界 `v_limit`；倍率较小会降低相应弧长区间的候选峰值速度，因此会同时影响 `Time` 和基于跃度的 `Vibration`。

旧版 `chromosome_final.txt` 的目标值对应旧速度规划器，不能与新版目标值直接比较；正式使用时应重新运行优化。

## 综合速度约束

所有速度规划都在弧长域 `S` 中进行。`generate_velocity_constraints.m` 分别保存下列约束：

```text
v_limit(s) = min(v_command, v_curvature, v_chord, v_axis,
                 v_geometric_jerk, v_process)
```

其中：

```text
v_command   = Vmax
v_curvature = sqrt(Acmax / max(kappa, epsilon_kappa))
v_axis      = min(Vxmax/|Tx|, Vymax/|Ty|, Vzmax/|Tz|)
v_geometric_jerk = min_axis((Jxyz_max/|q_sss|)^(1/3))
```

直线区的曲率接近零时，曲率约束为 `Inf`，不会错误地产生零速度。可选弦高约束为：

```text
v_chord = 2/Ts * sqrt(max(0, 2*delta_max/max(kappa,epsilon_kappa) - delta_max^2))
```

默认 `cfg.speed.enable_chord_constraint = false`。工艺速度 `v_process` 默认是 `Inf`，也可配置为标量或沿路径分布的向量。

约束数据分为三层，避免含义混淆：

- `v_limit`：指令、曲率、弦高、轴速度、几何跃度和工艺速度形成的物理综合上界。
- `v_candidate`：加入动态速度倍率节点后的候选上界。
- `v_planning_limit`：经过 XYZ 动态约束迭代修正后，实际交给规划器的上界。

`active_constraint` 标记最终起作用的约束，包括 `command`、`curvature`、`chord`、`axis`、`geometric_jerk`、`process`、`speed_scale` 和 `dynamic`。其中 `geometric_jerk` 在速度规划前根据解析 NURBS 三阶导数限制 `q_sss*v^3`，优先在尖角附近局部降速。

## LMSC 检测和 SPU 划分

`detect_lmsc_points.m` 的基础局部极小条件为：

```text
v_limit(i) <= v_limit(i-1)
v_limit(i) <  v_limit(i+1)
```

同时执行以下过滤：

- 起点和终点强制保留。
- 用户在 `cfg.speed.forced_limit_indices` 指定的限速点强制保留。
- 速度下降量或局部显著度必须超过配置阈值。
- 显著曲率峰值对应的速度谷值予以保留。
- 相邻候选点距离过近时，优先保留强制点；否则保留速度更低的点。
- 小幅数值波动不会生成大量 SPU。

动态约束迭代自动加入的强制 LMSC 使用动态最小弧长间距聚类。两个候选点之间所需的七次过渡时间由速度差、`Amax` 和 `Jmax` 共同计算，过渡距离为平均速度乘以该时间，再乘 `lmsc_transition_distance_factor` 安全系数；同时使用 `lmsc_min_sample_intervals` 与路径弧长采样间距形成分辨率下限。物理过渡距离和采样下限先取较大值，再乘独立的 `lmsc_merge_distance_factor`。当前该倍率为 `1.5`，因此最终合并范围比基础动态阈值扩大 50%；它只改变相邻内部 LMSC 是否合并，不修改七次过渡时间公式。速度差越大要求间距越大，路径采样越密时分辨率下限越小。同一动态邻域只保留速度上限最低（同速时曲率最大）的点，避免连续动态修正生成无法独立完成七次过渡的极短 SPU。起终点仅定义零速度边界，不会吞并最近的内部动态特征；边界与内部 LMSC 的可达性由后续七次扫描处理。起终点和用户在 `forced_limit_indices` 中明确配置的强制点属于受保护点，不会被自动删除。

`build_speed_planning_units.m` 将相邻两个 LMSC 之间定义为一个 SPU。每个 SPU 保存起止索引、弧长范围、完整局部速度上界、候选峰值、规划峰值和三段持续时间。

## 七次柔性变速模型

令 `tau=t/T`，七次速度过渡函数为：

```text
g(tau) = 35*tau^4 - 84*tau^5 + 70*tau^6 - 20*tau^7
v(t)   = v0 + (v1-v0)*g(tau)
a(t)   = (v1-v0)/T * 140*tau^3*(1-tau)^3
j(t)   = (v1-v0)/T^2 * 420*tau^2*(1-tau)^2*(1-2*tau)
```

位移为：

```text
s(t) = s0 + v0*t + (v1-v0)*T*
       (7*tau^5 - 14*tau^6 + 10*tau^7 - 2.5*tau^8)
```

由路径加速度和跃度上限得到最短过渡时间：

```text
T_A = 35/16 * abs(v1-v0) / Amax
T_J = sqrt(84*sqrt(5)/25 * abs(v1-v0) / Jmax)
T   = max(T_A, T_J, minimum_transition_time)
L7  = (v0+v1)/2 * T
```

零跃度基线在每个过渡段两端满足速度指定值，且加速度、跃度和 Snap 均为零。导出所选解时还可使用广义七次速度多项式：保持端点速度、零加速度和零 Snap，同时允许内部 LMSC 使用指定的非零共享跃度。相邻 SPU 在同一 LMSC 使用同一个跃度值，因此 `S/V/A/J` 仍连续；路径起点、终点以及变速段与匀速段的连接仍保持零跃度。

## 七次可达性扫描

`forward_backward_scan.m` 先用：

```text
Aeq7 = 16/35 * Amax
```

形成连续的加减速可达速度包络：

```text
v_reachable(i+1) = sqrt(v(i)^2 + 2*Aeq7*delta_S)
```

扫描不会在每个离散采样间隔上重新执行一次首末加速度、跃度为零的完整七次过渡，避免路径采样越密，计算出的可达峰值反而越低。随后 `velocity_planning_7th.m` 在每个完整 SPU 的可用距离上通过有界二分严格检查：

```text
L7(va, vb) <= available_distance
```

反向扫描从终点零速度开始，形成减速包络；正向扫描从起点零速度开始，形成加速包络。完整 SPU 求峰值时仍使用七次过渡的真实距离、加速度和跃度限制，生成后继续执行内部速度上界与 XYZ 动态约束检查。扫描结果不做可能破坏约束的后处理平滑或逐点裁剪。

连续包络生成后，`enforce_spu_boundary_reachability.m` 只在 LMSC/SPU 边界之间执行完整七次可达性扫描，保证相邻短 SPU 的边界速度在跃度限制下也能连接。该检查以整个 SPU 长度为可用距离，不会恢复逐采样点重启七次过渡造成的采样密度保守性。

这也是大曲率区会提前降速的原因：曲率峰值先降低未来位置的速度上界，反向扫描再把减速起点向前传播，而不是到曲率峰值才突然截断速度。

## SPU 内速度规划与检查

对于长度为 `L`、边界速度为 `vs/ve`、候选峰值为 `vp` 的 SPU：

```text
La = L7(vs, vp)
Ld = L7(vp, ve)
```

- 当 `La + Ld <= L` 时，生成七次加速、匀速、七次减速三段。
- 当 `La + Ld > L` 时，用有界二分求实际峰值 `vm`，使 `L7(vs,vm)+L7(vm,ve)=L`，不生成匀速段。

每个生成段至少使用 `cfg.speed.internal_constraint_check_samples` 个内部采样点；当段持续时间较长时会自动增加采样数，使检查间隔不大于 `Ts`。程序根据七次位移函数计算 `s(t)`，再以线性插值检查：

```text
v(tk) <= v_planning_limit(s(tk)) + tolerance
```

发现内部穿越时，程序降低该 SPU 峰值，必要时把最大超限位置加入新的 LMSC 并重新划分 SPU。程序不会使用 `v=min(v,v_limit)` 逐点裁剪，因为那会破坏加速度和跃度连续性。

## XYZ 运动学和动态约束修正

`xyz_motion.m` 在同一弧长参数和同一时间轴上使用链式法则。NURBS 关于参数的一、二、三阶导数由基函数解析计算，再严格转换成弧长导数；不再在包含极小节点间隔的非均匀弧长网格上连续执行三次 `gradient`。对于 `q(s)`（X、Y 或 Z）：

```text
q_dot   = q_s*v
q_ddot  = q_ss*v^2 + q_s*a
q_dddot = q_sss*v^3 + 3*q_ss*v*a + q_s*j
```

硬约束直接使用未削峰的真实运动量检查：

```text
abs(q_dot)   <= xyz_vmax
abs(q_ddot)  <= xyz_amax
abs(q_dddot) <= xyz_jmax
abs(d(q_dddot)/dt) <= xyz_smax
```

若某个时间点超限，`plan_speed_profile.m` 会：

1. 定位最大超限对应的弧长位置。
2. 根据速度、加速度、跃度和 Snap 超限比例估算局部速度缩放量。
3. 乘以 `cfg.speed.constraint_safety_factor` 安全系数。
4. 在该位置邻域平滑降低 `v_planning_limit`，并强制加入 LMSC。
5. 重新执行 SPU 划分、前后向扫描、七次规划和 XYZ 检查。

达到最大修正次数后仍不满足约束的染色体会标记为不可行并施加罚值；`evaluate_objective.m` 也会把发生数值异常的个别染色体转换为有限大罚值，避免中止整轮 NSGA-II。动力学判定使用路径速度、加速度、跃度以及 XYZ 分解后的速度、加速度、跃度；当 `cfg.constraints.enable_xyz_snap_limit = true` 时还包含 Snap 最大约束比值。任一启用的比值大于 1，或者连续性、弧长单调性等严格可行性检查失败，都会产生 `dynamic_penalty`。该罚值只用于 NSGA-II 排序；`raw_time`、`raw_error`、`raw_vibration` 及 Pareto 候选展示仍保留不加罚的原始值。

局部 LMSC/SPU 修正达到上限后若仍有少量超限，可选择整体拉伸已经生成的七次段时间。该解析缩放保持弧长轨迹和七次段结构不变，使速度、加速度、跃度和 Snap 分别按 `1/gamma`、`1/gamma^2`、`1/gamma^3` 和 `1/gamma^4` 降低，并在固定 `Ts` 时间轴上重新检查。`cfg.speed.enable_global_time_scaling = true` 开启这一可行性兜底；当前项目配置为 `true`。

初次缩放先使用带安全余量的因子找到可行上界；当 `cfg.speed.minimize_global_time_scale = true` 时，再在 `[1, 可行上界]` 内二分回缩，寻找满足全部已启用 V/A/J/Snap 约束的最小可行 `global_time_scale`。关闭 Snap 硬约束后，Snap 不再影响该缩放因子。`cfg.speed.global_time_scale_refinement_iterations` 给出基础细化次数（默认 12 次）；若初始区间较宽，程序会按区间宽度自动补足达到 `cfg.speed.global_time_scale_refinement_tolerance`（默认相对容差 `1e-4`）所需的少量迭代。因此缩放因子会尽可能接近 1，同时避免固定增加所有候选解的计算量；若原始轨迹确实超限，不会为追求等于 1 而放宽约束。

`cfg.smoothing.enable_xyz` 只生成可选的诊断平滑跃度 `filtered_jerk`。硬约束、CSV 导出和 `Vibration` 目标都使用未滤波运动量，滤波不能掩盖超限。`Vibration` 先对 `Jx/Jy/Jz` 分别计算 RMS，再按照 `cfg.objectives.vibration_axis_weights` 加权，并加入 Snap、Crackle 和可选机床共振频带目标；XYZ 跃度峰值保留为诊断量，但当前权重为 0，不进入目标函数。

## 全局时间样条二次光顺

LMSC/SPU 七次规划首先产生满足局部拓扑和动力学约束的参考 `S(t)`。当 `cfg.time_smoothing.enabled=true` 时，`global_time_spline_smoothing.m` 再拟合一个全局七次 B 样条时间律：

- 使用简单内部节点，内部连续性为 `C6`，避免逐 SPU 连接处的高阶导数抖动。
- 首尾连续 5 个控制系数相等，因此路径速度、加速度、跃度和 Snap 在全局起终点均为零。
- 控制系数投影为单调非降，保证 `S(t)` 不倒退、路径速度不为负。
- 对控制系数四阶、五阶差分正则化，分别抑制 Snap 和 Crackle 波动。
- 在固定 `Ts` 时间轴上重新检查路径与 XYZ 的 V/A/J，以及启用时的 Snap；必要时只增加总时长，不逐点裁剪。
- 只有平滑评分改善且全部约束满足时才采用；否则按 `fallback_to_septic=true` 自动保留原七次轨迹。

时间样条控制点数量按 LMSC 数量动态取值，默认范围 `64～128`。`maximum_duration_factor=1.5` 防止为了光顺把加工时间无限拉长。`apply_during_objective=true` 保证 NSGA-II 目标、Pareto 候选表和最终导出使用同一条轨迹，避免目标函数值与仿真结果不一致。

## 机床共振与 Snap 优化

`machine_resonance_metric.m` 对固定时间轴上的 XYZ 加速度和跃度执行单边 FFT。每个进给轴可配置多个二阶模态：

```matlab
cfg.resonance.mode_frequencies_hz = {[fx1 fx2], [fy1 fy2], [fz1]};
cfg.resonance.damping_ratios = {[zx1 zx2], [zy1 zy2], [zz1]};
cfg.resonance.modal_gains = {[gx1 gx2], [gy1 gy2], [gz1]};
cfg.resonance.enabled = true;
```

模态 FRF 对加速度频谱形成预测结构响应，并对固有频率附近的跃度能量提高权重。第三目标函数现在由下列部分组成：

```text
Vibration = XYZ jerk RMS
          + snap_weight * Ts * XYZ snap RMS
          + crackle_weight * Ts^2 * XYZ crackle RMS
          + resonance_objective
```

`response_rms_limits` 可对 X/Y/Z 的预测响应施加硬约束，超限时产生独立的 `resonance_penalty`。该罚值只参与 NSGA-II 排序；导出的 `raw_vibration` 仍是不加罚的目标值。速度节点属于染色体，因此启用共振目标后，NSGA-II 可通过局部调速改变激励频率和幅值，不依赖整体时间缩放。

当前示例配置启用共振评价，并给出 X/Y/Z 三轴 `[45]`、`[3]`、`[78] Hz` 的示例模态。正式加工前应替换为实测固有频率、阻尼比和模态增益；频率大于奈奎斯特频率 `1/(2*Ts)` 的模态不会参与计算，并记录为未解析模态。

Snap 按未滤波 XYZ 跃度的时间导数计算，Crackle 再由 Snap 的时间导数计算。`cfg.constraints.enable_xyz_snap_limit = true` 时，每个轴必须满足 `abs(snap) <= cfg.limits.xyz_smax`，任一轴超限会进入 `dynamic_penalty`；设为 `false` 后，Snap 和 Crackle 仍参与振动目标并继续计算、绘图和导出，但 Snap 不再参与局部限速、`global_time_scale`、可行性或罚值判定。当前三轴 Snap 上限均为 `7000000 mm/s^4`。

### 可选的旧版七次段端点共享跃度

当前默认关闭 `cfg.lmsc_jerk_boundary.enabled` 和 `evaluate_during_objective`，染色体末尾的兼容变量固定为 0；全局时间样条负责高阶连续光顺。只有在需要复现实验时才建议重新启用旧版共享端点跃度：启用后，每个通过轮廓预筛选的染色体使用末尾连续变量 `rho` 设置七次多项式段的端点共享跃度。对每个内部 LMSC 速度谷值，程序计算局部正跃度上限；对两个 LMSC 之间没有匀速平台、只有单一最高速度点的 SPU，程序计算该峰值的局部负跃度幅值上限。

内部 LMSC 使用 `+rho * local_maximum_jerk`，非匀速速度峰值使用 `-rho * local_maximum_peak_jerk`；同一物理连接点的左右七次段共享同一个值。若非零峰值跃度使两过渡段之间出现剩余距离，程序会共同拉伸两段直至填满 SPU，不允许从非零跃度直接跳到匀速段的零跃度。该 `rho` 轨迹重新检查局部速度上界、路径 V/A/J、XYZ V/A/J、可选 Snap、共振约束和整体时间缩放，其实际时间和振动直接进入 Pareto 目标。不可行的非零 `rho` 会被罚掉，不会自动回退到零。

默认 `enforce_vibration_increase_limit=false`，由振动第三目标、XYZ 动态硬约束和共振硬约束共同选择利用率；如需额外限制相对零跃度基线的振动增幅，可打开该开关，并通过 `maximum_vibration_increase` 设置比例。

浮点求解有时会在两个七次过渡段之间留下 `10^-11～10^-9 s` 量级的“匀速段”。它远小于固定采样周期，不是实际匀速平台。程序使用 `numerical_peak_plateau_time_tolerance` 将这类残余段合并成一个速度峰值，并在该峰值设置 `-rho * 局部最大可行跃度幅值`；真正超过阈值的匀速平台仍要求首末跃度为零，以保证与恒速段连续。

胜出利用率对应轨迹的实际加工时间和实际振动值直接进入 NSGA-II 三目标评价；轮廓误差目标保持不变，所有罚值仍只用于内部约束排序。路径起终点以及七次变速段与匀速段的连接必须保持零跃度，因为静止端点和匀速段的跃度都是零；这些位置的局部最大可行共享跃度按连续性定义为零。内部 LMSC 谷值和无匀速平台的单点速度峰值保持零加速度、零 Snap，并由左右相邻七次段共享同一个非零跃度。因此，“每段两端”都被赋值，但只有满足连续性条件的内部连接允许非零值。

## 固定时间采样

最终时间序列以 `cfg.interpolation.Ts` 为周期：

```text
0, Ts, 2*Ts, ..., floor(T_total/Ts)*Ts
```

若 `T_total` 不是 `Ts` 的整数倍，会追加精确终点，因此只有最后一个时间间隔可能短于 `Ts`。同一时间轴上计算并导出：

```text
S, path_velocity, path_acceleration, path_jerk,
X, Y, Z, Vx, Vy, Vz, Ax, Ay, Az, Jx, Jy, Jz
```

所有数组长度一致，末点满足 `S(end)=total_length`，首末路径速度、加速度和跃度均为零。

## 主要配置与单位

项目统一使用：

```text
长度 mm
时间 s
速度 mm/s
加速度 mm/s^2
跃度 mm/s^3
Snap mm/s^4
```

主要配置位于 `trajectory_config.m`：

```matlab
cfg.interpolation.Ts

cfg.limits.Vmax
cfg.limits.Amax
cfg.limits.Jmax
cfg.limits.Acmax
cfg.limits.xyz_vmax
cfg.limits.xyz_amax
cfg.limits.xyz_jmax
cfg.limits.xyz_smax        % 当前 [7000000 7000000 7000000] mm/s^4
cfg.constraints.enable_xyz_snap_limit  % Snap 硬约束开关

cfg.objectives.vibration_axis_weights
cfg.objectives.vibration_peak_weight
cfg.objectives.vibration_snap_weight
cfg.objectives.vibration_crackle_weight

cfg.resonance.enabled
cfg.resonance.enable_diagnostics
cfg.resonance.mode_frequencies_hz
cfg.resonance.damping_ratios
cfg.resonance.modal_gains
cfg.resonance.axis_weights
cfg.resonance.spectral_peak_gain
cfg.resonance.objective_weight
cfg.resonance.response_rms_limits
cfg.resonance.penalty_weight

cfg.pareto.max_time_ratio
cfg.evaluation.max_time_samples  % NSGA-II 候选固定周期采样数上限
cfg.evaluation.enable_contour_prefilter
cfg.evaluation.contour_prefilter_relative_margin

cfg.nurbs.eval_points
cfg.nurbs.corner_density_gain
cfg.nurbs.control_point_interpolation
cfg.nurbs.corner_knot_width_factor
cfg.nurbs.corner_knot_multiplicity
cfg.nurbs.force_unit_weights
cfg.nurbs.bidirectional_contour_error
cfg.nurbs.fit_samples_per_segment_max
cfg.nurbs.fit_regularization
cfg.nurbs.fit_third_difference_regularization
cfg.nurbs.fit_fourth_difference_regularization
cfg.dynamic_dimensions.min_control_points
cfg.dynamic_dimensions.max_control_points
cfg.dynamic_dimensions.control_point_step
cfg.dynamic_dimensions.min_speed_nodes
cfg.dynamic_dimensions.max_speed_nodes
cfg.dynamic_dimensions.control_points_per_speed_node
cfg.dynamic_dimensions.speed_nodes_per_feature
cfg.dynamic_dimensions.feature_angle_threshold_deg

cfg.speed.lmsc_velocity_threshold
cfg.speed.lmsc_prominence_threshold
cfg.speed.lmsc_transition_distance_factor
cfg.speed.lmsc_min_sample_intervals
cfg.speed.lmsc_merge_distance_factor
cfg.speed.constraint_correction_radius_samples
cfg.speed.minimum_transition_time
cfg.speed.velocity_tolerance
cfg.speed.distance_tolerance
cfg.speed.binary_search_max_iterations
cfg.speed.max_constraint_iterations
cfg.speed.enable_global_time_scaling
cfg.speed.max_global_time_scaling_iterations
cfg.speed.minimize_global_time_scale
cfg.speed.global_time_scale_refinement_iterations
cfg.speed.global_time_scale_refinement_tolerance
cfg.speed.constraint_safety_factor
cfg.speed.internal_constraint_check_samples
cfg.speed.enable_chord_constraint

cfg.time_smoothing.enabled
cfg.time_smoothing.apply_during_objective
cfg.time_smoothing.degree
cfg.time_smoothing.min_control_points
cfg.time_smoothing.max_control_points
cfg.time_smoothing.control_points_per_lmsc
cfg.time_smoothing.fit_samples
cfg.time_smoothing.fourth_difference_regularization
cfg.time_smoothing.fifth_difference_regularization
cfg.time_smoothing.max_scaling_iterations
cfg.time_smoothing.maximum_duration_factor
cfg.time_smoothing.minimum_smoothness_improvement
cfg.time_smoothing.fallback_to_septic

cfg.lmsc_jerk_boundary.enabled
cfg.lmsc_jerk_boundary.evaluate_during_objective
cfg.lmsc_jerk_boundary.local_bound_safety_factor
cfg.lmsc_jerk_boundary.maximum_vibration_increase
cfg.lmsc_jerk_boundary.enforce_vibration_increase_limit
cfg.lmsc_jerk_boundary.minimum_duration_fraction
cfg.lmsc_jerk_boundary.maximum_duration_factor
cfg.lmsc_jerk_boundary.time_search_samples
cfg.lmsc_jerk_boundary.time_refinement_iterations
cfg.lmsc_jerk_boundary.peak_search_samples
cfg.lmsc_jerk_boundary.numerical_peak_plateau_time_tolerance
cfg.bounds.jerk_utilization_lb
cfg.bounds.jerk_utilization_ub
```

## 输入文件

默认输入为：

```text
SR5NURBS/SR5NURBS/5stars.csv
```

CSV 至少包含三列 `X,Y,Z`。程序会移除非有限行和相邻重复点，坐标单位必须与上述配置一致。

## 运行优化

在 MATLAB 中执行：

```matlab
cd('SR5NURBS/SR5NURBS')
main_function
```

当前配置为 20 个个体、20 代；可通过 `cfg.nsga.pop` 和 `cfg.nsga.gen` 调整。完整优化包含大量染色体评价，运行时间明显长于单条轨迹测试。调试阶段建议先用 20～40 个个体、20～40 代验证流程，再恢复正式规模。建议修改参数或算法前先运行：

```matlab
test_machine_resonance       % 独立检查共振频率加权、响应硬约束和 Snap
test_speed_planning(false)   % 8 个合成速度规划用例
```

再运行完整测试：

```matlab
test_speed_planning(true)    % 另含 5stars 和临时 Pareto 导出
```

### 七次段端点共享跃度优化的使用方式

该功能改变了时间和振动目标函数的计算方式，因此启用后必须从头重新运行 NSGA-II，不能把旧的 `chromosome_iter.txt` 接着迭代：

```matlab
clear functions
rehash path
main_function
```

优化完成后按原命令导出 Pareto 解：

```matlab
export_selected_solution(1)
```

NSGA-II 每次完整染色体评价、Pareto 候选重算和所选解导出都读取染色体末尾的同一个连续利用率变量。程序先建立零跃度基线以计算各连接点的局部上限，再只评价该染色体指定的一个非零利用率，不再枚举固定候选列表。命令行和 CSV 会显示该利用率、时间变化以及 `global_time_scale` 的变化。若关闭此功能，将 `cfg.lmsc_jerk_boundary.enabled` 或 `cfg.lmsc_jerk_boundary.evaluate_during_objective` 设为 `false`。

导出阶段不会再执行第二次独立精修。`pareto_candidates.csv`、`pareto_front_selected.png` 和所选解导出都表示已经选好端点跃度利用率的同一条轨迹，避免候选目标值与仿真结果不一致。

### 验证状态

2026-08-27 使用 MATLAB R2026a 完成当前双光顺版本验证：

- `test_speed_planning(false)` 的 8 个合成用例全部通过。
- `test_speed_planning(true)` 的 10 个完整用例全部通过，包括当前 `5stars.csv` 和临时 Pareto 导出。
- 当前输入动态得到 24 个几何控制点、9 个速度节点，单位权重基线双向轮廓误差为 `0.0809704 mm`，满足 `0.15 mm` 上限。
- 当前参考解的全局时间样条使用 128 个控制点，被接受时的总时长倍率为 1，平滑评分改善约 `0.238%`。
- Pareto 原始目标重算和选中解导出都启用同一全局时间样条，测试确认导出目标与候选表一致。
- 回归还检查了 Crackle 指标、固定 `Ts`、局部速度上界、路径/XYZ 动力学约束、首末零状态、时间/弧长单调性及连续性。
- 当前正方形输入在 `0.15 mm` 误差和尖角窄限速条件下，参考解总时间约 `79.252 s`、`global_time_scale≈103.16`。若需显著缩短正方形加工时间，应优先放宽尖角轮廓误差、明确允许圆角半径，或采用“角点精确停机 + 分段 C4”模式，而不是增大速度上限。

2026-07-30 在 MATLAB R2026a 中运行旧版速度与振动指标时，10 个用例全部通过：

- 速度节点下界和上界用例确认动态速度节点会实际影响运行时间。
- 真实 `5stars.csv` 用例总时间为 `6.8005 s`，包含 7 个 LMSC 和 6 个 SPU。
- 临时 Pareto 候选解导出成功，所有规定 CSV 均已生成并检查。
- 所有用例均通过局部速度上界、路径加速度/跃度、XYZ 单轴动态约束、首末零状态、时间/弧长单调性和段间连续性检查。
- 多曲率峰、真实路径和 Pareto 候选解分别使用 `1.1171`、`1.1283` 和 `1.1782` 的整体时间缩放兜底；其余用例无需整体缩放。该缩放保持七次段连续，不进行逐点速度裁剪。

解析 NURBS 三阶导数、当前三轴 `200000 mm/s^3` 跃度上限、几何跃度局部限速及 RMS+峰值振动目标会改变上述历史数值。代码已加入解析导数、几何跃度限速和峰值目标的回归断言；修改后应先在已打开的 MATLAB 会话运行 `test_speed_planning(false)`，通过后再重新运行完整优化。

2026-08-02 新增的共振频谱、FRF 加权响应、响应硬约束和 Snap 代码已通过 MATLAB `mlint` 静态分析，并加入 40 Hz 合成模态回归断言。独立 MATLAB R2026a 仍在启动阶段发生访问冲突，尚未执行到项目代码，因此这些新增断言需要在已打开的 MATLAB 会话中运行确认。

2026-08-25 LMSC 正共享跃度、非匀速 SPU 峰值负共享跃度、整体缩放和导出代码已通过 MATLAB `mlint` 静态分析。独立数值验证中峰值左右 `V/A/J/Snap` 最大连接误差为 `2.91e-11`，拉伸后 SPU 距离误差为 `5.33e-15`。当前 MCP 无法附加到 MATLAB 会话，独立 MATLAB 又报 `File system inconsistency`，因此完整动态回归仍需在已打开的 MATLAB 会话执行。

同日已先把统一利用率搜索移入 `evaluate_trajectory`，随后又将其改为染色体末尾的一个连续决策变量。评价、候选缓存、导出与绘图链路已通过 `mlint`；由于上述 MATLAB 会话故障，仍需在可用会话中运行 `test_speed_planning(false)` 和 `test_speed_planning(true)` 完成动态确认。

## Pareto 解选择与导出

优化完成后会生成：

```text
chromosome_final.txt
output/pareto_candidates.csv
output/pareto_candidates_cache.mat
```

候选表默认只保留真实加工时间不超过最短候选时间 10 倍的解：

```text
Time <= cfg.pareto.max_time_ratio * min(Time)
cfg.pareto.max_time_ratio = 10
```

这里使用重新仿真得到的不含罚值时间。过滤仅影响 `pareto_candidates.csv` 和可导出的候选编号，不会删除 `chromosome_final.txt` 中的原始最终种群；设为 `Inf` 可关闭过滤。查看过滤后的 `pareto_candidates.csv` 后，按原接口导出指定解：

原始目标重算按内部时间目标从小到大处理。得到当前最短原始时间后，如果后续候选在零跃度基线或染色体端点跃度利用率评价中已能证明其时间超过 `cfg.pareto.max_time_ratio` 对应上限，程序会把它识别为超长候选并直接跳过，而不是误报为不可行跃度并中止候选重算。整体缩放的可行因子搜索只计算 XYZ 约束摘要，选定最小因子后才生成一次完整 XYZ 序列，从而降低 Pareto 批量重算的峰值内存。

候选表中的 `Time`、`Error` 和 `Vibration` 是该染色体端点共享跃度利用率对应的原始物理指标，不包含轮廓误差罚值、动态约束罚值或共振罚值；最后一列 `endpoint_jerk_utilization` 就是该染色体参与交叉和变异的连续变量。NSGA-II 的非支配等级使用同一轨迹的实际目标再叠加约束罚值，因此端点跃度会实际影响优化和选解结果。

最终 `pareto_front.csv`、`pareto_front_selected.png`、`plot_pareto.m` 以及优化结束时显示的 Pareto 图均使用上述原始目标值，不显示罚值。迭代过程图仍显示 NSGA-II 内部罚后目标以便观察约束淘汰过程，并在图标题中明确标为 `internal penalized objectives`。

```matlab
export_selected_solution(3)
```

首次生成候选表时，程序会逐项处理全部 `rank=1` 解，使用各自染色体保存的利用率重新评价并得到不含罚值的原始物理目标，然后写入 `pareto_candidates_cache.mat`。后续执行 `export_selected_solution` 时，如果最终种群、输入 CSV、配置和项目 MATLAB 源码均未变化，将直接复用候选缓存，只按相同规则重新评价被选中的一个解；不会重复评价其余 Pareto 候选。上述任一内容发生变化时缓存自动失效并刷新。

NSGA-II 保存最终种群后会先从 `chromosome_final.txt` 按相同精度重新载入，再建立候选缓存，避免内存浮点值与 ASCII 序列化值的微小差异导致主程序重复评价。命令行候选列表逐项显示原始数值，不使用 MATLAB 的统一科学计数缩放。

如需忽略有效缓存并强制重新计算全部候选，可执行：

```matlab
export_selected_solution(3, true)
```

如果只想刷新过滤后的候选表、暂不导出任何单个解，可执行：

```matlab
refresh_pareto_candidates
```

结果写入：

```text
output/solution_003/
```

也可在 `trajectory_config.m` 设置：

```matlab
cfg.selection.selected_solution_id = 3;
```

## 导出文件

保留的主要 CSV：

| 文件 | 内容 |
|---|---|
| `smooth_path.csv` | NURBS 光顺 XYZ 路径 |
| `position.csv` | 固定时间轴上的 `time,S,X,Y,Z` |
| `velocity.csv` | `time,S,Vx,Vy,Vz` |
| `acceleration.csv` | `time,S,Ax,Ay,Az` |
| `jerk.csv` | `time,S,Jx,Jy,Jz` |
| `jerk_components.csv` | XYZ 跃度的几何项、加速度耦合项、切向路径跃度项及合计值 |
| `snap.csv` | `time,S,Sx,Sy,Sz`，由未平滑 XYZ 跃度计算 |
| `crackle.csv` | `time,S,Cx,Cy,Cz`，由未平滑 Snap 的时间导数计算 |
| `curvature.csv` | `S,kappa,radius` |
| `path_motion.csv` | `time,S,path_velocity,path_acceleration,path_jerk,path_snap` |
| `selected_solution.csv` | 解编号、三个目标、本次动态确定的全部变量及胜出端点跃度利用率 |
| `exported_solution_objectives.csv` | Pareto 候选目标与按相同端点跃度规则重新评价所得目标的对照 |
| `constraint_status.csv` | 轮廓、动力学、Snap、共振和全局时间样条状态，约束利用率、迭代次数和连接跳变量 |
| `vibration_metrics.csv` | X/Y/Z 跃度、Snap、Crackle 指标，Snap 上限/利用率、共振频带跃度、预测响应及最终振动目标 |

新增诊断 CSV：

| 文件 | 内容 |
|---|---|
| `speed_constraints.csv` | 各速度约束（含 `v_geometric_jerk`）、`v_candidate`、`v_limit`、`v_planning_limit`、`v_planned` 和激活约束 |
| `lmsc_points.csv` | LMSC 编号、索引、弧长、速度上界、曲率和最终共享路径跃度 |
| `spu_info.csv` | SPU 边界、长度、候选/规划峰值、三段时间及首/峰/末共享跃度 |
| `lmsc_jerk_boundary_search.csv` | LMSC/非匀速峰值利用率是否可行、时间/振动变化和整体缩放变化 |
| `lmsc_jerk_boundary_candidates.csv` | 当前染色体利用率的一次评价结果、最终时间、振动和整体缩放 |
| `lmsc_jerk_boundary_values.csv` | 各 LMSC 的局部上限、请求跃度和缩放后实际跃度 |
| `spu_peak_jerk_boundary_values.csv` | 各 SPU 单点速度峰值的位置、局部幅值上限、请求负跃度和缩放后实际跃度 |
| `resonance_spectrum.csv` | XYZ 加速度/跃度功率谱、共振权重、加权跃度谱和 FRF 预测响应谱 |

图形包括原有轨迹、曲率、路径 V/A/J、XYZ V/A/J 和 Pareto 图，并新增：

```text
speed_constraints_vs_arc_length.png
xyz_snap_vs_time.png
resonance_spectrum.png
```

该图用统一弧长横轴显示各约束曲线、最终规划速度、LMSC 和 SPU 边界，可直接检查大曲率前是否开始降速，以及各区间由哪一种约束激活。`v curvature` 和 `v chord` 仅在低于 `v command`、能够实际限制指令速度的区间绘制；其余非激活部分隐藏，避免超大限速值压缩纵轴。`v planned` 最后绘制并使用加粗蓝线，保证与重合约束仍能清楚区分。

XYZ 跃度图在固定 `Ts` 样本之外额外插入七次段的精确连接时刻。零跃度连接使用连续几何三阶导数项；内部 LMSC 谷值和非匀速峰值的非零连接还加入 `q_s * path_jerk` 切向项。该显示方式使用左右段共享的精确边界跃度，避免固定采样未命中连接时刻而形成视觉假跳变，不改变导出的固定周期 CSV 或硬约束计算。

`xyz_snap_vs_time.png` 使用与硬约束相同的未滤波 Snap 数据，并绘制三轴正负 `xyz_smax` 边界。

图形统一使用 `exportgraphics` 写入 PNG。再次导出时会先关闭上一批带项目标记的图窗，避免反复运行造成图形资源累积；如果目标 PNG 正被图片查看器占用，程序会将新图保存为带时间戳的同名备用文件并给出警告，不会因此中断其余结果导出。

## 性能说明

针对“初始种群很久未完成”的主要瓶颈，NURBS 曲线求值已由递归计算全部基函数改为迭代计算每个参数处仅 `degree+1` 个非零基函数，并复用已经计算的光顺路径进行弧长重参数化，避免同一染色体重复求值 NURBS。该修改与原 Cox-de Boor 数学定义等价，不改变控制点、权重或节点编码。

另外保留：

- 输入 CSV 缓存。
- NSGA-II 候选评价在创建固定 `Ts` 时间轴之前检查 `cfg.evaluation.max_time_samples`。默认最多 `100000` 点；当前 `Ts = 0.002 s` 时对应约 `199.998 s`。超过上限的异常慢候选直接获得无效候选罚值，不再分配完整 XYZ、Snap 和共振频谱数组。该保护仅用于遗传优化和批量候选筛选，不限制 `export_selected_solution` 的最终高精度导出。
- 遗传优化默认先计算轮廓误差。误差超过约束上限 `2%` 以上的候选保留随超限量递增的轮廓罚值，但跳过速度、XYZ、Snap 和共振评价；约束边界附近的候选仍执行完整评价。该预筛选只作用于 `evaluate_objective`，不会改变 Pareto 原始目标重算或指定解导出。
- `chromosome_iter.txt` 默认每 10 代保存一次。
- Pareto 图默认每 20 代刷新一次。
- 二分搜索和动态修正均设最大迭代次数。
- 当前端点跃度利用率变量固定为 0，不执行旧版利用率内部搜索。全局时间样条直接参与每个通过轮廓预筛选候选的目标评价；它只求解一个最多 128 控制点的线性最小二乘问题，并在不能改善平滑度或超过 1.5 倍时长时立即回退。
- 每完成一代都会输出该代耗时、子代数量和不同染色体数量，便于区分“单代计算较慢”和“程序卡死”。
- 遗传算子寻找不同父代时设置了有限重选次数；种群收敛到重复染色体后，会继续执行，不会停在父代选择循环中。
- 默认变异概率为 `20%`，由 `cfg.nsga.mutation_probability` 配置。
- 每代种群和交叉父代池至少保留 `cfg.nsga.min_unique_chromosomes = 10` 个不同的动态维度染色体。环境选择会用父子联合种群中排名最好的缺失个体替换重复副本；若联合种群本身多样性不足，遗传算子会补充强制变异个体，必要时加入有界随机个体。
- 多样性补充会增加少量目标函数评价，因此单代子代数量和计算时间可能比修改前略高。

若优化被中断，可从 `chromosome_iter.txt` 中最近保存的完整代继续：

```matlab
clear functions
rehash path
resume_optimization
```

例如文件中已有 200 行且种群规模为 20，会从第 10 代之后继续。检查点每 10 代写入一次，因此中断点之后尚未保存的代需要重新计算。若希望从头运行，仍使用 `main_function`。

断点文件只适用于相同配置和相同速度规划算法。修改速度约束、目标函数或前后向可达性算法后，应使用 `main_function` 重新优化，不能继续使用旧 `chromosome_iter.txt`，否则同一种群中会混合新旧目标值。

如只需评估或导出已有解，不要重新运行完整优化，直接调用 `evaluate_trajectory` 或 `export_selected_solution`。

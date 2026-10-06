function cfg = trajectory_config(project_root)
% 集中管理三轴机床轨迹规划、NURBS 和 NSGA-II 参数。

if nargin < 1 || isempty(project_root)
    project_root = fileparts(mfilename('fullpath'));
end

cfg.project_root = project_root;
cfg.output_dir = fullfile(project_root, 'output');

cfg.files.input_csv = fullfile(project_root, '5stars.csv');
cfg.files.chromosome_initial = fullfile(project_root, 'chromosome_initial.txt');
cfg.files.chromosome_iter = fullfile(project_root, 'chromosome_iter.txt');
cfg.files.chromosome_final = fullfile(project_root, 'chromosome_final.txt');

cfg.nurbs.degree = 5;
cfg.nurbs.num_control_points = 20; % Replaced by dynamic selection below.
cfg.nurbs.eval_points = 2000;
cfg.nurbs.adaptive_control_point_distribution = true;
cfg.nurbs.corner_density_gain = 4.0;
cfg.nurbs.control_point_interpolation = 'linear';
cfg.nurbs.corner_knot_width_factor = 1.0;
% A degree-5 B-spline with simple interior knots is C4 continuous. Keep
% unit weights so chromosome weights cannot reintroduce rational oscillation.
cfg.nurbs.corner_knot_multiplicity = 1;
cfg.nurbs.force_unit_weights = true;
cfg.nurbs.bidirectional_contour_error = true;
cfg.nurbs.fit_samples_per_segment_max = 400;
cfg.nurbs.fit_regularization = 1e-10;
cfg.nurbs.fit_third_difference_regularization = 0.003;
cfg.nurbs.fit_fourth_difference_regularization = 0.001;

cfg.sampling.num_points = 600;

% Units used throughout the project: mm, s, mm/s, mm/s^2 and mm/s^3.
cfg.units.length = 'mm';
cfg.units.time = 's';

% Fixed interpolation period for all exported time-domain series.
cfg.interpolation.Ts = 0.002;

% XYZ 速度/加速度/跃度由离散轨迹数值微分得到。适当平滑切向和中间导数，
% 可避免局部曲率变化在二次微分后被放大成不真实的跃度尖峰。
cfg.smoothing.enable_xyz = true;
cfg.smoothing.window = 21;

cfg.limits.Vmax = 500;       % 路径最大速度
cfg.limits.Amax = 10000;       % 路径最大加速度
cfg.limits.Jmax = 200000;       % 路径最大跃度
cfg.limits.Acmax = 10000;      % 曲率法向加速度上限
cfg.limits.xyz_vmax = [500 500 500];
cfg.limits.xyz_amax = [10000 10000 10000];
cfg.limits.xyz_jmax = [200000 200000 200000];
cfg.limits.xyz_smax = [7000000 7000000 7000000]; % XYZ snap = d(jerk)/dt

% 最大轮廓误差约束。单位与 CSV 坐标一致；设为 Inf 可关闭约束。
% 当最大轮廓误差超过该阈值时，目标函数会增加罚值，避免超限解进入前沿。
cfg.constraints.max_contour_error = 0.15;
cfg.constraints.contour_error_penalty = 1e6;
cfg.constraints.invalid_candidate_penalty = 1e12;
cfg.constraints.max_chord_error = 0.01;
% true：xyz_smax 参与局部修正、整体时间缩放、可行性和罚值判定。
% false：仍计算、绘制和导出 Snap，但不把 xyz_smax 作为硬约束。
cfg.constraints.enable_xyz_snap_limit = false;

% Memory guard used only while NSGA-II evaluates candidate chromosomes.
% A candidate requiring more fixed-Ts samples is not competitive for the
% time objective and is penalized before allocating its full XYZ arrays.
% Selected-solution export is intentionally uncapped.
cfg.evaluation.max_time_samples = 100000;
cfg.evaluation.enable_contour_prefilter = true;
cfg.evaluation.contour_prefilter_relative_margin = 0.02;

% Vibration uses the weighted L2 norm of X/Y/Z jerk RMS. The XYZ jerk
% peak remains available as an exported diagnostic but is disabled in the
% objective by setting its weight to zero.
cfg.objectives.vibration_axis_weights = [1 1 1];
cfg.objectives.vibration_peak_weight = 0;
% Penalize rapid changes of XYZ jerk while keeping the term in jerk units:
% equivalent_snap = Ts * RMS(d(jerk)/dt).
cfg.objectives.vibration_snap_weight = 0.10;
% Penalize changes of Snap while retaining jerk units:
% equivalent_crackle = Ts^2 * RMS(d(Snap)/dt).
cfg.objectives.vibration_crackle_weight = 0.02;

% Machine-resonance model. Frequencies, damping ratios and modal gains are
% specified per X/Y/Z axis. Leave the mode lists empty until measured modal
% data are available; selected-solution exports still include FFT spectra.
cfg.resonance.enabled = true;
cfg.resonance.enable_diagnostics = true;
cfg.resonance.mode_frequencies_hz = {45, 3, 78};
cfg.resonance.damping_ratios = {0.03, 0.025, 0.02};
cfg.resonance.modal_gains = {1, 1.2, 0.8};
cfg.resonance.axis_weights = [1 1 1];
cfg.resonance.spectral_peak_gain = 10;
cfg.resonance.objective_weight = 1;
cfg.resonance.response_rms_limits = [Inf Inf Inf];
cfg.resonance.penalty_weight = 1e6;

% LMSC/SPU seventh-order speed-planning parameters.
cfg.speed.enable_lmsc = true;
cfg.speed.enable_spu = true;
cfg.speed.enable_chord_constraint = false;
cfg.speed.lmsc_velocity_threshold = 5.0;
cfg.speed.lmsc_prominence_threshold = 5.0;
cfg.speed.lmsc_transition_distance_factor = 1.05;
cfg.speed.lmsc_min_sample_intervals = 2.0;
cfg.speed.lmsc_merge_distance_factor = 1.5;
cfg.speed.constraint_correction_radius_samples = 3.0;
cfg.speed.minimum_transition_time = 1e-3;
cfg.speed.velocity_tolerance = 1e-6;
cfg.speed.distance_tolerance = 1e-8;
cfg.speed.binary_search_max_iterations = 60;
cfg.speed.max_constraint_iterations = 8;
cfg.speed.enable_global_time_scaling = true;
cfg.speed.max_global_time_scaling_iterations = 3;
cfg.speed.minimize_global_time_scale = true;
cfg.speed.global_time_scale_refinement_iterations = 12;
cfg.speed.global_time_scale_refinement_tolerance = 1e-4;
cfg.speed.constraint_safety_factor = 0.98;
cfg.speed.internal_constraint_check_samples = 65;
cfg.speed.epsilon_kappa = 1e-9;
cfg.speed.process_speed_limit = Inf;
cfg.speed.forced_limit_indices = [];

% Optional second smoothing layer. The existing LMSC/SPU seventh-order
% planner supplies a feasible reference law; a global degree-7 B-spline
% then fits S(t), enforces zero V/A/J/Snap at both ends, and is uniformly
% time-scaled until all active path and XYZ constraints are satisfied.
cfg.time_smoothing.enabled = true;
% Keep true so Pareto objective values and the exported simulation are
% evaluated from the same globally smoothed time law.
cfg.time_smoothing.apply_during_objective = true;
cfg.time_smoothing.degree = 7;
cfg.time_smoothing.min_control_points = 64;
cfg.time_smoothing.max_control_points = 128;
cfg.time_smoothing.control_points_per_lmsc = 16;
cfg.time_smoothing.fit_samples = 401;
cfg.time_smoothing.fourth_difference_regularization = 0.005;
cfg.time_smoothing.fifth_difference_regularization = 0.001;
cfg.time_smoothing.max_scaling_iterations = 4;
cfg.time_smoothing.maximum_duration_factor = 1.5;
cfg.time_smoothing.minimum_smoothness_improvement = 1e-4;
cfg.time_smoothing.fallback_to_septic = true;

% Shared endpoint-jerk optimization for every fully evaluated chromosome.
% Internal LMSC minima use positive jerk and direct transition-transition
% speed peaks use negative jerk. Global endpoints and transition/constant
% joins retain zero jerk so that the complete trajectory stays continuous.
% A local feasible magnitude is computed at each eligible join. One global
% utilization ratio is stored as the last continuous chromosome variable;
% its trajectory supplies the time and vibration objectives to NSGA-II.
cfg.lmsc_jerk_boundary.enabled = false;
cfg.lmsc_jerk_boundary.evaluate_during_objective = false;
cfg.lmsc_jerk_boundary.local_bound_safety_factor = 0.98;
cfg.lmsc_jerk_boundary.maximum_vibration_increase = 1000;
% false: let the vibration Pareto objective and machine hard limits select
% the utilization. true: additionally reject a candidate whose vibration
% exceeds the zero-jerk baseline by maximum_vibration_increase.
cfg.lmsc_jerk_boundary.enforce_vibration_increase_limit = false;
cfg.lmsc_jerk_boundary.minimum_duration_fraction = 0.4;
cfg.lmsc_jerk_boundary.maximum_duration_factor = 2.0;
cfg.lmsc_jerk_boundary.time_search_samples = 9;
cfg.lmsc_jerk_boundary.time_refinement_iterations = 6;
cfg.lmsc_jerk_boundary.peak_search_samples = 11;
% A constant-speed plateau shorter than this is numerical roundoff from
% the two transition distances, not a physical dwell. Treat its two
% transition boundaries as one shared nonzero-jerk speed peak.
cfg.lmsc_jerk_boundary.numerical_peak_plateau_time_tolerance = ...
    1e-6 * cfg.interpolation.Ts;

cfg.nsga.pop = 20;
cfg.nsga.gen = 20;
cfg.nsga.mutation_probability = 0.2;
cfg.nsga.min_unique_chromosomes = 10;
cfg.nsga.save_interval = 10;
cfg.nsga.plot_interval = 20;

% Determine the chromosome dimensions once from path complexity and the
% requested contour-error tolerance. They remain fixed during NSGA-II.
cfg.dynamic_dimensions.enabled = true;
cfg.dynamic_dimensions.min_control_points = 12;
cfg.dynamic_dimensions.max_control_points = 60;
cfg.dynamic_dimensions.control_point_step = 4;
cfg.dynamic_dimensions.probe_eval_points = 161;
cfg.dynamic_dimensions.baseline_error_fraction = 1.0;
cfg.dynamic_dimensions.fallback_control_points = 20;
cfg.dynamic_dimensions.min_speed_nodes = 7;
cfg.dynamic_dimensions.max_speed_nodes = 25;
cfg.dynamic_dimensions.control_points_per_speed_node = 3.0;
cfg.dynamic_dimensions.speed_nodes_per_feature = 2;
cfg.dynamic_dimensions.feature_angle_threshold_deg = 8.0;
cfg.dynamic_dimensions.fallback_speed_nodes = 15;
cfg = configure_dynamic_dimensions(cfg);

cfg.nsga.num_weights = cfg.nurbs.num_control_points;
cfg.nsga.num_jerk_utilization_variables = 1;
cfg.nsga.num_variables = cfg.nsga.num_weights + cfg.nsga.num_time_nodes + ...
    cfg.nsga.num_jerk_utilization_variables;

% 前沿解选择。先运行一次生成 output/pareto_candidates.csv，再把需要导出的
% solution_id 填到这里；保持 NaN 时只列出候选解，不导出单个解。
cfg.selection.selected_solution_id = NaN;

% Only expose Pareto candidates whose raw traversal time is within this
% multiple of the shortest raw traversal time. Inf disables the filter.
cfg.pareto.max_time_ratio = 10;

if cfg.nurbs.force_unit_weights
    cfg.bounds.weights_lb = ones(1, cfg.nsga.num_weights);
    cfg.bounds.weights_ub = ones(1, cfg.nsga.num_weights);
else
    cfg.bounds.weights_lb = 0.1 * ones(1, cfg.nsga.num_weights);
    cfg.bounds.weights_ub = 1.0 * ones(1, cfg.nsga.num_weights);
end
cfg.bounds.time_lb = 0.6 * ones(1, cfg.nsga.num_time_nodes);
cfg.bounds.time_ub = 1.8 * ones(1, cfg.nsga.num_time_nodes);
cfg.bounds.jerk_utilization_lb = 0;
cfg.bounds.jerk_utilization_ub = double(cfg.lmsc_jerk_boundary.enabled && ...
    cfg.lmsc_jerk_boundary.evaluate_during_objective);
cfg.nsga.reference_decision = [ones(1, cfg.nsga.num_weights), ...
    ones(1, cfg.nsga.num_time_nodes), 0];
end

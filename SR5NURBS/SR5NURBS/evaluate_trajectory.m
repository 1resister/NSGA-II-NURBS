function result = evaluate_trajectory( ...
        x, cfg, keep_series, maximum_raw_time, enable_contour_prefilter)
% 根据一条染色体完成完整 XYZ 轨迹评价。

if nargin < 3
    keep_series = false;
end
if nargin < 4
    maximum_raw_time = Inf;
end
if nargin < 5
    enable_contour_prefilter = false;
end
if ~keep_series
    maximum_raw_time = min(maximum_raw_time, configured_batch_time_limit(cfg));
end

raw_points = load_path_points(cfg.files.input_csv);
if numel(x) ~= cfg.nsga.num_variables
    error('evaluate_trajectory:DecisionDimensionMismatch', ...
        'Expected %d decision variables but received %d.', ...
        cfg.nsga.num_variables, numel(x));
end
weight_end = cfg.nsga.num_weights;
time_end = weight_end + cfg.nsga.num_time_nodes;
weights = x(1:weight_end);
time_scale = x(weight_end + 1:time_end);
jerk_utilization = x(time_end + 1);

[smooth_path, curve, u] = nurbs_global_smoothing(raw_points, weights, time_scale, cfg);
bidirectional_error = isfield(cfg.nurbs, ...
    'bidirectional_contour_error') && ...
    cfg.nurbs.bidirectional_contour_error;
err = nurbs_contour_error(raw_points, curve, u, bidirectional_error);
if enable_contour_prefilter && should_prefilter_contour(err, cfg)
    result = contour_prefilter_result(err, cfg);
    return
end
arc = arc_length_parameterization(curve, u, cfg.sampling.num_points, smooth_path);
geom = compute_curvature(arc.sample_points, arc.S, arc.derivatives);
constraints = generate_velocity_constraints(geom.kappa, geom.tangent, cfg, ...
    arc.S, time_scale, arc.derivatives);
planning_cfg = cfg;
if ~keep_series && isfield(planning_cfg, 'time_smoothing') && ...
        isfield(planning_cfg.time_smoothing, 'apply_during_objective') && ...
        ~planning_cfg.time_smoothing.apply_during_objective
    planning_cfg.time_smoothing.enabled = false;
end
[constraints, profile, motion, xyz, speed_diagnostics] = plan_speed_profile( ...
    arc.sample_points, arc.S, geom, constraints, time_scale, planning_cfg, ...
    arc.derivatives, maximum_raw_time);

vibration_metrics = xyz_vibration_metric( ...
    xyz.jerk, cfg, motion.time, xyz.acceleration, keep_series);
vibration = vibration_metrics.objective_value;

raw_time = motion.time(end);
raw_vibration = vibration;
[time_obj, error_obj, vibration_obj, constraint] = apply_contour_error_constraint( ...
    raw_time, err, raw_vibration, cfg);

dynamic_ratio = max([motion.validation.max_velocity_ratio, ...
    motion.validation.max_acceleration_ratio, motion.validation.max_jerk_ratio, ...
    xyz.max_velocity_ratio, xyz.max_acceleration_ratio, xyz.max_jerk_ratio, ...
    xyz.max_enforced_snap_ratio]);
dynamic_violation = max(dynamic_ratio - 1, 0);
if ~speed_diagnostics.feasible
    dynamic_violation = max(dynamic_violation, 1);
end
dynamic_penalty = cfg.constraints.contour_error_penalty * dynamic_violation;
time_obj = time_obj + dynamic_penalty;
error_obj = error_obj + dynamic_penalty;
vibration_obj = vibration_obj + dynamic_penalty;
constraint.dynamic_ratio = dynamic_ratio;
constraint.dynamic_violation = dynamic_violation;
constraint.dynamic_penalty = dynamic_penalty;
constraint.path_feasible = speed_diagnostics.path_feasible;
constraint.xyz_feasible = speed_diagnostics.xyz_feasible;
constraint.xyz_snap_limit_enabled = xyz.snap_limit_enabled;
constraint.xyz_snap_limit = cfg.limits.xyz_smax;
constraint.xyz_snap_peak = vibration_metrics.axis_snap_peak;
constraint.xyz_snap_axis_ratio = vibration_metrics.axis_snap_ratio;
constraint.xyz_snap_ratio = xyz.max_snap_ratio;
constraint.xyz_snap_violation = max(xyz.max_enforced_snap_ratio - 1, 0);
constraint.xyz_snap_feasible = xyz.max_enforced_snap_ratio <= ...
    1 + cfg.speed.velocity_tolerance;
constraint.dynamic_feasible = speed_diagnostics.feasible && ...
    dynamic_violation <= cfg.speed.velocity_tolerance;
constraint.feasible = constraint.feasible && constraint.dynamic_feasible;
constraint.iterations = speed_diagnostics.iterations;
constraint.global_time_scaling_enabled = ...
    speed_diagnostics.global_time_scaling_enabled;
constraint.global_time_scaling_iterations = ...
    speed_diagnostics.global_time_scaling_iterations;
constraint.global_time_scale = speed_diagnostics.global_time_scale;
constraint.global_time_scale_minimized = ...
    speed_diagnostics.global_time_scale_minimized;
constraint.global_time_scale_refinement_iterations = ...
    speed_diagnostics.global_time_scale_refinement_iterations;
constraint.global_time_scale_lower_bound = ...
    speed_diagnostics.global_time_scale_lower_bound;
constraint.time_smoothing_enabled = ...
    speed_diagnostics.time_smoothing_enabled;
constraint.time_smoothing_accepted = ...
    speed_diagnostics.time_smoothing_accepted;
constraint.time_smoothing_duration_scale = ...
    speed_diagnostics.time_smoothing_duration_scale;
constraint.time_smoothing_attempted_duration_scale = ...
    speed_diagnostics.time_smoothing.attempted_duration_scale;
constraint.time_smoothing_baseline_score = ...
    speed_diagnostics.time_smoothing.baseline_smoothness;
constraint.time_smoothing_final_score = ...
    speed_diagnostics.time_smoothing.final_smoothness;
constraint.time_smoothing_improvement = ...
    speed_diagnostics.time_smoothing.smoothness_improvement;
constraint.time_smoothing_reason = ...
    speed_diagnostics.time_smoothing.reason;

resonance = vibration_metrics.resonance;
resonance_violation = max(resonance.max_response_ratio - 1, 0);
resonance_penalty_weight = cfg.resonance.penalty_weight;
if ~isscalar(resonance_penalty_weight) || ...
        ~isfinite(resonance_penalty_weight) || resonance_penalty_weight < 0
    error('cfg.resonance.penalty_weight must be a finite nonnegative scalar.');
end
resonance_penalty = resonance_penalty_weight * resonance_violation;
time_obj = time_obj + resonance_penalty;
error_obj = error_obj + resonance_penalty;
vibration_obj = vibration_obj + resonance_penalty;
constraint.resonance_enabled = resonance.enabled;
constraint.resonance_configured = resonance.configured;
constraint.resonance_ratio = resonance.max_response_ratio;
constraint.resonance_violation = resonance_violation;
constraint.resonance_penalty = resonance_penalty;
constraint.resonance_feasible = resonance.feasible;
constraint.feasible = constraint.feasible && resonance.feasible;

result.objectives.time = time_obj;
result.objectives.error = error_obj;
result.objectives.vibration = vibration_obj;
result.raw_objectives.time = raw_time;
result.raw_objectives.error = err;
result.raw_objectives.vibration = raw_vibration;
result.raw_objectives.vibration_xyz = vibration_metrics.axis_rms;
result.vibration_metrics = vibration_metrics;
result.constraint = constraint;
result.contour_prefiltered = false;

% The shared-jerk decision evaluation needs the complete baseline topology
% even during a lightweight NSGA-II objective evaluation. Populate it once,
% evaluate the chromosome utilization, then remove the heavy arrays when
% the caller did not request time series.
result.raw_points = raw_points;
result.smooth_path = smooth_path;
result.curve = curve;
result.u = u;
result.arc = arc;
result.geom = geom;
result.constraints = constraints;
result.profile = profile;
result.lmsc = motion.lmsc;
result.spus = motion.spus;
result.motion = motion;
result.xyz = xyz;

if shared_jerk_objective_enabled(cfg) && constraint.violation <= 0
    [result, jerk_report] = refine_selected_lmsc_jerk_boundary( ...
        result, cfg, keep_series, jerk_utilization, maximum_raw_time);
    if ~jerk_report.decision_variable_feasible
        failure_reason = string(jerk_report.reason);
        if contains(failure_reason, ...
                "plan_speed_profile:MaximumTimeExceeded")
            error('plan_speed_profile:MaximumTimeExceeded', ...
                ['Endpoint-jerk utilization %.6g exceeded the active ' ...
                'maximum raw time %.9g s.'], ...
                jerk_utilization, maximum_raw_time);
        end
        error('evaluate_trajectory:InfeasibleJerkUtilization', ...
            'Endpoint-jerk utilization %.6g is infeasible: %s.', ...
            jerk_utilization, char(failure_reason));
    end
end

if ~keep_series
    result = remove_heavy_series(result);
end
end

function enabled = shared_jerk_objective_enabled(cfg)
enabled = isfield(cfg, 'lmsc_jerk_boundary') && ...
    isfield(cfg.lmsc_jerk_boundary, 'enabled') && ...
    logical(cfg.lmsc_jerk_boundary.enabled) && ...
    isfield(cfg.lmsc_jerk_boundary, 'evaluate_during_objective') && ...
    logical(cfg.lmsc_jerk_boundary.evaluate_during_objective);
end

function result = remove_heavy_series(result)
heavy_fields = {'raw_points','smooth_path','curve','u','arc','geom', ...
    'constraints','profile','lmsc','spus','motion','xyz'};
present = heavy_fields(isfield(result, heavy_fields));
if ~isempty(present)
    result = rmfield(result, present);
end
end

function prefilter = should_prefilter_contour(err, cfg)
prefilter = false;
if ~isfield(cfg, 'evaluation') || ...
        ~isfield(cfg.evaluation, 'enable_contour_prefilter') || ...
        ~logical(cfg.evaluation.enable_contour_prefilter)
    return
end
maximum_error = cfg.constraints.max_contour_error;
if ~isfinite(maximum_error)
    return
end
margin = cfg.evaluation.contour_prefilter_relative_margin;
if ~isscalar(margin) || ~isfinite(margin) || margin < 0
    error(['cfg.evaluation.contour_prefilter_relative_margin must be a ' ...
        'finite nonnegative scalar.']);
end
prefilter = err > maximum_error * (1 + margin);
end

function result = contour_prefilter_result(err, cfg)
% Preserve a smooth selection pressure toward the feasible contour boundary
% without evaluating an already dominated candidate's time-domain dynamics.
[time_obj, error_obj, vibration_obj, constraint] = ...
    apply_contour_error_constraint(0, err, 0, cfg);
result.objectives.time = time_obj;
result.objectives.error = error_obj;
result.objectives.vibration = vibration_obj;
result.raw_objectives.time = NaN;
result.raw_objectives.error = err;
result.raw_objectives.vibration = NaN;
result.raw_objectives.vibration_xyz = nan(1, 3);
result.constraint = constraint;
result.contour_prefiltered = true;
end

function maximum_raw_time = configured_batch_time_limit(cfg)
maximum_raw_time = Inf;
if ~isfield(cfg, 'evaluation') || ...
        ~isfield(cfg.evaluation, 'max_time_samples')
    return
end
maximum_samples = cfg.evaluation.max_time_samples;
if ~isscalar(maximum_samples) || isnan(maximum_samples) || ...
        maximum_samples < 2 || ...
        (~isinf(maximum_samples) && maximum_samples ~= floor(maximum_samples))
    error('cfg.evaluation.max_time_samples must be an integer >= 2 or Inf.');
end
if ~isinf(maximum_samples)
    maximum_raw_time = (maximum_samples - 1) * cfg.interpolation.Ts;
end
end

function points = load_path_points(input_csv)
% 缓存 CSV 离散点，避免 NSGA-II 目标函数评估时反复读盘。
persistent cached_file cached_points cached_stamp
info = dir(input_csv);
if isempty(info)
    error('找不到输入 CSV 文件：%s', input_csv);
end

if isempty(cached_points) || ~strcmp(cached_file, input_csv) || cached_stamp ~= info.datenum
    cached_file = input_csv;
    cached_stamp = info.datenum;
    cached_points = preprocess_points(readmatrix(input_csv));
end
points = cached_points;
end

function points = preprocess_points(points)
% 清理 CSV 点数据，保证格式为 N x 3 且去除相邻重复点。
points = points(:, 1:3);
points = points(all(isfinite(points), 2), :);
if size(points, 1) < 2
    error('输入 CSV 至少需要包含两个三维路径点。');
end
d = [true; vecnorm(diff(points), 2, 2) > 1e-10];
points = points(d, :);
end

function err = contour_error(raw_points, smooth_path) %#ok<DEFNU>
% 最大轮廓误差：每个原始离散点到 NURBS 光顺轨迹的最近距离最大值。
block_size = 1000;
dist_min = inf(size(raw_points, 1), 1);
for start_id = 1:block_size:size(smooth_path, 1)
    stop_id = min(start_id + block_size - 1, size(smooth_path, 1));
    block = smooth_path(start_id:stop_id, :);
    dx = raw_points(:, 1) - block(:, 1)';
    dy = raw_points(:, 2) - block(:, 2)';
    dz = raw_points(:, 3) - block(:, 3)';
    dist_min = min(dist_min, min(sqrt(dx.^2 + dy.^2 + dz.^2), [], 2));
end
err = max(dist_min);
end

function [time_obj, error_obj, vibration_obj, constraint] = apply_contour_error_constraint(time, err, vibration, cfg)
% 最大轮廓误差约束。
%
% NSGA-II 原框架没有显式约束接口，因此采用罚函数法：
% 当 err > max_contour_error 时，对三个目标同时增加罚值，使超限解
% 在非支配排序中处于明显劣势。
if ~isfield(cfg, 'constraints') || ~isfield(cfg.constraints, 'max_contour_error')
    max_error = Inf;
else
    max_error = cfg.constraints.max_contour_error;
end
if ~isfield(cfg, 'constraints') || ~isfield(cfg.constraints, 'contour_error_penalty')
    penalty_weight = 1e6;
else
    penalty_weight = cfg.constraints.contour_error_penalty;
end

violation = max(0, err - max_error);
if isinf(max_error)
    violation_ratio = 0;
else
    violation_ratio = violation / max(max_error, eps);
end
penalty = penalty_weight * violation_ratio;

time_obj = time + penalty;
error_obj = err + penalty;
vibration_obj = vibration + penalty;

constraint.max_contour_error = max_error;
constraint.contour_error = err;
constraint.violation = violation;
constraint.violation_ratio = violation_ratio;
constraint.penalty = penalty;
constraint.feasible = violation <= 0;
end

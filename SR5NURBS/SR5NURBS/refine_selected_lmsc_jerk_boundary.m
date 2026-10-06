function [result, report] = refine_selected_lmsc_jerk_boundary( ...
        result, cfg, keep_spectrum, utilization, maximum_time)
% Apply the chromosome's shared endpoint-jerk utilization.
%
% Each internal LMSC speed minimum shares one positive path jerk. Each SPU
% that has a single nonconstant speed maximum shares one negative peak
% jerk. Velocity, acceleration and jerk remain continuous at both types of
% boundary. Local magnitudes are limited by path Jmax and the XYZ jerk
% margin. The one global utilization is a continuous NSGA-II decision
% variable; this function evaluates that value once and returns its actual
% time and vibration objectives. It never silently falls back to zero when
% a nonzero decision is infeasible.

if nargin < 3
    keep_spectrum = true;
end
if nargin < 4 || isempty(utilization)
    error('A chromosome endpoint-jerk utilization is required.');
end
if nargin < 5 || isempty(maximum_time)
    maximum_time = Inf;
end

validate_configuration(cfg);
baseline = result;
internal_count = max(numel(result.lmsc.indices) - 2, 0);
report = initial_report(result, cfg, internal_count);
report.selected_utilization = utilization;

lower_bound = cfg.bounds.jerk_utilization_lb;
upper_bound = cfg.bounds.jerk_utilization_ub;
if ~isscalar(utilization) || ~isfinite(utilization) || ...
        utilization < lower_bound || utilization > upper_bound
    report.reason = "decision_variable_out_of_bounds";
    report.decision_variable_feasible = false;
    result.lmsc_jerk_boundary_search = report;
    return
end

if ~logical(cfg.lmsc_jerk_boundary.enabled)
    report.reason = "disabled";
    report.trials.utilization(1) = utilization;
    report.trials.reason(1) = report.reason;
    result.lmsc_jerk_boundary_search = report;
    return
end
[local_upper, active] = local_feasible_upper_bounds(result, cfg);
[peak_upper, peak_active, peak_S] = ...
    spu_peak_feasible_upper_bounds(result, cfg);
report.local_maximum_jerk = local_upper;
report.active_internal_lmsc = active;
report.local_maximum_peak_jerk = peak_upper;
report.active_spu_peak = peak_active;
report.spu_peak_S = peak_S;
jerk_tolerance = max(cfg.limits.Jmax * 1e-12, ...
    cfg.speed.velocity_tolerance);
if utilization <= jerk_tolerance
    report.reason = "zero_jerk_decision";
    report.trials.utilization(1) = utilization;
    report.trials.reason(1) = report.reason;
    result.lmsc_jerk_boundary_search = report;
    return
end
if ~any(active) && ~any(peak_active)
    report.reason = "no_local_jerk_margin";
    report.trials.utilization(1) = utilization;
    report.trials.reason(1) = report.reason;
    result.lmsc_jerk_boundary_search = report;
    return
end

requested_jerk = [0; utilization .* local_upper; 0];
requested_peak_jerk = -utilization .* peak_upper;
[final_result, final_ok, final_reason] = evaluate_candidate( ...
    baseline, requested_jerk, requested_peak_jerk, cfg, keep_spectrum, ...
    maximum_time, false);
report.evaluations = report.evaluations + 1;
if ~final_ok || ~final_result.constraint.feasible
    report.reason = "decision_variable_evaluation_failed:" + final_reason;
    report.decision_variable_feasible = false;
    report.rejected_evaluations = report.rejected_evaluations + 1;
    report.trials = failed_trial(utilization, final_reason);
    result.lmsc_jerk_boundary_search = report;
    return
end

vibration_ok = true;
if cfg.lmsc_jerk_boundary.enforce_vibration_increase_limit
    maximum_vibration = baseline.raw_objectives.vibration * ...
        (1 + cfg.lmsc_jerk_boundary.maximum_vibration_increase);
    vibration_ok = final_result.raw_objectives.vibration <= ...
        maximum_vibration + 10 * eps(max(maximum_vibration, 1));
end
if ~vibration_ok
    report.reason = "vibration_increase_limit";
    report.decision_variable_feasible = false;
    report.rejected_evaluations = report.rejected_evaluations + 1;
    report.trials = candidate_trial(utilization, final_result, false, ...
        report.reason);
    result.lmsc_jerk_boundary_search = report;
    return
end

report.applied = true;
report.reason = "decision_variable_accepted";
report.requested_unscaled_lmsc_jerk = requested_jerk;
report.applied_scaled_lmsc_jerk = final_result.motion.lmsc_jerk;
report.requested_unscaled_spu_peak_jerk = requested_peak_jerk;
report.applied_scaled_spu_peak_jerk = ...
    final_result.motion.spu_peak_jerk;
report.applied_spu_peak_S = actual_spu_peak_positions(final_result.motion);
report.refined_time = final_result.raw_objectives.time;
report.time_gain = report.baseline_time - report.refined_time;
report.time_gain_ratio = report.time_gain / max(report.baseline_time, eps);
report.refined_vibration = final_result.raw_objectives.vibration;
report.vibration_change_ratio = ...
    (report.refined_vibration - report.baseline_vibration) / ...
    max(abs(report.baseline_vibration), eps);
report.refined_global_time_scale = ...
    final_result.constraint.global_time_scale;
report.best_trial_time = report.refined_time;
report.best_trial_vibration = report.refined_vibration;
report.best_trial_global_time_scale = report.refined_global_time_scale;
report.trials = candidate_trial(utilization, final_result, true, ...
    report.reason);
result = final_result;
result.lmsc_jerk_boundary_search = report;
end

function trial = failed_trial(utilization, reason)
trial = table(utilization, false, false, NaN, NaN, NaN, string(reason), ...
    'VariableNames', {'utilization','feasible','acceptable','final_time', ...
    'vibration','global_time_scale','reason'});
end

function trial = candidate_trial(utilization, candidate, acceptable, reason)
trial = table(utilization, candidate.constraint.feasible, acceptable, ...
    candidate.raw_objectives.time, candidate.raw_objectives.vibration, ...
    candidate.constraint.global_time_scale, string(reason), ...
    'VariableNames', {'utilization','feasible','acceptable','final_time', ...
    'vibration','global_time_scale','reason'});
end

function report = initial_report(result, cfg, internal_count)
report.enabled = logical(cfg.lmsc_jerk_boundary.enabled);
report.evaluated_during_objective = ...
    logical(cfg.lmsc_jerk_boundary.evaluate_during_objective);
report.decision_variable_feasible = true;
report.applied = false;
report.reason = "not_run";
report.internal_lmsc_count = internal_count;
report.evaluations = 0;
report.rejected_evaluations = 0;
report.baseline_time = result.raw_objectives.time;
report.refined_time = result.raw_objectives.time;
report.time_gain = 0;
report.time_gain_ratio = 0;
report.baseline_vibration = result.raw_objectives.vibration;
report.refined_vibration = result.raw_objectives.vibration;
report.vibration_change_ratio = 0;
report.baseline_global_time_scale = result.constraint.global_time_scale;
report.refined_global_time_scale = result.constraint.global_time_scale;
report.local_maximum_jerk = zeros(internal_count, 1);
report.active_internal_lmsc = false(internal_count, 1);
spu_count = numel(result.spus);
report.local_maximum_peak_jerk = zeros(spu_count, 1);
report.active_spu_peak = false(spu_count, 1);
report.spu_peak_S = nan(spu_count, 1);
report.applied_spu_peak_S = nan(spu_count, 1);
report.selected_utilization = 0;
report.requested_unscaled_lmsc_jerk = ...
    zeros(internal_count + 2, 1);
report.applied_scaled_lmsc_jerk = zeros(internal_count + 2, 1);
report.requested_unscaled_spu_peak_jerk = zeros(spu_count, 1);
report.applied_scaled_spu_peak_jerk = zeros(spu_count, 1);
report.best_trial_time = result.raw_objectives.time;
report.best_trial_vibration = result.raw_objectives.vibration;
report.best_trial_global_time_scale = result.constraint.global_time_scale;
report.trials = table(0, result.constraint.feasible, ...
    result.constraint.feasible, result.raw_objectives.time, ...
    result.raw_objectives.vibration, result.constraint.global_time_scale, ...
    "zero_jerk_baseline", ...
    'VariableNames', {'utilization','feasible','acceptable','final_time', ...
    'vibration','global_time_scale','reason'});
end

function [upper_bound, active] = local_feasible_upper_bounds(result, cfg)
internal_ids = (2:numel(result.lmsc.indices)-1)';
count = numel(internal_ids);
upper_bound = zeros(count, 1);
active = false(count, 1);
if count == 0 || ~isfield(result.arc, 'derivatives') || ...
        ~isfield(result.arc.derivatives, 'q_s') || ...
        ~isfield(result.arc.derivatives, 'q_sss')
    return
end

lmsc_S = result.lmsc.S(internal_ids);
q_s = interp1(result.arc.S, result.arc.derivatives.q_s, ...
    lmsc_S, 'pchip', 'extrap');
q_sss = interp1(result.arc.S, result.arc.derivatives.q_sss, ...
    lmsc_S, 'pchip', 'extrap');
boundary_speed = result.profile.V(result.lmsc.indices(internal_ids));
reference_scale = max(result.constraint.global_time_scale, 1);
axis_limit = cfg.limits.xyz_jmax(:)' .* reference_scale^3;
velocity_tolerance = cfg.speed.velocity_tolerance;

for i = 1:count
    lmsc_id = internal_ids(i);
    left_peak = result.spus(lmsc_id - 1).candidate_peak_speed;
    right_peak = result.spus(lmsc_id).candidate_peak_speed;
    if left_peak <= boundary_speed(i) + velocity_tolerance || ...
            right_peak <= boundary_speed(i) + velocity_tolerance
        continue
    end

    upper_bound(i) = signed_local_jerk_upper_bound( ...
        q_s(i, :), q_sss(i, :), boundary_speed(i), 1, ...
        axis_limit, cfg);
    active(i) = upper_bound(i) > velocity_tolerance;
end
end

function [upper_bound, active, peak_S] = ...
        spu_peak_feasible_upper_bounds(result, cfg)
spu_count = numel(result.spus);
upper_bound = zeros(spu_count, 1);
active = false(spu_count, 1);
peak_S = nan(spu_count, 1);
if spu_count == 0 || ~isfield(result.arc, 'derivatives') || ...
        ~isfield(result.arc.derivatives, 'q_s') || ...
        ~isfield(result.arc.derivatives, 'q_sss')
    return
end

segments = result.motion.segments;
reference_scale = max(result.constraint.global_time_scale, 1);
axis_limit = cfg.limits.xyz_jmax(:)' .* reference_scale^3;
velocity_tolerance = cfg.speed.velocity_tolerance;
for spu_id = 1:spu_count
    segment_ids = find([segments.spu_id] == spu_id);
    [left, right, numerical_plateau, pair_ok] = ...
        peak_transition_pair(segments, segment_ids, cfg);
    if ~pair_ok
        continue
    end
    if left.type ~= "transition" || right.type ~= "transition" || ...
            left.v1 <= left.v0 + velocity_tolerance || ...
            right.v0 <= right.v1 + velocity_tolerance
        continue
    end
    join_tolerance = max(cfg.speed.distance_tolerance, ...
        10 * eps(max(abs([left.s1, right.s0, 1]))));
    if numerical_plateau
        join_tolerance = max(join_tolerance, ...
            abs(right.s0 - left.s1) + ...
            10 * eps(max(abs([left.s1, right.s0, 1]))));
    end
    if abs(left.s1 - right.s0) > join_tolerance
        continue
    end

    peak_S(spu_id) = 0.5 * (left.s1 + right.s0);
    peak_speed = 0.5 * (left.v1 + right.v0) * reference_scale;
    q_s = interp1(result.arc.S, result.arc.derivatives.q_s, ...
        peak_S(spu_id), 'pchip', 'extrap');
    q_sss = interp1(result.arc.S, result.arc.derivatives.q_sss, ...
        peak_S(spu_id), 'pchip', 'extrap');
    % A speed maximum requires negative path jerk. The optimized variable
    % is its nonnegative magnitude, hence direction = -1.
    upper_bound(spu_id) = signed_local_jerk_upper_bound( ...
        q_s, q_sss, peak_speed, -1, axis_limit, cfg);
    active(spu_id) = upper_bound(spu_id) > velocity_tolerance;
end
end

function [left, right, numerical_plateau, feasible] = ...
        peak_transition_pair(segments, segment_ids, cfg)
left = struct();
right = struct();
numerical_plateau = false;
feasible = false;
if numel(segment_ids) == 2
    left = segments(segment_ids(1));
    right = segments(segment_ids(2));
    feasible = left.type == "transition" && right.type == "transition";
    return
end
if numel(segment_ids) ~= 3
    return
end

left = segments(segment_ids(1));
plateau = segments(segment_ids(2));
right = segments(segment_ids(3));
plateau_time = max(plateau.t1 - plateau.t0, 0);
time_tolerance = ...
    cfg.lmsc_jerk_boundary.numerical_peak_plateau_time_tolerance;
numerical_plateau = left.type == "transition" && ...
    plateau.type == "constant" && right.type == "transition" && ...
    plateau_time <= time_tolerance;
feasible = numerical_plateau;
end

function upper = signed_local_jerk_upper_bound( ...
        q_s, q_sss, speed, direction, axis_limit, cfg)
lower = 0;
upper = cfg.limits.Jmax;
geometric_jerk = q_sss .* speed^3;
velocity_tolerance = cfg.speed.velocity_tolerance;
for axis_id = 1:3
    coefficient = direction * q_s(axis_id);
    if abs(coefficient) <= 1e-12
        if abs(geometric_jerk(axis_id)) > axis_limit(axis_id) * ...
                (1 + velocity_tolerance)
            upper = 0;
            return
        end
        continue
    end
    endpoint_a = (-axis_limit(axis_id) - ...
        geometric_jerk(axis_id)) / coefficient;
    endpoint_b = ( axis_limit(axis_id) - ...
        geometric_jerk(axis_id)) / coefficient;
    lower = max(lower, min(endpoint_a, endpoint_b));
    upper = min(upper, max(endpoint_a, endpoint_b));
end
% The zero-jerk baseline is feasible at the reference scale. Preserve zero
% in the interval against small interpolation roundoff.
lower = min(lower, 0);
if upper > max(lower, velocity_tolerance)
    upper = cfg.lmsc_jerk_boundary.local_bound_safety_factor * upper;
else
    upper = 0;
end
end

function peak_S = actual_spu_peak_positions(motion)
spu_count = numel(motion.spus);
peak_S = nan(spu_count, 1);
segments = motion.segments;
for spu_id = 1:spu_count
    if ~isfield(motion, 'spu_peak_jerk') || ...
            abs(motion.spu_peak_jerk(spu_id)) <= 0
        continue
    end
    segment_ids = find([segments.spu_id] == spu_id);
    for i = 1:numel(segment_ids)-1
        left = segments(segment_ids(i));
        right = segments(segment_ids(i + 1));
        if left.type == "transition" && right.type == "transition"
            peak_S(spu_id) = 0.5 * (left.s1 + right.s0);
            break
        end
    end
end
end

function [candidate, ok, reason] = evaluate_candidate( ...
        baseline, full_lmsc_jerk, spu_peak_jerk, cfg, keep_spectrum, ...
        maximum_time, require_time_improvement)
candidate = baseline;
ok = false;
try
    motion = velocity_planning_7th( ...
        baseline.arc.S, baseline.profile, [], cfg, baseline.lmsc, ...
        baseline.spus, baseline.geom.kappa, maximum_time, ...
        full_lmsc_jerk, spu_peak_jerk);
    [motion, xyz, scale_diagnostics, dynamically_feasible, ...
        dynamics_reason] = ...
        enforce_candidate_dynamics(motion, baseline, cfg, maximum_time, ...
        require_time_improvement);
    if ~dynamically_feasible
        reason = dynamics_reason;
        return
    end
    candidate = assemble_candidate(baseline, motion, xyz, ...
        scale_diagnostics, cfg, keep_spectrum);
    ok = true;
    reason = "evaluated";
catch exception
    if strlength(string(exception.identifier)) > 0
        reason = string(exception.identifier);
    else
        reason = "evaluation_exception";
    end
end
end

function [motion, xyz, diagnostics, feasible, reason] = ...
        enforce_candidate_dynamics( ...
        motion, baseline, cfg, maximum_time, require_time_improvement)
sample_points = baseline.arc.sample_points;
S = baseline.arc.S;
arc_derivatives = baseline.arc.derivatives;
v_limit = baseline.constraints.v_planning_limit;
unscaled_motion = motion;
xyz_summary = xyz_motion( ...
    sample_points, S, motion, cfg, arc_derivatives, true);
scale = 1;
scaling_iterations = 0;
reason = "dynamic_constraint_failure";

while (~motion.validation.feasible || ~xyz_summary.feasible) && ...
        logical(cfg.speed.enable_global_time_scaling) && ...
        scaling_iterations < cfg.speed.max_global_time_scaling_iterations
    required_scale = max([motion.validation.max_velocity_ratio, ...
        sqrt(motion.validation.max_acceleration_ratio), ...
        nthroot(motion.validation.max_jerk_ratio, 3), ...
        xyz_summary.max_velocity_ratio, ...
        sqrt(xyz_summary.max_acceleration_ratio), ...
        nthroot(xyz_summary.max_jerk_ratio, 3), ...
        nthroot(xyz_summary.max_enforced_snap_ratio, 4), 1]);
    increment = max(required_scale / cfg.speed.constraint_safety_factor, ...
        1 + cfg.speed.velocity_tolerance);
    motion = rescale_septic_motion( ...
        motion, increment, S, v_limit, cfg, maximum_time);
    xyz_summary = xyz_motion( ...
        sample_points, S, motion, cfg, arc_derivatives, true);
    scale = scale * increment;
    scaling_iterations = scaling_iterations + 1;
end

lower_scale = 1;
refinement_iterations = 0;
minimized = logical(cfg.speed.enable_global_time_scaling) && ...
    logical(cfg.speed.minimize_global_time_scale);
if minimized && scale > 1 && ...
        motion.validation.feasible && xyz_summary.feasible
    upper_scale = scale;
    best_motion = motion;
    tolerance = cfg.speed.global_time_scale_refinement_tolerance;
    configured_iterations = round( ...
        cfg.speed.global_time_scale_refinement_iterations);
    required_iterations = max(0, ceil(log2(max( ...
        (upper_scale - lower_scale) / tolerance, 1))));
    maximum_iterations = max(configured_iterations, required_iterations);
    for iteration = 1:maximum_iterations
        if upper_scale - lower_scale <= tolerance * max(upper_scale, 1)
            break
        end
        trial_scale = 0.5 * (lower_scale + upper_scale);
        trial_motion = rescale_septic_motion( ...
            unscaled_motion, trial_scale, S, v_limit, cfg, maximum_time);
        trial_xyz = xyz_motion( ...
            sample_points, S, trial_motion, cfg, arc_derivatives, true);
        refinement_iterations = iteration;
        if trial_motion.validation.feasible && trial_xyz.feasible
            upper_scale = trial_scale;
            best_motion = trial_motion;
        else
            lower_scale = trial_scale;
        end
    end
    motion = best_motion;
    scale = upper_scale;
    xyz_summary = xyz_motion( ...
        sample_points, S, motion, cfg, arc_derivatives, true);
end

feasible = motion.validation.feasible && xyz_summary.feasible;
if feasible
    baseline_time = baseline.raw_objectives.time;
    time_tolerance = 1e-12 * max(baseline_time, 1);
    if require_time_improvement && ...
            motion.time(end) >= baseline_time - time_tolerance
        xyz = struct();
        feasible = false;
        reason = "not_faster_than_zero_jerk_baseline";
    else
        xyz = xyz_motion(sample_points, S, motion, cfg, arc_derivatives);
        feasible = xyz.feasible;
        if feasible
            reason = "feasible";
        end
    end
else
    xyz = struct();
end
diagnostics.global_time_scale = scale;
diagnostics.global_time_scaling_iterations = scaling_iterations;
diagnostics.global_time_scale_refinement_iterations = ...
    refinement_iterations;
diagnostics.global_time_scale_lower_bound = lower_scale;
diagnostics.global_time_scale_minimized = minimized;
end

function candidate = assemble_candidate( ...
        baseline, motion, xyz, scale_diagnostics, cfg, keep_spectrum)
candidate = baseline;
candidate.motion = motion;
candidate.xyz = xyz;
candidate.lmsc = motion.lmsc;
candidate.spus = motion.spus;

[unique_S, unique_index] = unique(motion.S, 'stable');
if numel(unique_S) >= 2
    candidate.constraints.v_planned = interp1( ...
        unique_S, motion.path_velocity(unique_index), ...
        candidate.constraints.S, 'linear', 'extrap');
else
    candidate.constraints.v_planned = zeros(size(candidate.constraints.S));
end

vibration = xyz_vibration_metric( ...
    xyz.jerk, cfg, motion.time, xyz.acceleration, keep_spectrum);
candidate.vibration_metrics = vibration;
raw_time = motion.time(end);
raw_vibration = vibration.objective_value;
dynamic_ratio = max([motion.validation.max_velocity_ratio, ...
    motion.validation.max_acceleration_ratio, ...
    motion.validation.max_jerk_ratio, xyz.max_velocity_ratio, ...
    xyz.max_acceleration_ratio, xyz.max_jerk_ratio, ...
    xyz.max_enforced_snap_ratio]);
dynamic_violation = max(dynamic_ratio - 1, 0);
dynamic_penalty = cfg.constraints.contour_error_penalty * dynamic_violation;
resonance = vibration.resonance;
resonance_violation = max(resonance.max_response_ratio - 1, 0);
resonance_penalty = cfg.resonance.penalty_weight * resonance_violation;
contour_penalty = baseline.constraint.penalty;
total_penalty = contour_penalty + dynamic_penalty + resonance_penalty;

candidate.raw_objectives.time = raw_time;
candidate.raw_objectives.vibration = raw_vibration;
candidate.raw_objectives.vibration_xyz = vibration.axis_rms;
candidate.objectives.time = raw_time + total_penalty;
candidate.objectives.error = baseline.raw_objectives.error + total_penalty;
candidate.objectives.vibration = raw_vibration + total_penalty;

constraint = baseline.constraint;
constraint.dynamic_ratio = dynamic_ratio;
constraint.dynamic_violation = dynamic_violation;
constraint.dynamic_penalty = dynamic_penalty;
constraint.path_feasible = motion.validation.feasible;
constraint.xyz_feasible = xyz.feasible;
constraint.xyz_snap_limit_enabled = xyz.snap_limit_enabled;
constraint.xyz_snap_peak = vibration.axis_snap_peak;
constraint.xyz_snap_axis_ratio = vibration.axis_snap_ratio;
constraint.xyz_snap_ratio = xyz.max_snap_ratio;
constraint.xyz_snap_violation = max(xyz.max_enforced_snap_ratio - 1, 0);
constraint.xyz_snap_feasible = xyz.max_enforced_snap_ratio <= ...
    1 + cfg.speed.velocity_tolerance;
constraint.dynamic_feasible = motion.validation.feasible && xyz.feasible && ...
    dynamic_violation <= cfg.speed.velocity_tolerance;
constraint.global_time_scaling_enabled = ...
    logical(cfg.speed.enable_global_time_scaling);
constraint.global_time_scaling_iterations = ...
    scale_diagnostics.global_time_scaling_iterations;
constraint.global_time_scale = scale_diagnostics.global_time_scale;
constraint.global_time_scale_minimized = ...
    scale_diagnostics.global_time_scale_minimized;
constraint.global_time_scale_refinement_iterations = ...
    scale_diagnostics.global_time_scale_refinement_iterations;
constraint.global_time_scale_lower_bound = ...
    scale_diagnostics.global_time_scale_lower_bound;
constraint.resonance_enabled = resonance.enabled;
constraint.resonance_configured = resonance.configured;
constraint.resonance_ratio = resonance.max_response_ratio;
constraint.resonance_violation = resonance_violation;
constraint.resonance_penalty = resonance_penalty;
constraint.resonance_feasible = resonance.feasible;
constraint.feasible = baseline.constraint.violation <= 0 && ...
    constraint.dynamic_feasible && resonance.feasible;
candidate.constraint = constraint;
end

function validate_configuration(cfg)
settings = cfg.lmsc_jerk_boundary;
if ~isfield(settings, 'evaluate_during_objective') || ...
        ~isscalar(settings.evaluate_during_objective)
    error(['cfg.lmsc_jerk_boundary.evaluate_during_objective must be a ' ...
        'logical scalar.']);
end
if ~isfield(settings, 'enforce_vibration_increase_limit') || ...
        ~isscalar(settings.enforce_vibration_increase_limit)
    error(['cfg.lmsc_jerk_boundary.enforce_vibration_increase_limit ' ...
        'must be a logical scalar.']);
end
if ~isfield(cfg, 'bounds') || ...
        ~isfield(cfg.bounds, 'jerk_utilization_lb') || ...
        ~isfield(cfg.bounds, 'jerk_utilization_ub') || ...
        cfg.bounds.jerk_utilization_lb < 0 || ...
        cfg.bounds.jerk_utilization_ub > 1 || ...
        cfg.bounds.jerk_utilization_lb > cfg.bounds.jerk_utilization_ub
    error('Endpoint-jerk utilization bounds must lie within [0, 1].');
end
if ~isscalar(settings.local_bound_safety_factor) || ...
        ~isfinite(settings.local_bound_safety_factor) || ...
        settings.local_bound_safety_factor <= 0 || ...
        settings.local_bound_safety_factor > 1
    error('local_bound_safety_factor must be in (0, 1].');
end
nonnegative = {'maximum_vibration_increase'};
nonnegative{end + 1} = 'numerical_peak_plateau_time_tolerance';
for i = 1:numel(nonnegative)
    value = settings.(nonnegative{i});
    if ~isscalar(value) || ~isfinite(value) || value < 0
        error('cfg.lmsc_jerk_boundary.%s must be nonnegative.', ...
            nonnegative{i});
    end
end
positive_integers = {'time_search_samples', ...
    'time_refinement_iterations','peak_search_samples'};
for i = 1:numel(positive_integers)
    value = settings.(positive_integers{i});
    if ~isscalar(value) || ~isfinite(value) || value < 1 || ...
            value ~= floor(value)
        error('cfg.lmsc_jerk_boundary.%s must be a positive integer.', ...
            positive_integers{i});
    end
end
if ~isscalar(settings.minimum_duration_fraction) || ...
        ~isfinite(settings.minimum_duration_fraction) || ...
        settings.minimum_duration_fraction <= 0 || ...
        ~isscalar(settings.maximum_duration_factor) || ...
        ~isfinite(settings.maximum_duration_factor) || ...
        settings.maximum_duration_factor < settings.minimum_duration_fraction
    error('Invalid LMSC generalized-transition duration search range.');
end
end

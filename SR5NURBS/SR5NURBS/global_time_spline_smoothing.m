function [motion, xyz, diagnostics] = global_time_spline_smoothing( ...
        reference_motion, reference_xyz, sample_points, S_grid, ...
        v_limit, cfg, arc_derivatives, maximum_time)
% Fit one globally smooth S(t) B-spline and retain it only when beneficial.

diagnostics = default_diagnostics(reference_motion);
motion = reference_motion;
xyz = reference_xyz;
if ~time_smoothing_enabled(cfg)
    diagnostics.reason = "disabled";
    return
end
diagnostics.enabled = true;
diagnostics.attempted = true;

settings = cfg.time_smoothing;
degree = validate_integer(settings.degree, 7, 'degree');
if degree < 7
    diagnostics.reason = "degree_below_seven";
    return_or_error(settings, diagnostics.reason);
    return
end

minimum_control_points = max(validate_integer( ...
    settings.min_control_points, degree + 1, 'min_control_points'), 10);
maximum_control_points = max(validate_integer( ...
    settings.max_control_points, minimum_control_points, ...
    'max_control_points'), minimum_control_points);
lmsc_count = numel(reference_motion.lmsc.indices);
nctrl = round(minimum_control_points + ...
    settings.control_points_per_lmsc * lmsc_count);
nctrl = min(max(nctrl, minimum_control_points), maximum_control_points);
nctrl = max(nctrl, degree + 1);
if nctrl < 10
    diagnostics.reason = "insufficient_control_points";
    return_or_error(settings, diagnostics.reason);
    return
end
diagnostics.control_point_count = nctrl;

total_length = S_grid(end) - S_grid(1);
initial_duration = reference_motion.time(end);
if total_length <= cfg.speed.distance_tolerance || initial_duration <= 0
    diagnostics.reason = "degenerate_reference";
    return
end

knots = open_uniform_knots(nctrl, degree);
control = fit_monotone_time_law(reference_motion, total_length, ...
    nctrl, degree, knots, settings);
duration = initial_duration;
candidate_motion = reference_motion;
candidate_summary = [];

for iteration = 0:validate_integer(settings.max_scaling_iterations, ...
        0, 'max_scaling_iterations')
    diagnostics.attempted_duration_scale = duration / initial_duration;
    if exceeds_time_limit(duration, maximum_time) || ...
            exceeds_duration_factor(duration, initial_duration, settings)
        diagnostics.reason = "duration_limit_exceeded";
        return_or_error(settings, diagnostics.reason);
        return
    end
    candidate_motion = build_motion(reference_motion, control, degree, ...
        knots, duration, S_grid, v_limit, cfg);
    candidate_summary = xyz_motion(sample_points, S_grid, ...
        candidate_motion, cfg, arc_derivatives, true);
    diagnostics.scaling_iterations = iteration;
    if candidate_motion.validation.feasible && candidate_summary.feasible
        break
    end
    required_scale = required_time_scale(candidate_motion, ...
        candidate_summary);
    duration = duration * max( ...
        required_scale / cfg.speed.constraint_safety_factor, ...
        1 + cfg.speed.velocity_tolerance);
end

if ~candidate_motion.validation.feasible || ~candidate_summary.feasible
    diagnostics.reason = "dynamic_constraints_not_satisfied";
    return_or_error(settings, diagnostics.reason);
    return
end

candidate_xyz = xyz_motion(sample_points, S_grid, candidate_motion, ...
    cfg, arc_derivatives);
baseline_score = smoothness_score(reference_xyz, reference_motion.time, cfg);
candidate_score = smoothness_score(candidate_xyz, candidate_motion.time, cfg);
diagnostics.baseline_smoothness = baseline_score;
diagnostics.final_smoothness = candidate_score;
diagnostics.smoothness_improvement = ...
    (baseline_score - candidate_score) / max(baseline_score, eps);

minimum_improvement = settings.minimum_smoothness_improvement;
if ~isscalar(minimum_improvement) || ~isfinite(minimum_improvement) || ...
        minimum_improvement < 0
    error(['cfg.time_smoothing.minimum_smoothness_improvement must be a ' ...
        'finite nonnegative scalar.']);
end
if candidate_score > baseline_score * (1 - minimum_improvement) && ...
        baseline_score > eps
    diagnostics.reason = "smoothness_not_improved";
    return_or_error(settings, diagnostics.reason);
    return
end

candidate_motion.time_smoothing.accepted = true;
candidate_motion.time_smoothing.degree = degree;
candidate_motion.time_smoothing.control_points = control;
candidate_motion.time_smoothing.knots = knots;
candidate_motion.time_smoothing.duration_scale = duration / initial_duration;
motion = candidate_motion;
xyz = candidate_xyz;
diagnostics.accepted = true;
diagnostics.reason = "accepted";
diagnostics.final_time = motion.time(end);
diagnostics.duration_scale = duration / initial_duration;
end

function control = fit_monotone_time_law(reference_motion, total_length, ...
        nctrl, degree, knots, settings)
sample_count = min(numel(reference_motion.time), ...
    validate_integer(settings.fit_samples, 2, 'fit_samples'));
sample_id = unique(round(linspace(1, numel(reference_motion.time), ...
    sample_count)));
tau = reference_motion.time(sample_id) / reference_motion.time(end);
target = reference_motion.S(sample_id) - reference_motion.S(1);
basis = bspline_basis_matrix(tau, nctrl, degree, knots);

endpoint_order = 4;
fixed = [1: endpoint_order + 1, nctrl - endpoint_order:nctrl];
free = setdiff(1:nctrl, fixed, 'stable');
control = zeros(nctrl, 1);
control(nctrl - endpoint_order:nctrl) = total_length;

rows = basis(:, free);
right_hand_side = target - basis(:, fixed) * control(fixed);
[rows, right_hand_side] = append_difference_penalty(rows, ...
    right_hand_side, control, free, fixed, 4, ...
    settings.fourth_difference_regularization);
[rows, right_hand_side] = append_difference_penalty(rows, ...
    right_hand_side, control, free, fixed, 5, ...
    settings.fifth_difference_regularization);
if isempty(free)
    error('The global time spline has no free control points.');
end
control(free) = rows \ right_hand_side;

% A nondecreasing B-spline control polygon guarantees dS/dt >= 0.
control = min(max(control, 0), total_length);
control(1:end) = cummax(control);
control(1:endpoint_order + 1) = 0;
control(nctrl - endpoint_order:nctrl) = total_length;
end

function [rows, right_hand_side] = append_difference_penalty( ...
        rows, right_hand_side, control, free, fixed, order, weight)
if ~isscalar(weight) || ~isfinite(weight) || weight < 0
    error('Time-spline difference regularization must be nonnegative.');
end
if weight == 0 || numel(control) <= order
    return
end
difference = diff(eye(numel(control)), order, 1);
scale = sqrt(weight);
rows = [rows; scale * difference(:, free)];
right_hand_side = [right_hand_side; ...
    -scale * difference(:, fixed) * control(fixed)];
end

function motion = build_motion(reference_motion, control, degree, knots, ...
        duration, S_grid, v_limit, cfg)
time = fixed_time_axis(duration, cfg.interpolation.Ts);
tau = time / duration;
derivatives = evaluate_nurbs_derivatives(control, ones(size(control)), ...
    degree, knots, tau, 4);

path_S = S_grid(1) + reshape(derivatives(:, 1, 1), [], 1);
path_velocity = reshape(derivatives(:, 1, 2), [], 1) / duration;
path_acceleration = reshape(derivatives(:, 1, 3), [], 1) / duration^2;
path_jerk = reshape(derivatives(:, 1, 4), [], 1) / duration^3;
path_snap = reshape(derivatives(:, 1, 5), [], 1) / duration^4;
path_S(1) = S_grid(1);
path_S(end) = S_grid(end);
path_velocity([1 end]) = 0;
path_acceleration([1 end]) = 0;
path_jerk([1 end]) = 0;
path_snap([1 end]) = 0;

motion = reference_motion;
motion.time = time;
motion.S = path_S;
motion.path_velocity = path_velocity;
motion.path_acceleration = path_acceleration;
motion.path_jerk = path_jerk;
motion.path_snap = path_snap;
motion.constraint_scale = reference_motion.constraint_scale * ...
    duration / reference_motion.time(end);
motion.validation = validate_global_motion(path_S, path_velocity, ...
    path_acceleration, path_jerk, S_grid, v_limit, cfg);
motion.time_smoothing.accepted = false;
end

function validation = validate_global_motion(path_S, path_velocity, ...
        path_acceleration, path_jerk, S_grid, v_limit, cfg)
sampled_limit = interp1(S_grid, v_limit, path_S, 'linear', 'extrap');
velocity_excess = path_velocity - sampled_limit;
[maximum_violation, violation_id] = max(velocity_excess);
validation.maximum_violation = max(maximum_violation, 0);
validation.violation_S = NaN;
if validation.maximum_violation > 0
    validation.violation_S = path_S(violation_id);
end
validation.violation_spu_id = 0;
validation.max_velocity_ratio = max(path_velocity ./ ...
    max(sampled_limit, cfg.speed.velocity_tolerance));
validation.max_acceleration_ratio = max(abs(path_acceleration)) / ...
    max(cfg.limits.Amax, eps);
validation.max_jerk_ratio = max(abs(path_jerk)) / ...
    max(cfg.limits.Jmax, eps);
validation.minimum_velocity = min(path_velocity);
validation.max_join_jump = zeros(1, 4);
validation.continuous = true;
validation.monotone_S = all(diff(path_S) >= -cfg.speed.distance_tolerance);
validation.feasible = validation.maximum_violation <= ...
    cfg.speed.velocity_tolerance && ...
    validation.max_velocity_ratio <= 1 + cfg.speed.velocity_tolerance && ...
    validation.max_acceleration_ratio <= 1 + cfg.speed.velocity_tolerance && ...
    validation.max_jerk_ratio <= 1 + cfg.speed.velocity_tolerance && ...
    validation.minimum_velocity >= -cfg.speed.velocity_tolerance && ...
    validation.monotone_S;
end

function scale = required_time_scale(motion, xyz)
scale = max([motion.validation.max_velocity_ratio, ...
    sqrt(max(motion.validation.max_acceleration_ratio, 0)), ...
    nthroot(max(motion.validation.max_jerk_ratio, 0), 3), ...
    xyz.max_velocity_ratio, sqrt(max(xyz.max_acceleration_ratio, 0)), ...
    nthroot(max(xyz.max_jerk_ratio, 0), 3), ...
    nthroot(max(xyz.max_enforced_snap_ratio, 0), 4), 1]);
end

function score = smoothness_score(xyz, time, cfg)
weights = cfg.objectives.vibration_axis_weights(:)';
axis_jerk_rms = sqrt(mean(xyz.jerk.^2, 1));
score = norm(weights .* axis_jerk_rms, 2);
if numel(time) < 2
    return
end
sample_time = median(diff(time));
axis_snap_rms = zeros(1, 3);
axis_crackle_rms = zeros(1, 3);
for axis_id = 1:3
    snap_axis = gradient(xyz.jerk(:, axis_id), sample_time);
    crackle_axis = gradient(snap_axis, sample_time);
    axis_snap_rms(axis_id) = sqrt(mean(snap_axis.^2));
    axis_crackle_rms(axis_id) = sqrt(mean(crackle_axis.^2));
end
score = score + cfg.objectives.vibration_snap_weight * ...
    norm(weights .* (sample_time * axis_snap_rms), 2) + ...
    cfg.objectives.vibration_crackle_weight * ...
    norm(weights .* (sample_time^2 * axis_crackle_rms), 2);
end

function knots = open_uniform_knots(nctrl, degree)
internal_count = nctrl - degree - 1;
if internal_count > 0
    internal = (1:internal_count) / (internal_count + 1);
else
    internal = zeros(1, 0);
end
knots = [zeros(1, degree + 1), internal, ones(1, degree + 1)];
end

function time = fixed_time_axis(total_time, Ts)
tolerance = max(10 * eps(max(total_time, 1)), Ts * 1e-12);
time = (0:Ts:total_time)';
if total_time - time(end) > tolerance
    time(end + 1, 1) = total_time;
else
    time(end) = total_time;
end
if isscalar(time)
    time(end + 1, 1) = total_time;
end
end

function exceeds = exceeds_time_limit(duration, maximum_time)
exceeds = ~isinf(maximum_time) && duration > maximum_time + ...
    1e-12 * max([1, duration, maximum_time]);
end

function exceeds = exceeds_duration_factor(duration, initial_duration, settings)
maximum_factor = settings.maximum_duration_factor;
if ~isscalar(maximum_factor) || isnan(maximum_factor) || maximum_factor < 1
    error(['cfg.time_smoothing.maximum_duration_factor must be a scalar ' ...
        'greater than or equal to one.']);
end
exceeds = duration > initial_duration * maximum_factor * (1 + 1e-12);
end

function value = validate_integer(value, minimum, name)
if ~isscalar(value) || ~isfinite(value) || value < minimum || ...
        value ~= round(value)
    error('cfg.time_smoothing.%s must be an integer >= %d.', name, minimum);
end
end

function enabled = time_smoothing_enabled(cfg)
enabled = isfield(cfg, 'time_smoothing') && ...
    isfield(cfg.time_smoothing, 'enabled') && ...
    logical(cfg.time_smoothing.enabled);
end

function return_or_error(settings, reason)
if ~isfield(settings, 'fallback_to_septic') || ...
        ~logical(settings.fallback_to_septic)
    error('global_time_spline_smoothing:Rejected', ...
        'Global time spline rejected: %s.', char(reason));
end
end

function diagnostics = default_diagnostics(reference_motion)
diagnostics.enabled = false;
diagnostics.attempted = false;
diagnostics.accepted = false;
diagnostics.reason = "not_attempted";
diagnostics.control_point_count = 0;
diagnostics.scaling_iterations = 0;
diagnostics.initial_time = reference_motion.time(end);
diagnostics.final_time = reference_motion.time(end);
diagnostics.duration_scale = 1;
diagnostics.attempted_duration_scale = 1;
diagnostics.baseline_smoothness = NaN;
diagnostics.final_smoothness = NaN;
diagnostics.smoothness_improvement = NaN;
end

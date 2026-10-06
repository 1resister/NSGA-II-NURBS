function [constraints, profile, motion, xyz, diagnostics] = plan_speed_profile( ...
        sample_points, S, geom, constraints, speed_scale_nodes, cfg, ...
        arc_derivatives, maximum_time)
% Plan a dynamically feasible LMSC/SPU speed profile in the arc-length domain.

if nargin < 7
    arc_derivatives = [];
end
if nargin < 8 || isempty(maximum_time)
    maximum_time = Inf;
end
if ~isscalar(maximum_time) || isnan(maximum_time) || maximum_time <= 0
    error('maximum_time must be a positive scalar or Inf.');
end

S = S(:);
working_limit = constraints.v_candidate;
forced_indices = [];

for constraint_iteration = 1:cfg.speed.max_constraint_iterations
    lmsc = detect_lmsc_points(S, working_limit, geom.kappa, cfg, forced_indices);
    spus = build_speed_planning_units(S, working_limit, lmsc);
    profile = forward_backward_scan(S, working_limit, cfg);
    profile = enforce_spu_boundary_reachability(profile, lmsc, cfg);
    motion = velocity_planning_7th(S, profile, speed_scale_nodes, cfg, ...
        lmsc, spus, geom.kappa, maximum_time);
    enforce_candidate_time_limit(motion, maximum_time);
    xyz = xyz_motion(sample_points, S, motion, cfg, arc_derivatives);

    if motion.validation.feasible && xyz.feasible
        break
    end
    if constraint_iteration == cfg.speed.max_constraint_iterations
        break
    end

    [violation_S, required_scale] = locate_constraint_violation(motion, xyz, cfg);
    [~, path_index] = min(abs(S - violation_S));
    forced_indices = unique([forced_indices; path_index], 'stable');

    center_factor = min(max(cfg.speed.constraint_safety_factor / required_scale, 0), 1);
    sample_step = median(diff(S));
    radius = max(cfg.speed.constraint_correction_radius_samples * sample_step, ...
        cfg.speed.distance_tolerance);
    taper = max(1 - abs(S - S(path_index)) / radius, 0);
    correction = 1 - (1 - center_factor) .* taper;
    revised_limit = min(working_limit, max(working_limit .* correction, ...
        cfg.speed.velocity_tolerance));
    if max(abs(revised_limit - working_limit)) <= cfg.speed.velocity_tolerance
        break
    end
    working_limit = revised_limit;
end

global_time_scale = 1;
global_scaling_iterations = 0;
global_refinement_iterations = 0;
global_scale_lower_bound = 1;
enable_global_time_scaling = logical(cfg.speed.enable_global_time_scaling);
unscaled_motion = motion;
maximum_global_scale = maximum_time / max(unscaled_motion.time(end), eps);
xyz_is_compact = false;
while enable_global_time_scaling && ...
        (~motion.validation.feasible || ~xyz.feasible) && ...
        global_scaling_iterations < cfg.speed.max_global_time_scaling_iterations
    required_scale = max([motion.validation.max_velocity_ratio, ...
        sqrt(motion.validation.max_acceleration_ratio), ...
        nthroot(motion.validation.max_jerk_ratio, 3), ...
        xyz.max_velocity_ratio, sqrt(xyz.max_acceleration_ratio), ...
        nthroot(xyz.max_jerk_ratio, 3), ...
        nthroot(xyz.max_enforced_snap_ratio, 4), 1]);
    scale = max(required_scale / cfg.speed.constraint_safety_factor, ...
        1 + cfg.speed.velocity_tolerance);
    reached_time_cap = global_time_scale * scale > maximum_global_scale;
    if reached_time_cap
        scale = maximum_global_scale / global_time_scale;
        if scale <= 1 + cfg.speed.velocity_tolerance
            error('plan_speed_profile:MaximumTimeExceeded', ...
                ['The candidate cannot satisfy all dynamic constraints ' ...
                'within the configured evaluation-time limit of %.9g s.'], ...
                maximum_time);
        end
    end
    % Once global scaling starts, only scalar constraint summaries are
    % needed until the smallest feasible scale has been selected.
    clear xyz
    motion = rescale_septic_motion( ...
        motion, scale, S, working_limit, cfg, maximum_time);
    xyz = xyz_motion( ...
        sample_points, S, motion, cfg, arc_derivatives, true);
    xyz_is_compact = true;
    global_time_scale = global_time_scale * scale;
    global_scaling_iterations = global_scaling_iterations + 1;
    if reached_time_cap && ...
            (~motion.validation.feasible || ~xyz.feasible)
        error('plan_speed_profile:MaximumTimeExceeded', ...
            ['The candidate cannot satisfy all dynamic constraints ' ...
            'within the configured evaluation-time limit of %.9g s.'], ...
            maximum_time);
    end
end

minimize_global_time_scale = enable_global_time_scaling && ...
    logical(cfg.speed.minimize_global_time_scale);
if minimize_global_time_scale && global_time_scale > 1 && ...
        motion.validation.feasible && xyz.feasible
    clear xyz
    [motion, xyz, global_time_scale, global_scale_lower_bound, ...
        global_refinement_iterations] = refine_global_time_scale( ...
        unscaled_motion, motion, global_time_scale, sample_points, ...
        S, working_limit, cfg, arc_derivatives, maximum_time);
elseif xyz_is_compact
    xyz = xyz_motion(sample_points, S, motion, cfg, arc_derivatives);
end

[motion, xyz, time_smoothing_diagnostics] = ...
    global_time_spline_smoothing(motion, xyz, sample_points, S, ...
    working_limit, cfg, arc_derivatives, maximum_time);

constraints.v_dynamic = working_limit;
constraints.v_planning_limit = working_limit;
dynamic_mask = working_limit < constraints.v_candidate - cfg.speed.velocity_tolerance;
constraints.active_constraint(dynamic_mask) = "dynamic";
constraints.Vlimit = working_limit;

[unique_S, unique_id] = unique(motion.S, 'stable');
if numel(unique_S) < 2
    constraints.v_planned = zeros(size(S));
else
    constraints.v_planned = interp1(unique_S, motion.path_velocity(unique_id), ...
        S, 'linear', 'extrap');
end

diagnostics.iterations = constraint_iteration;
diagnostics.global_time_scaling_enabled = enable_global_time_scaling;
diagnostics.global_time_scaling_iterations = global_scaling_iterations;
diagnostics.global_time_scale = global_time_scale;
diagnostics.global_time_scale_minimized = minimize_global_time_scale;
diagnostics.global_time_scale_refinement_iterations = ...
    global_refinement_iterations;
diagnostics.global_time_scale_lower_bound = global_scale_lower_bound;
diagnostics.time_smoothing = time_smoothing_diagnostics;
diagnostics.time_smoothing_enabled = time_smoothing_diagnostics.enabled;
diagnostics.time_smoothing_accepted = time_smoothing_diagnostics.accepted;
diagnostics.time_smoothing_duration_scale = ...
    time_smoothing_diagnostics.duration_scale;
diagnostics.corrected = any(dynamic_mask);
diagnostics.path_feasible = motion.validation.feasible;
diagnostics.xyz_feasible = xyz.feasible;
diagnostics.feasible = diagnostics.path_feasible && diagnostics.xyz_feasible;
end

function [best_motion, best_xyz, upper_scale, lower_scale, iterations] = ...
        refine_global_time_scale(unscaled_motion, best_motion, ...
        upper_scale, sample_points, S, working_limit, cfg, arc_derivatives, ...
        maximum_time)
% Find the smallest feasible uniform time scale inside [1, feasible upper].
lower_scale = 1;
iterations = 0;
tolerance = cfg.speed.global_time_scale_refinement_tolerance;
if ~isscalar(tolerance) || ~isfinite(tolerance) || tolerance <= 0
    error(['cfg.speed.global_time_scale_refinement_tolerance must be a ' ...
        'positive finite scalar.']);
end
configured_iterations = cfg.speed.global_time_scale_refinement_iterations;
if ~isscalar(configured_iterations) || ~isfinite(configured_iterations) || ...
        configured_iterations < 0
    error(['cfg.speed.global_time_scale_refinement_iterations must be a ' ...
        'nonnegative finite scalar.']);
end
configured_iterations = round(configured_iterations);
max_iterations = configured_iterations;
if configured_iterations > 0
    % A fixed iteration count cannot guarantee the requested tolerance when
    % the first feasible scale is far from one. Extend only as much as the
    % initial bracket mathematically requires; this avoids both a false test
    % failure and an unnecessarily large fixed cost for ordinary cases.
    initial_width = max(upper_scale - lower_scale, 0);
    required_iterations = max(0, ...
        ceil(log2(max(initial_width / tolerance, 1))));
    max_iterations = max(configured_iterations, required_iterations);
end

for iteration = 1:max_iterations
    if upper_scale - lower_scale <= tolerance * max(upper_scale, 1)
        break
    end
    trial_scale = 0.5 * (lower_scale + upper_scale);
    trial_motion = rescale_septic_motion( ...
        unscaled_motion, trial_scale, S, working_limit, cfg, maximum_time);
    trial_xyz = xyz_motion( ...
        sample_points, S, trial_motion, cfg, arc_derivatives, true);
    iterations = iteration;
    if trial_motion.validation.feasible && trial_xyz.feasible
        upper_scale = trial_scale;
        best_motion = trial_motion;
    else
        lower_scale = trial_scale;
    end
end
best_xyz = xyz_motion( ...
    sample_points, S, best_motion, cfg, arc_derivatives);
end

function enforce_candidate_time_limit(motion, maximum_time)
if isinf(maximum_time)
    return
end
total_time = motion.time(end);
tolerance = 1e-12 * max([1, total_time, maximum_time]);
if total_time > maximum_time + tolerance
    error('plan_speed_profile:MaximumTimeExceeded', ...
        ['The candidate minimum unscaled traversal time %.9g s exceeds ' ...
        'the configured evaluation-time limit %.9g s.'], ...
        total_time, maximum_time);
end
end

function [violation_S, required_scale] = locate_constraint_violation(motion, xyz, cfg)
path_required_scale = max([motion.validation.max_velocity_ratio, ...
    sqrt(motion.validation.max_acceleration_ratio), ...
    nthroot(motion.validation.max_jerk_ratio, 3), 1]);
if motion.validation.maximum_violation > cfg.speed.velocity_tolerance && ...
        isfinite(motion.validation.violation_S)
    violation_S = motion.validation.violation_S;
    required_scale = path_required_scale;
    return
end

if ~xyz.feasible
    [~, time_index] = max(xyz.point_constraint_ratio);
    violation_S = motion.S(time_index);
    local_v = max(xyz.velocity_ratio(time_index, :));
    local_a = max(xyz.acceleration_ratio(time_index, :));
    local_j = max(xyz.jerk_ratio(time_index, :));
    local_s = 0;
    if xyz.snap_limit_enabled
        local_s = max(xyz.snap_ratio(time_index, :));
    end
    required_scale = max([local_v, sqrt(local_a), nthroot(local_j, 3), ...
        nthroot(local_s, 4), 1]);
    return
end

violation_S = motion.validation.violation_S;
required_scale = path_required_scale;
if ~isfinite(violation_S)
    path_ratio = max([motion.path_velocity / max(cfg.limits.Vmax, eps), ...
        abs(motion.path_acceleration) / max(cfg.limits.Amax, eps), ...
        abs(motion.path_jerk) / max(cfg.limits.Jmax, eps)], [], 2);
    [~, time_index] = max(path_ratio);
    violation_S = motion.S(time_index);
end
end

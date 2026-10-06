function motion = rescale_septic_motion( ...
        motion, scale, S, v_limit, cfg, maximum_time)
% Uniformly stretch an existing seventh-order motion without changing q(s).

scale = max(real(scale), 1);
if nargin < 6 || isempty(maximum_time)
    maximum_time = Inf;
end
segments = motion.segments;
for i = 1:numel(segments)
    segments(i).t0 = segments(i).t0 * scale;
    segments(i).t1 = segments(i).t1 * scale;
    segments(i).v0 = segments(i).v0 / scale;
    segments(i).v1 = segments(i).v1 / scale;
    if isfield(segments, 'j0')
        segments(i).j0 = segments(i).j0 / scale^3;
        segments(i).j1 = segments(i).j1 / scale^3;
    end
end

spus = motion.spus;
for i = 1:numel(spus)
    spus(i).start_speed_limit = spus(i).start_speed_limit / scale;
    spus(i).end_speed_limit = spus(i).end_speed_limit / scale;
    spus(i).planned_peak_speed = spus(i).planned_peak_speed / scale;
    spus(i).acceleration_time = spus(i).acceleration_time * scale;
    spus(i).constant_time = spus(i).constant_time * scale;
    spus(i).deceleration_time = spus(i).deceleration_time * scale;
    if isfield(spus, 'start_jerk')
        spus(i).start_jerk = spus(i).start_jerk / scale^3;
        if isfield(spus, 'peak_jerk')
            spus(i).peak_jerk = spus(i).peak_jerk / scale^3;
        end
        spus(i).end_jerk = spus(i).end_jerk / scale^3;
    end
end

total_time = segments(end).t1;
enforce_candidate_time_limit(total_time, maximum_time);
time = fixed_time_axis(total_time, cfg.interpolation.Ts);
[path_S, path_V, path_A, path_J] = evaluate_segments(segments, time);
path_S(1) = S(1);
path_S(end) = S(end);
path_V([1 end]) = 0;
path_A([1 end]) = 0;
path_J([1 end]) = 0;

motion.time = time;
motion.S = path_S;
motion.path_velocity = path_V;
motion.path_acceleration = path_A;
motion.path_jerk = path_J;
motion.segments = segments;
motion.spus = spus;
if isfield(motion, 'lmsc_jerk')
    motion.lmsc_jerk = motion.lmsc_jerk / scale^3;
end
if isfield(motion, 'spu_peak_jerk')
    motion.spu_peak_jerk = motion.spu_peak_jerk / scale^3;
end
motion.constraint_scale = motion.constraint_scale * scale;
motion.validation = validate_motion(segments, time, path_S, path_V, ...
    path_A, path_J, S, v_limit, cfg);
end

function enforce_candidate_time_limit(total_time, maximum_time)
if isinf(maximum_time)
    return
end
tolerance = 1e-12 * max([1, total_time, maximum_time]);
if total_time > maximum_time + tolerance
    error('plan_speed_profile:MaximumTimeExceeded', ...
        ['The globally scaled candidate time %.9g s would require more ' ...
        'than the configured fixed-Ts evaluation budget (limit %.9g s).'], ...
        total_time, maximum_time);
end
end

function validation = validate_motion(segments, time, path_S, path_V, ...
        path_A, path_J, S, v_limit, cfg)
validation.maximum_violation = 0;
validation.violation_S = NaN;
validation.violation_spu_id = 0;
validation.max_velocity_ratio = 0;
validation.max_acceleration_ratio = 0;
validation.max_jerk_ratio = 0;

for i = 1:numel(segments)
    segment = segments(i);
    duration = segment.t1 - segment.t0;
    check_samples = max(cfg.speed.internal_constraint_check_samples, ...
        ceil(duration / cfg.interpolation.Ts) + 1);
    local_time = linspace(0, duration, check_samples)';
    [sq, vq, aq, jq] = evaluate_one_segment(segment, local_time);
    local_limit = interp1(S, v_limit, sq, 'linear', 'extrap');
    [maximum, index] = max(vq - local_limit);
    if maximum > validation.maximum_violation
        validation.maximum_violation = maximum;
        validation.violation_S = sq(index);
        validation.violation_spu_id = segment.spu_id;
    end
    validation.max_velocity_ratio = max(validation.max_velocity_ratio, ...
        max(vq ./ max(local_limit, cfg.speed.velocity_tolerance)));
    validation.max_acceleration_ratio = max(validation.max_acceleration_ratio, ...
        max(abs(aq)) / max(cfg.limits.Amax, eps));
    validation.max_jerk_ratio = max(validation.max_jerk_ratio, ...
        max(abs(jq)) / max(cfg.limits.Jmax, eps));
end

sampled_limit = interp1(S, v_limit, path_S, 'linear', 'extrap');
[sampled_violation, sampled_index] = max(path_V - sampled_limit);
if sampled_violation > validation.maximum_violation
    validation.maximum_violation = sampled_violation;
    validation.violation_S = path_S(sampled_index);
    segment_index = find(time(sampled_index) >= [segments.t0] & ...
        time(sampled_index) <= [segments.t1], 1, 'last');
    if ~isempty(segment_index)
        validation.violation_spu_id = segments(segment_index).spu_id;
    end
end
validation.max_velocity_ratio = max(validation.max_velocity_ratio, ...
    max(path_V ./ max(sampled_limit, cfg.speed.velocity_tolerance)));
validation.max_acceleration_ratio = max(validation.max_acceleration_ratio, ...
    max(abs(path_A)) / max(cfg.limits.Amax, eps));
validation.max_jerk_ratio = max(validation.max_jerk_ratio, ...
    max(abs(path_J)) / max(cfg.limits.Jmax, eps));
validation.minimum_velocity = min(path_V);
validation.max_join_jump = segment_join_jumps(segments);
join_tolerance = [cfg.speed.distance_tolerance, cfg.speed.velocity_tolerance, ...
    max(cfg.limits.Amax * 1e-9, cfg.speed.velocity_tolerance), ...
    max(cfg.limits.Jmax * 1e-9, cfg.speed.velocity_tolerance)];
validation.continuous = all(validation.max_join_jump <= join_tolerance);
validation.monotone_S = all(diff(path_S) >= -cfg.speed.distance_tolerance);
validation.feasible = ...
    validation.maximum_violation <= cfg.speed.velocity_tolerance && ...
    validation.max_velocity_ratio <= 1 + cfg.speed.velocity_tolerance && ...
    validation.max_acceleration_ratio <= 1 + cfg.speed.velocity_tolerance && ...
    validation.max_jerk_ratio <= 1 + cfg.speed.velocity_tolerance && ...
    validation.continuous && validation.monotone_S;
end

function time = fixed_time_axis(total_time, Ts)
time_tolerance = max(10 * eps(max(total_time, 1)), Ts * 1e-12);
time = (0:Ts:total_time)';
if isempty(time)
    time = 0;
end
if total_time - time(end) > time_tolerance
    time(end + 1, 1) = total_time;
else
    time(end) = total_time;
end
if isscalar(time)
    time(end + 1, 1) = total_time;
end
end

function [S, V, A, J] = evaluate_segments(segments, time)
S = zeros(size(time));
V = zeros(size(time));
A = zeros(size(time));
J = zeros(size(time));
for i = 1:numel(segments)
    if i < numel(segments)
        selected = time >= segments(i).t0 & time < segments(i).t1;
    else
        selected = time >= segments(i).t0 & time <= segments(i).t1;
    end
    local_time = time(selected) - segments(i).t0;
    [S(selected), V(selected), A(selected), J(selected)] = ...
        evaluate_one_segment(segments(i), local_time);
end
end

function [S, V, A, J] = evaluate_one_segment(segment, local_time)
duration = segment.t1 - segment.t0;
if segment.type == "constant"
    S = segment.s0 + segment.v0 .* local_time;
    V = segment.v0 * ones(size(local_time));
    A = zeros(size(local_time));
    J = zeros(size(local_time));
else
    [S, V, A, J] = evaluate_septic_transition(local_time, ...
        segment.s0, segment.v0, segment.v1, duration, ...
        segment.j0, segment.j1);
end
end

function jumps = segment_join_jumps(segments)
jumps = zeros(1, 4);
for i = 1:numel(segments)-1
    left_time = segments(i).t1 - segments(i).t0;
    [Sl, Vl, Al, Jl] = evaluate_one_segment(segments(i), left_time);
    [Sr, Vr, Ar, Jr] = evaluate_one_segment(segments(i + 1), 0);
    jumps = max(jumps, abs([Sl-Sr, Vl-Vr, Al-Ar, Jl-Jr]));
end
end

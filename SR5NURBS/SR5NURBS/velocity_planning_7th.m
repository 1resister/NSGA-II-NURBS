function motion = velocity_planning_7th( ...
        S, profile, ~, cfg, lmsc, spus, kappa, maximum_time, ...
        lmsc_jerk, spu_peak_jerk)
% Generate an LMSC/SPU speed plan with seventh-order flexible transitions.

S = S(:);
if nargin < 8 || isempty(maximum_time)
    maximum_time = Inf;
end
if nargin < 5 || isempty(lmsc)
    if nargin < 7 || isempty(kappa)
        kappa = zeros(size(S));
    end
    lmsc = detect_lmsc_points(S, profile.Vlimit, kappa, cfg);
end
if nargin < 6 || isempty(spus)
    spus = build_speed_planning_units(S, profile.Vlimit, lmsc);
end
if nargin < 7 || isempty(kappa)
    kappa = zeros(size(S));
end
fixed_lmsc_jerk = nargin >= 9 && ~isempty(lmsc_jerk);
if fixed_lmsc_jerk
    lmsc_jerk = lmsc_jerk(:);
    if numel(lmsc_jerk) ~= numel(lmsc.indices) || ...
            any(~isfinite(lmsc_jerk))
        error('velocity_planning_7th:InvalidLMSCJerk', ...
            'lmsc_jerk must contain one finite value per LMSC point.');
    end
    jerk_tolerance = max(cfg.limits.Jmax * 1e-12, ...
        cfg.speed.velocity_tolerance);
    if abs(lmsc_jerk(1)) > jerk_tolerance || ...
            abs(lmsc_jerk(end)) > jerk_tolerance
        error('velocity_planning_7th:EndpointLMSCJerk', ...
            'The first and last LMSC jerk values must remain zero.');
    end
    lmsc_jerk([1 end]) = 0;
else
    lmsc_jerk = zeros(numel(lmsc.indices), 1);
end
fixed_spu_peak_jerk = nargin >= 10 && ~isempty(spu_peak_jerk);
if fixed_spu_peak_jerk
    spu_peak_jerk = spu_peak_jerk(:);
    if numel(spu_peak_jerk) ~= numel(spus) || ...
            any(~isfinite(spu_peak_jerk))
        error('velocity_planning_7th:InvalidSPUPeakJerk', ...
            'spu_peak_jerk must contain one finite value per SPU.');
    end
else
    spu_peak_jerk = zeros(numel(spus), 1);
end
fixed_boundary_jerk = fixed_lmsc_jerk || fixed_spu_peak_jerk;
profile = enforce_spu_boundary_reachability(profile, lmsc, cfg);
% speed_scale_nodes remains in the public interface for compatibility; its
% effect has already been included in profile.Vlimit by the caller.

extra_forced = [];
peak_caps = [spus.candidate_peak_speed]';
internal = struct('feasible', false, 'max_velocity_ratio', Inf, ...
    'max_acceleration_ratio', Inf, 'max_jerk_ratio', Inf, ...
    'maximum_violation', Inf, 'violation_S', NaN, 'violation_spu_id', 0);

for iteration = 1:cfg.speed.max_constraint_iterations
    [segments, spus, build_ok] = plan_spus( ...
        profile, spus, peak_caps, lmsc_jerk, spu_peak_jerk, cfg);
    if ~build_ok
        error('velocity_planning_7th:NoFeasibleSPU', ...
            'No feasible SPU trajectory could be constructed.');
    end
    enforce_candidate_time_limit(segments(end).t1, maximum_time);
    internal = validate_segments(segments, S, profile.Vlimit, cfg);
    if internal.feasible
        break
    end

    % A supplied LMSC-jerk vector fixes the topology so that each shared
    % boundary value keeps the same physical meaning. A violation is
    % handled by rejecting or uniformly scaling that candidate rather than
    % inserting another LMSC and changing the vector length.
    if fixed_boundary_jerk
        break
    end

    [~, new_index] = min(abs(S - internal.violation_S));
    if ~ismember(new_index, lmsc.indices) && new_index > 1 && new_index < numel(S)
        extra_forced(end + 1, 1) = new_index; %#ok<AGROW>
        lmsc = detect_lmsc_points(S, profile.Vlimit, kappa, cfg, extra_forced);
        profile = enforce_spu_boundary_reachability(profile, lmsc, cfg);
        spus = build_speed_planning_units(S, profile.Vlimit, lmsc);
        lmsc_jerk = zeros(numel(lmsc.indices), 1);
        spu_peak_jerk = zeros(numel(spus), 1);
        peak_caps = [spus.candidate_peak_speed]';
    elseif internal.violation_spu_id >= 1 && internal.violation_spu_id <= numel(spus)
        id = internal.violation_spu_id;
        local_limit = interp1(S, profile.Vlimit, internal.violation_S, 'linear');
        peak_caps(id) = min(peak_caps(id), ...
            cfg.speed.constraint_safety_factor * local_limit);
    else
        break
    end
end

if ~exist('segments', 'var') || isempty(segments)
    error('velocity_planning_7th:NoFeasibleSPU', ...
        'No feasible SPU trajectory could be constructed.');
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
motion.lmsc = lmsc;
motion.lmsc_jerk = lmsc_jerk;
motion.spu_peak_jerk = spu_peak_jerk;
motion.spus = spus;
motion.validation = internal;
motion.validation.minimum_velocity = min(path_V);
motion.validation.max_join_jump = segment_join_jumps(segments);
sampled_limit = interp1(S, profile.Vlimit, path_S, 'linear', 'extrap');
[sampled_violation, sampled_index] = max(path_V - sampled_limit);
if sampled_violation > motion.validation.maximum_violation
    motion.validation.maximum_violation = sampled_violation;
    motion.validation.violation_S = path_S(sampled_index);
    segment_index = find(time(sampled_index) >= [segments.t0] & ...
        time(sampled_index) <= [segments.t1], 1, 'last');
    if ~isempty(segment_index)
        motion.validation.violation_spu_id = segments(segment_index).spu_id;
    end
end
motion.validation.max_velocity_ratio = max(motion.validation.max_velocity_ratio, ...
    max(path_V ./ max(sampled_limit, cfg.speed.velocity_tolerance)));
motion.validation.max_acceleration_ratio = max(motion.validation.max_acceleration_ratio, ...
    max(abs(path_A)) / max(cfg.limits.Amax, eps));
motion.validation.max_jerk_ratio = max(motion.validation.max_jerk_ratio, ...
    max(abs(path_J)) / max(cfg.limits.Jmax, eps));
join_tolerance = [cfg.speed.distance_tolerance, cfg.speed.velocity_tolerance, ...
    max(cfg.limits.Amax * 1e-9, cfg.speed.velocity_tolerance), ...
    max(cfg.limits.Jmax * 1e-9, cfg.speed.velocity_tolerance)];
motion.validation.continuous = all(motion.validation.max_join_jump <= join_tolerance);
motion.validation.monotone_S = all(diff(path_S) >= -cfg.speed.distance_tolerance);
motion.validation.feasible = ...
    motion.validation.maximum_violation <= cfg.speed.velocity_tolerance && ...
    motion.validation.max_velocity_ratio <= 1 + cfg.speed.velocity_tolerance && ...
    motion.validation.max_acceleration_ratio <= 1 + cfg.speed.velocity_tolerance && ...
    motion.validation.max_jerk_ratio <= 1 + cfg.speed.velocity_tolerance && ...
    motion.validation.continuous && motion.validation.monotone_S;
motion.constraint_scale = 1;
end

function enforce_candidate_time_limit(total_time, maximum_time)
if isinf(maximum_time)
    return
end
tolerance = 1e-12 * max([1, total_time, maximum_time]);
if total_time > maximum_time + tolerance
    error('plan_speed_profile:MaximumTimeExceeded', ...
        ['The candidate traversal time %.9g s would require more than the ' ...
        'configured fixed-Ts evaluation budget (limit %.9g s).'], ...
        total_time, maximum_time);
end
end

function [segments, spus, feasible] = plan_spus( ...
        profile, spus, peak_caps, lmsc_jerk, spu_peak_jerk, cfg)
segments = empty_segments();
feasible = true;
t_cursor = 0;

for id = 1:numel(spus)
    i0 = spus(id).start_index;
    i1 = spus(id).end_index;
    L = spus(id).length;
    vs = max(profile.V(i0), 0);
    ve = max(profile.V(i1), 0);
    local_profile = profile.V(i0:i1);
    candidate_peak = min(max(local_profile), peak_caps(id));
    candidate_peak = max(candidate_peak, max(vs, ve));
    start_jerk = lmsc_jerk(id);
    end_jerk = lmsc_jerk(id + 1);
    peak_jerk = spu_peak_jerk(id);
    spus(id).start_jerk = start_jerk;
    spus(id).peak_jerk = peak_jerk;
    spus(id).end_jerk = end_jerk;

    if L <= cfg.speed.distance_tolerance
        spus(id).planned_peak_speed = max(vs, ve);
        continue
    end

    [peak, peak_ok] = solve_spu_peak( ...
        vs, ve, start_jerk, peak_jerk, end_jerk, ...
        candidate_peak, L, cfg);
    [La, Ta, acceleration_ok] = transition_metrics( ...
        vs, peak, start_jerk, peak_jerk, cfg);
    [Ld, Td, deceleration_ok] = transition_metrics( ...
        peak, ve, peak_jerk, end_jerk, cfg);
    if ~peak_ok || ~acceleration_ok || ~deceleration_ok || ...
            La + Ld > L + cfg.speed.distance_tolerance
        feasible = false;
        return
    end

    jerk_tolerance = max(cfg.limits.Jmax * 1e-12, ...
        cfg.speed.velocity_tolerance);
    if abs(peak_jerk) > jerk_tolerance && ...
            L - La - Ld > cfg.speed.distance_tolerance
        [La, Ta, Ld, Td, stretched_ok] = ...
            stretch_transition_pair_to_distance(vs, peak, ve, ...
            start_jerk, peak_jerk, end_jerk, Ta, Td, L, cfg);
        if ~stretched_ok
            feasible = false;
            return
        end
    end

    Lc = max(L - La - Ld, 0);
    if Lc > cfg.speed.distance_tolerance && peak <= cfg.speed.velocity_tolerance
        feasible = false;
        return
    end
    Tc = Lc / max(peak, cfg.speed.velocity_tolerance);

    s_cursor = spus(id).start_s;
    if Ta > 0
        segment = make_segment("transition", id, t_cursor, t_cursor + Ta, ...
            s_cursor, s_cursor + La, vs, peak, start_jerk, peak_jerk);
        segments(end + 1, 1) = segment; %#ok<AGROW>
        t_cursor = t_cursor + Ta;
        s_cursor = s_cursor + La;
    end
    if Lc > cfg.speed.distance_tolerance
        segment = make_segment("constant", id, t_cursor, t_cursor + Tc, ...
            s_cursor, s_cursor + Lc, peak, peak, 0, 0);
        segments(end + 1, 1) = segment; %#ok<AGROW>
        t_cursor = t_cursor + Tc;
        s_cursor = s_cursor + Lc;
    end
    if Td > 0
        segment = make_segment("transition", id, t_cursor, t_cursor + Td, ...
            s_cursor, spus(id).end_s, peak, ve, peak_jerk, end_jerk);
        segments(end + 1, 1) = segment; %#ok<AGROW>
        t_cursor = t_cursor + Td;
    end

    spus(id).candidate_peak_speed = candidate_peak;
    spus(id).planned_peak_speed = peak;
    spus(id).start_speed_limit = vs;
    spus(id).end_speed_limit = ve;
    spus(id).acceleration_time = Ta;
    spus(id).constant_time = Tc;
    spus(id).deceleration_time = Td;
end

if isempty(segments)
    feasible = false;
end
end

function [peak, feasible] = solve_spu_peak( ...
        vs, ve, start_jerk, peak_jerk, end_jerk, upper, L, cfg)
base = max(vs, ve);
upper = max(upper, base);
feasible = false;
[required_upper, upper_ok] = required_spu_distance( ...
    vs, ve, start_jerk, peak_jerk, end_jerk, upper, cfg);
if upper_ok && required_upper <= L + cfg.speed.distance_tolerance
    peak = upper;
    feasible = true;
    return
end

jerk_tolerance = max(cfg.limits.Jmax * 1e-12, ...
    cfg.speed.velocity_tolerance);
if abs(start_jerk) > jerk_tolerance || ...
        abs(peak_jerk) > jerk_tolerance || ...
        abs(end_jerk) > jerk_tolerance
    [peak, feasible] = solve_nonzero_jerk_peak( ...
        vs, ve, start_jerk, peak_jerk, end_jerk, ...
        base, upper, L, cfg);
    return
end

[required_base, base_ok] = required_spu_distance( ...
    vs, ve, start_jerk, peak_jerk, end_jerk, base, cfg);
if ~base_ok || required_base > L + cfg.speed.distance_tolerance
    peak = base;
    return
end

low = base;
high = upper;
for iter = 1:cfg.speed.binary_search_max_iterations
    middle = 0.5 * (low + high);
    [required, middle_ok] = required_spu_distance( ...
        vs, ve, start_jerk, peak_jerk, end_jerk, middle, cfg);
    if middle_ok && required <= L
        low = middle;
    else
        high = middle;
    end
    if high - low <= cfg.speed.velocity_tolerance || ...
            abs(required - L) <= cfg.speed.distance_tolerance
        break
    end
end
peak = low;
feasible = true;
end

function [peak, feasible] = solve_nonzero_jerk_peak( ...
        vs, ve, start_jerk, peak_jerk, end_jerk, ...
        base, upper, L, cfg)
samples = cfg.lmsc_jerk_boundary.peak_search_samples;
peak_grid = linspace(base, upper, samples);
valid = false(size(peak_grid));
distance = inf(size(peak_grid));
for i = 1:numel(peak_grid)
    [distance(i), transition_ok] = required_spu_distance( ...
        vs, ve, start_jerk, peak_jerk, end_jerk, peak_grid(i), cfg);
    valid(i) = transition_ok && ...
        distance(i) <= L + cfg.speed.distance_tolerance;
end
last_valid = find(valid, 1, 'last');
if isempty(last_valid)
    peak = base;
    feasible = false;
    return
end

peak = peak_grid(last_valid);
feasible = true;
if last_valid == numel(peak_grid)
    return
end
low = peak_grid(last_valid);
high = peak_grid(last_valid + 1);
for iteration = 1:cfg.speed.binary_search_max_iterations
    middle = 0.5 * (low + high);
    [required, middle_ok] = required_spu_distance( ...
        vs, ve, start_jerk, peak_jerk, end_jerk, middle, cfg);
    if middle_ok && required <= L
        low = middle;
    else
        high = middle;
    end
    if high - low <= cfg.speed.velocity_tolerance
        break
    end
end
peak = low;
end

function [distance, feasible] = required_spu_distance( ...
        vs, ve, start_jerk, peak_jerk, end_jerk, peak, cfg)
[La, ~, acceleration_ok] = transition_metrics( ...
    vs, peak, start_jerk, peak_jerk, cfg);
[Ld, ~, deceleration_ok] = transition_metrics( ...
    peak, ve, peak_jerk, end_jerk, cfg);
feasible = acceleration_ok && deceleration_ok;
if feasible
    distance = La + Ld;
else
    distance = Inf;
end
end

function [La, Ta, Ld, Td, feasible] = ...
        stretch_transition_pair_to_distance(vs, peak, ve, ...
        start_jerk, peak_jerk, end_jerk, Ta, Td, target_distance, cfg)
% A nonzero peak jerk cannot connect directly to a constant-speed segment.
% Stretch the two transitions by one common factor until they consume the
% complete SPU distance, preserving the shared jerk at their common peak.
La = Inf;
Ld = Inf;
feasible = false;
if Ta <= 0 || Td <= 0
    return
end

maximum_Ta = max(Ta, cfg.lmsc_jerk_boundary.maximum_duration_factor * ...
    septic_transition_time(vs, peak, cfg));
maximum_Td = max(Td, cfg.lmsc_jerk_boundary.maximum_duration_factor * ...
    septic_transition_time(peak, ve, cfg));
maximum_scale = min(maximum_Ta / Ta, maximum_Td / Td);
if maximum_scale <= 1 + cfg.speed.velocity_tolerance
    return
end

scale_grid = linspace(1, maximum_scale, ...
    max(cfg.lmsc_jerk_boundary.time_search_samples, 3));
lower_scale = 1;
[lower_distance, lower_ok] = transition_pair_distance( ...
    vs, peak, ve, start_jerk, peak_jerk, end_jerk, ...
    Ta, Td, lower_scale, cfg);
if ~lower_ok || lower_distance > ...
        target_distance + cfg.speed.distance_tolerance
    return
end

upper_scale = NaN;
for i = 2:numel(scale_grid)
    trial_scale = scale_grid(i);
    [trial_distance, trial_ok] = transition_pair_distance( ...
        vs, peak, ve, start_jerk, peak_jerk, end_jerk, ...
        Ta, Td, trial_scale, cfg);
    if ~trial_ok
        continue
    end
    if trial_distance >= target_distance
        upper_scale = trial_scale;
        break
    end
    lower_scale = trial_scale;
end
if ~isfinite(upper_scale)
    return
end

for iteration = 1:cfg.speed.binary_search_max_iterations
    middle_scale = 0.5 * (lower_scale + upper_scale);
    [middle_distance, middle_ok] = transition_pair_distance( ...
        vs, peak, ve, start_jerk, peak_jerk, end_jerk, ...
        Ta, Td, middle_scale, cfg);
    if ~middle_ok
        return
    end
    if abs(middle_distance - target_distance) <= ...
            cfg.speed.distance_tolerance
        upper_scale = middle_scale;
        break
    elseif middle_distance >= target_distance
        upper_scale = middle_scale;
    else
        lower_scale = middle_scale;
    end
end

Ta = Ta * upper_scale;
Td = Td * upper_scale;
[La, acceleration_ok] = transition_distance_at_duration( ...
    vs, peak, start_jerk, peak_jerk, Ta, cfg);
[Ld, deceleration_ok] = transition_distance_at_duration( ...
    peak, ve, peak_jerk, end_jerk, Td, cfg);
distance_error = abs(La + Ld - target_distance);
feasible = acceleration_ok && deceleration_ok && ...
    distance_error <= max(cfg.speed.distance_tolerance, ...
    10 * eps(max(target_distance, 1)));
end

function [distance, feasible] = transition_pair_distance( ...
        vs, peak, ve, start_jerk, peak_jerk, end_jerk, ...
        Ta, Td, scale, cfg)
[La, acceleration_ok] = transition_distance_at_duration( ...
    vs, peak, start_jerk, peak_jerk, Ta * scale, cfg);
[Ld, deceleration_ok] = transition_distance_at_duration( ...
    peak, ve, peak_jerk, end_jerk, Td * scale, cfg);
distance = La + Ld;
feasible = acceleration_ok && deceleration_ok && isfinite(distance);
end

function [distance, feasible] = transition_distance_at_duration( ...
        v0, v1, j0, j1, duration, cfg)
feasible = generalized_transition_feasible( ...
    v0, v1, duration, j0, j1, cfg);
if ~feasible
    distance = Inf;
    return
end
[Sq, ~, ~, ~] = evaluate_septic_transition( ...
    [0; duration], 0, v0, v1, duration, j0, j1);
distance = Sq(end);
feasible = isfinite(distance) && distance >= -cfg.speed.distance_tolerance;
end

function [distance, time, feasible] = transition_metrics( ...
        v0, v1, j0, j1, cfg)
jerk_tolerance = max(cfg.limits.Jmax * 1e-12, ...
    cfg.speed.velocity_tolerance);
if abs(v1 - v0) <= cfg.speed.velocity_tolerance
    feasible = abs(j0) <= jerk_tolerance && abs(j1) <= jerk_tolerance;
    if feasible
        distance = 0;
        time = 0;
    else
        distance = Inf;
        time = Inf;
    end
elseif abs(j0) <= jerk_tolerance && abs(j1) <= jerk_tolerance
    [distance, time] = septic_transition_distance(v0, v1, cfg);
    feasible = true;
else
    [distance, time, feasible] = minimum_generalized_transition( ...
        v0, v1, j0, j1, cfg);
end
end

function [distance, best_time, feasible] = minimum_generalized_transition( ...
        v0, v1, j0, j1, cfg)
baseline_time = septic_transition_time(v0, v1, cfg);
minimum_time = max(cfg.speed.minimum_transition_time, ...
    cfg.lmsc_jerk_boundary.minimum_duration_fraction * baseline_time);
maximum_time = max(minimum_time, ...
    cfg.lmsc_jerk_boundary.maximum_duration_factor * baseline_time);
time_grid = linspace(minimum_time, maximum_time, ...
    cfg.lmsc_jerk_boundary.time_search_samples);
feasible_grid = false(size(time_grid));
for i = 1:numel(time_grid)
    feasible_grid(i) = generalized_transition_feasible( ...
        v0, v1, time_grid(i), j0, j1, cfg);
end
first_feasible = find(feasible_grid, 1, 'first');
if isempty(first_feasible)
    distance = Inf;
    best_time = Inf;
    feasible = false;
    return
end

best_time = time_grid(first_feasible);
if first_feasible > 1
    low = time_grid(first_feasible - 1);
    high = best_time;
    for iteration = 1:cfg.lmsc_jerk_boundary.time_refinement_iterations
        middle = 0.5 * (low + high);
        if generalized_transition_feasible(v0, v1, middle, j0, j1, cfg)
            high = middle;
        else
            low = middle;
        end
    end
    best_time = high;
end
[Sq, ~, ~, ~] = evaluate_septic_transition( ...
    [0; best_time], 0, v0, v1, best_time, j0, j1);
distance = Sq(end);
feasible = isfinite(distance) && distance >= -cfg.speed.distance_tolerance;
end

function feasible = generalized_transition_feasible( ...
        v0, v1, duration, j0, j1, cfg)
check_samples = max(cfg.speed.internal_constraint_check_samples, ...
    ceil(duration / cfg.interpolation.Ts) + 1);
local_time = linspace(0, duration, check_samples)';
[Sq, Vq, Aq, Jq] = evaluate_septic_transition( ...
    local_time, 0, v0, v1, duration, j0, j1);
tolerance = cfg.speed.velocity_tolerance;
lower_velocity = min(v0, v1) - tolerance;
upper_velocity = max(v0, v1) + tolerance;
if v1 >= v0
    monotone_velocity = all(diff(Vq) >= -tolerance);
else
    monotone_velocity = all(diff(Vq) <= tolerance);
end
feasible = all(isfinite([Sq; Vq; Aq; Jq])) && ...
    all(diff(Sq) >= -cfg.speed.distance_tolerance) && ...
    all(Vq >= max(lower_velocity, -tolerance)) && ...
    all(Vq <= upper_velocity) && monotone_velocity && ...
    max(abs(Aq)) <= cfg.limits.Amax * (1 + tolerance) && ...
    max(abs(Jq)) <= cfg.limits.Jmax * (1 + tolerance);
end

function validation = validate_segments(segments, S, v_limit, cfg)
validation.feasible = true;
validation.max_velocity_ratio = 0;
validation.max_acceleration_ratio = 0;
validation.max_jerk_ratio = 0;
validation.maximum_violation = 0;
validation.violation_S = NaN;
validation.violation_spu_id = 0;

for i = 1:numel(segments)
    segment = segments(i);
    duration = segment.t1 - segment.t0;
    time_based_samples = ceil(duration / cfg.interpolation.Ts) + 1;
    check_samples = max(cfg.speed.internal_constraint_check_samples, ...
        time_based_samples);
    local_time = linspace(0, duration, check_samples)';
    [sq, vq, aq, jq] = evaluate_one_segment(segment, local_time);
    local_limit = interp1(S, v_limit, sq, 'linear', 'extrap');
    violation = vq - local_limit;
    [maximum, idx] = max(violation);
    velocity_ratio = max(vq ./ max(local_limit, cfg.speed.velocity_tolerance));
    validation.max_velocity_ratio = max(validation.max_velocity_ratio, velocity_ratio);
    validation.max_acceleration_ratio = max(validation.max_acceleration_ratio, ...
        max(abs(aq)) / max(cfg.limits.Amax, eps));
    validation.max_jerk_ratio = max(validation.max_jerk_ratio, ...
        max(abs(jq)) / max(cfg.limits.Jmax, eps));
    if maximum > validation.maximum_violation
        validation.maximum_violation = maximum;
        validation.violation_S = sq(idx);
        validation.violation_spu_id = segment.spu_id;
    end
end

validation.feasible = validation.maximum_violation <= cfg.speed.velocity_tolerance && ...
    validation.max_acceleration_ratio <= 1 + cfg.speed.velocity_tolerance && ...
    validation.max_jerk_ratio <= 1 + cfg.speed.velocity_tolerance;
end

function time = fixed_time_axis(total_time, Ts)
time_tolerance = max(10 * eps(max(total_time, 1)), Ts * 1e-12);
if total_time <= time_tolerance
    time = [0; total_time];
    return
end
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
        ids = time >= segments(i).t0 & time < segments(i).t1;
    else
        ids = time >= segments(i).t0 & time <= segments(i).t1;
    end
    local_time = time(ids) - segments(i).t0;
    [S(ids), V(ids), A(ids), J(ids)] = evaluate_one_segment(segments(i), local_time);
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

function segments = empty_segments()
template = make_segment("constant", 0, 0, 0, 0, 0, 0, 0, 0, 0);
segments = repmat(template, 0, 1);
end

function segment = make_segment( ...
        type, spu_id, t0, t1, s0, s1, v0, v1, j0, j1)
segment.type = type;
segment.spu_id = spu_id;
segment.t0 = t0;
segment.t1 = t1;
segment.s0 = s0;
segment.s1 = s1;
segment.v0 = v0;
segment.v1 = v1;
segment.j0 = j0;
segment.j1 = j1;
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

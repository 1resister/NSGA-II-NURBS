function xyz = xyz_motion( ...
        sample_points, S, motion, cfg, arc_derivatives, constraint_only)
% Cartesian motion from q(s) derivatives and the common path time law.

S = S(:);
if nargin < 5
    arc_derivatives = [];
end
if nargin < 6
    constraint_only = false;
end
if ~isempty(arc_derivatives)
    q_s = arc_derivatives.q_s;
    q_ss = arc_derivatives.q_ss;
    q_sss = arc_derivatives.q_sss;
else
    q_s = zeros(size(sample_points));
    q_ss = zeros(size(sample_points));
    q_sss = zeros(size(sample_points));
    for axis_id = 1:3
        q_s(:, axis_id) = gradient(sample_points(:, axis_id), S);
        q_ss(:, axis_id) = gradient(q_s(:, axis_id), S);
        q_sss(:, axis_id) = gradient(q_ss(:, axis_id), S);
    end
end

if constraint_only
    xyz = constraint_summary(q_s, q_ss, q_sss, S, motion, cfg);
    return
end

points_t = interp1(S, sample_points, motion.S, 'pchip');
q_s_t = interp1(S, q_s, motion.S, 'pchip', 'extrap');
q_ss_t = interp1(S, q_ss, motion.S, 'pchip', 'extrap');
q_sss_t = interp1(S, q_sss, motion.S, 'pchip', 'extrap');

v = motion.path_velocity;
a = motion.path_acceleration;
j = motion.path_jerk;
velocity = q_s_t .* v;
acceleration = q_ss_t .* (v.^2) + q_s_t .* a;
jerk = q_sss_t .* (v.^3) + 3 * q_ss_t .* (v .* a) + q_s_t .* j;
jerk_geometry = q_sss_t .* (v.^3);
jerk_coupling = 3 * q_ss_t .* (v .* a);
jerk_tangential = q_s_t .* j;

points_t(1, :) = sample_points(1, :);
points_t(end, :) = sample_points(end, :);
velocity([1 end], :) = 0;
acceleration([1 end], :) = 0;
jerk([1 end], :) = 0;
jerk_geometry([1 end], :) = 0;
jerk_coupling([1 end], :) = 0;
jerk_tangential([1 end], :) = 0;

snap = zeros(size(jerk));
if numel(motion.time) >= 2
    sample_time = median(diff(motion.time));
    for axis_id = 1:3
        snap(:, axis_id) = gradient(jerk(:, axis_id), sample_time);
    end
end

velocity_ratio = abs(velocity) ./ cfg.limits.xyz_vmax;
acceleration_ratio = abs(acceleration) ./ cfg.limits.xyz_amax;
jerk_ratio = abs(jerk) ./ cfg.limits.xyz_jmax;
snap_limit = cfg.limits.xyz_smax(:)';
if numel(snap_limit) ~= 3 || any(isnan(snap_limit)) || ...
        any(snap_limit <= 0)
    error('cfg.limits.xyz_smax must contain three positive values.');
end
snap_ratio = abs(snap) ./ snap_limit;
snap_limit_enabled = xyz_snap_limit_enabled(cfg);
% Reduce one N-by-3 array at a time. Concatenating these four matrices would
% allocate an avoidable N-by-12 temporary and can exhaust memory for a long
% time-scaled candidate.
point_ratio = max(velocity_ratio, acceleration_ratio);
point_ratio = max(point_ratio, jerk_ratio);
if snap_limit_enabled
    point_ratio = max(point_ratio, snap_ratio);
end
point_ratio = max(point_ratio, [], 2);

xyz.points = points_t;
xyz.velocity = velocity;
xyz.acceleration = acceleration;
xyz.jerk = jerk;
xyz.snap = snap;
xyz.jerk_geometry = jerk_geometry;
xyz.jerk_coupling = jerk_coupling;
xyz.jerk_tangential = jerk_tangential;
xyz.q_s = q_s_t;
xyz.q_ss = q_ss_t;
xyz.q_sss = q_sss_t;
xyz.velocity_ratio = velocity_ratio;
xyz.acceleration_ratio = acceleration_ratio;
xyz.jerk_ratio = jerk_ratio;
xyz.snap_ratio = snap_ratio;
xyz.snap_limit_enabled = snap_limit_enabled;
xyz.point_constraint_ratio = point_ratio;
xyz.max_velocity_ratio = max(velocity_ratio, [], 'all');
xyz.max_acceleration_ratio = max(acceleration_ratio, [], 'all');
xyz.max_jerk_ratio = max(jerk_ratio, [], 'all');
xyz.max_snap_ratio = max(snap_ratio, [], 'all');
xyz.max_enforced_snap_ratio = 0;
if snap_limit_enabled
    xyz.max_enforced_snap_ratio = xyz.max_snap_ratio;
end
xyz.feasible = max(point_ratio) <= 1 + cfg.speed.velocity_tolerance;

% Optional filtered series are diagnostics only and never used for hard
% constraints or the Vibration objective.
xyz.filtered_jerk = jerk;
if isfield(cfg, 'smoothing') && cfg.smoothing.enable_xyz
    win = min(cfg.smoothing.window, size(jerk, 1));
    if mod(win, 2) == 0
        win = win - 1;
    end
    if win >= 5
        xyz.filtered_jerk = smoothdata(jerk, 1, 'sgolay', win);
    end
end
end

function xyz = constraint_summary(q_s, q_ss, q_sss, S, motion, cfg)
% Evaluate only the scalar XYZ constraint summary used during scale search.
% Plotting/export series are deliberately omitted to keep bisection memory
% independent of the previously retained feasible trajectory.
q_s_t = interp1(S, q_s, motion.S, 'pchip', 'extrap');
v = motion.path_velocity;
a = motion.path_acceleration;
j = motion.path_jerk;

velocity = q_s_t .* v;
xyz.max_velocity_ratio = max( ...
    abs(velocity) ./ cfg.limits.xyz_vmax, [], 'all');
clear velocity

q_ss_t = interp1(S, q_ss, motion.S, 'pchip', 'extrap');
acceleration = q_ss_t .* (v.^2) + q_s_t .* a;
xyz.max_acceleration_ratio = max( ...
    abs(acceleration) ./ cfg.limits.xyz_amax, [], 'all');
clear acceleration

q_sss_t = interp1(S, q_sss, motion.S, 'pchip', 'extrap');
jerk = q_sss_t .* (v.^3) + 3 * q_ss_t .* (v .* a) + q_s_t .* j;
jerk([1 end], :) = 0;
xyz.max_jerk_ratio = max( ...
    abs(jerk) ./ cfg.limits.xyz_jmax, [], 'all');

snap_limit = cfg.limits.xyz_smax(:)';
if numel(snap_limit) ~= 3 || any(isnan(snap_limit)) || ...
        any(snap_limit <= 0)
    error('cfg.limits.xyz_smax must contain three positive values.');
end
xyz.max_snap_ratio = 0;
if numel(motion.time) >= 2
    sample_time = median(diff(motion.time));
    for axis_id = 1:3
        snap_axis = gradient(jerk(:, axis_id), sample_time);
        xyz.max_snap_ratio = max(xyz.max_snap_ratio, ...
            max(abs(snap_axis)) / snap_limit(axis_id));
    end
end

xyz.snap_limit_enabled = xyz_snap_limit_enabled(cfg);
xyz.max_enforced_snap_ratio = 0;
if xyz.snap_limit_enabled
    xyz.max_enforced_snap_ratio = xyz.max_snap_ratio;
end
maximum_ratio = max([xyz.max_velocity_ratio, ...
    xyz.max_acceleration_ratio, xyz.max_jerk_ratio, ...
    xyz.max_enforced_snap_ratio]);
xyz.feasible = maximum_ratio <= 1 + cfg.speed.velocity_tolerance;
end

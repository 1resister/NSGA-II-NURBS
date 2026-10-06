function constraints = generate_velocity_constraints( ...
        kappa, tangent, cfg, S, speed_scale_nodes, arc_derivatives)
% Build all pointwise speed-constraint components in the arc-length domain.

kappa = max(kappa(:), 0);
n = numel(kappa);
if nargin < 4 || isempty(S)
    S = linspace(0, 1, n)';
else
    S = S(:);
end
if nargin < 5
    speed_scale_nodes = [];
end
if nargin < 6
    arc_derivatives = [];
end

epsilon_kappa = cfg.speed.epsilon_kappa;
v_command = cfg.limits.Vmax * ones(n, 1);
v_curvature = sqrt(cfg.limits.Acmax ./ max(kappa, epsilon_kappa));
v_curvature(kappa <= epsilon_kappa) = Inf;

v_chord = inf(n, 1);
if cfg.speed.enable_chord_constraint && ...
        isfield(cfg, 'interpolation') && isfield(cfg.interpolation, 'Ts') && ...
        isfield(cfg.constraints, 'max_chord_error') && ...
        cfg.interpolation.Ts > 0 && cfg.constraints.max_chord_error > 0
    delta = cfg.constraints.max_chord_error;
    chord_term = 2 * delta ./ max(kappa, epsilon_kappa) - delta^2;
    v_chord = 2 / cfg.interpolation.Ts * sqrt(max(chord_term, 0));
    v_chord(kappa <= epsilon_kappa) = Inf;
end

tangent_limit = tangent;
tangent_limit = tangent_limit ./ max(vecnorm(tangent_limit, 2, 2), 1e-12);
v_axis_components = inf(n, 3);
for axis_id = 1:3
    v_axis_components(:, axis_id) = cfg.limits.xyz_vmax(axis_id) ./ ...
        max(abs(tangent_limit(:, axis_id)), 1e-12);
end
v_axis = min(v_axis_components, [], 2);

v_geometric_jerk_components = inf(n, 3);
if ~isempty(arc_derivatives) && isfield(arc_derivatives, 'q_sss')
    for axis_id = 1:3
        geometric_rate = abs(arc_derivatives.q_sss(:, axis_id));
        v_geometric_jerk_components(:, axis_id) = ...
            nthroot(cfg.limits.xyz_jmax(axis_id) ./ ...
            max(geometric_rate, 1e-12), 3);
        v_geometric_jerk_components(geometric_rate <= 1e-12, axis_id) = Inf;
    end
end
v_geometric_jerk = min(v_geometric_jerk_components, [], 2);

process_limit = cfg.speed.process_speed_limit;
if isscalar(process_limit)
    v_process = process_limit * ones(n, 1);
elseif numel(process_limit) == n
    v_process = process_limit(:);
else
    source_S = linspace(S(1), S(end), numel(process_limit))';
    v_process = interp1(source_S, process_limit(:), S, 'linear', 'extrap');
end
v_process = max(v_process, 0);

if isempty(speed_scale_nodes)
    speed_scale = ones(n, 1);
elseif isscalar(speed_scale_nodes)
    speed_scale = speed_scale_nodes * ones(n, 1);
else
    node_S = linspace(S(1), S(end), numel(speed_scale_nodes))';
    speed_scale = interp1(node_S, speed_scale_nodes(:), S, 'pchip', 'extrap');
end
speed_scale = max(speed_scale, 0);
v_scale = cfg.limits.Vmax * speed_scale;

physical_limits = [v_command, v_curvature, v_chord, v_axis, ...
    v_geometric_jerk, v_process];
[v_limit, active_id] = min(physical_limits, [], 2);
v_limit(~isfinite(v_limit)) = cfg.limits.Vmax;
v_limit = max(v_limit, cfg.speed.velocity_tolerance);
v_candidate = min(v_limit, v_scale);
v_candidate = max(v_candidate, cfg.speed.velocity_tolerance);

active_names = ["command", "curvature", "chord", "axis", ...
    "geometric_jerk", "process"];
base_active_constraint = reshape(active_names(active_id), [], 1);
active_constraint = base_active_constraint;
active_constraint(v_scale < v_limit - cfg.speed.velocity_tolerance) = "speed_scale";
constraints.S = S;
constraints.v_command = v_command;
constraints.v_curvature = v_curvature;
constraints.v_chord = v_chord;
constraints.v_axis = v_axis;
constraints.v_axis_components = v_axis_components;
constraints.v_geometric_jerk = v_geometric_jerk;
constraints.v_geometric_jerk_components = v_geometric_jerk_components;
constraints.v_process = v_process;
constraints.speed_scale = speed_scale;
constraints.v_scale = v_scale;
constraints.v_candidate = v_candidate;
constraints.v_limit = v_limit;
constraints.active_constraint_id = active_id;
constraints.base_active_constraint = base_active_constraint;
constraints.active_constraint = active_constraint;

% Legacy field names retained for existing callers.
constraints.Vlimit = constraints.v_candidate;
constraints.Vc = constraints.v_curvature;
end

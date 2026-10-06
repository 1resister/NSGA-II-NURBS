function [smooth_path, curve, u] = nurbs_global_smoothing(raw_points, weights, time_scale, cfg)
% Dynamic baselines use fixed high-continuity corner knots. Speed nodes do
% not alter the fitted geometry; the legacy generator is only a fallback.
% 五次 NURBS 全局光顺。
%
% 控制点由原始离散点按弦长参数自动重采样得到；权重来自 NSGA-II 染色体。
% 节点向量由时间节点倍率生成，使优化能够影响曲线参数分布。

p = cfg.nurbs.degree;
nctrl = cfg.nurbs.num_control_points;
if isfield(cfg.nurbs, 'base_control_points') && ...
        size(cfg.nurbs.base_control_points, 1) == nctrl && ...
        isfield(cfg.nurbs, 'base_knots')
    control_points = cfg.nurbs.base_control_points;
    knotvec = cfg.nurbs.base_knots;
else
    raw_s = cumulative_distance(raw_points);
    target_s = control_point_parameters(raw_points, raw_s, nctrl, cfg);
    method = cfg.nurbs.control_point_interpolation;
    control_points = interp1(raw_s, raw_points, target_s, method);
    knotvec = optimized_open_knot_vector(nctrl, p, time_scale);
end

if isfield(cfg.nurbs, 'force_unit_weights') && ...
        cfg.nurbs.force_unit_weights
    weights = ones(1, nctrl);
else
    weights = max(weights(:)', 1e-6);
end

curve.control_points = control_points;
curve.weights = weights;
curve.degree = p;
curve.knots = knotvec;

base_u = linspace(0, 1, cfg.nurbs.eval_points)';
unique_knots = unique(knotvec(:));
span_midpoints = 0.5 * (unique_knots(1:end-1) + unique_knots(2:end));
u = unique([base_u; unique_knots; span_midpoints])';
smooth_path = evaluate_nurbs_curve(control_points, weights, p, knotvec, u);
end

function target_s = control_point_parameters(raw_points, raw_s, nctrl, cfg)
if ~cfg.nurbs.adaptive_control_point_distribution || ...
        size(raw_points, 1) < 3
    target_s = linspace(0, raw_s(end), nctrl);
    return
end

segments = diff(raw_points, 1, 1);
segment_length = sqrt(sum(segments.^2, 2));
unit_tangent = segments ./ max(segment_length, eps);
cosine = sum(unit_tangent(1:end-1, :) .* unit_tangent(2:end, :), 2);
turning_angle = [0; acos(max(min(cosine, 1), -1)); 0];
vertex_density = 1 + cfg.nurbs.corner_density_gain * turning_angle / pi;
segment_density = 0.5 * (vertex_density(1:end-1) + vertex_density(2:end));
adaptive_coordinate = [0; cumsum(segment_length .* segment_density)];

if adaptive_coordinate(end) <= eps
    target_s = linspace(0, raw_s(end), nctrl);
else
    query = linspace(0, adaptive_coordinate(end), nctrl);
    target_s = interp1(adaptive_coordinate, raw_s, query, 'linear');
    target_s([1 end]) = raw_s([1 end]);
end
end

function s = cumulative_distance(points)
d = [0; vecnorm(diff(points), 2, 2)];
s = cumsum(d);
if s(end) < eps
    s = (0:size(points, 1)-1)';
end
end

function knotvec = optimized_open_knot_vector(nctrl, p, time_scale)
% 构造开放节点向量，内部节点由时间倍率的累计比例确定。
n_internal = nctrl - p - 1;
if n_internal <= 0
    knotvec = [zeros(1, p + 1) ones(1, p + 1)];
    return
end
t = abs(time_scale(:)');
if isempty(t) || sum(t) < eps
    internal = linspace(0, 1, n_internal + 2);
    internal = internal(2:end-1);
else
    tq = interp1(linspace(0, 1, numel(t)), t, linspace(0, 1, n_internal), 'linear', 'extrap');
    internal = cumsum(max(tq, 1e-6));
    internal = internal / (internal(end) + max(tq(end), 1e-6));
end
knotvec = [zeros(1, p + 1) internal ones(1, p + 1)];
end

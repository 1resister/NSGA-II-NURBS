function cfg = configure_dynamic_dimensions(cfg)
% Select fixed NSGA-II dimensions dynamically for the current input path.
%
% The dimensions are selected once per input/configuration and cached.  They
% remain fixed during an optimization run, which is required by NSGA-II.

defaults.control_points = cfg.dynamic_dimensions.fallback_control_points;
defaults.speed_nodes = cfg.dynamic_dimensions.fallback_speed_nodes;
defaults.baseline_error = Inf;
defaults.target_met = false;
defaults.significant_features = 0;
defaults.path_points = 0;
defaults.base_control_points = [];
defaults.base_knots = [];

if ~cfg.dynamic_dimensions.enabled || ~isfile(cfg.files.input_csv)
    cfg = apply_dimensions(cfg, defaults);
    return
end

info = dir(cfg.files.input_csv);
dimension_settings = cfg.dynamic_dimensions;
nurbs_settings = struct('degree', cfg.nurbs.degree, ...
    'corner_density_gain', cfg.nurbs.corner_density_gain, ...
    'corner_knot_width_factor', cfg.nurbs.corner_knot_width_factor, ...
    'corner_knot_multiplicity', cfg.nurbs.corner_knot_multiplicity, ...
    'force_unit_weights', cfg.nurbs.force_unit_weights, ...
    'bidirectional_contour_error', cfg.nurbs.bidirectional_contour_error, ...
    'fit_samples_per_segment_max', cfg.nurbs.fit_samples_per_segment_max, ...
    'fit_regularization', cfg.nurbs.fit_regularization, ...
    'fit_third_difference_regularization', ...
        cfg.nurbs.fit_third_difference_regularization, ...
    'fit_fourth_difference_regularization', ...
        cfg.nurbs.fit_fourth_difference_regularization, ...
    'control_point_interpolation', cfg.nurbs.control_point_interpolation);
cache_key = sprintf('%s|%.15g|%d|%.15g|%s|%s', ...
    cfg.files.input_csv, info.datenum, info.bytes, ...
    cfg.constraints.max_contour_error, jsonencode(dimension_settings), ...
    jsonencode(nurbs_settings));

persistent cached_key cached_dimensions
if ~isempty(cached_key) && strcmp(cached_key, cache_key)
    cfg = apply_dimensions(cfg, cached_dimensions);
    return
end

raw_points = readmatrix(cfg.files.input_csv);
if size(raw_points, 2) < 3
    cfg = apply_dimensions(cfg, defaults);
    return
end
raw_points = raw_points(:, 1:3);
raw_points = raw_points(all(isfinite(raw_points), 2), :);
if size(raw_points, 1) >= 2
    segment_length = sqrt(sum(diff(raw_points, 1, 1).^2, 2));
    raw_points = raw_points([true; segment_length > 1e-10], :);
end
if size(raw_points, 1) < 2
    cfg = apply_dimensions(cfg, defaults);
    return
end

feature_count = count_significant_features(raw_points, ...
    cfg.dynamic_dimensions.feature_angle_threshold_deg);
minimum_count = max([cfg.nurbs.degree + 1, ...
    cfg.dynamic_dimensions.min_control_points, ...
    cfg.nurbs.degree + 1 + ...
        cfg.nurbs.corner_knot_multiplicity * feature_count]);
maximum_count = max(cfg.nurbs.degree + 1, ...
    cfg.dynamic_dimensions.max_control_points);
minimum_count = min(minimum_count, maximum_count);
step = max(1, round(cfg.dynamic_dimensions.control_point_step));
candidates = unique([minimum_count:step:maximum_count, maximum_count]);
target_error = cfg.constraints.max_contour_error * ...
    cfg.dynamic_dimensions.baseline_error_fraction;

selected_count = maximum_count;
selected_error = Inf;
selected_curve = struct('control_points', [], 'knots', []);
for candidate = candidates
    [selected_curve, selected_error] = ...
        fit_adaptive_nurbs_path(raw_points, candidate, cfg);
    selected_count = candidate;
    if selected_error <= target_error
        break
    end
end

nodes_from_control_points = ceil(selected_count / ...
    cfg.dynamic_dimensions.control_points_per_speed_node);
nodes_from_features = cfg.dynamic_dimensions.speed_nodes_per_feature * ...
    feature_count + 3;
speed_nodes = max([cfg.dynamic_dimensions.min_speed_nodes, ...
    nodes_from_control_points, nodes_from_features]);
speed_nodes = min(speed_nodes, cfg.dynamic_dimensions.max_speed_nodes);

dimensions.control_points = selected_count;
dimensions.speed_nodes = speed_nodes;
dimensions.baseline_error = selected_error;
dimensions.target_met = selected_error <= target_error;
dimensions.significant_features = feature_count;
dimensions.path_points = size(raw_points, 1);
dimensions.base_control_points = selected_curve.control_points;
dimensions.base_knots = selected_curve.knots;

cached_key = cache_key;
cached_dimensions = dimensions;
cfg = apply_dimensions(cfg, dimensions);
end

function cfg = apply_dimensions(cfg, dimensions)
cfg.nurbs.num_control_points = dimensions.control_points;
cfg.nsga.num_time_nodes = dimensions.speed_nodes;
cfg.dynamic_dimensions.selected_control_points = dimensions.control_points;
cfg.dynamic_dimensions.selected_speed_nodes = dimensions.speed_nodes;
cfg.dynamic_dimensions.baseline_error = dimensions.baseline_error;
cfg.dynamic_dimensions.baseline_target_met = dimensions.target_met;
cfg.dynamic_dimensions.significant_features = dimensions.significant_features;
cfg.dynamic_dimensions.path_points = dimensions.path_points;
if ~isempty(dimensions.base_control_points)
    cfg.nurbs.base_control_points = dimensions.base_control_points;
    cfg.nurbs.base_knots = dimensions.base_knots;
end
end

function feature_count = count_significant_features(points, threshold_deg)
if size(points, 1) < 3
    feature_count = 0;
    return
end
segments = diff(points, 1, 1);
lengths = sqrt(sum(segments.^2, 2));
unit_tangent = segments ./ max(lengths, eps);
cosine = sum(unit_tangent(1:end-1, :) .* unit_tangent(2:end, :), 2);
angles_deg = acos(max(min(cosine, 1), -1)) * 180 / pi;
feature_count = sum(angles_deg >= threshold_deg);
end

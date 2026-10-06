function [curve, fit_error, diagnostics] = ...
        fit_adaptive_nurbs_path(raw_points, nctrl, cfg)
% Fit a high-continuity degree-p B-spline baseline to a polyline.

degree = cfg.nurbs.degree;
[vertex_u, segment_length] = chord_parameters(raw_points);
corner_ids = significant_corners(raw_points, ...
    cfg.dynamic_dimensions.feature_angle_threshold_deg);
target_error = max(cfg.constraints.max_contour_error, eps);
total_length = sum(segment_length);
normalized_width = cfg.nurbs.corner_knot_width_factor * ...
    target_error / max(total_length, eps);
knots = clustered_corner_knots(nctrl, degree, vertex_u, corner_ids, ...
    normalized_width, cfg.nurbs.corner_knot_multiplicity);

[fit_u, fit_points] = polyline_fit_samples(raw_points, vertex_u, ...
    segment_length, target_error, cfg.nurbs.fit_samples_per_segment_max);
basis = bspline_basis_matrix(fit_u, nctrl, degree, knots);

control_points = zeros(nctrl, 3);
control_points(1, :) = raw_points(1, :);
control_points(end, :) = raw_points(end, :);
interior_basis = basis(:, 2:end-1);
right_hand_side = fit_points - ...
    basis(:, 1) * control_points(1, :) - ...
    basis(:, end) * control_points(end, :);
regularization = cfg.nurbs.fit_regularization;
regularization_rows = sqrt(max(regularization, 0)) * eye(nctrl - 2);
third_difference = diff(eye(nctrl), 3, 1);
fairness_weight = cfg.nurbs.fit_third_difference_regularization;
fairness_rows = sqrt(max(fairness_weight, 0)) * ...
    third_difference(:, 2:end-1);
fairness_endpoint = sqrt(max(fairness_weight, 0)) * ...
    (third_difference(:, 1) * control_points(1, :) + ...
    third_difference(:, end) * control_points(end, :));
fourth_difference = diff(eye(nctrl), 4, 1);
fourth_weight = cfg.nurbs.fit_fourth_difference_regularization;
fourth_rows = sqrt(max(fourth_weight, 0)) * ...
    fourth_difference(:, 2:end-1);
fourth_endpoint = sqrt(max(fourth_weight, 0)) * ...
    (fourth_difference(:, 1) * control_points(1, :) + ...
    fourth_difference(:, end) * control_points(end, :));
control_points(2:end-1, :) = ...
    [interior_basis; regularization_rows; fairness_rows; fourth_rows] \ ...
    [right_hand_side; zeros(nctrl - 2, 3); ...
        -fairness_endpoint; -fourth_endpoint];

curve.control_points = control_points;
curve.weights = ones(1, nctrl);
curve.degree = degree;
curve.knots = knots;

unique_knots = unique(knots(:));
span_midpoints = 0.5 * (unique_knots(1:end-1) + unique_knots(2:end));
coarse_u = unique([linspace(0, 1, ...
    cfg.dynamic_dimensions.probe_eval_points)'; unique_knots; span_midpoints]);
[vertex_error, ~, ~, curve_error] = nurbs_contour_error( ...
    raw_points, curve, coarse_u, cfg.nurbs.bidirectional_contour_error);
matched_points = basis * control_points;
matched_error = max(sqrt(sum((matched_points - fit_points).^2, 2)));
fit_error = max([vertex_error, curve_error, matched_error]);

diagnostics.vertex_error = vertex_error;
diagnostics.curve_to_polyline_error = curve_error;
diagnostics.matched_polyline_error = matched_error;
diagnostics.corner_count = numel(corner_ids);
diagnostics.corner_knot_multiplicity = cfg.nurbs.corner_knot_multiplicity;
diagnostics.continuity_order = degree - cfg.nurbs.corner_knot_multiplicity;
diagnostics.fit_sample_count = numel(fit_u);
end

function [vertex_u, segment_length] = chord_parameters(points)
segment_length = sqrt(sum(diff(points, 1, 1).^2, 2));
distance = [0; cumsum(segment_length)];
vertex_u = distance / max(distance(end), eps);
end

function corner_ids = significant_corners(points, threshold_deg)
if size(points, 1) < 3
    corner_ids = zeros(0, 1);
    return
end
segments = diff(points, 1, 1);
lengths = sqrt(sum(segments.^2, 2));
unit_tangent = segments ./ max(lengths, eps);
cosine = sum(unit_tangent(1:end-1, :) .* unit_tangent(2:end, :), 2);
angle = acos(max(min(cosine, 1), -1)) * 180 / pi;
corner_ids = find(angle >= threshold_deg) + 1;
end

function knots = clustered_corner_knots(nctrl, degree, vertex_u, ...
        corner_ids, normalized_width, corner_multiplicity)
n_internal = nctrl - degree - 1;
if n_internal <= 0
    knots = [zeros(1, degree + 1), ones(1, degree + 1)];
    return
end

corner_u = vertex_u(corner_ids);
internal = zeros(1, n_internal);
count = 0;
if ~isscalar(corner_multiplicity) || ~isfinite(corner_multiplicity) || ...
        corner_multiplicity < 1 || corner_multiplicity > degree || ...
        corner_multiplicity ~= floor(corner_multiplicity)
    error('corner_knot_multiplicity must be an integer in [1, degree].');
end
for repetition = 1:corner_multiplicity
    for corner_id = 1:numel(corner_u)
        if count >= n_internal
            break
        end
        count = count + 1;
        internal(count) = corner_u(corner_id);
    end
end

ring = 1;
attempts = 0;
while count < n_internal && ~isempty(corner_u) && attempts < 1000
    for corner_id = 1:numel(corner_u)
        center = corner_u(corner_id);
        if corner_id == 1
            lower_bound = 0;
        else
            lower_bound = 0.5 * (corner_u(corner_id - 1) + center);
        end
        if corner_id == numel(corner_u)
            upper_bound = 1;
        else
            upper_bound = 0.5 * (center + corner_u(corner_id + 1));
        end
        maximum_offset = 0.45 * min(center - lower_bound, ...
            upper_bound - center);
        offset = min(ring * max(normalized_width, 1e-8), maximum_offset);
        for direction = [-1, 1]
            if count >= n_internal
                break
            end
            candidate = center + direction * offset;
            if candidate > 0 && candidate < 1 && ...
                    ~any(abs(internal(1:count) - candidate) <= 1e-12)
                count = count + 1;
                internal(count) = candidate;
            end
        end
    end
    ring = ring + 1;
    attempts = attempts + 1;
end

if count < n_internal
    uniform = linspace(0, 1, 4 * n_internal + 2);
    uniform = uniform(2:end-1);
    for candidate = uniform
        if count >= n_internal
            break
        end
        if ~any(abs(internal(1:count) - candidate) <= 1e-12)
            count = count + 1;
            internal(count) = candidate;
        end
    end
end
internal = sort(internal(1:n_internal));
knots = [zeros(1, degree + 1), internal, ones(1, degree + 1)];
end

function [parameters, points] = polyline_fit_samples(raw_points, vertex_u, ...
        segment_length, target_error, maximum_per_segment)
parameters = zeros(0, 1);
points = zeros(0, 3);
spacing = max(0.5 * target_error, eps);
for segment_id = 1:numel(segment_length)
    count = max(20, ceil(segment_length(segment_id) / spacing));
    count = min(count, maximum_per_segment);
    fraction = (0:count-1)' / count;
    parameters = [parameters; vertex_u(segment_id) + ...
        fraction * (vertex_u(segment_id + 1) - vertex_u(segment_id))]; %#ok<AGROW>
    points = [points; raw_points(segment_id, :) + ...
        fraction .* (raw_points(segment_id + 1, :) - ...
        raw_points(segment_id, :))]; %#ok<AGROW>
end
parameters(end + 1, 1) = 1;
points(end + 1, :) = raw_points(end, :);
end

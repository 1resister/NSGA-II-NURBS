function [err, point_error, projected_u, curve_error, ...
        curve_point_error] = nurbs_contour_error( ...
        raw_points, curve, coarse_u, bidirectional)
% Maximum contour error, optionally measured in both directions.

% The original path points are first matched to a coarse parameter grid.
% Each match is then refined with a bounded one-dimensional minimization,
% avoiding dependence on the number of plotted curve samples.

if nargin < 3 || isempty(coarse_u)
    coarse_u = linspace(0, 1, max(101, 4 * numel(curve.weights)))';
else
    coarse_u = unique(max(min(coarse_u(:), 1), 0));
end
if nargin < 4 || isempty(bidirectional)
    bidirectional = false;
end
if numel(coarse_u) < 3
    coarse_u = linspace(0, 1, 101)';
end

coarse_points = evaluate_nurbs_curve(curve.control_points, curve.weights, ...
    curve.degree, curve.knots, coarse_u);
n = size(raw_points, 1);
point_error = zeros(n, 1);
projected_u = zeros(n, 1);
options = optimset('Display', 'off', 'TolX', 1e-10);

for point_id = 1:n
    point = raw_points(point_id, :);
    squared_distance = sum((coarse_points - point).^2, 2);
    [coarse_minimum, nearest_id] = min(squared_distance);
    left_id = max(nearest_id - 1, 1);
    right_id = min(nearest_id + 1, numel(coarse_u));
    lower = coarse_u(left_id);
    upper = coarse_u(right_id);

    if upper - lower <= eps
        best_u = coarse_u(nearest_id);
        best_squared_distance = coarse_minimum;
    else
        objective = @(parameter) point_distance_squared( ...
            parameter, point, curve);
        [best_u, best_squared_distance] = fminbnd( ...
            objective, lower, upper, options);
        if coarse_minimum < best_squared_distance
            best_u = coarse_u(nearest_id);
            best_squared_distance = coarse_minimum;
        end
    end

    projected_u(point_id) = best_u;
    point_error(point_id) = sqrt(max(best_squared_distance, 0));
end

curve_error = 0;
curve_point_error = zeros(0, 1);
if bidirectional
    curve_point_error = point_to_polyline_distance(coarse_points, raw_points);
    curve_error = max(curve_point_error);
end
err = max([point_error; curve_error]);
end

function value = point_distance_squared(parameter, point, curve)
curve_point = evaluate_nurbs_curve(curve.control_points, curve.weights, ...
    curve.degree, curve.knots, parameter);
difference = curve_point - point;
value = sum(difference.^2);
end

function distance = point_to_polyline_distance(points, polyline)
% Distance from each sampled spline point to the original input polyline.
if size(polyline, 1) == 1
    distance = vecnorm(points - polyline, 2, 2);
    return
end

minimum_squared_distance = inf(size(points, 1), 1);
for segment_id = 1:size(polyline, 1) - 1
    start_point = polyline(segment_id, :);
    segment = polyline(segment_id + 1, :) - start_point;
    segment_squared_length = sum(segment.^2);
    if segment_squared_length <= eps
        squared_distance = sum((points - start_point).^2, 2);
    else
        parameter = ((points - start_point) * segment') / ...
            segment_squared_length;
        parameter = min(max(parameter, 0), 1);
        projection = start_point + parameter .* segment;
        squared_distance = sum((points - projection).^2, 2);
    end
    minimum_squared_distance = min(minimum_squared_distance, ...
        squared_distance);
end
distance = sqrt(max(minimum_squared_distance, 0));
end

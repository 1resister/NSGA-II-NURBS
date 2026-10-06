function derivatives = evaluate_nurbs_derivatives( ...
        control_points, weights, degree, knots, parameters, maximum_order)
% Evaluate a NURBS curve and its parameter derivatives analytically.
if nargin < 6
    maximum_order = 3;
end
maximum_order = min(maximum_order, degree);
parameters = parameters(:);
weights = weights(:);
if numel(weights) ~= size(control_points, 1)
    error('The number of NURBS weights must equal the control-point count.');
end

homogeneous_control = [control_points .* weights, weights];
homogeneous_derivatives = zeros(numel(parameters), ...
    size(homogeneous_control, 2), maximum_order + 1);
derivative_control = homogeneous_control;
derivative_degree = degree;
derivative_knots = knots(:)';

for order = 0:maximum_order
    basis = bspline_basis_matrix(parameters, size(derivative_control, 1), ...
        derivative_degree, derivative_knots);
    homogeneous_derivatives(:, :, order + 1) = ...
        basis * derivative_control;
    if order < maximum_order
        derivative_control = differentiate_control_polygon( ...
            derivative_control, derivative_degree, derivative_knots);
        derivative_degree = derivative_degree - 1;
        derivative_knots = derivative_knots(2:end-1);
    end
end

dimension = size(control_points, 2);
derivatives = zeros(numel(parameters), dimension, maximum_order + 1);
weight_derivatives = homogeneous_derivatives(:, end, :);
base_weight = max(weight_derivatives(:, 1, 1), eps);
for order = 0:maximum_order
    value = homogeneous_derivatives(:, 1:dimension, order + 1);
    for lower_order = 1:order
        value = value - nchoosek(order, lower_order) .* ...
            weight_derivatives(:, 1, lower_order + 1) .* ...
            derivatives(:, :, order - lower_order + 1);
    end
    derivatives(:, :, order + 1) = value ./ base_weight;
end
end

function derivative_control = differentiate_control_polygon( ...
        control, degree, knots)
count = size(control, 1) - 1;
derivative_control = zeros(count, size(control, 2));
for control_id = 1:count
    denominator = knots(control_id + degree + 1) - ...
        knots(control_id + 1);
    if abs(denominator) > eps
        derivative_control(control_id, :) = degree * ...
            (control(control_id + 1, :) - control(control_id, :)) / ...
            denominator;
    end
end
end

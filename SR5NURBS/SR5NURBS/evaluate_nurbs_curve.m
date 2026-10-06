function points = evaluate_nurbs_curve(control_points, weights, degree, knots, u)
% Evaluate a NURBS curve with the iterative local Cox-de Boor algorithm.
% This is mathematically equivalent to the previous recursive evaluator but
% evaluates only the degree+1 nonzero basis functions at each parameter.

nctrl = size(control_points, 1);
weights = weights(:)';
u = u(:);
points = zeros(numel(u), size(control_points, 2));

for sample_id = 1:numel(u)
    parameter = min(max(u(sample_id), knots(1)), knots(end));
    if parameter >= knots(end)
        points(sample_id, :) = control_points(end, :);
        continue
    end

    span = find(parameter >= knots, 1, 'last') - 1; % zero-based span
    span = min(max(span, degree), nctrl - 1);
    basis = local_basis_functions(span, parameter, degree, knots);
    control_ids = (span - degree + 1):(span + 1);
    weighted_basis = basis .* weights(control_ids);
    denominator = sum(weighted_basis);
    if denominator <= eps
        points(sample_id, :) = control_points(control_ids(1), :);
    else
        points(sample_id, :) = ...
            (weighted_basis * control_points(control_ids, :)) / denominator;
    end
end

points(1, :) = control_points(1, :);
points(end, :) = control_points(end, :);
end

function basis = local_basis_functions(span, parameter, degree, knots)
basis = zeros(1, degree + 1);
basis(1) = 1;
left = zeros(1, degree);
right = zeros(1, degree);

for j = 1:degree
    left(j) = parameter - knots(span + 2 - j);
    right(j) = knots(span + 1 + j) - parameter;
    saved = 0;
    for r = 0:j-1
        denominator = right(r + 1) + left(j - r);
        if abs(denominator) <= eps
            temporary = 0;
        else
            temporary = basis(r + 1) / denominator;
        end
        basis(r + 1) = saved + right(r + 1) * temporary;
        saved = left(j - r) * temporary;
    end
    basis(j + 1) = saved;
end
end

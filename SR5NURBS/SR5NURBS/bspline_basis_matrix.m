function matrix = bspline_basis_matrix(parameters, nctrl, degree, knots)
% Evaluate all nonzero B-spline basis functions on a parameter vector.

parameters = parameters(:);
matrix = zeros(numel(parameters), nctrl);
for sample_id = 1:numel(parameters)
    parameter = min(max(parameters(sample_id), knots(degree + 1)), ...
        knots(nctrl + 1));
    if parameter >= knots(nctrl + 1)
        span = nctrl - 1;
    else
        span = find_knot_span(nctrl - 1, degree, parameter, knots);
    end
    basis = local_basis_functions(span, parameter, degree, knots);
    control_ids = (span - degree + 1):(span + 1);
    matrix(sample_id, control_ids) = basis;
end
end

function span = find_knot_span(n, degree, parameter, knots)
low = degree;
high = n + 1;
middle = floor((low + high) / 2);
while parameter < knots(middle + 1) || parameter >= knots(middle + 2)
    if parameter < knots(middle + 1)
        high = middle;
    else
        low = middle;
    end
    middle = floor((low + high) / 2);
end
span = middle;
end

function basis = local_basis_functions(span, parameter, degree, knots)
basis = zeros(1, degree + 1);
left = zeros(1, degree + 1);
right = zeros(1, degree + 1);
basis(1) = 1;
for order = 1:degree
    left(order + 1) = parameter - knots(span - order + 2);
    right(order + 1) = knots(span + order + 1) - parameter;
    saved = 0;
    for index = 0:order-1
        denominator = right(index + 2) + left(order - index + 1);
        if abs(denominator) <= eps
            value = 0;
        else
            value = basis(index + 1) / denominator;
        end
        basis(index + 1) = saved + right(index + 2) * value;
        saved = left(order - index + 1) * value;
    end
    basis(order + 1) = saved;
end
end

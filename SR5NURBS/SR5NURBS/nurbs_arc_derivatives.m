function arc_derivatives = nurbs_arc_derivatives(curve, parameters)
% Convert analytic NURBS parameter derivatives to arc-length derivatives.
parameter_derivatives = evaluate_nurbs_derivatives( ...
    curve.control_points, curve.weights, curve.degree, curve.knots, ...
    parameters, 3);

points = parameter_derivatives(:, :, 1);
r1 = parameter_derivatives(:, :, 2);
r2 = parameter_derivatives(:, :, 3);
r3 = parameter_derivatives(:, :, 4);

speed_u = max(vecnorm(r1, 2, 2), 1e-12);
r1_dot_r2 = sum(r1 .* r2, 2);
speed_u_first = r1_dot_r2 ./ speed_u;
speed_u_second = (sum(r2.^2, 2) + sum(r1 .* r3, 2)) ./ speed_u - ...
    r1_dot_r2.^2 ./ speed_u.^3;

q_s = r1 ./ speed_u;
q_ss = r2 ./ speed_u.^2 - ...
    r1 .* speed_u_first ./ speed_u.^3;
q_sss = r3 ./ speed_u.^3 - ...
    3 * r2 .* speed_u_first ./ speed_u.^4 - ...
    r1 .* speed_u_second ./ speed_u.^4 + ...
    3 * r1 .* speed_u_first.^2 ./ speed_u.^5;

arc_derivatives.points = points;
arc_derivatives.q_s = q_s;
arc_derivatives.q_ss = q_ss;
arc_derivatives.q_sss = q_sss;
arc_derivatives.parameter_derivatives = parameter_derivatives;
end

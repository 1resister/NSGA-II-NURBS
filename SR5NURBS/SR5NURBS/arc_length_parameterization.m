function arc = arc_length_parameterization(curve, u, num_samples, precomputed_points)
% 弧长重参数化：建立 u->S 和 S->u，并输出等弧长采样点。

if nargin >= 4 && ~isempty(precomputed_points)
    points = precomputed_points;
else
    points = evaluate_nurbs_curve(curve.control_points, curve.weights, curve.degree, curve.knots, u);
end
d = [0; vecnorm(diff(points), 2, 2)];
S_raw = cumsum(d);
[S_unique, ia] = unique(S_raw, 'stable');
u_unique = u(ia);

base_S = linspace(0, S_unique(end), num_samples)';
unique_knots = unique(curve.knots(:));
span_midpoints = 0.5 * (unique_knots(1:end-1) + unique_knots(2:end));
critical_u = unique([unique_knots; span_midpoints]);
critical_u = critical_u(critical_u >= u_unique(1) & ...
    critical_u <= u_unique(end));
critical_S = interp1(u_unique, S_unique, critical_u, 'linear');
S = unique([base_S; critical_S]);
u_sample = interp1(S_unique, u_unique, S, 'linear', 'extrap');
derivatives = nurbs_arc_derivatives(curve, u_sample);
sample_points = derivatives.points;

arc.S = S;
arc.ds = median(diff(S));
arc.u = u_sample;
arc.sample_points = sample_points;
arc.derivatives = derivatives;
arc.total_length = S(end);
end

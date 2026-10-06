function geom = compute_curvature(points, S, arc_derivatives)
% 基于弧长参数计算切向、法向和曲率。
%
% 曲率公式：kappa = |r' x r''| / |r'|^3。

if nargin >= 3 && ~isempty(arc_derivatives)
    r1 = arc_derivatives.q_s;
    r2 = arc_derivatives.q_ss;
    speed = max(vecnorm(r1, 2, 2), 1e-9);
    tangent = r1 ./ speed;
    cross_term = cross(r1, r2, 2);
    kappa = vecnorm(cross_term, 2, 2) ./ speed.^3;
    normal_raw = r2;
else
    r1 = zeros(size(points));
    r2 = zeros(size(points));
    for j = 1:3
        r1(:, j) = gradient(points(:, j), S);
        r2(:, j) = gradient(r1(:, j), S);
    end
    speed = max(vecnorm(r1, 2, 2), 1e-9);
    tangent = r1 ./ speed;
    cross_term = cross(r1, r2, 2);
    kappa = vecnorm(cross_term, 2, 2) ./ speed.^3;
    normal_raw = zeros(size(points));
    for j = 1:3
        normal_raw(:, j) = gradient(tangent(:, j), S);
    end
end
kappa(~isfinite(kappa)) = 0;
normal_norm = vecnorm(normal_raw, 2, 2);
normal = normal_raw ./ max(normal_norm, 1e-9);

geom.kappa = kappa;
geom.radius = 1 ./ max(kappa, 1e-9);
geom.tangent = tangent;
geom.normal = normal;
end

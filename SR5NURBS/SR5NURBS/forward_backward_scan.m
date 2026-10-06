function profile = forward_backward_scan(S, Vlimit, cfg)
% Continuous acceleration envelope used before full-SPU seventh-order planning.
%
% Do not solve a complete zero-acceleration/zero-jerk transition over every
% sampling interval. Doing so restarts the seventh-order transition at each
% point and makes the reachable speed depend strongly on sampling density.
% Aeq7 gives the exact acceleration-limited distance envelope for the
% seventh-order transition; jerk and full transition distance are checked
% later over each complete SPU by velocity_planning_7th.

S = S(:);
Vlimit = max(real(Vlimit(:)), 0);
n = numel(S);
Aeq7 = (16 / 35) * cfg.limits.Amax;

v_back = min(Vlimit, cfg.limits.Vmax);
v_back(end) = 0;
for i = n-1:-1:1
    ds = max(S(i + 1) - S(i), 0);
    approximate = sqrt(max(v_back(i + 1)^2 + 2 * Aeq7 * ds, 0));
    v_back(i) = min([Vlimit(i), cfg.limits.Vmax, approximate]);
end

v_forward = zeros(n, 1);
v_forward(1) = 0;
for i = 1:n-1
    ds = max(S(i + 1) - S(i), 0);
    approximate = sqrt(max(v_forward(i)^2 + 2 * Aeq7 * ds, 0));
    v_forward(i + 1) = min([v_back(i + 1), Vlimit(i + 1), approximate]);
end

v_forward(1) = 0;
v_forward(end) = 0;
profile.S = S;
profile.V = min(v_forward, Vlimit);
profile.Vlimit = Vlimit;
profile.v_backward = v_back;
profile.v_forward = v_forward;
profile.Aeq7 = Aeq7;
end

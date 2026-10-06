function profile = enforce_spu_boundary_reachability(profile, lmsc, cfg)
%ENFORCE_SPU_BOUNDARY_REACHABILITY Make LMSC boundary speeds septic-feasible.
%
% The continuous pointwise envelope determines candidate peak speeds. LMSC
% boundaries additionally require zero-acceleration/zero-jerk seventh-order
% transitions over each complete SPU. Applying the exact reachability solve
% only at these boundaries avoids the sampling-density conservatism caused
% by restarting a full transition at every path sample.

indices = unique(lmsc.indices(:), 'stable');
if isempty(indices)
    return
end
if indices(1) ~= 1
    indices = [1; indices];
end
if indices(end) ~= numel(profile.S)
    indices(end + 1, 1) = numel(profile.S);
end

boundary_S = profile.S(indices);
boundary_limit = min(profile.V(indices), profile.Vlimit(indices));
boundary_limit = min(boundary_limit, cfg.limits.Vmax);
boundary_limit = max(boundary_limit, 0);
n = numel(indices);

backward = boundary_limit;
backward(end) = 0;
for i = n-1:-1:1
    distance = max(boundary_S(i + 1) - boundary_S(i), 0);
    upper = boundary_limit(i);
    if upper > backward(i + 1)
        backward(i) = solve_reachable_speed_7th( ...
            backward(i + 1), distance, upper, cfg);
    else
        backward(i) = upper;
    end
end

forward = zeros(n, 1);
forward(1) = 0;
for i = 1:n-1
    distance = max(boundary_S(i + 1) - boundary_S(i), 0);
    upper = min(backward(i + 1), boundary_limit(i + 1));
    if upper > forward(i)
        forward(i + 1) = solve_reachable_speed_7th( ...
            forward(i), distance, upper, cfg);
    else
        forward(i + 1) = upper;
    end
end
forward([1 end]) = 0;
profile.V(indices) = min(profile.V(indices), forward);
profile.lmsc_boundary_indices = indices;
profile.lmsc_boundary_velocity = profile.V(indices);
end

function lmsc = detect_lmsc_points(S, v_limit, kappa, cfg, extra_forced_indices)
% Detect and filter local minimum speed-constraint points (LMSC).

S = S(:);
v_limit = v_limit(:);
kappa = kappa(:);
n = numel(S);
if nargin < 5
    extra_forced_indices = [];
end

endpoint_forced = false(n, 1);
endpoint_forced([1 end]) = true;
configured_forced = false(n, 1);
dynamic_forced = false(n, 1);
configured = cfg.speed.forced_limit_indices(:);
configured = configured(configured >= 1 & configured <= n);
extra_forced_indices = extra_forced_indices(:);
extra_forced_indices = extra_forced_indices(extra_forced_indices >= 1 & extra_forced_indices <= n);
configured_forced(configured) = true;
dynamic_forced(extra_forced_indices) = true;
forced = endpoint_forced | configured_forced | dynamic_forced;
protected_forced = endpoint_forced | configured_forced;

candidate = forced;
window = max(2, round(n * 0.01));
for i = 2:n-1
    is_minimum = v_limit(i) <= v_limit(i - 1) && v_limit(i) < v_limit(i + 1);
    left = max(v_limit(max(1, i-window):i-1));
    right = max(v_limit(i+1:min(n, i+window)));
    prominence = min(left - v_limit(i), right - v_limit(i));
    local_drop = max(v_limit([i-1 i+1])) - v_limit(i);

    curvature_peak = kappa(i) >= kappa(i - 1) && kappa(i) > kappa(i + 1);
    curvature_prominence = kappa(i) - min(kappa([i-1 i+1]));
    keep_curvature_peak = curvature_peak && curvature_prominence > 0 && ...
        prominence >= 0.5 * cfg.speed.lmsc_prominence_threshold;

    if is_minimum && (local_drop >= cfg.speed.lmsc_velocity_threshold || ...
            prominence >= cfg.speed.lmsc_prominence_threshold || keep_curvature_peak)
        candidate(i) = true;
    end
end

indices = find(candidate);
kept = indices(1);
for k = 2:numel(indices)
    current = indices(k);
    previous = kept(end);
    previous_speed = v_limit(previous);
    current_speed = v_limit(current);
    if previous == 1 || previous == n
        previous_speed = 0;
    end
    if current == 1 || current == n
        current_speed = 0;
    end

    % Endpoints define zero-speed boundary conditions but must not absorb
    % the first or last interior geometric/dynamic speed feature. Boundary
    % reachability is handled later by the seventh-order scan.
    if endpoint_forced(previous) || endpoint_forced(current)
        kept(end + 1, 1) = current; %#ok<AGROW>
        continue
    end

    minimum_spacing = dynamic_lmsc_spacing( ...
        previous_speed, current_speed, S, cfg);
    if S(current) - S(previous) >= minimum_spacing || ...
            (protected_forced(current) && protected_forced(previous))
        kept(end + 1, 1) = current; %#ok<AGROW>
        continue
    end

    if protected_forced(current) && ~protected_forced(previous)
        kept(end) = current;
    elseif ~protected_forced(current) && protected_forced(previous)
        continue
    elseif v_limit(current) < v_limit(previous) || ...
            (v_limit(current) == v_limit(previous) && kappa(current) > kappa(previous))
        kept(end) = current;
    end
end

if kept(1) ~= 1
    kept = [1; kept];
end
if kept(end) ~= n
    kept(end + 1, 1) = n;
end
kept = unique(kept, 'stable');

lmsc.indices = kept;
lmsc.index = kept;
lmsc.S = S(kept);
lmsc.velocity = v_limit(kept);
lmsc.curvature = kappa(kept);
lmsc.is_forced = forced(kept);
lmsc.is_configured = configured_forced(kept);
lmsc.is_dynamic = dynamic_forced(kept);
lmsc.is_protected = protected_forced(kept);
end

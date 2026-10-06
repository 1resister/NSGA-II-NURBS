function [minimum_spacing, transition_distance, sample_floor, transition_time] = ...
        dynamic_lmsc_spacing(v0, v1, S, cfg)
% Compute a physics- and resolution-based minimum spacing for two LMSCs.

% Two distinct LMSC boundary states need enough arc length for a seventh-
% order transition with zero acceleration and jerk at both ends.  A small
% resolution floor also prevents several samples of the same local feature
% from being treated as independent speed-planning units.

v0 = sanitize_speed(v0, cfg);
v1 = sanitize_speed(v1, cfg);

if abs(v1 - v0) <= cfg.speed.velocity_tolerance
    transition_distance = 0;
    transition_time = 0;
else
    [transition_distance, transition_time] = ...
        septic_transition_distance(v0, v1, cfg);
end

S = S(:);
steps = diff(S);
steps = steps(isfinite(steps) & steps > cfg.speed.distance_tolerance);
if isempty(steps)
    sample_step = 0;
else
    sample_step = median(steps);
end
sample_floor = cfg.speed.lmsc_min_sample_intervals * sample_step;

base_spacing = max( ...
    cfg.speed.lmsc_transition_distance_factor * transition_distance, ...
    sample_floor);
minimum_spacing = cfg.speed.lmsc_merge_distance_factor * base_spacing;
end

function value = sanitize_speed(value, cfg)
value = value(1);
if ~isfinite(value)
    value = cfg.limits.Vmax;
end
value = min(max(value, 0), cfg.limits.Vmax);
end

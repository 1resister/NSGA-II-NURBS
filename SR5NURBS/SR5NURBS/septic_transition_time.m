function [T, T_A, T_J] = septic_transition_time(v0, v1, cfg)
% Minimum seventh-order transition time satisfying path A/J limits.

delta_v = abs(v1 - v0);
T_A = (35 / 16) * delta_v / max(cfg.limits.Amax, eps);
T_J = sqrt((84 * sqrt(5) / 25) * delta_v / max(cfg.limits.Jmax, eps));
T = max([T_A, T_J, cfg.speed.minimum_transition_time]);
end

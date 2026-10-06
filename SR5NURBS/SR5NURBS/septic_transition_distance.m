function [distance, T] = septic_transition_distance(v0, v1, cfg)
% Distance travelled by a minimum-time seventh-order speed transition.

T = septic_transition_time(v0, v1, cfg);
distance = 0.5 * (max(v0, 0) + max(v1, 0)) * T;
end

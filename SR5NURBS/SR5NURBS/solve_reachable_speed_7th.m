function reachable = solve_reachable_speed_7th(start_speed, available_distance, upper_speed, cfg)
% Maximum upper speed reachable with a seventh-order transition.

start_speed = max(real(start_speed), 0);
upper_speed = max(real(upper_speed), 0);
available_distance = max(real(available_distance), 0);

if upper_speed <= start_speed || available_distance <= cfg.speed.distance_tolerance
    reachable = min(start_speed, upper_speed);
    return
end

if septic_transition_distance(start_speed, upper_speed, cfg) <= ...
        available_distance + cfg.speed.distance_tolerance
    reachable = upper_speed;
    return
end

low = start_speed;
high = upper_speed;
for iter = 1:cfg.speed.binary_search_max_iterations
    middle = 0.5 * (low + high);
    required = septic_transition_distance(start_speed, middle, cfg);
    if required <= available_distance
        low = middle;
    else
        high = middle;
    end
    if high - low <= cfg.speed.velocity_tolerance || ...
            abs(required - available_distance) <= cfg.speed.distance_tolerance
        break
    end
end
reachable = max(low, 0);
end

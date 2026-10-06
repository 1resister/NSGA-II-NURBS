function [keep, minimum_time, maximum_time] = ...
        filter_pareto_by_time(raw_objectives, maximum_ratio)
% Keep Pareto candidates within a multiple of the shortest raw time.
if isempty(raw_objectives)
    keep = false(0, 1);
    minimum_time = NaN;
    maximum_time = NaN;
    return
end
if size(raw_objectives, 2) < 1 || ...
        any(~isfinite(raw_objectives(:, 1))) || ...
        any(raw_objectives(:, 1) < 0)
    error('Pareto candidate times must be finite and nonnegative.');
end
if ~isscalar(maximum_ratio) || isnan(maximum_ratio) || maximum_ratio < 1
    error('The Pareto maximum time ratio must be a scalar greater than or equal to one.');
end

minimum_time = min(raw_objectives(:, 1));
if isinf(maximum_ratio)
    maximum_time = Inf;
    keep = true(size(raw_objectives, 1), 1);
else
    maximum_time = maximum_ratio * minimum_time;
    tolerance = 1e-12 * max(1, maximum_time);
    keep = raw_objectives(:, 1) <= maximum_time + tolerance;
end
end

function display_pareto_candidates(candidates)
% Display raw Pareto objectives without MATLAB's shared exponent scaling.

if isempty(candidates)
    disp('No Pareto candidates are available.');
    return
end
if size(candidates, 2) < 4
    error('Pareto candidates must contain solution_id, Time, Error and Vibration.');
end

fprintf(['Raw Pareto objectives with optimized endpoint jerk ' ...
    '(penalties excluded):\n']);
fprintf('%-12s %-16s %-16s %-16s %-16s\n', ...
    'solution_id', 'Time', 'Error', 'Vibration', 'jerk_utilization');
for i = 1:size(candidates, 1)
    fprintf('%-12d %-16.9g %-16.9g %-16.9g %-16.9g\n', ...
        round(candidates(i, 1)), candidates(i, 2), ...
        candidates(i, 3), candidates(i, 4), candidates(i, end));
end
end

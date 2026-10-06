function resume_optimization
%RESUME_OPTIMIZATION Resume NSGA-II from chromosome_iter.txt.

close all;
clc;

project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
cfg = trajectory_config(project_root);
fprintf('Current dynamic dimensions: %d control points, %d speed nodes.\n', ...
    cfg.nurbs.num_control_points, cfg.nsga.num_time_nodes);

pop = cfg.nsga.pop;
gen = cfg.nsga.gen;
M = 3;
V = cfg.nsga.num_variables;
min_range = [cfg.bounds.weights_lb cfg.bounds.time_lb ...
    cfg.bounds.jerk_utilization_lb];
max_range = [cfg.bounds.weights_ub cfg.bounds.time_ub ...
    cfg.bounds.jerk_utilization_ub];

NSGA2(pop, gen, M, V, min_range, max_range, true);

if isfile(cfg.files.chromosome_final)
    chromosome = load(cfg.files.chromosome_final);
    candidates = list_pareto_candidates(chromosome, cfg);
    disp('Pareto candidates were written to output/pareto_candidates.csv.');
    display_pareto_candidates(candidates);

    if ~isnan(cfg.selection.selected_solution_id)
        selected_id = cfg.selection.selected_solution_id;
        row = find(candidates(:, 1) == selected_id, 1);
        if isempty(row)
            error('selected_solution_id=%d is not on the Pareto front.', ...
                selected_id);
        end
        chromosome_row = candidates(row, 5);
        best_x = chromosome(chromosome_row, 1:V);
        result = evaluate_trajectory(best_x, cfg, true);
        export_results(result, chromosome, candidates, selected_id, cfg);
        visualization(result, chromosome, candidates, selected_id, cfg);
    else
        disp(['Set cfg.selection.selected_solution_id in ' ...
            'trajectory_config.m to export one solution.']);
    end
end
end

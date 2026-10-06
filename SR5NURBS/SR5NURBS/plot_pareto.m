

clear; close all; clc

project_root = fileparts(mfilename('fullpath'));
cfg = trajectory_config(project_root);
V = cfg.nsga.num_variables;

xf = load('chromosome_final.txt');
if size(xf, 2) < V + 3
    error(['Saved chromosomes do not match the current dynamic ' ...
        'configuration. Run main_function to generate a new population.']);
end
candidates = list_pareto_candidates(xf, cfg);

figure(1); clf
plot3(candidates(:, 2), candidates(:, 3), candidates(:, 4), ...
    'o', 'MarkerFaceColor', [0.2 0.45 0.8])
xlabel('Time [s]')
ylabel('Contour error [mm]')
zlabel('Vibration')
title('Pareto front: raw objectives with optimized endpoint jerk')
legend('optimized raw-objective Pareto candidates')
grid on
axis tight

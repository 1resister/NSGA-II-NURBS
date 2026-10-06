

clear; close all; clc

project_root = fileparts(mfilename('fullpath'));
cfg = trajectory_config(project_root);
V = cfg.nsga.num_variables;

x = load('chromosome_final.txt');

n = min(6, size(x, 1));
objective_columns = V + (1:3);
if size(x, 2) < objective_columns(end)
    error(['Saved chromosomes do not match the current dynamic ' ...
        'configuration. Run main_function to generate a new population.']);
end

figure(1)
clf
figure(1)
subplot(1, 3, 1)
hold on
plot(x(n, objective_columns(1)), x(n, objective_columns(2)), ...
    'o', 'MarkerFaceColor', 'b', 'MarkerSize', 10)
plot(x(:, objective_columns(1)), x(:, objective_columns(2)), 'o')
hold off
legend('selected solution')
xlabel('f_1: traveling time [sec]')
ylabel('f_2: mean acceleration [rad/s^2]')
grid on
axis tight
subplot(1, 3, 2)
hold on
plot(x(n, objective_columns(1)), x(n, objective_columns(3)), ...
    'o', 'MarkerFaceColor', 'b', 'MarkerSize', 10)
plot(x(:, objective_columns(1)), x(:, objective_columns(3)), 'o')
hold off
% legend('initial populations', '100th generation', '500th generation', 'final generation')
xlabel('f_1: traveling time [sec]')
ylabel('f_3: mean jerk [rad/s^3]')
grid on
axis tight
subplot(1, 3, 3)
hold on
plot(x(n, objective_columns(2)), x(n, objective_columns(3)), ...
    'o', 'MarkerFaceColor', 'b', 'MarkerSize', 10)
plot(x(:, objective_columns(2)), x(:, objective_columns(3)), 'o')
hold off
% legend('initial populations', '100th generation', '500th generation', 'final generation')
xlabel('f_2: mean acceleration [rad/s^2]')
ylabel('f_3: mean jerk [rad/s^3]')
grid on
axis tight

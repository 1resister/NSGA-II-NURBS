function visualization(result, ~, candidates, selected_id, cfg)
% 自动生成所选解图形：轨迹、曲率-弧长、路径运动量-时间、XYZ三方向运动量-时间。

if ~exist(cfg.output_dir, 'dir')
    mkdir(cfg.output_dir);
end
solution_dir = fullfile(cfg.output_dir, sprintf('solution_%03d', selected_id));
if ~exist(solution_dir, 'dir')
    mkdir(solution_dir);
end

% Repeated exports used to leave every old figure open. Besides consuming
% graphics resources, that could make the PNG backend fail on a later plot.
old_figures = findall(groot, 'Type', 'figure', ...
    'Tag', 'SR5NURBSVisualization');
if ~isempty(old_figures)
    close(old_figures);
end

fig = new_visualization_figure('XYZ Trajectory');
plot3(result.raw_points(:,1), result.raw_points(:,2), result.raw_points(:,3), 'ro--'); hold on
plot3(result.smooth_path(:,1), result.smooth_path(:,2), result.smooth_path(:,3), 'b-', 'LineWidth', 1.5);
grid on; axis equal; xlabel('X'); ylabel('Y'); zlabel('Z');
legend('原始离散点', 'NURBS光顺轨迹');
export_figure(fig, fullfile(solution_dir, 'trajectory.png'));

plot_series(result.arc.S, result.geom.kappa, 'S', 'kappa', fullfile(solution_dir, 'curvature_vs_arc_length.png'));
plot_speed_constraints(result, fullfile(solution_dir, 'speed_constraints_vs_arc_length.png'));
plot_series(result.motion.time, result.motion.path_velocity, 'time', 'V', fullfile(solution_dir, 'velocity_vs_time.png'));
plot_series(result.motion.time, result.motion.path_acceleration, 'time', 'A', fullfile(solution_dir, 'acceleration_vs_time.png'));
plot_series(result.motion.time, result.motion.path_jerk, 'time', 'J', fullfile(solution_dir, 'jerk_vs_time.png'));
plot_xyz_series(result.motion.time, result.xyz.velocity, 'time', 'velocity', ...
    {'Vx','Vy','Vz'}, fullfile(solution_dir, 'xyz_velocity_vs_time.png'));
plot_xyz_series(result.motion.time, result.xyz.acceleration, 'time', 'acceleration', ...
    {'Ax','Ay','Az'}, fullfile(solution_dir, 'xyz_acceleration_vs_time.png'));
plot_xyz_jerk_with_boundaries(result, ...
    fullfile(solution_dir, 'xyz_jerk_vs_time.png'));
plot_xyz_snap_with_limits(result, cfg, ...
    fullfile(solution_dir, 'xyz_snap_vs_time.png'));
plot_resonance_spectrum(result, cfg, ...
    fullfile(solution_dir, 'resonance_spectrum.png'));

fig = new_visualization_figure('Pareto Front');
h_front = plot3(candidates(:,2), candidates(:,3), candidates(:,4), 'o'); hold on
selected = candidates(candidates(:, 1) == selected_id, :);
h_selected = plot3(selected(2), selected(3), selected(4), 'rp', ...
    'MarkerSize', 14, 'MarkerFaceColor', 'r');
grid on; xlabel('Time'); ylabel('Error'); zlabel('Vibration');
title('Pareto front: raw objectives with optimized endpoint jerk');
legend([h_front,h_selected], ...
    '端点跃度利用率优化后的 Pareto 候选解', '当前选择');
export_figure(fig, fullfile(solution_dir, 'pareto_front_selected.png'));
end

function plot_speed_constraints(result, filename)
c = result.constraints;
fig = new_visualization_figure('Speed constraints');
hold on
h_command = plot(c.S, c.v_command, '--', 'LineWidth', 1.0);
h_curvature = plot(c.S, constraint_limit_for_plot(c.v_curvature, c.v_command), ...
    'LineWidth', 1.0);
h_chord = plot(c.S, constraint_limit_for_plot(c.v_chord, c.v_command), ...
    'LineWidth', 1.0);
h_axis = plot(c.S, finite_for_plot(c.v_axis), 'LineWidth', 1.0);
h_geometric_jerk = plot(c.S, ...
    constraint_limit_for_plot(c.v_geometric_jerk, c.v_command), ...
    'LineWidth', 1.0);
h_limit = plot(c.S, c.v_limit, 'k-', 'LineWidth', 1.4);
h_planning = plot(c.S, c.v_planning_limit, 'g-', 'LineWidth', 2.6);
for i = 1:numel(result.spus)
    xline(result.spus(i).start_s, ':', ...
        'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
end
h_boundary = xline(result.spus(end).end_s, ':', ...
    'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
% Draw the planned speed last so it stays visible above overlapping limits.
h_planned = plot(c.S, c.v_planned, 'b-', 'LineWidth', 2.6);
h_lmsc = plot(result.lmsc.S, result.lmsc.velocity, 'ro', ...
    'MarkerFaceColor', 'r');
grid on
xlabel('S'); ylabel('velocity');
legend([h_command,h_curvature,h_chord,h_axis,h_geometric_jerk,h_limit, ...
    h_planning,h_planned,h_lmsc,h_boundary], ...
    'v command','v curvature','v chord','v axis','v geometric jerk','v limit', ...
    'v planning limit','v planned','LMSC','SPU boundary');
export_figure(fig, filename);
end

function y = finite_for_plot(y)
y(~isfinite(y)) = NaN;
end

function y = constraint_limit_for_plot(v_constraint, v_command)
% Hide a constraint wherever it is above the commanded speed. Such values
% are inactive and can otherwise compress the useful vertical plot range.
y = finite_for_plot(v_constraint);
tolerance = 1e-9 .* max(1, abs(v_command));
active = y < v_command - tolerance;
y(~active) = NaN;
end

function plot_series(x, y, xl, yl, filename)
fig = new_visualization_figure(yl);
plot(x, y, 'LineWidth', 1.2);
grid on; xlabel(xl); ylabel(yl);
export_figure(fig, filename);
end

function plot_xyz_series(x, y, xl, yl, names, filename)
fig = new_visualization_figure(yl);
plot(x, y(:, 1), 'LineWidth', 1.2); hold on
plot(x, y(:, 2), 'LineWidth', 1.2);
plot(x, y(:, 3), 'LineWidth', 1.2);
grid on; xlabel(xl); ylabel(yl);
legend(names{:});
export_figure(fig, filename);
end

function plot_xyz_snap_with_limits(result, cfg, filename)
colors = lines(3);
fig = new_visualization_figure('XYZ snap');
hold on
for axis_id = 1:3
    plot(result.motion.time, result.xyz.snap(:, axis_id), ...
        'Color', colors(axis_id, :), 'LineWidth', 1.2);
end
limits = cfg.limits.xyz_smax(:)';
for axis_id = 1:3
    if isfinite(limits(axis_id))
        yline(limits(axis_id), '--', 'Color', colors(axis_id, :), ...
            'HandleVisibility', 'off');
        yline(-limits(axis_id), '--', 'Color', colors(axis_id, :), ...
            'HandleVisibility', 'off');
    end
end
grid on; xlabel('time'); ylabel('snap');
legend('Sx','Sy','Sz');
export_figure(fig, filename);
end

function plot_resonance_spectrum(result, cfg, filename)
resonance = result.vibration_metrics.resonance;
if isempty(resonance.frequency_hz)
    return
end

frequency = resonance.frequency_hz;
colors = lines(3);
fig = new_visualization_figure('Machine resonance spectrum');
subplot(2, 1, 1);
hold on
for axis_id = 1:3
    semilogy(frequency, max(resonance.jerk_power(:, axis_id), realmin), ...
        'Color', colors(axis_id, :), 'LineWidth', 1.0);
end
plot_mode_lines(cfg, colors);
grid on
xlabel('frequency (Hz)'); ylabel('jerk power');
legend('X','Y','Z');
title('XYZ command jerk spectrum');

subplot(2, 1, 2);
hold on
if resonance.configured
    spectrum = resonance.response_power;
    y_label = 'predicted response power';
    plot_title = 'FRF-weighted acceleration response';
else
    spectrum = resonance.acceleration_power;
    y_label = 'acceleration power';
    plot_title = 'XYZ acceleration spectrum (no resonance modes configured)';
end
for axis_id = 1:3
    semilogy(frequency, max(spectrum(:, axis_id), realmin), ...
        'Color', colors(axis_id, :), 'LineWidth', 1.0);
end
plot_mode_lines(cfg, colors);
grid on
xlabel('frequency (Hz)'); ylabel(y_label);
legend('X','Y','Z');
title(plot_title);
export_figure(fig, filename);
end

function plot_mode_lines(cfg, colors)
if ~isfield(cfg, 'resonance') || ...
        ~isfield(cfg.resonance, 'mode_frequencies_hz')
    return
end
modes = cfg.resonance.mode_frequencies_hz;
for axis_id = 1:3
    if iscell(modes)
        frequencies = modes{axis_id};
    else
        frequencies = modes;
    end
    for mode_id = 1:numel(frequencies)
        xline(frequencies(mode_id), ':', 'Color', colors(axis_id, :), ...
            'HandleVisibility', 'off');
    end
end
end

function fig = new_visualization_figure(name)
fig = figure('Name', name, 'Tag', 'SR5NURBSVisualization');
end

function export_figure(fig, filename)
% Use exportgraphics instead of saveas/print and preserve a new image even
% when an existing target PNG is temporarily locked by another program.
drawnow;
try
    exportgraphics(fig, filename, 'Resolution', 150);
catch primary_error
    [folder, name, extension] = fileparts(filename);
    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
    fallback = fullfile(folder, sprintf('%s_%s%s', ...
        name, timestamp, extension));
    try
        exportgraphics(fig, fallback, 'Resolution', 150);
        warning('visualization:AlternateImageFile', ...
            ['Could not overwrite %s (%s). The new plot was saved as ' ...
            '%s instead. Close any program displaying the old PNG before ' ...
            'exporting again.'], filename, primary_error.message, fallback);
    catch fallback_error
        warning('visualization:ImageExportFailed', ...
            'Could not export %s: %s | fallback failed: %s', ...
            filename, primary_error.message, fallback_error.message);
    end
end
end

function plot_xyz_jerk_with_boundaries(result, filename)
plot_time = result.motion.time;
plot_jerk = result.xyz.jerk;
segments = result.motion.segments;
globally_smoothed = isfield(result.motion, 'time_smoothing') && ...
    isfield(result.motion.time_smoothing, 'accepted') && ...
    result.motion.time_smoothing.accepted;

if numel(segments) > 1 && ~globally_smoothed
    join_time = [segments(1:end-1).t1]';
    join_S = [segments(1:end-1).s1]';
    join_velocity = [segments(1:end-1).v1]';
    join_path_jerk = [segments(1:end-1).j1]';
    [unique_S, unique_index] = unique(result.motion.S, 'stable');
    if numel(unique_S) >= 2
        q_sss_join = interp1(unique_S, result.xyz.q_sss(unique_index, :), ...
            join_S, 'pchip', 'extrap');
        q_s_join = interp1(unique_S, result.xyz.q_s(unique_index, :), ...
            join_S, 'pchip', 'extrap');
        join_jerk = q_sss_join .* (join_velocity.^3) + ...
            q_s_join .* join_path_jerk;

        time_tolerance = max(10 * eps(max(plot_time(end), 1)), ...
            result.motion.time(end) * 1e-12);
        for i = 1:numel(join_time)
            coincident = abs(plot_time - join_time(i)) <= time_tolerance;
            plot_time(coincident) = [];
            plot_jerk(coincident, :) = [];
        end
        plot_time = [plot_time; join_time];
        plot_jerk = [plot_jerk; join_jerk];
        [plot_time, order] = sort(plot_time);
        plot_jerk = plot_jerk(order, :);
    end
end

plot_xyz_series(plot_time, plot_jerk, 'time', 'jerk', ...
    {'Jx','Jy','Jz'}, filename);
end

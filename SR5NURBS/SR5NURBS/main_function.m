%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 三轴机床 NURBS 轨迹规划与 NSGA-II 多目标优化主入口
%
% 流程：
% CSV 离散点 -> NURBS 全局光顺 -> 弧长重参数化 -> 曲率/约束 ->
% 七次多项式速度规划 -> XYZ 速度/加速度/跃度 -> NSGA-II Pareto 优化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc

project_root = fileparts(mfilename('fullpath'));
addpath(project_root);

cfg = trajectory_config(project_root);
fprintf(['Dynamic dimensions: %d NURBS control points, %d speed nodes, ' ...
    '%d variables. Baseline contour error: %.6g mm (target %.6g mm).\n'], ...
    cfg.nurbs.num_control_points, cfg.nsga.num_time_nodes, ...
    cfg.nsga.num_variables, cfg.dynamic_dimensions.baseline_error, ...
    cfg.constraints.max_contour_error);
if ~cfg.dynamic_dimensions.baseline_target_met
    warning('main_function:DynamicFitTargetNotMet', ...
        ['The maximum dynamic control-point count was reached before the ' ...
        'unit-weight baseline met the contour-error target. NSGA-II may ' ...
        'still improve it through the NURBS weights.']);
end

mode_count = 0;
highest_mode_hz = 0;
configured_modes = cfg.resonance.mode_frequencies_hz;
for axis_id = 1:3
    if iscell(configured_modes)
        axis_modes = configured_modes{axis_id};
    else
        axis_modes = configured_modes;
    end
    mode_count = mode_count + numel(axis_modes);
    if ~isempty(axis_modes)
        highest_mode_hz = max(highest_mode_hz, max(axis_modes));
    end
end
fprintf('Resonance optimization: enabled=%d, configured modes=%d, Nyquist=%.3f Hz.\n', ...
    cfg.resonance.enabled, mode_count, 1/(2*cfg.interpolation.Ts));
fprintf(['XYZ snap hard limit: enabled=%d, limits=[%.6g %.6g %.6g] ' ...
    'mm/s^4.\n'], xyz_snap_limit_enabled(cfg), cfg.limits.xyz_smax);
fprintf(['Shared endpoint-jerk decision variable: enabled=%d, ' ...
    'bounds=[%.3g, %.3g], relative-vibration cap enabled=%d.\n'], ...
    cfg.lmsc_jerk_boundary.enabled && ...
        cfg.lmsc_jerk_boundary.evaluate_during_objective, ...
    cfg.bounds.jerk_utilization_lb, cfg.bounds.jerk_utilization_ub, ...
    cfg.lmsc_jerk_boundary.enforce_vibration_increase_limit);
if isfinite(cfg.evaluation.max_time_samples)
    optimization_time_limit = ...
        (cfg.evaluation.max_time_samples - 1) * cfg.interpolation.Ts;
    fprintf(['NSGA-II candidate memory guard: at most %d fixed-Ts samples ' ...
        '(%.6g s); selected-solution export remains uncapped.\n'], ...
        cfg.evaluation.max_time_samples, optimization_time_limit);
end
if cfg.evaluation.enable_contour_prefilter
    fprintf(['Contour prefilter: candidates above %.3g%% of the configured ' ...
        'error limit are penalized before time-domain evaluation.\n'], ...
        100 * (1 + cfg.evaluation.contour_prefilter_relative_margin));
end
if cfg.resonance.enabled && mode_count == 0
    warning('main_function:ResonanceModesMissing', ...
        ['Resonance optimization is enabled but no machine modes are ' ...
        'configured. The resonance objective and hard limit are inactive.']);
end
if highest_mode_hz >= 1/(2*cfg.interpolation.Ts)
    warning('main_function:ResonanceModeAboveNyquist', ...
        ['At least one configured machine mode is at or above the Nyquist ' ...
        'frequency and will be excluded from resonance evaluation.']);
end

%% NSGA-II 参数。算法框架保持原项目接口，仅替换变量含义和目标函数。
pop = cfg.nsga.pop;
gen = cfg.nsga.gen;
M = 3;
V = cfg.nsga.num_variables;

min_range = [cfg.bounds.weights_lb cfg.bounds.time_lb ...
    cfg.bounds.jerk_utilization_lb];
max_range = [cfg.bounds.weights_ub cfg.bounds.time_ub ...
    cfg.bounds.jerk_utilization_ub];

   NSGA2(pop, gen, M, V, min_range, max_range);

%% 输出最终种群 Pareto 前沿解，按配置选择解后导出轨迹、CSV 和图形。
if isfile(cfg.files.chromosome_final)
    chromosome = load(cfg.files.chromosome_final);
    candidates = list_pareto_candidates(chromosome, cfg);
    disp('Pareto 前沿候选解已写入 output/pareto_candidates.csv。');
    display_pareto_candidates(candidates);

    if ~isnan(cfg.selection.selected_solution_id)
        selected_id = cfg.selection.selected_solution_id;
        row = find(candidates(:, 1) == selected_id, 1);
        if isempty(row)
            error('selected_solution_id=%d 不在 Pareto 前沿候选解中。', selected_id);
        end
        chromosome_row = candidates(row, 5);
        best_x = chromosome(chromosome_row, 1 : V);
        result = evaluate_trajectory(best_x, cfg, true);
        if isfield(result, 'lmsc_jerk_boundary_search') && ...
                result.lmsc_jerk_boundary_search.applied
            jerk_report = result.lmsc_jerk_boundary_search;
            fprintf(['LMSC/peak shared jerk decision applied: ' ...
                'utilization %.3g, time improvement %.3f%%.\n'], ...
                jerk_report.selected_utilization, ...
                100 * jerk_report.time_gain_ratio);
        elseif isfield(result, 'lmsc_jerk_boundary_search')
            jerk_report = result.lmsc_jerk_boundary_search;
            fprintf('LMSC/peak shared jerk kept at zero: %s.\n', ...
                char(jerk_report.reason));
        else
            fprintf('Shared endpoint-jerk objective search is disabled.\n');
        end
        export_results(result, chromosome, candidates, selected_id, cfg);
        visualization(result, chromosome, candidates, selected_id, cfg);
    else
        disp('请在 trajectory_config.m 中设置 cfg.selection.selected_solution_id 后重新运行，以导出指定解。');
    end
end

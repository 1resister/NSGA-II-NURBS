function export_selected_solution(selected_solution_id, force_candidate_refresh)
% 根据已有 chromosome_final.txt 导出指定 Pareto 前沿解。
%
% 用法：
%   export_selected_solution(3)
%   export_selected_solution(3, true)  % 强制刷新全部候选目标值
%
% 如果不传入编号，则读取 trajectory_config.m 中的
% cfg.selection.selected_solution_id。

project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
cfg = trajectory_config(project_root);

if nargin >= 1
    cfg.selection.selected_solution_id = selected_solution_id;
end
if nargin < 2
    force_candidate_refresh = false;
end

if ~isfile(cfg.files.chromosome_final)
    error('未找到 chromosome_final.txt，请先运行 main_function 完成 NSGA-II 优化。');
end
if isnan(cfg.selection.selected_solution_id)
    error('请提供 selected_solution_id，例如 export_selected_solution(3)。');
end

chromosome = load(cfg.files.chromosome_final);
[candidates, cache_used] = list_pareto_candidates( ...
    chromosome, cfg, force_candidate_refresh);
selected_id = cfg.selection.selected_solution_id;
row = find(candidates(:, 1) == selected_id, 1);
if isempty(row)
    error('selected_solution_id=%d 不在 Pareto 前沿候选解中。', selected_id);
end

V = cfg.nsga.num_variables;
chromosome_row = candidates(row, 5);
x = chromosome(chromosome_row, 1:V);
result = evaluate_trajectory(x, cfg, true);
if isfield(result, 'lmsc_jerk_boundary_search')
    print_lmsc_jerk_search(result.lmsc_jerk_boundary_search);
else
    fprintf(['共享端点跃度利用率决策未进入目标评价；' ...
        '当前导出使用零跃度基线。\n']);
end
export_results(result, chromosome, candidates, selected_id, cfg);
visualization(result, chromosome, candidates, selected_id, cfg);

fprintf('已导出 Pareto 前沿解 %d 到 %s\n', selected_id, ...
    fullfile(cfg.output_dir, sprintf('solution_%03d', selected_id)));
if cache_used
    fprintf(['候选目标函数来自有效缓存，本次只按同一端点跃度' ...
        '目标评价规则重新评价了所选解。\n']);
end
end

function print_lmsc_jerk_search(report)
if report.applied
    fprintf(['LMSC/非匀速峰值共享跃度已采用：利用率 %.3g，' ...
        '时间 %.9g -> %.9g s ' ...
        '(时间改善 %.3f%%)，global_time_scale %.6g -> %.6g。\n'], ...
        report.selected_utilization, report.baseline_time, ...
        report.refined_time, 100 * report.time_gain_ratio, ...
        report.baseline_global_time_scale, ...
        report.refined_global_time_scale);
else
    fprintf(['当前利用率没有产生非零共享跃度，使用零跃度轨迹：' ...
        '%s。\n'], ...
        char(report.reason));
end
end

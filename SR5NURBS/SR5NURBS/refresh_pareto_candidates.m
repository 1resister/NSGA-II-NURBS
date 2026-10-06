function candidates = refresh_pareto_candidates()
% Force regeneration of the filtered raw-objective Pareto candidate table.
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
cfg = trajectory_config(project_root);
if ~isfile(cfg.files.chromosome_final)
    error('未找到 chromosome_final.txt，请先运行 main_function 完成 NSGA-II 优化。');
end

chromosome = load(cfg.files.chromosome_final);
candidates = list_pareto_candidates(chromosome, cfg, true);
fprintf('已刷新过滤后的 Pareto 候选表：%s\n', ...
    fullfile(cfg.output_dir, 'pareto_candidates.csv'));
display_pareto_candidates(candidates);
end

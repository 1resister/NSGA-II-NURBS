function [candidates, cache_used] = list_pareto_candidates( ...
        chromosome, cfg, force_refresh)
% 从最终种群提取 rank=1 的 Pareto 前沿候选解，并写出候选表。
%
% candidates 列含义：
% solution_id, 原始 Time, 原始 Error, 原始 Vibration, chromosome_row,
% weights..., time_nodes..., endpoint_jerk_utilization。NSGA-II 排名仍使用
% chromosome 中的罚后目标，
% 候选表则按染色体保存的共享端点跃度利用率重新评价轨迹，只展示不含
% 罚值的实际物理指标。

if nargin < 3
    force_refresh = false;
end
cache_used = false;

if ~exist(cfg.output_dir, 'dir')
    mkdir(cfg.output_dir);
end

V = cfg.nsga.num_variables;
M = 3;
rank_col = V + M + 1;
column_count = size(chromosome, 2);
valid_column_counts = [V + M, V + M + 2];
if ~ismember(column_count, valid_column_counts)
    error('list_pareto_candidates:DynamicDimensionMismatch', ...
        ['Saved chromosome has %d columns, but the current configuration ' ...
        'requires %d objective-only columns or %d ranked columns ' ...
        '(%d variables and %d objectives). Run ' ...
        'main_function to generate a new population; do not reuse the ' ...
        'old checkpoint after dynamic dimensions change.'], ...
        column_count, V + M, V + M + 2, V, M);
end

if column_count == V + M + 2
    front_rows = find(chromosome(:, rank_col) == 1);
else
    front_rows = (1:size(chromosome, 1))';
end

front = chromosome(front_rows, :);
[~, order] = sortrows(front(:, V + 1 : V + M), [1 2 3]);
front = front(order, :);
front_rows = front_rows(order);

cache_file = fullfile(cfg.output_dir, 'pareto_candidates_cache.mat');
signature = candidate_cache_signature(chromosome, cfg);
if ~force_refresh && isfile(cache_file)
    try
        saved = load(cache_file, 'candidate_cache');
        if isfield(saved, 'candidate_cache') && ...
                isfield(saved.candidate_cache, 'signature') && ...
                isfield(saved.candidate_cache, 'candidates') && ...
                isequaln(saved.candidate_cache.signature, signature)
            candidates = saved.candidate_cache.candidates;
            cache_used = true;
            write_candidate_table(candidates, cfg);
            fprintf(['Using cached raw objective values for %d Pareto ' ...
                'candidate(s).\n'], size(candidates, 1));
            return
        end
    catch cache_error
        warning('list_pareto_candidates:InvalidCache', ...
            'Ignoring invalid Pareto candidate cache: %s', ...
            cache_error.message);
    end
end

raw_objectives = zeros(size(front, 1), M);
evaluated = false(size(front, 1), 1);
shortest_raw_time = Inf;
screened_count = 0;
if ~isempty(front)
    fprintf('Re-evaluating %d Pareto candidate(s) for raw objective values.\n', ...
        size(front, 1));
end
for i = 1:size(front, 1)
    maximum_raw_time = Inf;
    if isfinite(cfg.pareto.max_time_ratio) && isfinite(shortest_raw_time)
        maximum_raw_time = cfg.pareto.max_time_ratio * shortest_raw_time;
    end
    try
        result = evaluate_trajectory( ...
            front(i, 1:V), cfg, false, maximum_raw_time);
    catch evaluation_error
        if is_maximum_time_exception(evaluation_error)
            screened_count = screened_count + 1;
            continue
        end
        rethrow(evaluation_error)
    end
    raw_objectives(i, :) = [result.raw_objectives.time, ...
        result.raw_objectives.error, result.raw_objectives.vibration];
    evaluated(i) = true;
    shortest_raw_time = min(shortest_raw_time, raw_objectives(i, 1));
    clear result
    if i == 1 || mod(i, 10) == 0 || i == size(front, 1)
        fprintf('  Raw candidate %d/%d evaluated; shortest time %.6g s.\n', ...
            i, size(front, 1), shortest_raw_time);
    end
end

if ~any(evaluated)
    error('No Pareto candidate could be evaluated for raw objective values.');
end

unfiltered_count = size(front, 1);
raw_objectives = raw_objectives(evaluated, :);
front_rows = front_rows(evaluated);
front = front(evaluated, :);
[keep, minimum_time, maximum_time] = filter_pareto_by_time( ...
    raw_objectives, cfg.pareto.max_time_ratio);
raw_objectives = raw_objectives(keep, :);
front_rows = front_rows(keep);
front = front(keep, :);
solution_id = (1:size(front, 1))';
candidates = [solution_id raw_objectives front_rows front(:, 1:V)];
if unfiltered_count > size(front, 1)
    fprintf(['Pareto time filter kept %d/%d candidate(s): shortest raw ' ...
        'time %.6g s, maximum allowed time %.6g s (ratio %.6g).\n'], ...
        size(front, 1), unfiltered_count, minimum_time, maximum_time, ...
        cfg.pareto.max_time_ratio);
end
if screened_count > 0
    fprintf(['Skipped %d candidate(s) before allocating their full XYZ ' ...
        'series because their minimum feasible time exceeded the active ' ...
        'Pareto limit.\n'], screened_count);
end

write_candidate_table(candidates, cfg);
candidate_cache.signature = signature;
candidate_cache.candidates = candidates;
save(cache_file, 'candidate_cache', '-v7');
end

function exceeded = is_maximum_time_exception(exception)
% The endpoint-jerk evaluator used to wrap the underlying time-screening
% identifier. Accept both forms so serialized populations remain
% exportable across that implementation change.
exceeded = strcmp(exception.identifier, ...
    'plan_speed_profile:MaximumTimeExceeded') || ...
    (strcmp(exception.identifier, ...
        'evaluate_trajectory:InfeasibleJerkUtilization') && ...
    contains(string(exception.message), ...
        'plan_speed_profile:MaximumTimeExceeded'));
end

function write_candidate_table(candidates, cfg)
variable_count = 6 + cfg.nsga.num_weights + cfg.nsga.num_time_nodes;
var_names = strings(1, variable_count);
var_names(1:5) = ...
    ["solution_id","Time","Error","Vibration","chromosome_row"];
for i = 1:cfg.nsga.num_weights
    var_names(5 + i) = "weight_" + string(i);
end
for i = 1:cfg.nsga.num_time_nodes
    var_names(5 + cfg.nsga.num_weights + i) = "time_node_" + string(i);
end
var_names(end) = "endpoint_jerk_utilization";

T = array2table(candidates, 'VariableNames', cellstr(var_names));
writetable(T, fullfile(cfg.output_dir, 'pareto_candidates.csv'));
end

function signature = candidate_cache_signature(chromosome, cfg)
% Any changed population, input, configuration or MATLAB source invalidates
% the cache. Selection and output paths do not affect objective values.
signature.schema_version = 4;
signature.chromosome = chromosome;

evaluation_cfg = cfg;
if isfield(evaluation_cfg, 'selection')
    evaluation_cfg = rmfield(evaluation_cfg, 'selection');
end
if isfield(evaluation_cfg, 'output_dir')
    evaluation_cfg = rmfield(evaluation_cfg, 'output_dir');
end
signature.config = evaluation_cfg;

if isfield(cfg, 'files') && isfield(cfg.files, 'input_csv') && ...
        isfile(cfg.files.input_csv)
    signature.input_csv = fileread(cfg.files.input_csv);
else
    signature.input_csv = '';
end

source_files = dir(fullfile(cfg.project_root, '*.m'));
[~, order] = sort(lower(string({source_files.name})));
source_files = source_files(order);
signature.source_names = string({source_files.name});
signature.source_bytes = [source_files.bytes];
signature.source_datenum = [source_files.datenum];
end

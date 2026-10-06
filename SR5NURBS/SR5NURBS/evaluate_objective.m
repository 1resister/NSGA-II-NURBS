function obj = evaluate_objective(x, M, V)
% 评价 NSGA-II 染色体的三目标函数。
%
% 输入：
%   x(1:num_weights) NURBS weight variables.
%   The next num_time_nodes variables are arc-length speed-scale nodes.
%   The final variable is the shared endpoint-jerk utilization in [0, 1].
%
% 输出：
%   obj = [总加工时间, 最大轮廓误差, 振动指标]
%
% 本函数不再调用 SR5、DH 或逆运动学，所有目标均直接基于 XYZ 轨迹计算。

if V ~= numel(x)
    x = x(1:V);
end

cfg = trajectory_config(fileparts(mfilename('fullpath')));
try
    maximum_raw_time = optimization_time_limit(cfg);
    result = evaluate_trajectory(x, cfg, false, maximum_raw_time, true);
    obj = [result.objectives.time, result.objectives.error, ...
        result.objectives.vibration];
    if any(~isfinite(obj))
        error('evaluate_objective:NonfiniteObjective', ...
            'The trajectory evaluation returned a nonfinite objective.');
    end
catch ME
    obj = cfg.constraints.invalid_candidate_penalty * ones(1, 3);
    warn_once(ME);
end

if length(obj) ~= M
    error('目标函数数量与 NSGA-II 配置不一致。');
end

function maximum_raw_time = optimization_time_limit(cfg)
if ~isfield(cfg, 'evaluation') || ...
        ~isfield(cfg.evaluation, 'max_time_samples')
    maximum_raw_time = Inf;
    return
end
maximum_samples = cfg.evaluation.max_time_samples;
if ~isscalar(maximum_samples) || isnan(maximum_samples) || ...
        maximum_samples < 2 || ...
        (~isinf(maximum_samples) && maximum_samples ~= floor(maximum_samples))
    error('cfg.evaluation.max_time_samples must be an integer >= 2 or Inf.');
end
if isinf(maximum_samples)
    maximum_raw_time = Inf;
else
    maximum_raw_time = (maximum_samples - 1) * cfg.interpolation.Ts;
end
end

function warn_once(ME)
% Keep a bad chromosome from stopping NSGA-II without flooding the console.
persistent warned_keys
if isempty(warned_keys)
    warned_keys = {};
end
key = ME.identifier;
if isempty(key)
    key = ME.message;
end
if ~any(strcmp(warned_keys, key))
    if strcmp(ME.identifier, 'plan_speed_profile:MaximumTimeExceeded')
        warning('evaluate_objective:CandidateTimeLimit', ...
            ['Candidate was penalized before allocating an oversized ' ...
            'time-domain trajectory: %s'], ME.message);
        warned_keys{end + 1} = key;
        return
    end
    warning('evaluate_objective:InvalidCandidate', ...
        'Candidate was penalized after trajectory evaluation failed: %s', ME.message);
    warned_keys{end + 1} = key;
end
end
end

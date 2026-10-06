function spus = build_speed_planning_units(S, v_limit, lmsc)
% Build speed-planning units (SPU) between adjacent LMSC points.

S = S(:);
v_limit = v_limit(:);
indices = lmsc.indices(:);
n_spu = max(numel(indices) - 1, 0);
template = struct('id', 0, 'start_index', 0, 'end_index', 0, ...
    'start_s', 0, 'end_s', 0, 'length', 0, ...
    'start_speed_limit', 0, 'end_speed_limit', 0, ...
    'local_S', [], 'local_speed_limit', [], ...
    'candidate_peak_speed', 0, 'planned_peak_speed', 0, ...
    'acceleration_time', 0, 'constant_time', 0, 'deceleration_time', 0, ...
    'start_jerk', 0, 'peak_jerk', 0, 'end_jerk', 0);
spus = repmat(template, n_spu, 1);

for i = 1:n_spu
    i0 = indices(i);
    i1 = indices(i + 1);
    local_ids = i0:i1;
    local_limit = v_limit(local_ids);
    spus(i).id = i;
    spus(i).start_index = i0;
    spus(i).end_index = i1;
    spus(i).start_s = S(i0);
    spus(i).end_s = S(i1);
    spus(i).length = max(S(i1) - S(i0), 0);
    spus(i).start_speed_limit = v_limit(i0);
    spus(i).end_speed_limit = v_limit(i1);
    spus(i).local_S = S(local_ids);
    spus(i).local_speed_limit = local_limit;
    spus(i).candidate_peak_speed = max(local_limit);
    spus(i).planned_peak_speed = spus(i).candidate_peak_speed;
end
end

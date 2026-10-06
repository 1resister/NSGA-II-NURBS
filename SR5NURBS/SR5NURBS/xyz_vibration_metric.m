function metrics = xyz_vibration_metric( ...
        xyz_jerk, cfg, time, xyz_acceleration, keep_spectrum)
% Compute the vibration objective from the decomposed X/Y/Z jerk signals.
if nargin < 3
    time = [];
end
if nargin < 4
    xyz_acceleration = [];
end
if nargin < 5
    keep_spectrum = false;
end
if size(xyz_jerk, 2) ~= 3 || isempty(xyz_jerk)
    error('xyz_jerk must be a nonempty N-by-3 matrix.');
end
if any(~isfinite(xyz_jerk), 'all')
    error('xyz_jerk contains a nonfinite value.');
end

weights = [1 1 1];
if isfield(cfg, 'objectives') && ...
        isfield(cfg.objectives, 'vibration_axis_weights')
    weights = cfg.objectives.vibration_axis_weights;
end
weights = weights(:)';
if numel(weights) ~= 3 || any(~isfinite(weights)) || any(weights < 0) || ...
        ~any(weights > 0)
    error('cfg.objectives.vibration_axis_weights must contain three finite, nonnegative values with at least one positive value.');
end

axis_rms = sqrt(mean(xyz_jerk.^2, 1));
weighted_axis_rms = weights .* axis_rms;
combined_rms = norm(weighted_axis_rms, 2);
axis_peak = max(abs(xyz_jerk), [], 1);
weighted_point_magnitude = vecnorm(xyz_jerk .* weights, 2, 2);
combined_peak = max(weighted_point_magnitude);
peak_weight = 0;
if isfield(cfg, 'objectives') && ...
        isfield(cfg.objectives, 'vibration_peak_weight')
    peak_weight = cfg.objectives.vibration_peak_weight;
end
if ~isscalar(peak_weight) || ~isfinite(peak_weight) || peak_weight < 0
    error('cfg.objectives.vibration_peak_weight must be a finite nonnegative scalar.');
end
if keep_spectrum
    snap = zeros(size(xyz_jerk));
    crackle = zeros(size(xyz_jerk));
else
    snap = zeros(0, 3);
    crackle = zeros(0, 3);
end
axis_snap_rms = zeros(1, 3);
axis_snap_peak = zeros(1, 3);
axis_crackle_rms = zeros(1, 3);
axis_crackle_peak = zeros(1, 3);
equivalent_snap = 0;
equivalent_crackle = 0;
if ~isempty(time)
    time = time(:);
    if numel(time) ~= size(xyz_jerk, 1) || ...
            any(~isfinite(time)) || any(diff(time) <= 0)
        error('time must be finite, strictly increasing and match xyz_jerk.');
    end
    if numel(time) >= 2
        sample_time = median(diff(time));
        for axis_id = 1:3
            snap_axis = gradient(xyz_jerk(:, axis_id), sample_time);
            crackle_axis = gradient(snap_axis, sample_time);
            axis_snap_rms(axis_id) = sqrt(mean(snap_axis.^2));
            axis_snap_peak(axis_id) = max(abs(snap_axis));
            axis_crackle_rms(axis_id) = sqrt(mean(crackle_axis.^2));
            axis_crackle_peak(axis_id) = max(abs(crackle_axis));
            if keep_spectrum
                snap(:, axis_id) = snap_axis;
                crackle(:, axis_id) = crackle_axis;
            end
        end
        equivalent_snap = norm(weights .* (sample_time * axis_snap_rms), 2);
        equivalent_crackle = norm( ...
            weights .* (sample_time^2 * axis_crackle_rms), 2);
    end
end

axis_snap_limit = inf(1, 3);
if isfield(cfg, 'limits') && isfield(cfg.limits, 'xyz_smax')
    axis_snap_limit = cfg.limits.xyz_smax(:)';
end
if numel(axis_snap_limit) ~= 3 || any(isnan(axis_snap_limit)) || ...
        any(axis_snap_limit <= 0)
    error('cfg.limits.xyz_smax must contain three positive values.');
end
axis_snap_ratio = axis_snap_peak ./ axis_snap_limit;
snap_limit_enabled = xyz_snap_limit_enabled(cfg);

snap_weight = 0;
if isfield(cfg, 'objectives') && ...
        isfield(cfg.objectives, 'vibration_snap_weight')
    snap_weight = cfg.objectives.vibration_snap_weight;
end
if ~isscalar(snap_weight) || ~isfinite(snap_weight) || snap_weight < 0
    error('cfg.objectives.vibration_snap_weight must be a finite nonnegative scalar.');
end

crackle_weight = 0;
if isfield(cfg, 'objectives') && ...
        isfield(cfg.objectives, 'vibration_crackle_weight')
    crackle_weight = cfg.objectives.vibration_crackle_weight;
end
if ~isscalar(crackle_weight) || ~isfinite(crackle_weight) || ...
        crackle_weight < 0
    error(['cfg.objectives.vibration_crackle_weight must be a finite ' ...
        'nonnegative scalar.']);
end

resonance = machine_resonance_metric( ...
    xyz_acceleration, xyz_jerk, time, cfg, keep_spectrum);
objective_value = combined_rms + peak_weight * combined_peak + ...
    snap_weight * equivalent_snap + ...
    crackle_weight * equivalent_crackle + resonance.objective_value;

axis_names = ["X" "Y" "Z"];
if combined_rms > 0
    contribution_ratio = weighted_axis_rms.^2 / combined_rms^2;
    [~, dominant_axis_id] = max(weighted_axis_rms);
    dominant_axis = axis_names(dominant_axis_id);
else
    contribution_ratio = zeros(1, 3);
    dominant_axis = "none";
end

metrics.axis_rms = axis_rms;
metrics.axis_weights = weights;
metrics.weighted_axis_rms = weighted_axis_rms;
metrics.combined_rms = combined_rms;
metrics.combined_peak = combined_peak;
metrics.peak_weight = peak_weight;
metrics.snap = snap;
metrics.axis_snap_rms = axis_snap_rms;
metrics.axis_snap_peak = axis_snap_peak;
metrics.axis_snap_limit = axis_snap_limit;
metrics.axis_snap_ratio = axis_snap_ratio;
metrics.max_snap_ratio = max(axis_snap_ratio);
metrics.snap_limit_enabled = snap_limit_enabled;
metrics.snap_feasible = ~snap_limit_enabled || metrics.max_snap_ratio <= ...
    1 + cfg.speed.velocity_tolerance;
metrics.combined_equivalent_snap = equivalent_snap;
metrics.snap_weight = snap_weight;
metrics.crackle = crackle;
metrics.axis_crackle_rms = axis_crackle_rms;
metrics.axis_crackle_peak = axis_crackle_peak;
metrics.combined_equivalent_crackle = equivalent_crackle;
metrics.crackle_weight = crackle_weight;
metrics.resonance = resonance;
metrics.objective_value = objective_value;
metrics.axis_peak = axis_peak;
metrics.contribution_ratio = contribution_ratio;
metrics.dominant_axis = dominant_axis;
end

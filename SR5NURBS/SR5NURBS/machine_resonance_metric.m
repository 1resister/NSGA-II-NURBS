function metrics = machine_resonance_metric( ...
        acceleration, jerk, time, cfg, keep_spectrum)
%MACHINE_RESONANCE_METRIC Frequency-weighted feed-axis vibration metric.
%
% The configured second-order modal transfer functions weight XYZ command
% spectra around measured machine resonances. No resonance frequency is
% invented: empty per-axis mode lists leave the resonance objective and hard
% constraint inactive while selected-solution exports can still retain the
% unweighted acceleration and jerk spectra.

if nargin < 5
    keep_spectrum = false;
end

enabled = resonance_option(cfg, 'enabled', false);
diagnostics_enabled = resonance_option(cfg, 'enable_diagnostics', true);
if ~isscalar(enabled) || ~isscalar(diagnostics_enabled)
    error('cfg.resonance enable flags must be scalar values.');
end
axis_weights = resonance_vector(cfg, 'axis_weights', [1 1 1]);
response_limits = resonance_vector( ...
    cfg, 'response_rms_limits', [Inf Inf Inf]);
if any(~isfinite(axis_weights)) || any(axis_weights < 0) || ...
        ~any(axis_weights > 0)
    error('cfg.resonance.axis_weights must be nonnegative with at least one positive value.');
end
if any(response_limits <= 0) || any(isnan(response_limits))
    error('cfg.resonance.response_rms_limits must contain positive values or Inf.');
end

metrics = empty_metrics(enabled, axis_weights, response_limits);
if isempty(time) || isempty(acceleration) || isempty(jerk)
    return
end
time = time(:);
if size(acceleration, 2) ~= 3 || size(jerk, 2) ~= 3 || ...
        size(acceleration, 1) ~= numel(time) || ...
        size(jerk, 1) ~= numel(time)
    error('Acceleration, jerk and time must have compatible N-by-3/N-by-1 sizes.');
end
if any(~isfinite(acceleration), 'all') || any(~isfinite(jerk), 'all') || ...
        any(~isfinite(time)) || any(diff(time) <= 0)
    error('Resonance inputs must be finite and time must be strictly increasing.');
end
if numel(time) < 4
    return
end

sample_time = median(diff(time));
sample_rate = 1 / sample_time;
nyquist = sample_rate / 2;
metrics.sample_time = sample_time;
metrics.sample_rate_hz = sample_rate;
metrics.nyquist_hz = nyquist;

if ~enabled && ~(keep_spectrum && diagnostics_enabled)
    return
end

[frequency, acceleration_power] = one_sided_power(acceleration, sample_time);
[~, jerk_power] = one_sided_power(jerk, sample_time);
n_frequency = numel(frequency);
spectral_weight = ones(n_frequency, 3);
response_power = zeros(n_frequency, 3);
axis_resonant_jerk_rms = zeros(1, 3);
axis_response_rms = zeros(1, 3);
configured_axes = false(1, 3);
unresolved_mode_count = zeros(1, 3);
spectral_peak_gain = resonance_option(cfg, 'spectral_peak_gain', 10);
objective_weight = resonance_option(cfg, 'objective_weight', 1);
if ~isscalar(spectral_peak_gain) || ~isfinite(spectral_peak_gain) || ...
        spectral_peak_gain < 0
    error('cfg.resonance.spectral_peak_gain must be finite and nonnegative.');
end
if ~isscalar(objective_weight) || ~isfinite(objective_weight) || ...
        objective_weight < 0
    error('cfg.resonance.objective_weight must be finite and nonnegative.');
end

for axis_id = 1:3
    frequencies = axis_mode_values(cfg, 'mode_frequencies_hz', axis_id);
    damping = expand_mode_values(axis_mode_values( ...
        cfg, 'damping_ratios', axis_id), numel(frequencies), 0.02);
    gains = expand_mode_values(axis_mode_values( ...
        cfg, 'modal_gains', axis_id), numel(frequencies), 1);
    if any(~isfinite(frequencies)) || any(frequencies <= 0)
        error('Configured resonance frequencies must be finite and positive.');
    end
    if any(~isfinite(damping)) || any(damping <= 0) || any(damping >= 1)
        error('Configured resonance damping ratios must lie strictly between 0 and 1.');
    end
    if any(~isfinite(gains)) || any(gains < 0)
        error('Configured resonance modal gains must be finite and nonnegative.');
    end

    resolved = frequencies < nyquist;
    unresolved_mode_count(axis_id) = sum(~resolved);
    frequencies = frequencies(resolved);
    damping = damping(resolved);
    gains = gains(resolved);
    configured_axes(axis_id) = ~isempty(frequencies);
    transfer = complex(zeros(n_frequency, 1));

    for mode_id = 1:numel(frequencies)
        frequency_ratio = frequency / frequencies(mode_id);
        denominator = 1 - frequency_ratio.^2 + ...
            1i * 2 * damping(mode_id) * frequency_ratio;
        modal_transfer = gains(mode_id) ./ denominator;
        transfer = transfer + modal_transfer;

        % Unit-height Lorentzian-like resonance emphasis. Modal gain enters
        % quadratically because the objective integrates spectral energy.
        normalized_peak = (2 * damping(mode_id))^2 ./ ...
            ((1 - frequency_ratio.^2).^2 + ...
            (2 * damping(mode_id) * frequency_ratio).^2);
        spectral_weight(:, axis_id) = spectral_weight(:, axis_id) + ...
            spectral_peak_gain * gains(mode_id)^2 * normalized_peak;
    end

    excess_jerk_power = max(spectral_weight(:, axis_id) - 1, 0) .* ...
        jerk_power(:, axis_id);
    axis_resonant_jerk_rms(axis_id) = sqrt(sum(excess_jerk_power));
    response_power(:, axis_id) = abs(transfer).^2 .* ...
        acceleration_power(:, axis_id);
    axis_response_rms(axis_id) = sqrt(sum(response_power(:, axis_id)));
end

configured = any(configured_axes);
response_ratio = zeros(1, 3);
finite_limit = isfinite(response_limits);
response_ratio(finite_limit) = axis_response_rms(finite_limit) ./ ...
    response_limits(finite_limit);
max_response_ratio = max(response_ratio);
tolerance = 1e-6;
if isfield(cfg, 'speed') && isfield(cfg.speed, 'velocity_tolerance')
    tolerance = cfg.speed.velocity_tolerance;
end

weighted_resonant_rms = axis_weights .* axis_resonant_jerk_rms;
objective_value = 0;
if enabled && configured
    objective_value = objective_weight * norm(weighted_resonant_rms, 2);
end

weighted_power = max(spectral_weight - 1, 0) .* jerk_power;
[dominant_value, linear_id] = max(weighted_power, [], 'all', 'linear');
dominant_axis = "none";
dominant_frequency = NaN;
if dominant_value > 0
    [frequency_id, axis_id] = ind2sub(size(weighted_power), linear_id);
    axis_names = ["X" "Y" "Z"];
    dominant_axis = axis_names(axis_id);
    dominant_frequency = frequency(frequency_id);
end

metrics.configured = configured;
metrics.configured_axes = configured_axes;
metrics.unresolved_mode_count = unresolved_mode_count;
metrics.axis_resonant_jerk_rms = axis_resonant_jerk_rms;
metrics.weighted_axis_resonant_jerk_rms = weighted_resonant_rms;
metrics.axis_response_rms = axis_response_rms;
metrics.response_ratio = response_ratio;
metrics.max_response_ratio = max_response_ratio;
metrics.feasible = ~enabled || ~configured || ...
    max_response_ratio <= 1 + tolerance;
metrics.objective_weight = objective_weight;
metrics.objective_value = objective_value;
metrics.dominant_axis = dominant_axis;
metrics.dominant_frequency_hz = dominant_frequency;

if keep_spectrum && diagnostics_enabled
    metrics.frequency_hz = frequency;
    metrics.acceleration_power = acceleration_power;
    metrics.jerk_power = jerk_power;
    metrics.spectral_weight = spectral_weight;
    metrics.weighted_jerk_power = spectral_weight .* jerk_power;
    metrics.response_power = response_power;
end
end

function metrics = empty_metrics(enabled, axis_weights, response_limits)
metrics.enabled = logical(enabled);
metrics.configured = false;
metrics.configured_axes = false(1, 3);
metrics.axis_weights = axis_weights;
metrics.response_rms_limits = response_limits;
metrics.sample_time = NaN;
metrics.sample_rate_hz = NaN;
metrics.nyquist_hz = NaN;
metrics.unresolved_mode_count = zeros(1, 3);
metrics.axis_resonant_jerk_rms = zeros(1, 3);
metrics.weighted_axis_resonant_jerk_rms = zeros(1, 3);
metrics.axis_response_rms = zeros(1, 3);
metrics.response_ratio = zeros(1, 3);
metrics.max_response_ratio = 0;
metrics.feasible = true;
metrics.objective_weight = 0;
metrics.objective_value = 0;
metrics.dominant_axis = "none";
metrics.dominant_frequency_hz = NaN;
metrics.frequency_hz = [];
metrics.acceleration_power = zeros(0, 3);
metrics.jerk_power = zeros(0, 3);
metrics.spectral_weight = zeros(0, 3);
metrics.weighted_jerk_power = zeros(0, 3);
metrics.response_power = zeros(0, 3);
end

function value = resonance_option(cfg, name, fallback)
value = fallback;
if isfield(cfg, 'resonance') && isfield(cfg.resonance, name)
    value = cfg.resonance.(name);
end
end

function value = resonance_vector(cfg, name, fallback)
value = resonance_option(cfg, name, fallback);
value = value(:)';
if numel(value) ~= 3 || any(~isfinite(value) & ~isinf(value))
    error('cfg.resonance.%s must contain three finite values or Inf.', name);
end
end

function values = axis_mode_values(cfg, name, axis_id)
values = [];
if ~isfield(cfg, 'resonance') || ~isfield(cfg.resonance, name)
    return
end
configured = cfg.resonance.(name);
if iscell(configured)
    if numel(configured) ~= 3
        error('cfg.resonance.%s must be a three-element cell array.', name);
    end
    values = configured{axis_id};
else
    values = configured;
end
values = values(:)';
end

function values = expand_mode_values(values, count, fallback)
if count == 0
    values = zeros(1, 0);
elseif isempty(values)
    values = fallback * ones(1, count);
elseif isscalar(values)
    values = values * ones(1, count);
elseif numel(values) ~= count
    error('Each resonance mode must have one matching parameter value.');
end
end

function [frequency, power] = one_sided_power(signal, sample_time)
n = size(signal, 1);
centered = signal - mean(signal, 1);
window = 0.5 - 0.5 * cos(2 * pi * (0:n-1)' / max(n-1, 1));
window_power = max(mean(window.^2), eps);
spectrum = fft(centered .* window, [], 1);
two_sided = abs(spectrum).^2 / (n^2 * window_power);
n_one_sided = floor(n / 2) + 1;
power = two_sided(1:n_one_sided, :);
if rem(n, 2) == 0
    if n_one_sided > 2
        power(2:end-1, :) = 2 * power(2:end-1, :);
    end
elseif n_one_sided > 1
    power(2:end, :) = 2 * power(2:end, :);
end
frequency = (0:n_one_sided-1)' / (n * sample_time);
end

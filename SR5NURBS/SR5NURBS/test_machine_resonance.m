function report = test_machine_resonance(cfg)
%TEST_MACHINE_RESONANCE Focused regression test for resonance and Snap metrics.
if nargin < 1 || isempty(cfg)
    project_root = fileparts(mfilename('fullpath'));
    cfg = trajectory_config(project_root);
end

test_cfg = cfg;
test_cfg.resonance.enabled = true;
test_cfg.resonance.mode_frequencies_hz = {40, [], []};
test_cfg.resonance.damping_ratios = {0.03, [], []};
test_cfg.resonance.modal_gains = {1, [], []};
test_cfg.resonance.response_rms_limits = [Inf Inf Inf];

time = (0:cfg.interpolation.Ts:2)';
at_mode = sin(2*pi*40*time);
away_from_mode = sin(2*pi*15*time);
acceleration_at = [100*at_mode, zeros(numel(time), 2)];
jerk_at = [1000*at_mode, zeros(numel(time), 2)];
acceleration_away = [100*away_from_mode, zeros(numel(time), 2)];
jerk_away = [1000*away_from_mode, zeros(numel(time), 2)];

at_metrics = xyz_vibration_metric( ...
    jerk_at, test_cfg, time, acceleration_at, true);
away_metrics = xyz_vibration_metric( ...
    jerk_away, test_cfg, time, acceleration_away, false);
assert(at_metrics.resonance.configured && ...
        ~isempty(at_metrics.resonance.frequency_hz), ...
    'Configured resonance modes did not produce spectrum diagnostics.');
assert(at_metrics.resonance.objective_value > ...
        away_metrics.resonance.objective_value, ...
    'The resonance objective does not emphasize excitation at a modal frequency.');
assert(all(isfinite(at_metrics.axis_snap_rms)) && ...
        max(at_metrics.axis_snap_rms) > 0, ...
    'XYZ snap diagnostics are missing or nonfinite.');

response = at_metrics.resonance.axis_response_rms(1);
assert(response > 0 && isfinite(response), ...
    'The configured modal response was not evaluated.');
test_cfg.resonance.response_rms_limits = [0.5*response Inf Inf];
limited = xyz_vibration_metric(jerk_at, test_cfg, time, acceleration_at, false);
assert(~limited.resonance.feasible && ...
        limited.resonance.max_response_ratio >= 2 - 1e-10, ...
    'The resonance response hard limit was not enforced.');

report.resonance_frequency_hz = 40;
report.off_resonance_frequency_hz = 15;
report.resonance_objective = at_metrics.resonance.objective_value;
report.off_resonance_objective = away_metrics.resonance.objective_value;
report.hard_limit_ratio = limited.resonance.max_response_ratio;
report.snap_rms = at_metrics.axis_snap_rms;
fprintf(['Machine-resonance test passed: resonance objective %.6g, ' ...
    'off-resonance %.6g, hard-limit ratio %.6g.\n'], ...
    report.resonance_objective, report.off_resonance_objective, ...
    report.hard_limit_ratio);
end

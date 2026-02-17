%% validate_and_analyze_single.m
%  Step-by-step validation and analysis of ONE experiment from ttl_1DRF.
%  Run this script section-by-section (Ctrl+Enter per section in MATLAB)
%  to inspect your data before committing to the full pipeline.
%
%  BEFORE RUNNING: make sure the nested_RF_stimulus/src folder is on your
%  MATLAB path. The easiest way:
%     addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'))

%% ========================================================================
%  0.  CONFIGURATION — edit these two lines, then run everything below
%  ========================================================================
DATA_ROOT   = '/Users/reiserm/Documents/ttl_1DRF';
EXP_NAME    = '2025_11_10_10_17';           % <-- pick any folder name
exp_folder  = fullfile(DATA_ROOT, EXP_NAME);

%% ========================================================================
%  1.  VALIDATE: check that all required files are present
%  ========================================================================
fprintf('\n=== FILE VALIDATION ===\n');

% currentExp.mat
f = fullfile(exp_folder, 'currentExp.mat');
assert(isfile(f), 'MISSING: currentExp.mat');
fprintf('  [OK]  currentExp.mat\n');

% Log file
log_dir = fullfile(exp_folder, 'Log Files');
log_files = dir(fullfile(log_dir, 'G4_TDMS*.mat'));
assert(~isempty(log_files), 'MISSING: G4_TDMS_Log*.mat in Log Files/');
fprintf('  [OK]  %s  (%.0f MB)\n', log_files(1).name, log_files(1).bytes/1e6);

% Params
param_dir = fullfile(exp_folder, 'params');
param_files = dir(fullfile(param_dir, '*.mat'));
assert(~isempty(param_files), 'MISSING: param .mat files');
fprintf('  [OK]  params/ contains %d file(s)\n', numel(param_files));
for k = 1:numel(param_files)
    fprintf('        - %s\n', param_files(k).name);
end

% Functions
func_dir = fullfile(exp_folder, 'Functions');
func_files = dir(fullfile(func_dir, '0001*.mat'));
assert(~isempty(func_files), 'MISSING: Function 0001*.mat');
fprintf('  [OK]  Functions/ lead file: %s\n', func_files(1).name);

fprintf('\nAll required files present.\n');

%% ========================================================================
%  2.  LOAD DATA — see what we are working with
%  ========================================================================
fprintf('\n=== LOADING DATA ===\n');

% Save and restore working directory (load_protocol2_data uses cd)
orig_dir = pwd;
cleanup = onCleanup(@() cd(orig_dir));

[date_str, time_str, Log, params, pfnparam] = load_protocol2_data(exp_folder);

% Restore working directory immediately
cd(orig_dir);

fprintf('  Experiment:  %s @ %s\n', date_str, time_str);
fprintf('  ON/OFF:      %s\n', params.on_off);
fprintf('  Center [x,y]: [%d, %d]\n', params.x, params.y);

% Load metadata from currentExp
ce = load(fullfile(exp_folder, 'currentExp.mat'), 'metadata', 'pattern_order', 'func_order');
metadata = ce.metadata;
fprintf('  Strain:      %s\n', metadata.Strain);
fprintf('  Age:         %s\n', metadata.Age);
fprintf('  Side:        %s\n', metadata.Side);
fprintf('  Frame:       %d\n', metadata.Frame);

% Extract voltage and frame traces
f_data   = Log.ADC.Volts(1, :);       % frame position signal
v_data   = Log.ADC.Volts(2, :) * 10;  % voltage (mV)
median_v = median(v_data);

fprintf('  Recording:   %.1f s  (%d samples @ 10 kHz)\n', numel(v_data)/10000, numel(v_data));
fprintf('  Median Vm:   %.1f mV\n', median_v);
fprintf('  Vm range:    [%.1f, %.1f] mV  (2nd–98th pctile)\n', ...
    prctile(v_data, 2), prctile(v_data, 98));

%% ========================================================================
%  3.  QUALITY CHECK — visual inspection of the full recording
%  ========================================================================
fprintf('\n=== QUALITY CHECK FIGURE ===\n');

figure('Name', 'Quality Check — Full Recording', 'Position', [50 200 1600 500]);
tiledlayout(2,1, 'TileSpacing', 'compact');

% Frame signal
ax1 = nexttile;
plot(f_data, 'Color', [0.3 0.3 0.8]);
ylabel('Frame #');
title(sprintf('%s — %s — %s — %s', EXP_NAME, metadata.Strain, params.on_off, date_str), 'Interpreter', 'none');
xlim([0 numel(f_data)]);

% Voltage trace
ax2 = nexttile;
plot(v_data, 'Color', [0.2 0.2 0.2]);
hold on;
yline(median_v, 'r-', sprintf('median = %.1f mV', median_v), 'LineWidth', 1.5);
ylabel('Vm (mV)');
xlabel('Sample # (10 kHz)');
xlim([0 numel(v_data)]);

linkaxes([ax1, ax2], 'x');

fprintf('  → Inspect this figure. You should see:\n');
fprintf('    - Frame signal: clear repeating block structure (3 reps)\n');
fprintf('    - Voltage: stable baseline without large drift or dropouts\n');
fprintf('  A "good" experiment has Vm range < 30 mV and no sudden jumps.\n');

%% ========================================================================
%  4.  RECORDING STATISTICS — quick quality metrics
%  ========================================================================
fprintf('\n=== RECORDING QUALITY METRICS ===\n');

% Slow drift: variance of 2s moving mean
filtered_v = movmean(v_data, 20000);
drift_var = var(filtered_v);
fprintf('  Slow drift variance:    %.2f mV²\n', drift_var);

% Noise: std of high-pass residual
residual = v_data - filtered_v;
noise_std = std(residual);
fprintf('  High-freq noise std:    %.2f mV\n', noise_std);

% Check for the 3-second gaps that define stimulus blocks
zero_mask = f_data == 0;
d = diff([0 zero_mask 0]);
start_idx = find(d == 1);
end_idx = find(d == -1) - 1;
gap_len = end_idx - start_idx + 1;
n_long_gaps = sum(gap_len >= 30000);
fprintf('  Long gaps (>=3s):       %d found\n', n_long_gaps);

if n_long_gaps < 17
    warning('Expected at least 17 long gaps for 3 full reps. Got %d. Recording may be incomplete.', n_long_gaps);
else
    fprintf('  → Structural integrity looks good (>= 17 gaps = 3 complete reps).\n');
end

%% ========================================================================
%  5.  RUN ANALYSIS — bar sweep direction selectivity
%  ========================================================================
fprintf('\n=== BAR SWEEP ANALYSIS ===\n');

% Parse bar data
bar_data = parse_bar_data(f_data, v_data);
fprintf('  Parsed %d bar conditions (expect 48 = 16 dir × 3 speeds)\n', size(bar_data, 1));

% Use plot_order to rearrange from paired to sequential
plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

% Compute bar sweep responses (slow bars only, rows 1-16)
sweep_opts.baseline_range = [1000 9000];
sweep_opts.stim_trim_end  = 7000;
sweep_opts.percentile     = 98;
max_v = compute_bar_sweep_responses(bar_data, plot_order, sweep_opts);

fprintf('  Max depolarisation per direction (slow bars):\n');
fprintf('    range: [%.1f, %.1f] mV\n', min(max_v), max(max_v));
fprintf('    mean:  %.1f mV\n', mean(max_v));

% Vector sum for PD
angles_rad = deg2rad(linspace(0, 360, 17));
angles_rad = angles_rad(1:16)';
max_v_polar = [max_v; max_v(1)]; % wrap for polar plot

[d_aligned, ord, magnitude, angle_rad, fwhm, cv, thetahat, kappa] = ...
    find_PD_and_order_idx(max_v_polar, median_v);

fprintf('\n  === Direction Selectivity Metrics (slow bars) ===\n');
fprintf('  Preferred direction:  %.1f°\n', rad2deg(angle_rad));
fprintf('  Vector sum magnitude: %.3f  (0=no selectivity, 1=perfect)\n', magnitude);
fprintf('  FWHM:                 %.0f°\n', fwhm);
fprintf('  Circular variance:    %.3f  (0=sharp, 1=broad)\n', cv);
fprintf('  Von Mises kappa:      %.2f  (higher=sharper)\n', kappa);

[sym_ratio, DSI, DSI_pdnd, vector_sum] = compute_bar_response_metrics(d_aligned);
fprintf('  DSI (vector sum):     %.3f\n', DSI);
fprintf('  DSI (PD/ND):          %.3f\n', DSI_pdnd);
fprintf('  Symmetry ratio:       %.3f  (1=symmetric)\n', sym_ratio);

%% ========================================================================
%  6.  POLAR PLOT — direction tuning (manual, so you can see it step by step)
%  ========================================================================
fprintf('\n=== POLAR PLOT ===\n');

angles_plot = linspace(0, 2*pi, 17);

figure('Name', 'Direction Selectivity — Slow Bars');
polarplot(angles_plot, max_v_polar, '-o', ...
    'Color', [0.2 0.4 0.7], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'w');
hold on;

% Add PD arrow
arrow_len = max(max_v_polar) * 0.8;
polarplot([angle_rad angle_rad], [0 arrow_len], '-', ...
    'Color', [0.8 0.2 0.2], 'LineWidth', 3);

title(sprintf('PD = %.0f°  |  DSI = %.2f  |  %s', ...
    rad2deg(angle_rad), DSI, metadata.Strain), 'Interpreter', 'none');

fprintf('  → Check: does the arrow point in a clear preferred direction?\n');
fprintf('    DSI > 0.3 typically indicates good direction selectivity.\n');

%% ========================================================================
%  7.  RUN THE FULL ANALYSIS PIPELINE (newer, cleaner version)
%  ========================================================================
fprintf('\n=== FULL PIPELINE: analyze_single_experiment ===\n');
fprintf('  This will generate 4 publication-quality figures and save PDFs.\n\n');

% Restore directory before calling (it uses cd internally with onCleanup)
cd(orig_dir);

opts = struct();
opts.save_figs = true;
opts.save_dir  = fullfile(exp_folder, 'analysis_output');
% opts.flash_ylim = [-15 35];   % Uncomment and adjust if needed

analyze_single_experiment(exp_folder, opts);

fprintf('\n  Figures saved to: %s\n', opts.save_dir);
fprintf('  Generated:\n');
fprintf('    1. slow_bar_sweep_polar.pdf       — polar timeseries with PD arrow\n');
fprintf('    2. bar_flash_all_orientations.pdf  — 8×11 bar flash heatmap\n');
fprintf('    3. bar_flash_PD_ND_axis.pdf        — 1×11 flashes along PD-ND\n');
fprintf('    4. bar_flash_orthogonal_axis.pdf   — 1×11 flashes along orthogonal\n');

%% ========================================================================
%  8.  SUMMARY — print a verdict
%  ========================================================================
fprintf('\n');
fprintf('  ══════════════════════════════════════════════════\n');
fprintf('  EXPERIMENT SUMMARY: %s\n', EXP_NAME);
fprintf('  ══════════════════════════════════════════════════\n');
fprintf('  Strain:         %s\n', metadata.Strain);
fprintf('  ON/OFF:         %s\n', params.on_off);
fprintf('  Median Vm:      %.1f mV\n', median_v);
fprintf('  Drift var:      %.2f mV²\n', drift_var);
fprintf('  Noise std:      %.2f mV\n', noise_std);
fprintf('  PD:             %.0f°\n', rad2deg(angle_rad));
fprintf('  DSI (vector):   %.3f\n', DSI);
fprintf('  DSI (PD/ND):    %.3f\n', DSI_pdnd);
fprintf('  FWHM:           %.0f°\n', fwhm);
fprintf('  Circ. variance: %.3f\n', cv);
fprintf('  Symmetry:       %.3f\n', sym_ratio);
fprintf('  ──────────────────────────────────────────────────\n');
if DSI > 0.3 && n_long_gaps >= 17
    fprintf('  VERDICT:  ✓ Looks like a good experiment\n');
elseif n_long_gaps < 17
    fprintf('  VERDICT:  ✗ Recording may be incomplete\n');
else
    fprintf('  VERDICT:  ? Weak direction selectivity (DSI=%.2f), inspect figures\n', DSI);
end
fprintf('  ══════════════════════════════════════════════════\n\n');

% ESTIMATE_RF_CENTER_FROM_SWEEPS  Estimate RF center position from PD/ND
%   bar sweep temporal offsets, using cross-correlation.
%
%   For each cell, the PD and ND bar sweeps traverse the same spatial axis
%   in opposite directions. The temporal offset between the two response
%   peaks encodes how far the RF center is from the midpoint of the sweep
%   path. A bar at 28 dps and 10 kHz sampling gives 0.0028 deg/sample.
%
%   Algorithm:
%     1. Extract PD (direction 5) and ND (direction 13) sweep traces
%        (mean across 3 reps, baseline-subtracted, PD-aligned via VecNN)
%     2. Cross-correlate PD and ND traces
%     3. Lag at max xcorr = 2x temporal offset of RF center from midpoint
%     4. Convert lag to bar flash position units (2.5 deg/position)
%     5. Compare with M2 (argmax) and M5 (FWHM centroid)
%
%   Uses only the late dataset (25 cells) since M2/M5 require bar flashes.
%
%   Usage:
%     run('scripts/estimate_rf_center_from_sweeps.m')

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root   = '/Users/reiserm/Documents/ttl_1DRF';
preview_dir = fullfile(data_root, 'figure_previews');
if ~isfolder(preview_dir), mkdir(preview_dir); end

plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

% Analysis parameters
BAR_SPEED_DPS    = 28;         % degrees per second
SAMPLING_RATE    = 10000;      % Hz
DEG_PER_SAMPLE   = BAR_SPEED_DPS / SAMPLING_RATE;  % 0.0028 deg/sample
POS_SPACING_DEG  = 2.5;        % degrees per bar flash position
BL_START         = 1000;       % baseline window (samples)
BL_END           = 9000;
STIM_START       = 9000;       % stimulus window start
STIM_TRIM_END    = 7000;       % samples to trim from end

% --- Sweep geometry (from generate_bar_pos_fns.m, px_crop=30) ---
% Position function sawtooth ramps frame index from 11 to 62.
% Each frame = 1.25 pixels on the arena.  The nominal speeds (28, 56, 168
% "dps") are approximately pixels/second.  speed × duration ≈ 65 for all
% three speeds — that is the total sweep span in "speed-units" (d-units).
SWEEP_FR_LOW   = 11;
SWEEP_FR_HIGH  = 62;
SWEEP_NFRAMES  = SWEEP_FR_HIGH - SWEEP_FR_LOW;          % 51 frame intervals
SWEEP_DUR_SEC  = containers.Map({28, 56, 168}, ...
                                {2.328, 1.175, 0.400});  % from pos-fn params

% Bar flash positions (from generate_bar_flash_stimulus_xy.m):
% centre_frames ≈ 33.4 (avg for PD-relevant patterns 1-8)
% Positions 1–11 at frames cf-10, cf-8, ..., cf, ..., cf+10  (step 2)
% Position 1 is at frame (cf-10), position 6 at frame cf, position 11 at
% frame (cf+10).  In d-units the spacing = 2 * SWEEP_SPAN / SWEEP_NFRAMES.
CF_AVG          = 33.4;   % average centre_frame for ON patterns 1-8
FLASH_POS1_FRAME = CF_AVG - 10;   % frame for position 1 (= 23.4)

% Gaussian matched-filter widths to sweep (FWHM in seconds)
GAUSS_FWHM_SEC = [0.3, 0.4, 0.5, 0.6, 0.7];  % sweep these widths

% PD-aligned direction indices (after circular shift so PD = position 5)
PD_DIR_IDX = 5;    % PD = 90 deg in aligned frame
ND_DIR_IDX = 13;   % ND = 270 deg in aligned frame

% LUT
lut_path = fullfile(fileparts(mfilename('fullpath')), ...
    '..', 'src', 'analysis', 'protocol2', 'bar_lut.mat');
S_lut = load(lut_path, 'Tbl');
Tbl = S_lut.Tbl;

%% Load batch results for M2/M5 and PD directions
fprintf('Loading late batch results...\n');
S_batch = load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
batch_results = S_batch.results;
fprintf('  %d cells in batch\n', numel(batch_results));

%% Step 1: Extract raw bar sweep traces for all 25 late cells
fprintf('\n=== Extracting PD/ND bar sweep traces ===\n');

d_folders = dir(data_root);
d_folders = d_folders([d_folders.isdir]);
d_folders = d_folders(~startsWith({d_folders.name}, '.'));
valid = false(numel(d_folders), 1);
for i = 1:numel(d_folders)
    valid(i) = isfile(fullfile(data_root, d_folders(i).name, 'currentExp.mat'));
end
d_folders = d_folders(valid);
n_exp = numel(d_folders);
fprintf('Found %d experiment folders\n', n_exp);

% Storage
cell_data = struct([]);
cell_idx = 0;

for exp_idx = 1:n_exp
    folder_name = d_folders(exp_idx).name;
    exp_folder  = fullfile(data_root, folder_name);
    fprintf('[%02d/%02d] %s ... ', exp_idx, n_exp, folder_name);

    try
        orig_dir = pwd;
        cleanup  = onCleanup(@() cd(orig_dir));

        % Load raw data
        [date_str, ~, Log, ~, ~] = load_protocol2_data(exp_folder);
        f_data = Log.ADC.Volts(1, :);
        v_data = Log.ADC.Volts(2, :) * 10;

        ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
            'pattern_order', 'func_order', 'metadata');

        % Parse bar sweeps
        bar_data = parse_bar_data(f_data, v_data);

        % LUT directions
        [lut_directions, ~, ~, ~] = ...
            verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, plot_order);

        % Find matching batch entry
        batch_entry = find_batch_entry(batch_results, folder_name);
        if isempty(batch_entry)
            fprintf('SKIP (no batch match)\n');
            continue;
        end

        % Extract all 3 speeds: rows 1-16 = 28dps, 17-32 = 56dps, 33-48 = 168dps
        n_bar_rows = size(bar_data, 1);
        speed_offsets = [0, 16, 32];  % row offsets for 28, 56, 168 dps
        speeds_dps   = [28, 56, 168];
        speed_labels = {'28dps', '56dps', '168dps'};

        cell_idx = cell_idx + 1;
        cell_data(cell_idx).folder     = folder_name;
        cell_data(cell_idx).date_str   = date_str;
        cell_data(cell_idx).is_on      = batch_entry.is_on;
        cell_data(cell_idx).is_ttl     = batch_entry.is_ttl;
        cell_data(cell_idx).pd_dir     = batch_entry.pd_direction;

        for si = 1:numel(speeds_dps)
            offset = speed_offsets(si);
            if offset + 16 > n_bar_rows
                % Speed not available (fewer than expected rows)
                continue;
            end
            % PD-align traces for this speed (row_offset shifts bar_data rows)
            [traces_spd, ~, ~] = extract_and_align_traces( ...
                bar_data, plot_order, lut_directions, batch_entry.pd_direction, offset);

            % Baseline-subtract all 16 directions
            bl_end_spd = min(BL_END, 9000);  % same baseline for all speeds
            for di = 1:16
                tr = traces_spd{di};
                if isempty(tr), continue; end
                bl = mean(tr(BL_START:min(bl_end_spd, numel(tr))));
                traces_spd{di} = tr - bl;
            end

            % Store per-speed PD and ND traces
            field_pd  = sprintf('pd_trace_%s', speed_labels{si});
            field_nd  = sprintf('nd_trace_%s', speed_labels{si});
            field_all = sprintf('all_traces_%s', speed_labels{si});
            cell_data(cell_idx).(field_pd)  = traces_spd{PD_DIR_IDX};
            cell_data(cell_idx).(field_nd)  = traces_spd{ND_DIR_IDX};
            cell_data(cell_idx).(field_all) = traces_spd;
        end

        % Backward compat: keep .pd_trace/.nd_trace as 28dps
        cell_data(cell_idx).pd_trace   = cell_data(cell_idx).pd_trace_28dps;
        cell_data(cell_idx).nd_trace   = cell_data(cell_idx).nd_trace_28dps;
        cell_data(cell_idx).all_traces = cell_data(cell_idx).all_traces_28dps;

        % M2 and M5 from batch results
        cell_data(cell_idx).m2_peak_pos = batch_entry.peak_pos;
        if isfield(batch_entry, 'centroid_m5')
            cell_data(cell_idx).m5_centroid = batch_entry.centroid_m5;
        else
            cell_data(cell_idx).m5_centroid = NaN;
        end

        fprintf('OK (PD=%.0f deg, M2=%d, M5=%.2f)\n', ...
            batch_entry.pd_direction, batch_entry.peak_pos, ...
            cell_data(cell_idx).m5_centroid);

    catch ME
        fprintf('ERROR: %s\n', ME.message);
        continue;
    end
end

n_cells = numel(cell_data);
fprintf('\nLoaded %d cells\n\n', n_cells);

%% Step 2: Gaussian matched-filter RF center estimation
%  For each cell: convolve PD and ND traces independently with a Gaussian
%  bump template, find peak of the convolution output in the stimulus
%  window, compute temporal offset, convert to position units.
%
%  Sweep multiple Gaussian widths to find which best matches M2.

fprintf('=== Gaussian matched-filter RF center estimation ===\n');
fprintf('  Sweeping FWHM: %s seconds\n\n', mat2str(GAUSS_FWHM_SEC));

n_widths = numel(GAUSS_FWHM_SEC);

% Pre-build Gaussian kernels
gauss_kernels = cell(1, n_widths);
for wi = 1:n_widths
    fwhm_samples = GAUSS_FWHM_SEC(wi) * SAMPLING_RATE;
    sigma = fwhm_samples / (2 * sqrt(2 * log(2)));
    % Kernel: ±3σ (captures 99.7% of area)
    half_len = ceil(3 * sigma);
    x = (-half_len:half_len)';
    kernel = exp(-x.^2 / (2 * sigma^2));
    kernel = kernel / sum(kernel);  % normalize to unit area
    gauss_kernels{wi} = kernel;
end

% Storage: one gauss_center per width per cell
gauss_centers_all = NaN(n_cells, n_widths);   % position estimate
gauss_pd_peak_t   = NaN(n_cells, n_widths);   % PD peak time (samples)
gauss_nd_peak_t   = NaN(n_cells, n_widths);   % ND peak time (samples)
gauss_lags         = NaN(n_cells, n_widths);   % raw lag (samples)

for ci = 1:n_cells
    pd_tr = cell_data(ci).pd_trace;
    nd_tr = cell_data(ci).nd_trace;

    if isempty(pd_tr) || isempty(nd_tr), continue; end

    % Trim to common length
    min_len = min(numel(pd_tr), numel(nd_tr));
    pd_tr = pd_tr(1:min_len);
    nd_tr = nd_tr(1:min_len);

    stim_end = max(1, min_len - STIM_TRIM_END);
    if STIM_START >= stim_end, continue; end

    for wi = 1:n_widths
        kernel = gauss_kernels{wi};

        % Convolve and find peak in stimulus window
        pd_conv = conv(pd_tr, kernel, 'same');
        nd_conv = conv(nd_tr, kernel, 'same');

        [~, pd_pk_idx] = max(pd_conv(STIM_START:stim_end));
        pd_pk_time = STIM_START + pd_pk_idx - 1;

        [~, nd_pk_idx] = max(nd_conv(STIM_START:stim_end));
        nd_pk_time = STIM_START + nd_pk_idx - 1;

        lag = pd_pk_time - nd_pk_time;  % positive = PD peaks later

        % Convert lag to position units
        offset_deg = lag * DEG_PER_SAMPLE / 2;
        offset_pos = offset_deg / POS_SPACING_DEG;
        center_est = 6 + offset_pos;

        gauss_centers_all(ci, wi) = center_est;
        gauss_pd_peak_t(ci, wi)   = pd_pk_time;
        gauss_nd_peak_t(ci, wi)   = nd_pk_time;
        gauss_lags(ci, wi)         = lag;
    end
end

%% Step 2b: Sign calibration and width selection
%  For each width, check if sign needs flipping, then compute corr with M2.

m2_vals    = [cell_data.m2_peak_pos]';
m5_vals    = [cell_data.m5_centroid]';

fprintf('--- Width sweep results (corr with M2) ---\n');
fprintf('  %8s  %8s  %8s  %8s  %8s  %8s\n', ...
    'FWHM(s)', 'r(+)', 'r(flip)', 'best_r', 'mean_d', 'std_d');

width_corrs = NaN(1, n_widths);
width_means = NaN(1, n_widths);
width_stds  = NaN(1, n_widths);
gauss_centers_cal = gauss_centers_all;  % calibrated (sign-corrected)

for wi = 1:n_widths
    gc = gauss_centers_all(:, wi);
    valid = ~isnan(gc) & ~isnan(m2_vals);

    r_pos  = corr(gc(valid), m2_vals(valid));
    r_flip = corr(12 - gc(valid), m2_vals(valid));

    if r_flip > r_pos
        gauss_centers_cal(:, wi) = 12 - gc;
        best_r = r_flip;
        gc_cal = 12 - gc;
    else
        best_r = r_pos;
        gc_cal = gc;
    end

    diff_m2 = gc_cal(valid) - m2_vals(valid);
    width_corrs(wi) = best_r;
    width_means(wi) = mean(diff_m2);
    width_stds(wi)  = std(diff_m2);

    fprintf('  %8.2f  %+7.3f  %+7.3f  %+7.3f  %+7.2f  %7.2f\n', ...
        GAUSS_FWHM_SEC(wi), r_pos, r_flip, best_r, ...
        width_means(wi), width_stds(wi));
end

% Select best width
[best_corr, best_wi] = max(width_corrs);
best_fwhm = GAUSS_FWHM_SEC(best_wi);
fprintf('\n  >> Best width: FWHM = %.2f s (r = %.3f with M2)\n', best_fwhm, best_corr);

% Store best-width results in cell_data
for ci = 1:n_cells
    cell_data(ci).gauss_center = gauss_centers_cal(ci, best_wi);
    cell_data(ci).gauss_lag    = gauss_lags(ci, best_wi);
    cell_data(ci).gauss_pd_t   = gauss_pd_peak_t(ci, best_wi);
    cell_data(ci).gauss_nd_t   = gauss_nd_peak_t(ci, best_wi);
end

gauss_vals = gauss_centers_cal(:, best_wi);
valid_mask = ~isnan(gauss_vals) & ~isnan(m2_vals);

% Color scheme for scatter plots (used across all figures)
group_clrs = struct( ...
    'on_ctrl',  [0.3 0.3 0.8], ...
    'on_ttl',   [1 0.4 0.4], ...
    'off_ctrl', [0.3 0.7 0.3], ...
    'off_ttl',  [1 0.6 0.3]);
get_clr = @(c) ifelse(c.is_on, ...
    ifelse(c.is_ttl, group_clrs.on_ttl, group_clrs.on_ctrl), ...
    ifelse(c.is_ttl, group_clrs.off_ttl, group_clrs.off_ctrl));

%% Step 2c: Time-reversal + rising-edge comparison
%  The ND bar sweeps opposite to PD. Time-reversing the ND trace within
%  the stimulus window creates a "virtual PD sweep" — the response shape
%  should now resemble PD's temporal profile.
%
%  Three alignment features compared:
%    A) Gaussian peak — peak of Gaussian-convolved trace (from Step 2)
%    B) Time-reversed xcorr — xcorr between PD and flipped-ND
%    C) Rising edge — steepest ascent (max derivative of convolved trace)
%
%  Position estimation: average of PD and flipped-ND feature times,
%  calibrated to M2 via linear offset.

fprintf('\n=== Time-reversal comparison ===\n');

best_kernel = gauss_kernels{best_wi};
SR_ms = SAMPLING_RATE / 1000;  % samples per ms
POS_PER_SAMPLE = DEG_PER_SAMPLE / POS_SPACING_DEG;  % 0.00112 pos/sample

% Storage
pd_peak_t_stim   = NaN(n_cells, 1);  % PD Gauss peak (stim-relative samples)
flip_nd_peak_t   = NaN(n_cells, 1);  % flipped-ND Gauss peak (stim-relative)
pd_rise_t        = NaN(n_cells, 1);  % PD rising edge (stim-relative)
flip_nd_rise_t   = NaN(n_cells, 1);  % flipped-ND rising edge (stim-relative)
xcorr_flip_lag   = NaN(n_cells, 1);  % xcorr lag (PD vs flipped-ND)
xcorr_flip_rval  = NaN(n_cells, 1);  % xcorr peak correlation
pd_conv_stim_all = cell(n_cells, 1); % for diagnostic plotting
nd_flip_conv_all = cell(n_cells, 1);
stim_window_len  = NaN(n_cells, 1);  % stim window length per cell

for ci = 1:n_cells
    pd_tr = cell_data(ci).pd_trace;
    nd_tr = cell_data(ci).nd_trace;
    if isempty(pd_tr) || isempty(nd_tr), continue; end

    min_len = min(numel(pd_tr), numel(nd_tr));
    stim_end = max(1, min_len - STIM_TRIM_END);
    if STIM_START >= stim_end, continue; end

    % Gaussian-convolve full traces
    pd_conv = conv(pd_tr(1:min_len), best_kernel, 'same');
    nd_conv = conv(nd_tr(1:min_len), best_kernel, 'same');

    % Extract stimulus window
    pd_stim = pd_conv(STIM_START:stim_end);
    nd_stim = nd_conv(STIM_START:stim_end);
    n_stim = numel(pd_stim);
    stim_window_len(ci) = n_stim;

    % Time-reverse ND within stimulus window
    nd_flipped = flip(nd_stim);

    % Store for plotting
    pd_conv_stim_all{ci} = pd_stim;
    nd_flip_conv_all{ci} = nd_flipped;

    % --- Feature A: Gaussian peaks (stim-relative) ---
    [~, pd_pk_idx] = max(pd_stim);
    pd_peak_t_stim(ci) = pd_pk_idx;

    [~, fnd_pk_idx] = max(nd_flipped);
    flip_nd_peak_t(ci) = fnd_pk_idx;

    % --- Feature B: Xcorr between PD and flipped-ND ---
    max_lag_xc = round(n_stim / 3);
    [xc, lags_xc] = xcorr(pd_stim, nd_flipped, max_lag_xc, 'coeff');
    [xc_max, xc_max_idx] = max(xc);
    xcorr_flip_lag(ci) = lags_xc(xc_max_idx);
    xcorr_flip_rval(ci) = xc_max;

    % --- Feature C: Rising edge (max positive derivative) ---
    pd_deriv = diff(pd_stim);
    fnd_deriv = diff(nd_flipped);

    [~, pd_rise_idx] = max(pd_deriv);
    pd_rise_t(ci) = pd_rise_idx;

    [~, fnd_rise_idx] = max(fnd_deriv);
    flip_nd_rise_t(ci) = fnd_rise_idx;
end

% Timing differences in ms
peak_diff_ms = (pd_peak_t_stim - flip_nd_peak_t) / SR_ms;
rise_diff_ms = (pd_rise_t - flip_nd_rise_t) / SR_ms;
xcorr_lag_ms = xcorr_flip_lag / SR_ms;

% Position estimates: average of PD and flipped-ND feature times
peak_avg_stim = (pd_peak_t_stim + flip_nd_peak_t) / 2;
rise_avg_stim = (pd_rise_t + flip_nd_rise_t) / 2;

% Calibrate to M2 via offset (slope fixed at POS_PER_SAMPLE)
valid_cal = ~isnan(peak_avg_stim) & ~isnan(m2_vals);

raw_peak_pos = peak_avg_stim * POS_PER_SAMPLE;
peak_offset = mean(m2_vals(valid_cal) - raw_peak_pos(valid_cal));
center_B = raw_peak_pos + peak_offset;

raw_rise_pos = rise_avg_stim * POS_PER_SAMPLE;
valid_rise = valid_cal & ~isnan(rise_avg_stim);
rise_offset = mean(m2_vals(valid_rise) - raw_rise_pos(valid_rise));
center_C = raw_rise_pos + rise_offset;

% Also compute position using only PD peak time (single-direction estimate)
raw_pd_only = pd_peak_t_stim * POS_PER_SAMPLE;
pd_only_offset = mean(m2_vals(valid_cal) - raw_pd_only(valid_cal));
center_PD_only = raw_pd_only + pd_only_offset;

% Correlations with M2
v_b = valid_cal;
v_c = valid_cal & ~isnan(center_C);
r_A = corr(gauss_vals(valid_mask), m2_vals(valid_mask));
r_B = corr(center_B(v_b), m2_vals(v_b));
r_C = corr(center_C(v_c), m2_vals(v_c));
r_PD = corr(center_PD_only(v_b), m2_vals(v_b));

fprintf('\n--- Position estimation methods vs M2 ---\n');
fprintf('  Method A  (Gauss lag PD-ND):    r=%.3f, mean diff=%+.2f\n', ...
    r_A, mean(gauss_vals(valid_mask) - m2_vals(valid_mask)));
fprintf('  Method B  (TR avg peak):        r=%.3f, mean diff=%+.2f\n', ...
    r_B, mean(center_B(v_b) - m2_vals(v_b)));
fprintf('  Method C  (TR avg rise edge):   r=%.3f, mean diff=%+.2f\n', ...
    r_C, mean(center_C(v_c) - m2_vals(v_c)));
fprintf('  Method PD (PD peak only):       r=%.3f, mean diff=%+.2f\n', ...
    r_PD, mean(center_PD_only(v_b) - m2_vals(v_b)));

fprintf('\n--- Timing differences (PD vs flipped-ND, ms) ---\n');
fprintf('  Peak diff:   mean=%+.1f, std=%.1f\n', nanmean(peak_diff_ms), nanstd(peak_diff_ms));
fprintf('  Rise diff:   mean=%+.1f, std=%.1f\n', nanmean(rise_diff_ms), nanstd(rise_diff_ms));
fprintf('  Xcorr lag:   mean=%+.1f, std=%.1f\n', nanmean(xcorr_lag_ms), nanstd(xcorr_lag_ms));
fprintf('  Xcorr qual:  mean=%.3f, std=%.3f\n', nanmean(xcorr_flip_rval), nanstd(xcorr_flip_rval));

% Store in cell_data
for ci = 1:n_cells
    cell_data(ci).center_B     = center_B(ci);
    cell_data(ci).center_C     = center_C(ci);
    cell_data(ci).center_PD    = center_PD_only(ci);
    cell_data(ci).peak_diff_ms = peak_diff_ms(ci);
    cell_data(ci).rise_diff_ms = rise_diff_ms(ci);
end

%% Step 2d: Intermediate figure — PD vs flipped-ND alignment per cell
%  Shows Gaussian-convolved PD (blue) and time-reversed ND (orange)
%  overlaid in the stimulus window, with peak (triangles) and rising-edge
%  (squares) markers.

fprintf('\n--- Generating time-reversal diagnostic figures ---\n');

fig_tr = figure('Position', [50 50 1600 1000]);
tl_tr = tiledlayout(5, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

for ci = 1:min(n_cells, 25)
    ax = nexttile;
    hold(ax, 'on');

    pd_stim = pd_conv_stim_all{ci};
    nd_flip = nd_flip_conv_all{ci};
    if isempty(pd_stim) || isempty(nd_flip), continue; end

    n_stim = numel(pd_stim);
    t_ms = (1:n_stim) / SR_ms;

    % Traces
    plot(ax, t_ms, pd_stim, 'b-', 'LineWidth', 1.2);
    plot(ax, t_ms, nd_flip, '-', 'Color', [0.9 0.4 0.1], 'LineWidth', 1.2);

    % Peak markers (filled triangles)
    pd_pk = pd_peak_t_stim(ci);
    fnd_pk = flip_nd_peak_t(ci);
    if ~isnan(pd_pk) && pd_pk <= n_stim
        plot(ax, pd_pk/SR_ms, pd_stim(pd_pk), 'bv', ...
            'MarkerSize', 7, 'MarkerFaceColor', 'b');
    end
    if ~isnan(fnd_pk) && fnd_pk <= n_stim
        plot(ax, fnd_pk/SR_ms, nd_flip(fnd_pk), 'v', 'MarkerSize', 7, ...
            'MarkerFaceColor', [0.9 0.4 0.1], 'MarkerEdgeColor', [0.9 0.4 0.1]);
    end

    % Rising-edge markers (filled squares)
    pd_r = pd_rise_t(ci);
    fnd_r = flip_nd_rise_t(ci);
    if ~isnan(pd_r) && pd_r <= n_stim
        plot(ax, pd_r/SR_ms, pd_stim(pd_r), 'bs', ...
            'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 1]);
    end
    if ~isnan(fnd_r) && fnd_r <= n_stim
        plot(ax, fnd_r/SR_ms, nd_flip(fnd_r), 's', 'MarkerSize', 6, ...
            'MarkerFaceColor', [1 0.7 0.4], 'MarkerEdgeColor', [0.9 0.4 0.1]);
    end

    % Vertical line connecting PD and fND peaks to show timing gap
    if ~isnan(pd_pk) && ~isnan(fnd_pk) && pd_pk <= n_stim && fnd_pk <= n_stim
        y_mid = mean([pd_stim(pd_pk), nd_flip(fnd_pk)]);
        plot(ax, [pd_pk fnd_pk]/SR_ms, [y_mid y_mid], 'k-', 'LineWidth', 0.8);
    end

    on_str  = ifelse(cell_data(ci).is_on, 'ON', 'OFF');
    ttl_str = ifelse(cell_data(ci).is_ttl, 'TTL', 'ctrl');
    short_id = cell_data(ci).folder;
    if numel(short_id) > 10, short_id = short_id(end-9:end); end

    title(ax, sprintf('%s %s %s\n\\DeltaPk=%+.0fms \\DeltaRise=%+.0fms M2=%d', ...
        short_id, on_str, ttl_str, ...
        peak_diff_ms(ci), rise_diff_ms(ci), cell_data(ci).m2_peak_pos), ...
        'FontSize', 6, 'Interpreter', 'tex');

    set(ax, 'FontSize', 6);
    xlabel(ax, 'ms from stim start');
    if ci == 1
        legend(ax, {'PD', 'flip-ND', 'PD pk', 'fND pk', 'PD rise', 'fND rise'}, ...
            'FontSize', 4, 'Location', 'northeast');
    end
end

title(tl_tr, sprintf('PD (blue) vs time-reversed ND (orange)  |  Gauss FWHM=%.2fs  |  v=peak  s=rise edge', best_fwhm));

exportgraphics(fig_tr, fullfile(preview_dir, 'rf_center_time_reversal_alignment.png'), ...
    'Resolution', 200);
fprintf('Saved: rf_center_time_reversal_alignment.png\n');

%% Step 2e: Timing comparison summary
fig_ts = figure('Position', [100 100 1400 400]);

% Panel 1: Peak timing diff histogram
subplot(1, 4, 1);
histogram(peak_diff_ms(~isnan(peak_diff_ms)), -300:25:300, 'FaceColor', [0.3 0.3 0.8]);
hold on; xline(0, 'k--'); xline(nanmean(peak_diff_ms), 'r-', 'LineWidth', 1.5);
xlabel('PD peak - fND peak (ms)');
ylabel('Count');
title(sprintf('Peak timing diff\nmean=%+.0f ms', nanmean(peak_diff_ms)));

% Panel 2: Rising-edge timing diff histogram
subplot(1, 4, 2);
histogram(rise_diff_ms(~isnan(rise_diff_ms)), -300:25:300, 'FaceColor', [0.3 0.7 0.3]);
hold on; xline(0, 'k--'); xline(nanmean(rise_diff_ms), 'r-', 'LineWidth', 1.5);
xlabel('PD rise - fND rise (ms)');
ylabel('Count');
title(sprintf('Rise-edge timing diff\nmean=%+.0f ms', nanmean(rise_diff_ms)));

% Panel 3: Peak diff vs rise diff scatter
subplot(1, 4, 3);
v_both = ~isnan(peak_diff_ms) & ~isnan(rise_diff_ms);
hold on;
for ci = 1:n_cells
    if ~v_both(ci), continue; end
    clr = get_clr(cell_data(ci));
    plot(peak_diff_ms(ci), rise_diff_ms(ci), 'o', 'MarkerSize', 7, ...
        'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end
r_pr = corr(peak_diff_ms(v_both), rise_diff_ms(v_both));
plot([-300 300], [-300 300], 'k--', 'LineWidth', 0.5);
xlabel('Peak timing diff (ms)');
ylabel('Rise timing diff (ms)');
title(sprintf('Peak vs rise diff (r=%.2f)', r_pr));
axis equal; grid on;
xl = xlim; yl = ylim; lim = [min(xl(1),yl(1)) max(xl(2),yl(2))];
xlim(lim); ylim(lim);

% Panel 4: Xcorr quality histogram
subplot(1, 4, 4);
histogram(xcorr_flip_rval(~isnan(xcorr_flip_rval)), 0:0.05:1, 'FaceColor', [0.7 0.3 0.7]);
hold on; xline(nanmean(xcorr_flip_rval), 'r-', 'LineWidth', 1.5);
xlabel('Xcorr peak (PD vs flip-ND)');
ylabel('Count');
title(sprintf('Xcorr quality\nmean=%.2f', nanmean(xcorr_flip_rval)));

sgtitle('PD vs time-reversed ND: timing feature comparison');

exportgraphics(fig_ts, fullfile(preview_dir, 'rf_center_timing_comparison.png'), ...
    'Resolution', 200);
fprintf('Saved: rf_center_timing_comparison.png\n');

%% Step 2f: Position estimation — 4 methods vs M2
fig_pos = figure('Position', [100 100 1600 400]);

method_names = {'A: Gauss lag', 'B: TR avg peak', 'C: TR avg rise', 'PD peak only'};
method_vals  = {gauss_vals, center_B, center_C, center_PD_only};
method_r     = [r_A, r_B, r_C, r_PD];
method_valid = {valid_mask, v_b, v_c, v_b};

for mi = 1:4
    subplot(1, 4, mi);
    hold on;
    mv = method_vals{mi};
    mvalid = method_valid{mi};
    for ci = 1:n_cells
        if ~mvalid(ci), continue; end
        clr = get_clr(cell_data(ci));
        plot(m2_vals(ci), mv(ci), 'o', 'MarkerSize', 7, ...
            'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end
    plot([0 12], [0 12], 'k--', 'LineWidth', 0.5);
    xlabel('M2 (bar flash peak)');
    ylabel(method_names{mi});
    title(sprintf('%s\nr=%.3f', method_names{mi}, method_r(mi)));
    xlim([0.5 11.5]); ylim([0.5 11.5]); axis square; grid on;
end

sgtitle('RF center from bar sweeps: 4 methods vs M2 (bar flashes)');

exportgraphics(fig_pos, fullfile(preview_dir, 'rf_center_4methods_vs_m2.png'), ...
    'Resolution', 200);
fprintf('Saved: rf_center_4methods_vs_m2.png\n');

%% Step 2g: Stimulus position vs time — calibration check
%  The bar sweeps at BAR_SPEED_DPS deg/s. Within the stimulus window,
%  the bar traverses a spatial path. This plot shows how sample number
%  in the stimulus window maps to bar flash position index, confirming
%  the time→position conversion used above.
%
%  PD bar: starts at one edge, sweeps to the other at 28 dps.
%  Assuming the 11 bar flash positions span the sweep range symmetrically,
%  position 1 = near sweep start, position 11 = near sweep end.
%  Time from STIM_START: t_pos(k) = (k - offset) * POS_SPACING_DEG / DEG_PER_SAMPLE

fprintf('\n--- Stimulus position vs time calibration ---\n');

fig_cal = figure('Position', [100 100 900 500]);

% Left panel: time→position mapping (theoretical)
subplot(1, 2, 1);
hold on;

% Theoretical curve: position = t * POS_PER_SAMPLE + offset
% We calibrated offset from M2 in Step 2c
t_range = 0:100:12000;  % samples from STIM_START
pos_from_t = t_range * POS_PER_SAMPLE + peak_offset;

plot(t_range / SR_ms, pos_from_t, 'k-', 'LineWidth', 1.5);

% Overlay PD peak positions from all cells
for ci = 1:n_cells
    if isnan(pd_peak_t_stim(ci)) || isnan(m2_vals(ci)), continue; end
    clr = get_clr(cell_data(ci));
    % PD peak time → expected position (from calibration)
    plot(pd_peak_t_stim(ci) / SR_ms, m2_vals(ci), 'o', 'MarkerSize', 7, ...
        'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

% Also show flipped-ND peaks
for ci = 1:n_cells
    if isnan(flip_nd_peak_t(ci)) || isnan(m2_vals(ci)), continue; end
    clr = get_clr(cell_data(ci));
    plot(flip_nd_peak_t(ci) / SR_ms, m2_vals(ci), 's', 'MarkerSize', 6, ...
        'MarkerFaceColor', clr, 'MarkerEdgeColor', clr, 'LineWidth', 0.5);
end

xlabel('Time from stim start (ms)');
ylabel('Bar flash position index (M2)');
title(sprintf('Stimulus position vs time\nslope=%.4f pos/sample, offset=%.2f', ...
    POS_PER_SAMPLE, peak_offset));
yticks(1:11);
ylim([0.5 11.5]); grid on;
legend({'Theoretical', 'PD peaks (o)', 'flip-ND peaks (s)'}, 'Location', 'northwest');

% Right panel: residuals — how far is each method from the calibration line?
subplot(1, 2, 2);
hold on;

% PD peak residuals: M2 - (pd_peak_t * POS_PER_SAMPLE + peak_offset)
pd_resid = m2_vals - (pd_peak_t_stim * POS_PER_SAMPLE + peak_offset);
fnd_resid = m2_vals - (flip_nd_peak_t * POS_PER_SAMPLE + peak_offset);

for ci = 1:n_cells
    if isnan(pd_resid(ci)), continue; end
    clr = get_clr(cell_data(ci));
    plot(ci - 0.15, pd_resid(ci), 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    if ~isnan(fnd_resid(ci))
        plot(ci + 0.15, fnd_resid(ci), 's', 'MarkerSize', 6, ...
            'MarkerFaceColor', clr, 'MarkerEdgeColor', clr, 'LineWidth', 0.5);
    end
end

yline(0, 'k--');
xlabel('Cell index');
ylabel('Residual (M2 - predicted position)');
title(sprintf('Calibration residuals\nPD peak: std=%.2f  |  flip-ND: std=%.2f', ...
    nanstd(pd_resid), nanstd(fnd_resid)));
grid on;
legend({'PD peak (o)', 'flip-ND peak (s)'}, 'Location', 'northeast');

sgtitle(sprintf('Bar sweep → position calibration (%.1f dps, %.1f deg/pos)', ...
    BAR_SPEED_DPS, POS_SPACING_DEG));

exportgraphics(fig_cal, fullfile(preview_dir, 'rf_center_stim_position_vs_time.png'), ...
    'Resolution', 200);
fprintf('Saved: rf_center_stim_position_vs_time.png\n');

%% Step 2h: Multi-speed RF center estimation (28, 56, 168 dps)
%  For each speed, detect Gaussian peak in PD and ND traces, convert
%  temporal peak position to spatial position, and combine across speeds.
%
%  Key idea: same RF center produces different peak times at different
%  speeds, but all should map to the same spatial position. Combining
%  estimates across speeds should reduce noise.
%
%  Also runs time-reversal (PD + flipped-ND average) per speed.

fprintf('\n=== Multi-speed RF center estimation ===\n');

speeds_dps   = [28, 56, 168];
speed_labels = {'28dps', '56dps', '168dps'};
speed_colors = [0 0.4 0.8; 0.8 0.4 0; 0.6 0 0.6];  % blue, orange, purple
n_speeds = numel(speeds_dps);

% Per-speed stim duration: use actual sweep durations from position function
% (NOT the old 37.5-deg approximation, which under-estimated by ~44%)
PADDING_SAMPLES = 9000;

% For each speed, find optimal Gaussian width and estimate position
% Storage: n_cells × n_speeds
gauss_center_by_speed = NaN(n_cells, n_speeds);
tr_peak_center_by_speed = NaN(n_cells, n_speeds);
tr_rise_center_by_speed = NaN(n_cells, n_speeds);
pd_peak_time_by_speed = NaN(n_cells, n_speeds);  % raw peak times (samples)
nd_peak_time_by_speed = NaN(n_cells, n_speeds);
pd_peak_amp_by_speed  = NaN(n_cells, n_speeds);  % peak amplitudes (mV)
nd_peak_amp_by_speed  = NaN(n_cells, n_speeds);

for si = 1:n_speeds
    spd = speeds_dps(si);
    lbl = speed_labels{si};
    deg_per_samp = spd / SAMPLING_RATE;
    pos_per_samp = deg_per_samp / POS_SPACING_DEG;

    % Actual stimulus duration from position function parameters
    stim_dur = round(SWEEP_DUR_SEC(spd) * SAMPLING_RATE);
    stim_start = PADDING_SAMPLES;
    stim_end = stim_start + stim_dur;

    % Gaussian FWHM: scale with speed (broader bumps at slower speeds)
    % Try several widths scaled to expected bump duration
    bump_dur_s = stim_dur / SAMPLING_RATE;  % expected stim duration in seconds
    fwhm_candidates = bump_dur_s * [0.15, 0.2, 0.25, 0.3, 0.4];
    fwhm_candidates = max(fwhm_candidates, 0.05);  % minimum 50ms

    % Build Gaussian kernels for this speed
    kernels_spd = cell(numel(fwhm_candidates), 1);
    for fi = 1:numel(fwhm_candidates)
        sigma = fwhm_candidates(fi) / 2.355 * SAMPLING_RATE;
        half_w = round(3 * sigma);
        x = (-half_w:half_w)';
        kernels_spd{fi} = exp(-x.^2 / (2*sigma^2));
        kernels_spd{fi} = kernels_spd{fi} / sum(kernels_spd{fi});
    end

    % Find best FWHM for this speed (highest corr with M2)
    best_r_spd = -Inf;
    best_fi_spd = 1;

    for fi = 1:numel(fwhm_candidates)
        kern = kernels_spd{fi};
        centers_tmp = NaN(n_cells, 1);

        for ci = 1:n_cells
            pd_field = sprintf('pd_trace_%s', lbl);
            nd_field = sprintf('nd_trace_%s', lbl);
            if ~isfield(cell_data, pd_field) || isempty(cell_data(ci).(pd_field))
                continue;
            end
            pd = cell_data(ci).(pd_field);
            nd = cell_data(ci).(nd_field);

            % Convolve
            pd_conv = conv(pd, kern, 'same');
            nd_conv = conv(nd, kern, 'same');

            % Find peaks in stimulus window (use per-trace length)
            kern_half = round(numel(kern)/2);
            win_st = max(1, stim_start - kern_half);
            pd_win_en = min(numel(pd_conv), stim_end + kern_half);
            nd_win_en = min(numel(nd_conv), stim_end + kern_half);
            [~, pd_pk] = max(pd_conv(win_st:pd_win_en));
            [~, nd_pk] = max(nd_conv(win_st:nd_win_en));
            pd_pk = pd_pk + win_st - 1;
            nd_pk = nd_pk + win_st - 1;

            % Lag-based center (same as Step 2a)
            lag = pd_pk - nd_pk;
            offset_deg = lag * deg_per_samp / 2;
            offset_pos = offset_deg / POS_SPACING_DEG;
            centers_tmp(ci) = 6 + offset_pos;
        end

        % Determine sign (check both orientations)
        valid_tmp = ~isnan(centers_tmp) & ~isnan(m2_vals);
        if sum(valid_tmp) < 3, continue; end
        r_pos = corr(centers_tmp(valid_tmp), m2_vals(valid_tmp));
        r_neg = corr(11 - centers_tmp(valid_tmp) + 1, m2_vals(valid_tmp));
        if abs(r_neg) > abs(r_pos)
            centers_tmp = 12 - centers_tmp;
            r_val = r_neg;
        else
            r_val = r_pos;
        end
        if r_val > best_r_spd
            best_r_spd = r_val;
            best_fi_spd = fi;
        end
    end

    % Run with best FWHM for this speed
    kern = kernels_spd{best_fi_spd};
    fprintf('\n  %s: best FWHM = %.3fs (r=%.3f with M2)\n', lbl, ...
        fwhm_candidates(best_fi_spd), best_r_spd);

    for ci = 1:n_cells
        pd_field = sprintf('pd_trace_%s', lbl);
        nd_field = sprintf('nd_trace_%s', lbl);
        if ~isfield(cell_data, pd_field) || isempty(cell_data(ci).(pd_field))
            continue;
        end
        pd = cell_data(ci).(pd_field);
        nd = cell_data(ci).(nd_field);

        % Convolve
        pd_conv = conv(pd, kern, 'same');
        nd_conv = conv(nd, kern, 'same');

        % Stimulus window peaks (per-trace length)
        kern_half = round(numel(kern)/2);
        win_st = max(1, stim_start - kern_half);
        pd_win_en = min(numel(pd_conv), stim_end + kern_half);
        nd_win_en = min(numel(nd_conv), stim_end + kern_half);
        [pd_pk_val, pd_pk] = max(pd_conv(win_st:pd_win_en));
        [nd_pk_val, nd_pk] = max(nd_conv(win_st:nd_win_en));
        pd_pk = pd_pk + win_st - 1;
        nd_pk = nd_pk + win_st - 1;

        pd_peak_time_by_speed(ci, si) = pd_pk;
        nd_peak_time_by_speed(ci, si) = nd_pk;
        pd_peak_amp_by_speed(ci, si)  = pd_pk_val;
        nd_peak_amp_by_speed(ci, si)  = nd_pk_val;

        % --- Method A: PD-ND lag ---
        lag = pd_pk - nd_pk;
        offset_deg = lag * deg_per_samp / 2;
        offset_pos = offset_deg / POS_SPACING_DEG;
        gauss_center_by_speed(ci, si) = 6 + offset_pos;

        % --- Method B: Time-reversal average peak ---
        stim_st = stim_start;
        stim_en = min(stim_end, min(numel(pd_conv), numel(nd_conv)));
        pd_stim = pd_conv(stim_st:stim_en);
        nd_stim = nd_conv(stim_st:stim_en);
        nd_flip = flip(nd_stim);
        [~, pd_pk_stim] = max(pd_stim);
        [~, fnd_pk_stim] = max(nd_flip);
        avg_t = (pd_pk_stim + fnd_pk_stim) / 2;
        tr_peak_center_by_speed(ci, si) = avg_t * pos_per_samp;  % raw, uncalibrated

        % --- Method C: Time-reversal average rise edge ---
        pd_deriv = diff(pd_stim);
        fnd_deriv = diff(nd_flip);
        [~, pd_rise] = max(pd_deriv);
        [~, fnd_rise] = max(fnd_deriv);
        avg_rise = (pd_rise + fnd_rise) / 2;
        tr_rise_center_by_speed(ci, si) = avg_rise * pos_per_samp;  % raw, uncalibrated
    end

    % Sign-calibrate the lag-based estimate
    valid_spd = ~isnan(gauss_center_by_speed(:,si)) & ~isnan(m2_vals);
    if sum(valid_spd) >= 3
        r_pos = corr(gauss_center_by_speed(valid_spd,si), m2_vals(valid_spd));
        r_neg = corr(12 - gauss_center_by_speed(valid_spd,si), m2_vals(valid_spd));
        if abs(r_neg) > abs(r_pos)
            gauss_center_by_speed(:,si) = 12 - gauss_center_by_speed(:,si);
        end
    end

    % Calibrate TR methods (offset to match M2 mean)
    for method_col = 1:2
        if method_col == 1
            raw = tr_peak_center_by_speed(:,si);
        else
            raw = tr_rise_center_by_speed(:,si);
        end
        v_cal = ~isnan(raw) & ~isnan(m2_vals);
        if sum(v_cal) >= 3
            cal_offset = mean(m2_vals(v_cal) - raw(v_cal));
            if method_col == 1
                tr_peak_center_by_speed(:,si) = raw + cal_offset;
            else
                tr_rise_center_by_speed(:,si) = raw + cal_offset;
            end
        end
    end
end

% --- Combined multi-speed estimate ---
% Average position estimate across speeds (each independently calibrated)
gauss_combined = nanmean(gauss_center_by_speed, 2);
tr_peak_combined = nanmean(tr_peak_center_by_speed, 2);
tr_rise_combined = nanmean(tr_rise_center_by_speed, 2);

% Print per-speed and combined correlations
fprintf('\n--- Per-speed and combined correlations with M2 ---\n');
fprintf('  %-12s  %8s  %8s  %8s\n', 'Speed', 'Lag', 'TR-Peak', 'TR-Rise');
for si = 1:n_speeds
    v_s = ~isnan(gauss_center_by_speed(:,si)) & ~isnan(m2_vals);
    r_lag = corr(gauss_center_by_speed(v_s,si), m2_vals(v_s));
    v_p = ~isnan(tr_peak_center_by_speed(:,si)) & ~isnan(m2_vals);
    r_trp = corr(tr_peak_center_by_speed(v_p,si), m2_vals(v_p));
    v_r = ~isnan(tr_rise_center_by_speed(:,si)) & ~isnan(m2_vals);
    r_trr = corr(tr_rise_center_by_speed(v_r,si), m2_vals(v_r));
    fprintf('  %-12s  %+8.3f  %+8.3f  %+8.3f\n', speed_labels{si}, r_lag, r_trp, r_trr);
end
v_c = ~isnan(gauss_combined) & ~isnan(m2_vals);
r_comb_lag = corr(gauss_combined(v_c), m2_vals(v_c));
v_cp = ~isnan(tr_peak_combined) & ~isnan(m2_vals);
r_comb_trp = corr(tr_peak_combined(v_cp), m2_vals(v_cp));
v_cr = ~isnan(tr_rise_combined) & ~isnan(m2_vals);
r_comb_trr = corr(tr_rise_combined(v_cr), m2_vals(v_cr));
fprintf('  %-12s  %+8.3f  %+8.3f  %+8.3f\n', 'Combined', r_comb_lag, r_comb_trp, r_comb_trr);

% --- Figure: Multi-speed comparison ---
fig_ms = figure('Position', [100 100 1600 500]);

% Panel 1: Per-speed lag-based scatter vs M2
subplot(1, 3, 1); hold on;
for si = 1:n_speeds
    v_s = ~isnan(gauss_center_by_speed(:,si)) & ~isnan(m2_vals);
    plot(m2_vals(v_s), gauss_center_by_speed(v_s,si), 'o', ...
        'MarkerSize', 6, 'MarkerFaceColor', speed_colors(si,:), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', speed_labels{si});
end
plot([1 11], [1 11], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
xlabel('M2 (bar flash peak)'); ylabel('Gauss lag estimate');
title('PD-ND lag method by speed');
legend('Location', 'northwest'); grid on;
xlim([0 12]); ylim([0 12]);

% Panel 2: Per-speed TR-peak scatter vs M2
subplot(1, 3, 2); hold on;
for si = 1:n_speeds
    v_s = ~isnan(tr_peak_center_by_speed(:,si)) & ~isnan(m2_vals);
    plot(m2_vals(v_s), tr_peak_center_by_speed(v_s,si), 'o', ...
        'MarkerSize', 6, 'MarkerFaceColor', speed_colors(si,:), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', speed_labels{si});
end
plot([1 11], [1 11], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
xlabel('M2 (bar flash peak)'); ylabel('TR avg peak estimate');
title('Time-reversal peak by speed');
legend('Location', 'northwest'); grid on;
xlim([0 12]); ylim([0 12]);

% Panel 3: Combined estimates (all 3 methods) vs M2
subplot(1, 3, 3); hold on;
plot(m2_vals(v_c), gauss_combined(v_c), 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.3 0.3 0.8], 'MarkerEdgeColor', 'k', ...
    'DisplayName', sprintf('Lag combined (r=%.2f)', r_comb_lag));
plot(m2_vals(v_cp), tr_peak_combined(v_cp), 's', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.8 0.3 0.3], 'MarkerEdgeColor', 'k', ...
    'DisplayName', sprintf('TR-peak combined (r=%.2f)', r_comb_trp));
plot(m2_vals(v_cr), tr_rise_combined(v_cr), '^', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.3 0.7 0.3], 'MarkerEdgeColor', 'k', ...
    'DisplayName', sprintf('TR-rise combined (r=%.2f)', r_comb_trr));
plot([1 11], [1 11], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
xlabel('M2 (bar flash peak)'); ylabel('Combined estimate');
title('3-speed combined estimates vs M2');
legend('Location', 'northwest'); grid on;
xlim([0 12]); ylim([0 12]);

sgtitle('Multi-speed RF center estimation (28 + 56 + 168 dps)');
exportgraphics(fig_ms, fullfile(preview_dir, 'rf_center_multispeed_vs_m2.png'), ...
    'Resolution', 200);
fprintf('Saved: rf_center_multispeed_vs_m2.png\n');

%% === Step 3: Speed-regression — extract latency + RF center ===
% Model: pd_peak_deg(speed) = RF_deg + latency_pd * speed
%        nd_peak_deg(speed) = (span - RF_deg) + latency_nd * speed
% Fit each cell separately using PD and ND peaks across 3 speeds.
% intercept → RF position (degrees), slope → response latency (seconds)

fprintf('\n=== Speed-regression: latency + RF center extraction ===\n');

% Per-speed sweep span in d-units (speed × duration)
sweep_span_by_speed = NaN(1, n_speeds);
for si = 1:n_speeds
    sweep_span_by_speed(si) = speeds_dps(si) * SWEEP_DUR_SEC(speeds_dps(si));
end
fprintf('  Sweep spans (d-units): %.1f  %.1f  %.1f\n', sweep_span_by_speed);

% Flash grid in d-units:
%   FLASH_OFFSET = d-value of position 1 (frame cf-10 from sweep start)
%   FLASH_SPACING = d-value increment per position (2 frames)
% Use average sweep span for the conversion factor (d-units per frame):
SWEEP_SPAN_AVG = mean(sweep_span_by_speed);  % ~66
D_PER_FRAME = SWEEP_SPAN_AVG / SWEEP_NFRAMES;           % ~1.27
FLASH_OFFSET  = (FLASH_POS1_FRAME - SWEEP_FR_LOW) * D_PER_FRAME;  % ~15.8
FLASH_SPACING = 2 * D_PER_FRAME;                                   % ~2.55
fprintf('  FLASH_OFFSET = %.2f d-units  (pos 1 from sweep start)\n', FLASH_OFFSET);
fprintf('  FLASH_SPACING = %.3f d-units  (per position step)\n', FLASH_SPACING);
fprintf('  Position 6 check: %.2f d-units  (sweep midpoint = %.2f)\n', ...
    FLASH_OFFSET + 5*FLASH_SPACING, SWEEP_SPAN_AVG/2);

% Convert peak sample indices to d-units from sweep start
pd_peak_deg_by_speed = NaN(n_cells, n_speeds);
nd_peak_deg_by_speed = NaN(n_cells, n_speeds);
% Also pre-convert ND to PD-start coordinates for cleaner regression
nd_peak_deg_pdcoord  = NaN(n_cells, n_speeds);

for si = 1:n_speeds
    spd = speeds_dps(si);
    deg_per_samp = spd / SAMPLING_RATE;
    span_si = sweep_span_by_speed(si);
    for ci = 1:n_cells
        if ~isnan(pd_peak_time_by_speed(ci, si))
            pd_peak_deg_by_speed(ci, si) = ...
                (pd_peak_time_by_speed(ci, si) - PADDING_SAMPLES) * deg_per_samp;
        end
        if ~isnan(nd_peak_time_by_speed(ci, si))
            nd_d = (nd_peak_time_by_speed(ci, si) - PADDING_SAMPLES) * deg_per_samp;
            nd_peak_deg_by_speed(ci, si) = nd_d;
            % Convert to PD-start coordinate: RF is at (span - nd_d) from PD start
            nd_peak_deg_pdcoord(ci, si) = span_si - nd_d;
        end
    end
end

% Per-cell linear regression: peak_deg = RF_d + latency * speed
% PD: directly gives RF_d from intercept
% ND: pre-converted to PD coords, so intercept also gives RF_d directly
rf_deg_from_pd   = NaN(n_cells, 1);  % intercept of PD fit (d-units from PD start)
rf_deg_from_nd   = NaN(n_cells, 1);  % intercept of ND fit (d-units from PD start)
latency_pd_ms    = NaN(n_cells, 1);  % slope of PD fit (ms)
latency_nd_ms    = NaN(n_cells, 1);  % slope of ND fit (ms)
rf_deg_combined  = NaN(n_cells, 1);  % average of PD and ND RF estimates
latency_avg_ms   = NaN(n_cells, 1);  % average latency

speed_vec = speeds_dps(:);  % [28; 56; 168]

for ci = 1:n_cells
    % PD fit: pd_peak_deg = RF_d + latency_pd * speed
    pd_vals = pd_peak_deg_by_speed(ci, :)';
    valid_pd = ~isnan(pd_vals);
    if sum(valid_pd) >= 2
        X_pd = [ones(sum(valid_pd), 1), speed_vec(valid_pd)];
        b_pd = X_pd \ pd_vals(valid_pd);
        rf_deg_from_pd(ci) = b_pd(1);           % intercept = RF position
        latency_pd_ms(ci)  = b_pd(2) * 1000;    % slope → ms
    end

    % ND fit (in PD-start coordinates):
    % nd_pdcoord = RF_d + latency_nd * speed  (note: latency_nd < 0 is OK)
    nd_vals = nd_peak_deg_pdcoord(ci, :)';
    valid_nd = ~isnan(nd_vals);
    if sum(valid_nd) >= 2
        X_nd = [ones(sum(valid_nd), 1), speed_vec(valid_nd)];
        b_nd = X_nd \ nd_vals(valid_nd);
        rf_deg_from_nd(ci) = b_nd(1);           % intercept = RF position (PD coords)
        latency_nd_ms(ci)  = b_nd(2) * 1000;
    end

    % Combined RF estimate (average of PD-derived and ND-derived)
    estimates = [rf_deg_from_pd(ci), rf_deg_from_nd(ci)];
    estimates = estimates(~isnan(estimates));
    if ~isempty(estimates)
        rf_deg_combined(ci) = mean(estimates);
    end

    % Average latency
    lats = [latency_pd_ms(ci), latency_nd_ms(ci)];
    lats = lats(~isnan(lats));
    if ~isempty(lats)
        latency_avg_ms(ci) = mean(lats);
    end
end

% Convert RF center from d-units to bar flash position units
% Position 1 starts at FLASH_OFFSET d-units from sweep start;
% each subsequent position is FLASH_SPACING d-units apart.
rf_pos_from_pd  = (rf_deg_from_pd  - FLASH_OFFSET) / FLASH_SPACING + 1;
rf_pos_from_nd  = (rf_deg_from_nd  - FLASH_OFFSET) / FLASH_SPACING + 1;
rf_pos_combined = (rf_deg_combined - FLASH_OFFSET) / FLASH_SPACING + 1;

% Correlations with M2
v_pd  = ~isnan(rf_pos_from_pd) & ~isnan(m2_vals);
v_nd  = ~isnan(rf_pos_from_nd) & ~isnan(m2_vals);
v_cmb = ~isnan(rf_pos_combined) & ~isnan(m2_vals);

r_pd  = corr(rf_pos_from_pd(v_pd), m2_vals(v_pd));
r_nd  = corr(rf_pos_from_nd(v_nd), m2_vals(v_nd));
r_cmb = corr(rf_pos_combined(v_cmb), m2_vals(v_cmb));

% Note: sign-flip logic removed — with correct spatial calibration (ND
% pre-converted to PD coords, FLASH_OFFSET/FLASH_SPACING), PD and ND
% should both give positions in the 1–11 range without needing flips.
% If correlations are negative, that indicates a deeper problem.

fprintf('\n--- Speed-regression RF center vs M2 ---\n');
fprintf('  PD-derived:  r = %.3f  (n=%d)\n', r_pd, sum(v_pd));
fprintf('  ND-derived:  r = %.3f  (n=%d)\n', r_nd, sum(v_nd));
fprintf('  Combined:    r = %.3f  (n=%d)\n', r_cmb, sum(v_cmb));

fprintf('\n--- Response latencies (ms) ---\n');
fprintf('  %-20s  %8s  %8s  %8s\n', '', 'PD', 'ND', 'PD-ND');
v_both = ~isnan(latency_pd_ms) & ~isnan(latency_nd_ms);
fprintf('  %-20s  %+7.1f  %+7.1f  %+7.1f\n', 'Mean', ...
    nanmean(latency_pd_ms), nanmean(latency_nd_ms), ...
    nanmean(latency_pd_ms(v_both) - latency_nd_ms(v_both)));
fprintf('  %-20s  %7.1f  %7.1f  %7.1f\n', 'Std', ...
    nanstd(latency_pd_ms), nanstd(latency_nd_ms), ...
    nanstd(latency_pd_ms(v_both) - latency_nd_ms(v_both)));
fprintf('  %-20s  %+7.1f  %+7.1f\n', 'Median', ...
    nanmedian(latency_pd_ms), nanmedian(latency_nd_ms));

% Per-cell table
fprintf('\n--- Per-cell speed regression results ---\n');
fprintf('  %-22s  %4s  %4s  %6s  %6s  %6s  %7s  %7s  %7s\n', ...
    'Folder', 'ON', 'TTL', 'M2', 'RF_pd', 'RF_nd', 'Lat_PD', 'Lat_ND', 'RF_cmb');
fprintf('  %s\n', repmat('-', 1, 90));
for ci = 1:n_cells
    on_str  = ifelse(cell_data(ci).is_on, 'ON', 'OFF');
    ttl_str = ifelse(cell_data(ci).is_ttl, 'TTL', 'ctrl');
    fprintf('  %-22s  %4s  %4s  %5d  %6.1f  %6.1f  %+6.1f  %+6.1f  %6.1f\n', ...
        cell_data(ci).folder, on_str, ttl_str, cell_data(ci).m2_peak_pos, ...
        rf_pos_from_pd(ci), rf_pos_from_nd(ci), ...
        latency_pd_ms(ci), latency_nd_ms(ci), rf_pos_combined(ci));
end

% --- Peak amplitude diagnostics ---
fprintf('\n--- PD vs ND peak amplitudes (convolved, mV) ---\n');
fprintf('  %-10s  %8s  %8s  %8s  %8s  %8s  %8s\n', ...
    'Speed', 'PD_mean', 'PD_std', 'ND_mean', 'ND_std', 'PD/ND', 'PD>ND%');
for si = 1:n_speeds
    pd_amps = pd_peak_amp_by_speed(:, si);
    nd_amps = nd_peak_amp_by_speed(:, si);
    v = ~isnan(pd_amps) & ~isnan(nd_amps);
    ratio = pd_amps(v) ./ nd_amps(v);
    pct_pd_bigger = 100 * mean(pd_amps(v) > nd_amps(v));
    fprintf('  %-10s  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %7.0f%%\n', ...
        speed_labels{si}, nanmean(pd_amps), nanstd(pd_amps), ...
        nanmean(nd_amps), nanstd(nd_amps), nanmedian(ratio), pct_pd_bigger);
end

% Per-cell amplitude table
fprintf('\n--- Per-cell PD/ND amplitude ratio by speed ---\n');
fprintf('  %-22s  %6s  %6s  %6s  %6s  %6s  %6s\n', ...
    'Folder', 'PD_28', 'ND_28', 'PD_56', 'ND_56', 'PD_168', 'ND_168');
fprintf('  %s\n', repmat('-', 1, 65));
for ci = 1:n_cells
    fprintf('  %-22s', cell_data(ci).folder);
    for si = 1:n_speeds
        fprintf('  %6.2f', pd_peak_amp_by_speed(ci, si));
        fprintf('  %6.2f', nd_peak_amp_by_speed(ci, si));
    end
    fprintf('\n');
end

%% === Step 3b: Rise-edge speed regression ===
% Use max-derivative (steepest rising edge) instead of peak to extract latency.
% Rise-edge timing is less affected by DS integration time and more robust
% on weak ND traces.

fprintf('\n=== Rise-edge speed regression ===\n');

pd_rise_time_by_speed = NaN(n_cells, n_speeds);
nd_rise_time_by_speed = NaN(n_cells, n_speeds);
pd_rise_deg_by_speed  = NaN(n_cells, n_speeds);
nd_rise_deg_by_speed  = NaN(n_cells, n_speeds);
nd_rise_deg_pdcoord   = NaN(n_cells, n_speeds);  % ND in PD-start coordinates

for si = 1:n_speeds
    spd = speeds_dps(si);
    lbl = speed_labels{si};
    deg_per_samp = spd / SAMPLING_RATE;

    stim_dur = round(SWEEP_DUR_SEC(spd) * SAMPLING_RATE);
    span_si_r = spd * SWEEP_DUR_SEC(spd);  % sweep span for this speed (d-units)
    stim_start_s = PADDING_SAMPLES;
    stim_end_s = stim_start_s + stim_dur;

    % Use same best kernel from peak-based analysis
    % (rebuild it — same FWHM selection approach)
    bump_dur_s = stim_dur / SAMPLING_RATE;
    fwhm_candidates = bump_dur_s * [0.15, 0.2, 0.25, 0.3, 0.4];
    fwhm_candidates = max(fwhm_candidates, 0.05);

    % Use the FWHM that gave best peak-based correlation
    % (recompute with a moderate width — 20% of stim duration)
    fwhm_rise = bump_dur_s * 0.2;
    fwhm_rise = max(fwhm_rise, 0.05);
    sigma_r = fwhm_rise / 2.355 * SAMPLING_RATE;
    half_w_r = round(3 * sigma_r);
    x_r = (-half_w_r:half_w_r)';
    kern_r = exp(-x_r.^2 / (2*sigma_r^2));
    kern_r = kern_r / sum(kern_r);

    for ci = 1:n_cells
        pd_field = sprintf('pd_trace_%s', lbl);
        nd_field = sprintf('nd_trace_%s', lbl);
        if ~isfield(cell_data, pd_field) || isempty(cell_data(ci).(pd_field))
            continue;
        end
        pd = cell_data(ci).(pd_field);
        nd = cell_data(ci).(nd_field);

        pd_conv = conv(pd, kern_r, 'same');
        nd_conv = conv(nd, kern_r, 'same');

        % Derivative of convolved trace
        pd_deriv = diff(pd_conv);
        nd_deriv = diff(nd_conv);

        % Search for max derivative in stimulus window
        kern_half = round(numel(kern_r)/2);
        win_st = max(1, stim_start_s - kern_half);
        pd_win_en = min(numel(pd_deriv), stim_end_s + kern_half);
        nd_win_en = min(numel(nd_deriv), stim_end_s + kern_half);

        [~, pd_rise_idx] = max(pd_deriv(win_st:pd_win_en));
        [~, nd_rise_idx] = max(nd_deriv(win_st:nd_win_en));
        pd_rise_idx = pd_rise_idx + win_st - 1;
        nd_rise_idx = nd_rise_idx + win_st - 1;

        pd_rise_time_by_speed(ci, si) = pd_rise_idx;
        nd_rise_time_by_speed(ci, si) = nd_rise_idx;
        pd_rise_deg_by_speed(ci, si) = (pd_rise_idx - PADDING_SAMPLES) * deg_per_samp;
        nd_d_r = (nd_rise_idx - PADDING_SAMPLES) * deg_per_samp;
        nd_rise_deg_by_speed(ci, si) = nd_d_r;
        % Pre-convert ND to PD-start coordinates
        nd_rise_deg_pdcoord(ci, si) = span_si_r - nd_d_r;
    end
end

% Per-cell rise-edge regression: rise_deg = RF_deg + latency * speed
rf_deg_rise_pd   = NaN(n_cells, 1);
rf_deg_rise_nd   = NaN(n_cells, 1);
latency_rise_pd_ms = NaN(n_cells, 1);
latency_rise_nd_ms = NaN(n_cells, 1);

for ci = 1:n_cells
    % PD rise-edge fit
    pd_vals = pd_rise_deg_by_speed(ci, :)';
    valid_pd = ~isnan(pd_vals);
    if sum(valid_pd) >= 2
        X_pd = [ones(sum(valid_pd), 1), speed_vec(valid_pd)];
        b_pd = X_pd \ pd_vals(valid_pd);
        rf_deg_rise_pd(ci) = b_pd(1);
        latency_rise_pd_ms(ci) = b_pd(2) * 1000;
    end

    % ND rise-edge fit (in PD-start coordinates — same as Step 3)
    nd_vals = nd_rise_deg_pdcoord(ci, :)';
    valid_nd = ~isnan(nd_vals);
    if sum(valid_nd) >= 2
        X_nd = [ones(sum(valid_nd), 1), speed_vec(valid_nd)];
        b_nd = X_nd \ nd_vals(valid_nd);
        rf_deg_rise_nd(ci) = b_nd(1);   % intercept = RF position (PD coords)
        latency_rise_nd_ms(ci) = b_nd(2) * 1000;
    end
end

% Convert to position units (same calibration as Step 3)
rf_pos_rise_pd  = (rf_deg_rise_pd  - FLASH_OFFSET) / FLASH_SPACING + 1;
rf_pos_rise_nd  = (rf_deg_rise_nd  - FLASH_OFFSET) / FLASH_SPACING + 1;

% Combined rise-edge RF estimate
rf_pos_rise_combined = NaN(n_cells, 1);
for ci = 1:n_cells
    est = [rf_pos_rise_pd(ci), rf_pos_rise_nd(ci)];
    est = est(~isnan(est));
    if ~isempty(est), rf_pos_rise_combined(ci) = mean(est); end
end

% Correlations
v_rpd = ~isnan(rf_pos_rise_pd) & ~isnan(m2_vals);
v_rnd = ~isnan(rf_pos_rise_nd) & ~isnan(m2_vals);
v_rc  = ~isnan(rf_pos_rise_combined) & ~isnan(m2_vals);
r_rise_pd  = corr(rf_pos_rise_pd(v_rpd), m2_vals(v_rpd));
r_rise_nd  = corr(rf_pos_rise_nd(v_rnd), m2_vals(v_rnd));
r_rise_cmb = corr(rf_pos_rise_combined(v_rc), m2_vals(v_rc));

fprintf('\n--- Rise-edge speed-regression RF center vs M2 ---\n');
fprintf('  PD-derived:  r = %.3f  (n=%d)\n', r_rise_pd, sum(v_rpd));
fprintf('  ND-derived:  r = %.3f  (n=%d)\n', r_rise_nd, sum(v_rnd));
fprintf('  Combined:    r = %.3f  (n=%d)\n', r_rise_cmb, sum(v_rc));

fprintf('\n--- Rise-edge latencies (ms) ---\n');
fprintf('  %-20s  %8s  %8s  %8s\n', '', 'PD', 'ND', 'PD-ND');
v_both_r = ~isnan(latency_rise_pd_ms) & ~isnan(latency_rise_nd_ms);
fprintf('  %-20s  %+7.1f  %+7.1f  %+7.1f\n', 'Mean', ...
    nanmean(latency_rise_pd_ms), nanmean(latency_rise_nd_ms), ...
    nanmean(latency_rise_pd_ms(v_both_r) - latency_rise_nd_ms(v_both_r)));
fprintf('  %-20s  %7.1f  %7.1f  %7.1f\n', 'Std', ...
    nanstd(latency_rise_pd_ms), nanstd(latency_rise_nd_ms), ...
    nanstd(latency_rise_pd_ms(v_both_r) - latency_rise_nd_ms(v_both_r)));
fprintf('  %-20s  %+7.1f  %+7.1f\n', 'Median', ...
    nanmedian(latency_rise_pd_ms), nanmedian(latency_rise_nd_ms));

% Comparison: peak-based vs rise-edge latencies
fprintf('\n--- Peak-based vs Rise-edge latency comparison ---\n');
fprintf('  %-20s  %8s  %8s  %8s  %8s\n', '', 'Peak_PD', 'Rise_PD', 'Peak_ND', 'Rise_ND');
fprintf('  %-20s  %+7.1f  %+7.1f  %+7.1f  %+7.1f\n', 'Mean (ms)', ...
    nanmean(latency_pd_ms), nanmean(latency_rise_pd_ms), ...
    nanmean(latency_nd_ms), nanmean(latency_rise_nd_ms));
fprintf('  %-20s  %7.1f  %7.1f  %7.1f  %7.1f\n', 'Std (ms)', ...
    nanstd(latency_pd_ms), nanstd(latency_rise_pd_ms), ...
    nanstd(latency_nd_ms), nanstd(latency_rise_nd_ms));

% --- Figure: Speed regression diagnostics (expanded to 3x3) ---
fig_sr = figure('Position', [50 50 1800 1600]);

% Panel 1: RF position estimates vs M2 (peak-based)
subplot(3, 3, 1); hold on;
plot(m2_vals(v_pd), rf_pos_from_pd(v_pd), 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0 0.4 0.8], 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
    'DisplayName', sprintf('PD-derived (r=%.2f)', r_pd));
plot(m2_vals(v_nd), rf_pos_from_nd(v_nd), 's', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
    'DisplayName', sprintf('ND-derived (r=%.2f)', r_nd));
plot([1 11], [1 11], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
xlabel('M2 (bar flash peak)'); ylabel('Speed-regression RF estimate');
title('Peak-based: RF center vs M2');
legend('Location', 'northwest'); grid on;
xlim([0 12]); ylim([0 12]);

% Panel 2: Combined RF estimate vs M2
subplot(3, 3, 2); hold on;
plot(m2_vals(v_cmb), rf_pos_combined(v_cmb), 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.4 0.2 0.6], 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
plot([1 11], [1 11], 'k--', 'LineWidth', 0.5);
xlabel('M2 (bar flash peak)'); ylabel('Speed-regress combined');
title(sprintf('Peak-based combined vs M2  r=%.2f', r_cmb));
grid on; xlim([0 12]); ylim([0 12]);

% Panel 3: PD latency vs ND latency (peak-based)
subplot(3, 3, 3); hold on;
v_lat = ~isnan(latency_pd_ms) & ~isnan(latency_nd_ms);
for ci = find(v_lat)'
    if cell_data(ci).is_on
        mc = [0.3 0.6 1]; % blue for ON
    else
        mc = [1 0.5 0.3]; % orange for OFF
    end
    if cell_data(ci).is_ttl
        mk = 's'; % square for TTL
    else
        mk = 'o'; % circle for ctrl
    end
    plot(latency_pd_ms(ci), latency_nd_ms(ci), mk, 'MarkerSize', 8, ...
        'MarkerFaceColor', mc, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end
ax_lim = [min([latency_pd_ms; latency_nd_ms])-5, max([latency_pd_ms; latency_nd_ms])+5];
plot(ax_lim, ax_lim, 'k--', 'LineWidth', 0.5);
xlabel('PD latency (ms)'); ylabel('ND latency (ms)');
title('Peak-based: PD vs ND latency');
grid on; axis equal;
% Legend entries
plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.3 0.6 1], 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'ON ctrl');
plot(NaN, NaN, 's', 'MarkerFaceColor', [0.3 0.6 1], 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'ON TTL');
plot(NaN, NaN, 'o', 'MarkerFaceColor', [1 0.5 0.3], 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'OFF ctrl');
plot(NaN, NaN, 's', 'MarkerFaceColor', [1 0.5 0.3], 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'OFF TTL');
legend('Location', 'northwest');

% Panel 4: Rise-edge RF position estimates vs M2
subplot(3, 3, 4); hold on;
plot(m2_vals(v_rpd), rf_pos_rise_pd(v_rpd), 'o', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0 0.4 0.8], 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
    'DisplayName', sprintf('PD rise (r=%.2f)', r_rise_pd));
plot(m2_vals(v_rnd), rf_pos_rise_nd(v_rnd), 's', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
    'DisplayName', sprintf('ND rise (r=%.2f)', r_rise_nd));
plot([1 11], [1 11], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
xlabel('M2 (bar flash peak)'); ylabel('Rise-edge RF estimate');
title('Rise-edge: RF center vs M2');
legend('Location', 'northwest'); grid on;
xlim([0 12]); ylim([0 12]);

% Panel 5: Rise-edge PD vs ND latency
subplot(3, 3, 5); hold on;
v_lat_r = ~isnan(latency_rise_pd_ms) & ~isnan(latency_rise_nd_ms);
for ci = find(v_lat_r)'
    if cell_data(ci).is_on, mc = [0.3 0.6 1]; else, mc = [1 0.5 0.3]; end
    if cell_data(ci).is_ttl, mk = 's'; else, mk = 'o'; end
    plot(latency_rise_pd_ms(ci), latency_rise_nd_ms(ci), mk, 'MarkerSize', 8, ...
        'MarkerFaceColor', mc, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end
all_rise_lat = [latency_rise_pd_ms; latency_rise_nd_ms];
ax_lim_r = [nanmin(all_rise_lat)-5, nanmax(all_rise_lat)+5];
plot(ax_lim_r, ax_lim_r, 'k--', 'LineWidth', 0.5);
xlabel('PD rise latency (ms)'); ylabel('ND rise latency (ms)');
title('Rise-edge: PD vs ND latency');
grid on; axis equal;
plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.3 0.6 1], 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'ON ctrl');
plot(NaN, NaN, 's', 'MarkerFaceColor', [0.3 0.6 1], 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'ON TTL');
plot(NaN, NaN, 'o', 'MarkerFaceColor', [1 0.5 0.3], 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'OFF ctrl');
plot(NaN, NaN, 's', 'MarkerFaceColor', [1 0.5 0.3], 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'OFF TTL');
legend('Location', 'northwest');

% Panel 6: Latency histograms (peak vs rise comparison)
subplot(3, 3, 6); hold on;
edges = -20:10:200;
edges = -20:10:200;
histogram(latency_pd_ms(~isnan(latency_pd_ms)), edges, ...
    'FaceColor', [0 0.4 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'w', ...
    'DisplayName', sprintf('Peak PD %.0f\\pm%.0fms', nanmean(latency_pd_ms), nanstd(latency_pd_ms)));
histogram(latency_nd_ms(~isnan(latency_nd_ms)), edges, ...
    'FaceColor', [0.8 0.2 0.2], 'FaceAlpha', 0.4, 'EdgeColor', 'w', ...
    'DisplayName', sprintf('Peak ND %.0f\\pm%.0fms', nanmean(latency_nd_ms), nanstd(latency_nd_ms)));
histogram(latency_rise_pd_ms(~isnan(latency_rise_pd_ms)), edges, ...
    'FaceColor', [0 0.4 0.8], 'FaceAlpha', 0.3, 'EdgeColor', [0 0.4 0.8], 'LineStyle', '--', ...
    'DisplayName', sprintf('Rise PD %.0f\\pm%.0fms', nanmean(latency_rise_pd_ms), nanstd(latency_rise_pd_ms)));
histogram(latency_rise_nd_ms(~isnan(latency_rise_nd_ms)), edges, ...
    'FaceColor', [0.8 0.2 0.2], 'FaceAlpha', 0.3, 'EdgeColor', [0.8 0.2 0.2], 'LineStyle', '--', ...
    'DisplayName', sprintf('Rise ND %.0f\\pm%.0fms', nanmean(latency_rise_nd_ms), nanstd(latency_rise_nd_ms)));
xlabel('Latency (ms)'); ylabel('Count');
title('Peak vs rise-edge latencies'); legend('Location', 'northeast', 'FontSize', 6);

% Panel 7: PD peak position vs speed (per-cell lines)
subplot(3, 3, 7); hold on;
for ci = 1:n_cells
    pd_vals = pd_peak_deg_by_speed(ci, :);
    valid = ~isnan(pd_vals);
    if sum(valid) < 2, continue; end
    if cell_data(ci).is_ttl
        lc = [0.8 0.3 0.3 0.3]; % red-ish, transparent
    else
        lc = [0.3 0.3 0.8 0.3]; % blue-ish, transparent
    end
    plot(speed_vec(valid), pd_vals(valid), '-o', 'Color', lc, ...
        'MarkerSize', 4, 'MarkerFaceColor', lc(1:3), 'LineWidth', 0.8);
    % Plot regression line
    if ~isnan(rf_deg_from_pd(ci)) && ~isnan(latency_pd_ms(ci))
        spd_range = [0 180];
        fit_line = rf_deg_from_pd(ci) + (latency_pd_ms(ci)/1000) * spd_range;
        plot(spd_range, fit_line, '-', 'Color', lc(1:3), 'LineWidth', 0.3);
    end
end
xlabel('Speed (deg/s)'); ylabel('PD peak position (deg)');
title('PD peak vs speed (slope = latency)');
xlim([-5 180]); grid on;
plot([0 180], [SWEEP_SPAN_AVG/2 SWEEP_SPAN_AVG/2], 'k:', 'LineWidth', 0.5);  % sweep center

% Panel 8: ND peak position vs speed (per-cell lines)
subplot(3, 3, 8); hold on;
for ci = 1:n_cells
    nd_vals = nd_peak_deg_by_speed(ci, :);
    valid = ~isnan(nd_vals);
    if sum(valid) < 2, continue; end
    if cell_data(ci).is_ttl
        lc = [0.8 0.3 0.3 0.3];
    else
        lc = [0.3 0.3 0.8 0.3];
    end
    plot(speed_vec(valid), nd_vals(valid), '-o', 'Color', lc, ...
        'MarkerSize', 4, 'MarkerFaceColor', lc(1:3), 'LineWidth', 0.8);
    % Plot regression line
    if ~isnan(rf_deg_from_nd(ci)) && ~isnan(latency_nd_ms(ci))
        spd_range = [0 180];
        rf_nd_intercept = SWEEP_SPAN_AVG - rf_deg_from_nd(ci);  % ND intercept (ND-start coords)
        fit_line = rf_nd_intercept + (latency_nd_ms(ci)/1000) * spd_range;
        plot(spd_range, fit_line, '-', 'Color', lc(1:3), 'LineWidth', 0.3);
    end
end
xlabel('Speed (deg/s)'); ylabel('ND peak position (deg)');
title('ND peak vs speed (slope = latency)');
xlim([-5 180]); grid on;
plot([0 180], [SWEEP_SPAN_AVG/2 SWEEP_SPAN_AVG/2], 'k:', 'LineWidth', 0.5);

% Panel 9: PD/ND amplitude ratio by speed
subplot(3, 3, 9); hold on;
amp_ratio_by_speed = pd_peak_amp_by_speed ./ nd_peak_amp_by_speed;
positions = [1 2 3];
for si = 1:n_speeds
    ratios = amp_ratio_by_speed(:, si);
    valid_r = ~isnan(ratios) & ~isinf(ratios);
    jitter = 0.15 * (rand(sum(valid_r), 1) - 0.5);
    % Color by ON/OFF
    for ci = find(valid_r)'
        if cell_data(ci).is_on, mc = [0.3 0.6 1]; else, mc = [1 0.5 0.3]; end
        if cell_data(ci).is_ttl, mk = 's'; else, mk = 'o'; end
        plot(si + 0.15*(rand-0.5), amp_ratio_by_speed(ci, si), mk, ...
            'MarkerSize', 6, 'MarkerFaceColor', mc, 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    end
    % Mean bar
    plot([si-0.3 si+0.3], [nanmedian(ratios) nanmedian(ratios)], 'k-', 'LineWidth', 2);
end
plot([0.5 3.5], [1 1], 'k:', 'LineWidth', 0.5);  % unity line
set(gca, 'XTick', 1:3, 'XTickLabel', speed_labels);
xlabel('Speed'); ylabel('PD/ND amplitude ratio');
title('PD vs ND peak amplitude');
grid on; ylim([0 max(amp_ratio_by_speed(:))*1.1]);

sgtitle('Speed-regression: peak vs rise-edge latency + RF center');
exportgraphics(fig_sr, fullfile(preview_dir, 'rf_center_speed_regression.png'), ...
    'Resolution', 200);
fprintf('Saved: rf_center_speed_regression.png\n');

% --- Figure: PD vs ND time-reversal alignment diagnostic ---
% Select 8 example cells spread across M2 range
[~, sort_by_m2] = sort(m2_vals);
n_examples = min(8, n_cells);
step = max(1, floor(n_cells / n_examples));
example_idx = sort_by_m2(1:step:end);
example_idx = example_idx(1:min(n_examples, numel(example_idx)));

fig_trdiag = figure('Position', [50 50 1800 1800]);
tl_trdiag = tiledlayout(n_examples, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for ei = 1:numel(example_idx)
    ci = example_idx(ei);
    on_str  = ifelse(cell_data(ci).is_on, 'ON', 'OFF');
    ttl_str = ifelse(cell_data(ci).is_ttl, 'TTL', 'ctrl');
    cell_lbl = sprintf('%s %s/%s M2=%d', cell_data(ci).folder(end-4:end), ...
        on_str, ttl_str, cell_data(ci).m2_peak_pos);

    for si = 1:n_speeds
        nexttile;
        hold on;

        pd_field = sprintf('pd_trace_%s', speed_labels{si});
        nd_field = sprintf('nd_trace_%s', speed_labels{si});
        if ~isfield(cell_data, pd_field) || isempty(cell_data(ci).(pd_field))
            text(0.5, 0.5, 'no data', 'Units', 'normalized', 'HorizontalAlignment', 'center');
            continue;
        end

        spd = speeds_dps(si);
        deg_per_s = spd / SAMPLING_RATE;

        % Build speed-appropriate Gaussian kernel
        stim_dur_samp = round(SWEEP_DUR_SEC(spd) * SAMPLING_RATE);
        stim_dur_s = stim_dur_samp / SAMPLING_RATE;
        fwhm_s = stim_dur_s * 0.25;  % 25% of stim duration
        sigma_s = fwhm_s / 2.355 * SAMPLING_RATE;
        half_w_s = round(3 * sigma_s);
        x_s = (-half_w_s:half_w_s)';
        kern_s = exp(-x_s.^2 / (2*sigma_s^2));
        kern_s = kern_s / sum(kern_s);

        pd_raw = cell_data(ci).(pd_field);
        nd_raw = cell_data(ci).(nd_field);

        % Convolve
        pd_conv_s = conv(pd_raw, kern_s, 'same');
        nd_conv_s = conv(nd_raw, kern_s, 'same');

        % Stimulus window
        stim_st = PADDING_SAMPLES;
        stim_en = stim_st + stim_dur_samp;
        stim_en_pd = min(stim_en, numel(pd_conv_s));
        stim_en_nd = min(stim_en, numel(nd_conv_s));
        stim_range_pd = stim_st:stim_en_pd;
        stim_range_nd = stim_st:stim_en_nd;

        % Time-reverse ND within stimulus window
        nd_stim = nd_conv_s(stim_range_nd);
        nd_flip = flip(nd_stim);

        % x-axis in degrees from sweep start
        t_deg_pd = ((1:numel(pd_conv_s)) - PADDING_SAMPLES) * deg_per_s;
        t_deg_stim = ((stim_range_pd) - PADDING_SAMPLES) * deg_per_s;

        % Plot raw traces (thin, transparent)
        plot(t_deg_pd(1:numel(pd_raw)), pd_raw, '-', 'Color', [0.3 0.3 1 0.2], 'LineWidth', 0.3);
        if numel(nd_raw) <= numel(t_deg_pd)
            t_deg_nd = ((1:numel(nd_raw)) - PADDING_SAMPLES) * deg_per_s;
        else
            t_deg_nd = ((1:numel(nd_raw)) - PADDING_SAMPLES) * deg_per_s;
        end
        plot(t_deg_nd, nd_raw, '-', 'Color', [1 0.3 0.3 0.2], 'LineWidth', 0.3);

        % Plot convolved PD (blue) and ND (red)
        plot(t_deg_pd(1:numel(pd_conv_s)), pd_conv_s, '-', 'Color', [0 0 0.8], 'LineWidth', 1.5);
        plot(t_deg_nd(1:numel(nd_conv_s)), nd_conv_s, '-', 'Color', [0.8 0 0], 'LineWidth', 1.5);

        % Plot flipped-ND in stimulus window (green dashed)
        % nd_flip aligns to the same spatial window as PD stim
        n_flip = numel(nd_flip);
        t_deg_flip = ((stim_st:stim_st+n_flip-1) - PADDING_SAMPLES) * deg_per_s;
        plot(t_deg_flip, nd_flip, '--', 'Color', [0 0.6 0], 'LineWidth', 1.5);

        % Mark peaks on convolved traces
        [pd_pk_val, pd_pk_i] = max(pd_conv_s(stim_range_pd));
        pd_pk_deg = t_deg_stim(pd_pk_i);
        plot(pd_pk_deg, pd_pk_val, 'v', 'MarkerSize', 8, 'Color', [0 0 0.8], ...
            'MarkerFaceColor', [0 0 0.8]);

        [nd_pk_val, nd_pk_i] = max(nd_conv_s(stim_range_nd));
        nd_pk_deg = ((stim_range_nd(nd_pk_i)) - PADDING_SAMPLES) * deg_per_s;
        plot(nd_pk_deg, nd_pk_val, 'v', 'MarkerSize', 8, 'Color', [0.8 0 0], ...
            'MarkerFaceColor', [0.8 0 0]);

        [fnd_pk_val, fnd_pk_i] = max(nd_flip);
        fnd_pk_deg = t_deg_flip(fnd_pk_i);
        plot(fnd_pk_deg, fnd_pk_val, 'v', 'MarkerSize', 8, 'Color', [0 0.6 0], ...
            'MarkerFaceColor', [0 0.6 0]);

        % Time-reversal average peak position (green vertical line)
        avg_pk_deg = (pd_pk_deg + fnd_pk_deg) / 2;
        yl = ylim;
        plot([avg_pk_deg avg_pk_deg], yl, '-', 'Color', [0 0.6 0 0.5], 'LineWidth', 1);

        % PD-ND lag midpoint (gray vertical line)
        mid_deg = (pd_pk_deg + nd_pk_deg) / 2;
        plot([mid_deg mid_deg], yl, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);

        xlim([-5 SWEEP_SPAN_AVG+5]);

        % Speed label on top row
        if ei == 1
            title(sprintf('%s  (FWHM=%.0fms)', speed_labels{si}, fwhm_s*1000), 'FontSize', 9);
        end
        % Cell label on first column
        if si == 1
            ylabel({cell_lbl, 'mV'}, 'FontSize', 7);
        end
    end
end

% Add legend at bottom
lg_ax = nexttile(tl_trdiag, n_examples*3, [1 1]);  % reuse last tile area
% Manual legend using annotation
annotation(fig_trdiag, 'textbox', [0.02 0.01 0.96 0.02], ...
    'String', 'Blue=PD  Red=ND  Green dashed=flipped-ND  |  v=peaks  Green line=TR avg  Gray dotted=PD-ND midpoint', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9);

sgtitle('PD vs ND time-reversal alignment by speed', 'FontSize', 13);
exportgraphics(fig_trdiag, fullfile(preview_dir, 'rf_center_multispeed_pd_nd_alignment.png'), ...
    'Resolution', 200);
fprintf('Saved: rf_center_multispeed_pd_nd_alignment.png\n');

%% Step 3: Summary table
fprintf('\n=== RF Center Comparison — Gauss FWHM=%.2fs (all %d cells) ===\n', ...
    best_fwhm, n_cells);
fprintf('%-20s  %3s  %4s  %4s  %5s  %6s  %5s  %6s  %6s\n', ...
    'Folder', 'ON', 'TTL', 'M2', 'M5', 'Gauss', 'G-M2', 'PD_t', 'ND_t');
fprintf('%s\n', repmat('-', 1, 85));

for ci = 1:n_cells
    on_str  = ifelse(cell_data(ci).is_on, 'ON', 'OFF');
    ttl_str = ifelse(cell_data(ci).is_ttl, 'TTL', 'ctrl');
    fprintf('%-20s  %3s  %4s  %4d  %5.2f  %5.2f  %+5.2f  %6d  %6d\n', ...
        cell_data(ci).folder, on_str, ttl_str, ...
        cell_data(ci).m2_peak_pos, cell_data(ci).m5_centroid, ...
        cell_data(ci).gauss_center, ...
        cell_data(ci).gauss_center - cell_data(ci).m2_peak_pos, ...
        cell_data(ci).gauss_pd_t, cell_data(ci).gauss_nd_t);
end

% Summary statistics
v = valid_mask;
v5 = v & ~isnan(m5_vals);

diff_gauss_m2 = gauss_vals - m2_vals;
diff_gauss_m5 = gauss_vals - m5_vals;

fprintf('\n--- Summary statistics (best width: FWHM=%.2fs) ---\n', best_fwhm);
fprintf('  Gauss vs M2:  mean diff = %+.2f pos, std = %.2f, corr = %.3f\n', ...
    mean(diff_gauss_m2(v)), std(diff_gauss_m2(v)), corr(gauss_vals(v), m2_vals(v)));
fprintf('  Gauss vs M5:  mean diff = %+.2f pos, std = %.2f, corr = %.3f\n', ...
    mean(diff_gauss_m5(v5)), std(diff_gauss_m5(v5)), corr(gauss_vals(v5), m5_vals(v5)));
fprintf('  M2 vs M5:     corr = %.3f\n', corr(m2_vals(v5), m5_vals(v5)));

% Per-width comparison table
fprintf('\n--- All widths vs M2 ---\n');
fprintf('  %8s  %8s  %8s  %8s\n', 'FWHM(s)', 'corr', 'mean_d', 'std_d');
for wi = 1:n_widths
    fprintf('  %8.2f  %+7.3f  %+7.2f  %7.2f', ...
        GAUSS_FWHM_SEC(wi), width_corrs(wi), width_means(wi), width_stds(wi));
    if wi == best_wi, fprintf('  << best'); end
    fprintf('\n');
end

%% Step 4: Comparison figures

% --- Figure 1: Scatter comparison (Gauss vs M2, M5, and width sweep) ---
fig1 = figure('Position', [100 100 1400 500]);

% Panel A: Gauss center vs M2
subplot(1, 3, 1);
hold on;
for ci = 1:n_cells
    if ~valid_mask(ci), continue; end
    clr = get_clr(cell_data(ci));
    plot(m2_vals(ci), gauss_vals(ci), 'o', 'MarkerSize', 7, ...
        'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end
plot([1 11], [1 11], 'k--', 'LineWidth', 0.5);
xlabel('M2 (peak position)');
ylabel(sprintf('Gauss center (FWHM=%.1fs)', best_fwhm));
title(sprintf('Gauss vs M2  (r=%.2f)', corr(gauss_vals(v), m2_vals(v))));
xlim([0.5 11.5]); ylim([0.5 11.5]);
axis square; grid on;

% Panel B: Gauss center vs M5
subplot(1, 3, 2);
hold on;
for ci = 1:n_cells
    if ~v5(ci), continue; end
    clr = get_clr(cell_data(ci));
    plot(m5_vals(ci), gauss_vals(ci), 'o', 'MarkerSize', 7, ...
        'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end
plot([1 11], [1 11], 'k--', 'LineWidth', 0.5);
xlabel('M5 (FWHM centroid)');
ylabel(sprintf('Gauss center (FWHM=%.1fs)', best_fwhm));
title(sprintf('Gauss vs M5  (r=%.2f)', corr(gauss_vals(v5), m5_vals(v5))));
xlim([0.5 11.5]); ylim([0.5 11.5]);
axis square; grid on;

% Panel C: Width sweep — corr with M2 as a function of FWHM
subplot(1, 3, 3);
hold on;
plot(GAUSS_FWHM_SEC, width_corrs, 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
plot(best_fwhm, best_corr, 'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
xlabel('Gaussian FWHM (s)');
ylabel('Correlation with M2');
title('Width sweep: optimal FWHM');
grid on;
ylim([min(width_corrs) - 0.1, min(1, max(width_corrs) + 0.1)]);

% Legend on panel A
subplot(1, 3, 1);
legend_entries = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
legend_colors  = {group_clrs.on_ctrl, group_clrs.on_ttl, group_clrs.off_ctrl, group_clrs.off_ttl};
for k = 1:4
    plot(NaN, NaN, 'o', 'MarkerSize', 7, ...
        'MarkerFaceColor', legend_colors{k}, 'MarkerEdgeColor', 'k', ...
        'DisplayName', legend_entries{k});
end
legend('Location', 'southeast', 'FontSize', 7);

sgtitle(sprintf('RF Center: Gaussian matched-filter (bar sweeps) vs M2/M5 (bar flashes)'));

exportgraphics(fig1, fullfile(preview_dir, 'rf_center_gauss_vs_m2_m5.png'), ...
    'Resolution', 200);
fprintf('\nSaved: rf_center_gauss_vs_m2_m5.png\n');

% --- Figure 2: Difference histograms ---
fig2 = figure('Position', [100 100 900 400]);

subplot(1, 2, 1);
histogram(diff_gauss_m2(v), -5.5:0.5:5.5, 'FaceColor', [0.4 0.4 0.8]);
hold on;
xline(0, 'k--');
xline(mean(diff_gauss_m2(v)), 'r-', 'LineWidth', 1.5);
xlabel('Gauss center - M2 (positions)');
ylabel('Count');
title(sprintf('Gauss vs M2 offset (mean=%+.2f)', mean(diff_gauss_m2(v))));

subplot(1, 2, 2);
histogram(diff_gauss_m5(v5), -5.5:0.5:5.5, 'FaceColor', [0.4 0.8 0.4]);
hold on;
xline(0, 'k--');
xline(mean(diff_gauss_m5(v5)), 'r-', 'LineWidth', 1.5);
xlabel('Gauss center - M5 (positions)');
ylabel('Count');
title(sprintf('Gauss vs M5 offset (mean=%+.2f)', mean(diff_gauss_m5(v5))));

sgtitle(sprintf('RF Center Differences (Gauss FWHM=%.2fs)', best_fwhm));

exportgraphics(fig2, fullfile(preview_dir, 'rf_center_gauss_diff_histograms.png'), ...
    'Resolution', 200);
fprintf('Saved: rf_center_gauss_diff_histograms.png\n');

% --- Figure 3: Per-cell diagnostics — PD/ND traces with Gaussian convolution ---
%  Show raw traces + convolution output + peak markers
best_kernel = gauss_kernels{best_wi};
n_diag = min(n_cells, 12);
fig3 = figure('Position', [50 50 1400 900]);
tl3 = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for ci = 1:n_diag
    ax = nexttile;
    hold(ax, 'on');

    pd_tr = cell_data(ci).pd_trace;
    nd_tr = cell_data(ci).nd_trace;
    if isempty(pd_tr) || isempty(nd_tr), continue; end

    min_len = min(numel(pd_tr), numel(nd_tr));
    t_ms = (1:min_len) / (SAMPLING_RATE / 1000);

    % Raw traces (light)
    plot(ax, t_ms, pd_tr(1:min_len), 'b', 'LineWidth', 0.4, 'Color', [0.6 0.6 1]);
    plot(ax, t_ms, nd_tr(1:min_len), 'r', 'LineWidth', 0.4, 'Color', [1 0.6 0.6]);

    % Gaussian-convolved traces (bold)
    pd_conv = conv(pd_tr(1:min_len), best_kernel, 'same');
    nd_conv = conv(nd_tr(1:min_len), best_kernel, 'same');
    plot(ax, t_ms, pd_conv, 'b', 'LineWidth', 1.5);
    plot(ax, t_ms, nd_conv, 'r', 'LineWidth', 1.5);

    % Peak markers
    pd_t = cell_data(ci).gauss_pd_t;
    nd_t = cell_data(ci).gauss_nd_t;
    if ~isnan(pd_t) && pd_t <= min_len
        plot(ax, pd_t / 10, pd_conv(pd_t), 'bv', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    end
    if ~isnan(nd_t) && nd_t <= min_len
        plot(ax, nd_t / 10, nd_conv(nd_t), 'rv', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end

    % Stimulus window
    stim_start_ms = STIM_START / 10;
    stim_end_ms   = (min_len - STIM_TRIM_END) / 10;
    xline(ax, stim_start_ms, 'g--', 'LineWidth', 0.5);
    xline(ax, stim_end_ms, 'g--', 'LineWidth', 0.5);

    on_str  = ifelse(cell_data(ci).is_on, 'ON', 'OFF');
    ttl_str = ifelse(cell_data(ci).is_ttl, 'TTL', 'ctrl');

    title(ax, sprintf('%s (%s %s)\nG=%.1f M2=%d M5=%.1f', ...
        cell_data(ci).folder, on_str, ttl_str, ...
        cell_data(ci).gauss_center, cell_data(ci).m2_peak_pos, ...
        cell_data(ci).m5_centroid), 'FontSize', 7, 'Interpreter', 'none');

    if ci == 1
        legend(ax, {'PD raw', 'ND raw', 'PD conv', 'ND conv'}, ...
            'FontSize', 5, 'Location', 'northeast');
    end
    ylabel(ax, 'mV');
    if ci > 8, xlabel(ax, 'ms'); end
    set(ax, 'FontSize', 7);
end

title(tl3, sprintf('Gaussian matched-filter (FWHM=%.2fs): PD vs ND peak detection', best_fwhm));

exportgraphics(fig3, fullfile(preview_dir, 'rf_center_gauss_diagnostics.png'), ...
    'Resolution', 200);
fprintf('Saved: rf_center_gauss_diagnostics.png\n');

% --- Figure 4: Remaining cells ---
if n_cells > 12
    n_remaining = n_cells - 12;
    fig4 = figure('Position', [50 50 1400 900]);
    tl4 = tiledlayout(3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

    for k = 1:min(n_remaining, 15)
        ci = 12 + k;
        ax = nexttile;
        hold(ax, 'on');

        pd_tr = cell_data(ci).pd_trace;
        nd_tr = cell_data(ci).nd_trace;
        if isempty(pd_tr) || isempty(nd_tr), continue; end

        min_len = min(numel(pd_tr), numel(nd_tr));
        t_ms = (1:min_len) / (SAMPLING_RATE / 1000);

        plot(ax, t_ms, pd_tr(1:min_len), 'b', 'LineWidth', 0.4, 'Color', [0.6 0.6 1]);
        plot(ax, t_ms, nd_tr(1:min_len), 'r', 'LineWidth', 0.4, 'Color', [1 0.6 0.6]);

        pd_conv = conv(pd_tr(1:min_len), best_kernel, 'same');
        nd_conv = conv(nd_tr(1:min_len), best_kernel, 'same');
        plot(ax, t_ms, pd_conv, 'b', 'LineWidth', 1.5);
        plot(ax, t_ms, nd_conv, 'r', 'LineWidth', 1.5);

        pd_t = cell_data(ci).gauss_pd_t;
        nd_t = cell_data(ci).gauss_nd_t;
        if ~isnan(pd_t) && pd_t <= min_len
            plot(ax, pd_t / 10, pd_conv(pd_t), 'bv', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        end
        if ~isnan(nd_t) && nd_t <= min_len
            plot(ax, nd_t / 10, nd_conv(nd_t), 'rv', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        end

        stim_start_ms = STIM_START / 10;
        stim_end_ms   = (min_len - STIM_TRIM_END) / 10;
        xline(ax, stim_start_ms, 'g--', 'LineWidth', 0.5);
        xline(ax, stim_end_ms, 'g--', 'LineWidth', 0.5);

        on_str  = ifelse(cell_data(ci).is_on, 'ON', 'OFF');
        ttl_str = ifelse(cell_data(ci).is_ttl, 'TTL', 'ctrl');

        title(ax, sprintf('%s (%s %s)\nG=%.1f M2=%d M5=%.1f', ...
            cell_data(ci).folder, on_str, ttl_str, ...
            cell_data(ci).gauss_center, cell_data(ci).m2_peak_pos, ...
            cell_data(ci).m5_centroid), 'FontSize', 7, 'Interpreter', 'none');
        set(ax, 'FontSize', 7);
    end

    title(tl4, sprintf('Gaussian matched-filter (FWHM=%.2fs): continued', best_fwhm));

    exportgraphics(fig4, fullfile(preview_dir, 'rf_center_gauss_diagnostics_2.png'), ...
        'Resolution', 200);
    fprintf('Saved: rf_center_gauss_diagnostics_2.png\n');
end

fprintf('\n=== Done ===\n');


%% ===== LOCAL FUNCTIONS =====

function entry = find_batch_entry(results, folder_name)
% FIND_BATCH_ENTRY  Find batch results entry matching a folder name.
    entry = [];
    for k = 1:numel(results)
        r = results(k);
        % Match by folder name (late dataset) or folder path (early)
        if isfield(r, 'folder')
            if strcmp(r.folder, folder_name) || endsWith(r.folder, folder_name)
                entry = r;
                return;
            end
        end
    end
end


function [traces_16, pd_shift, lut_dirs_ordered] = extract_and_align_traces( ...
    bar_data, plot_order, lut_directions, pd_direction, row_offset)
% EXTRACT_AND_ALIGN_TRACES  Extract 16 mean traces and PD-align.
%   Position 5 = PD (90 deg in aligned frame).
%   row_offset: optional offset added to plot_order for bar_data row access
%               (used for 56 dps = +16, 168 dps = +32). Default = 0.

    if nargin < 5, row_offset = 0; end
    n_reps = size(bar_data, 2) - 1;

    % Extract mean traces in subplot order
    traces_subplot_order = cell(16, 1);
    for si = 1:16
        data_row = plot_order(si) + row_offset;
        traces_subplot_order{si} = bar_data{data_row, n_reps + 1};
    end

    % Sort by direction angle (always use original plot_order for LUT)
    lut_dirs_subplot = lut_directions(plot_order);
    [sorted_angles, sort_idx] = sort(lut_dirs_subplot(:));
    traces_sorted = traces_subplot_order(sort_idx);

    % Find PD position and shift to slot 5
    angle_diffs = abs(mod(sorted_angles - pd_direction + 180, 360) - 180);
    [~, pd_sorted_idx] = min(angle_diffs);
    pd_shift = 5 - pd_sorted_idx;
    traces_16 = circshift(traces_sorted, pd_shift);
    lut_dirs_ordered = circshift(sorted_angles, pd_shift);
end


function s = ifelse(cond, a, b)
% IFELSE  Inline conditional for string selection.
    if cond, s = a; else, s = b; end
end

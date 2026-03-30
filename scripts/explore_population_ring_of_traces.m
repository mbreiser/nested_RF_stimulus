% EXPLORE_POPULATION_RING_OF_TRACES  Phase 0A: Population ring-of-traces
%   (ctrl vs TTL, xcorr-refined alignment, mean±SEM).
%
%   Re-loads raw bar sweep data from all 48 experiments (25 late + 23 early),
%   extracts 16 mean bar sweep traces per cell, PD-aligns by circular row
%   shift, applies dark-bar polarity correction, baseline-subtracts, then
%   temporally aligns using peak detection + xcorr refinement.
%
%   Outputs:
%     1. Population ring-of-traces PNGs (2): ON + OFF, xcorr-refined mean±SEM
%     2. Representative pair PNGs (2): ON + OFF, selected single cells
%     3. Candidate PNGs (~20): 5 per group for representative selection
%     4. Alignment diagnostic figure (1): per-cell traces color-coded
%     5. xcorr corrections diagnostic table
%     6. time_to_max data saved for RF position estimation
%
%   Usage:
%     run('scripts/explore_population_ring_of_traces.m')

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root    = '/Users/reiserm/Documents/ttl_1DRF';
preview_dir  = fullfile(data_root, 'figure_previews');
if ~isfolder(preview_dir), mkdir(preview_dir); end

plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

% Alignment parameters
SMOOTH_WIN     = 25;       % boxcar filter width for peak detection
BL_START       = 1000;     % baseline window start (samples)
BL_END         = 9000;     % baseline window end (samples)
STIM_START     = 9000;     % stimulus window start
STIM_TRIM_END  = 7000;     % samples to trim from end of stim
DISPLAY_HALF   = 14000;    % ±samples around peak for display (±1400ms)
DS_FACTOR      = 10;       % downsample factor for display

% LUT
lut_path = fullfile(fileparts(mfilename('fullpath')), ...
    '..', 'src', 'analysis', 'protocol2', 'bar_lut.mat');
S_lut = load(lut_path, 'Tbl');
Tbl = S_lut.Tbl;

% Load batch results for validated PD directions and polar tuning data
fprintf('Loading batch results for validated PD directions...\n');
S_late_batch = load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
late_batch = S_late_batch.results;
S_early_batch = load(fullfile(data_root, 'pre-bar-flash', 'population_results', ...
    'batch_results_pre_bf.mat'), 'results');
early_batch = S_early_batch.results;
fprintf('  Late batch: %d cells, Early batch: %d cells\n', ...
    numel(late_batch), numel(early_batch));

%% Step 1: Collect all experiment paths and metadata
fprintf('=== Phase 0A: Population Ring-of-Traces Exploration ===\n\n');

% --- Late dataset (25 cells) ---
late_root = data_root;
d_late = dir(late_root);
d_late = d_late([d_late.isdir]);
d_late = d_late(~startsWith({d_late.name}, '.'));
late_valid = false(numel(d_late), 1);
for i = 1:numel(d_late)
    late_valid(i) = isfile(fullfile(late_root, d_late(i).name, 'currentExp.mat'));
end
d_late = d_late(late_valid);

% --- Early dataset (23 cells) ---
early_root = fullfile(data_root, 'pre-bar-flash');
early_list = discover_early_experiments(early_root);

n_late  = numel(d_late);
n_early = numel(early_list);
n_total = n_late + n_early;
fprintf('Found %d late + %d early = %d total experiments\n\n', ...
    n_late, n_early, n_total);

%% Step 2: Extract raw bar sweep traces for all cells
%   For each cell, extract:
%     - 16 mean bar sweep traces (one per direction)
%     - PD-aligned ordering (circular shift so PD → position 5 = 90°)
%     - Metadata: is_on, is_ttl, date_str

all_cells = struct([]);
cell_idx = 0;

% --- Process late dataset ---
fprintf('--- Loading late dataset (bar sweeps from 3-speed protocol) ---\n');
for exp_idx = 1:n_late
    folder_name = d_late(exp_idx).name;
    exp_folder = fullfile(late_root, folder_name);
    fprintf('[%d/%d] %s ...', exp_idx, n_late, folder_name);

    try
        orig_dir = pwd;
        cleanup = onCleanup(@() cd(orig_dir));

        [date_str, ~, Log, ~, ~] = load_protocol2_data(exp_folder);
        f_data = Log.ADC.Volts(1, :);
        v_data = Log.ADC.Volts(2, :) * 10;

        ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
            'pattern_order', 'func_order', 'metadata');

        % Parse bar sweeps (3-speed protocol → 32×4, rows 1-16 = slow)
        bar_data = parse_bar_data(f_data, v_data);

        % Classify
        is_on  = ce.metadata.Frame > 129;
        is_ttl = contains(ce.metadata.Strain, 'ttl');

        % LUT directions
        [lut_directions, ~, ~, ~] = ...
            verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, plot_order);

        % Look up validated PD direction from batch results
        batch_entry = find_batch_entry(late_batch, folder_name);
        if isempty(batch_entry)
            fprintf(' SKIP: no batch match\n');
            continue;
        end

        % PD-align using validated batch PD direction
        [traces_16, pd_shift, lut_dirs_ordered] = extract_and_align_traces( ...
            bar_data, plot_order, lut_directions, batch_entry.pd_direction);

        cell_idx = cell_idx + 1;
        all_cells(cell_idx).traces_aligned = traces_16;  % 16 cell array of traces
        all_cells(cell_idx).pd_shift       = pd_shift;
        all_cells(cell_idx).is_on          = is_on;
        all_cells(cell_idx).is_ttl         = is_ttl;
        all_cells(cell_idx).date_str       = [date_str '_' folder_name(end-4:end)];
        all_cells(cell_idx).batch          = 'late';
        all_cells(cell_idx).lut_dirs       = lut_dirs_ordered;
        all_cells(cell_idx).max_v_aligned  = batch_entry.max_v_aligned;
        all_cells(cell_idx).dsi            = batch_entry.dsi_vector;

        group = classify_group(is_on, is_ttl);
        fprintf(' OK [%s] (batch PD=%.1f°)\n', group, batch_entry.pd_direction);

    catch ME
        fprintf(' ERROR: %s\n', ME.message);
        continue;
    end
end

% --- Process early dataset ---
fprintf('\n--- Loading early dataset (bar sweeps from 2-speed protocol) ---\n');
for exp_idx = 1:n_early
    ei = early_list(exp_idx);
    fprintf('[%d/%d] %s (%s/%s) ...', exp_idx, n_early, ...
        ei.date_str, ei.treatment, ei.cell_type);

    try
        orig_dir = pwd;
        cleanup = onCleanup(@() cd(orig_dir));

        [date_str, ~, Log, ~, ~] = load_protocol2_data(ei.folder);
        f_data = Log.ADC.Volts(1, :);
        v_data = Log.ADC.Volts(2, :) * 10;

        ce = load(fullfile(ei.folder, 'currentExp.mat'), ...
            'pattern_order', 'func_order', 'metadata');

        % Parse bar sweeps (2-speed protocol)
        bar_data = parse_bar_data_pre_bf(f_data, v_data);

        % Classify by directory path (not metadata)
        is_on  = strcmpi(ei.cell_type, 'ON');
        is_ttl = strcmpi(ei.treatment, 'ttl');
        is_off = ~is_on;

        % Apply dark-bar polarity correction for pre-Oct-15 OFF cells
        [bar_data, was_corrected] = correct_off_polarity_swap( ...
            bar_data, ce.pattern_order, ce.func_order, ei.date_str, is_off);

        % LUT directions
        [lut_directions, ~, ~, ~] = ...
            verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, plot_order);

        % Look up validated PD direction from early batch results
        batch_entry = find_batch_entry(early_batch, ei.folder);
        if isempty(batch_entry)
            fprintf(' SKIP: no batch match\n');
            continue;
        end

        % PD-align using validated batch PD direction
        [traces_16, pd_shift, lut_dirs_ordered] = extract_and_align_traces( ...
            bar_data, plot_order, lut_directions, batch_entry.pd_direction);

        cell_idx = cell_idx + 1;
        all_cells(cell_idx).traces_aligned = traces_16;
        all_cells(cell_idx).pd_shift       = pd_shift;
        all_cells(cell_idx).is_on          = is_on;
        all_cells(cell_idx).is_ttl         = is_ttl;
        all_cells(cell_idx).date_str       = ei.date_str;
        all_cells(cell_idx).batch          = 'early';
        all_cells(cell_idx).lut_dirs       = lut_dirs_ordered;
        all_cells(cell_idx).max_v_aligned  = batch_entry.max_v_aligned;
        all_cells(cell_idx).dsi            = batch_entry.dsi_vector;

        corr_tag = '';
        if was_corrected, corr_tag = ' *corrected*'; end
        group = classify_group(is_on, is_ttl);
        fprintf(' OK [%s]%s (batch PD=%.1f°)\n', group, corr_tag, batch_entry.pd_direction);

    catch ME
        fprintf(' ERROR: %s\n', ME.message);
        continue;
    end
end

n_cells = numel(all_cells);
fprintf('\nSuccessfully loaded %d / %d cells\n\n', n_cells, n_total);

%% Step 3: Baseline-subtract all traces
%   Baseline = mean of samples BL_START:BL_END (same as compute_bar_sweep_responses)

fprintf('Baseline-subtracting traces...\n');
for ci = 1:n_cells
    for di = 1:16
        tr = all_cells(ci).traces_aligned{di};
        if isempty(tr), continue; end
        bl = mean(tr(BL_START:min(BL_END, numel(tr))));
        all_cells(ci).traces_aligned{di} = tr - bl;
    end
end

%% Step 3b: Recompute 99.5th percentile peak amplitudes for polar plots
%   Uses same baseline (BL_START:BL_END) but 99.5th percentile (vs 98th
%   in compute_bar_sweep_responses) for consistency with rest of pipeline.
%   Stored in same format as max_v_aligned: 16×2 [angle, peak_amp].

fprintf('Recomputing 99.5th percentile peak amplitudes for polar plots...\n');
POLAR_PERCENTILE = 99.5;

for ci = 1:n_cells
    % Get PD-aligned angles from batch results
    angles = all_cells(ci).max_v_aligned(:, 1);

    peak_amps = NaN(16, 1);
    for di = 1:16
        tr = all_cells(ci).traces_aligned{di};  % already baseline-subtracted
        if isempty(tr), continue; end

        stim_end = max(1, numel(tr) - STIM_TRIM_END);
        if STIM_START > stim_end, continue; end

        d_stim = tr(STIM_START:stim_end);
        peak_amps(di) = prctile(d_stim, POLAR_PERCENTILE);
    end

    all_cells(ci).peak_amps_995 = [angles, peak_amps];
end
fprintf('  Done — stored as all_cells(ci).peak_amps_995\n');

%% Step 4: Temporal alignment via peak detection
%   For each cell × direction:
%     1. Smooth mean trace with 25-sample boxcar
%     2. Find max in stimulus window on smoothed trace
%     3. If max < 1 SD of baseline → flag as "no valid alignment"
%     4. Shift trace so max_time → common reference sample

fprintf('Temporally aligning traces (peak detection with %d-sample smooth)...\n', ...
    SMOOTH_WIN);

% Reference sample: where aligned peaks will map to
REF_SAMPLE = DISPLAY_HALF + 1;  % center of display window
DISPLAY_LEN = 2 * DISPLAY_HALF + 1;

% Storage for alignment data
time_to_max   = NaN(n_cells, 16);  % max_time in original trace
align_valid   = false(n_cells, 16);
aligned_traces = cell(n_cells, 16);  % NaN-padded, display-windowed traces

for ci = 1:n_cells
    % First pass: find valid alignment times for all 16 directions
    for di = 1:16
        tr = all_cells(ci).traces_aligned{di};
        if isempty(tr), continue; end

        % Baseline stats from raw trace
        bl_region = tr(BL_START:min(BL_END, numel(tr)));
        bl_std = std(bl_region);

        % Smooth for peak detection
        tr_smooth = movmean(tr, SMOOTH_WIN);

        % Stimulus window
        stim_end = max(1, numel(tr_smooth) - STIM_TRIM_END);
        stim_win = tr_smooth(STIM_START:stim_end);

        if isempty(stim_win), continue; end

        [max_val, max_idx] = max(stim_win);
        max_time = STIM_START + max_idx - 1;

        time_to_max(ci, di) = max_time;

        % Threshold: max must exceed 1 SD above baseline (which is ~0 after subtraction)
        if max_val > bl_std
            align_valid(ci, di) = true;
        end
    end

    % Second pass: borrow shifts for invalid directions from nearest valid neighbor
    for di = 1:16
        if align_valid(ci, di)
            shift_time = time_to_max(ci, di);
        else
            % Find nearest valid neighbor (circular)
            shift_time = borrow_nearest_shift(time_to_max(ci, :), align_valid(ci, :), di);
            if isnan(shift_time)
                % No valid neighbors at all — skip this cell×direction
                aligned_traces{ci, di} = NaN(DISPLAY_LEN, 1);
                continue;
            end
        end

        tr = all_cells(ci).traces_aligned{di};
        if isempty(tr)
            aligned_traces{ci, di} = NaN(DISPLAY_LEN, 1);
            continue;
        end

        % Shift trace so shift_time maps to REF_SAMPLE
        aligned_traces{ci, di} = shift_and_window(tr, shift_time, ...
            REF_SAMPLE, DISPLAY_LEN);
    end
end

% Report alignment statistics
n_valid = sum(align_valid(:));
n_total_dirs = n_cells * 16;
fprintf('  Valid alignments: %d / %d (%.1f%%)\n', ...
    n_valid, n_total_dirs, 100*n_valid/n_total_dirs);
fprintf('  Flagged (borrowed shift): %d\n', sum(~align_valid(:)));

%% Step 5: Downsample for display
fprintf('Downsampling by factor %d for display...\n', DS_FACTOR);

aligned_traces_ds = cell(size(aligned_traces));
for ci = 1:n_cells
    for di = 1:16
        tr = aligned_traces{ci, di};
        if isempty(tr) || all(isnan(tr))
            aligned_traces_ds{ci, di} = downsample(tr, DS_FACTOR);
        else
            aligned_traces_ds{ci, di} = downsample(tr, DS_FACTOR);
        end
    end
end

display_len_ds = ceil(DISPLAY_LEN / DS_FACTOR);
ref_sample_ds  = ceil(REF_SAMPLE / DS_FACTOR);

%% Step 6: Group cells
on_ctrl_idx  = find([all_cells.is_on] & ~[all_cells.is_ttl]);
on_ttl_idx   = find([all_cells.is_on] &  [all_cells.is_ttl]);
off_ctrl_idx = find(~[all_cells.is_on] & ~[all_cells.is_ttl]);
off_ttl_idx  = find(~[all_cells.is_on] &  [all_cells.is_ttl]);

fprintf('\nGroup sizes:\n');
fprintf('  ON  ctrl: %d, TTL: %d\n', numel(on_ctrl_idx), numel(on_ttl_idx));
fprintf('  OFF ctrl: %d, TTL: %d\n', numel(off_ctrl_idx), numel(off_ttl_idx));

%% Step 7: Compute population statistics per direction per group
%   Computes mean±SEM for timeseries and polar data.
fprintf('\nComputing population averages (mean±SEM)...\n');

% PD-aligned angles (16 directions, PD at position 5 = 90°)
pd_aligned_angles = (0:15)' * 22.5;  % 0, 22.5, ..., 337.5

groups = struct( ...
    'name',     {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'}, ...
    'indices',  {on_ctrl_idx, on_ttl_idx, off_ctrl_idx, off_ttl_idx}, ...
    'color',    {[0 0 0], [1 0 0], [0 0 0], [1 0 0]}, ...
    'linestyle',{'-', '-', '-', '-'} ...
);

% Compute per-direction statistics for each group
for g = 1:4
    idx = groups(g).indices;
    n_g = numel(idx);
    groups(g).n = n_g;

    % Initialize storage for timeseries stats
    groups(g).mean_traces   = cell(16, 1);
    groups(g).sem_traces    = cell(16, 1);

    % Polar: collect 99.5th percentile peak amplitudes
    polar_amps = NaN(16, n_g);
    for k = 1:n_g
        pa = all_cells(idx(k)).peak_amps_995;
        if ~isempty(pa) && size(pa, 2) == 2
            polar_amps(:, k) = pa(:, 2);
        end
    end

    % Polar statistics (for central polar plot)
    groups(g).polar_mean   = mean(polar_amps, 2, 'omitnan');
    groups(g).polar_sem    = std(polar_amps, 0, 2, 'omitnan') / sqrt(n_g);

    for di = 1:16
        % Stack all cells for this direction
        traces_mat = NaN(display_len_ds, n_g);
        for k = 1:n_g
            tr = aligned_traces_ds{idx(k), di};
            len = min(numel(tr), display_len_ds);
            traces_mat(1:len, k) = tr(1:len);
        end

        % Mean ± SEM
        groups(g).mean_traces{di}   = mean(traces_mat, 2, 'omitnan');
        groups(g).sem_traces{di}    = std(traces_mat, 0, 2, 'omitnan') / sqrt(n_g);
    end
end

%% Step 7b: Cross-correlation refined temporal alignment (population only)
%   Two-pass alignment:
%     Pass 1 (done above): Peak-based alignment — all traces centered on peak sample
%     Pass 2 (here):       xcorr refinement — use peak-aligned group mean as template,
%                          then refine each cell's shift by maximizing cross-correlation
%                          with the template. Search window: ±6000 samples (±600ms).
%   This reduces the "pinching" artifact where peak alignment forces all traces
%   to converge at a single sample.

fprintf('\n--- Cross-correlation refinement of temporal alignment ---\n');

XCORR_MAX_LAG  = 6000;   % ±600ms search window (at 10 kHz)
XCORR_MAX_LAG_DS = ceil(XCORR_MAX_LAG / DS_FACTOR);  % in downsampled units

% Storage for xcorr-refined traces (same size as aligned_traces_ds)
xcorr_aligned_traces_ds = aligned_traces_ds;  % start from peak-aligned

% Diagnostic table: cell_date | group | dir_idx | dir_angle | peak_shift |
%                   xcorr_correction | total_shift | xcorr_max_r
diag_rows = {};

for g = 1:4
    idx = groups(g).indices;
    n_g = numel(idx);
    if n_g < 2, continue; end  % need at least 2 cells to make a template

    for di = 1:16
        % Step 1: Compute group mean template from peak-aligned traces
        traces_mat = NaN(display_len_ds, n_g);
        for k = 1:n_g
            tr = aligned_traces_ds{idx(k), di};
            len = min(numel(tr), display_len_ds);
            traces_mat(1:len, k) = tr(1:len);
        end
        template = mean(traces_mat, 2, 'omitnan');

        % Replace NaN in template with 0 for xcorr
        template_clean = template;
        template_clean(isnan(template_clean)) = 0;

        if all(template_clean == 0), continue; end

        % Step 2: For each cell, refine alignment via xcorr
        for k = 1:n_g
            ci = idx(k);
            tr = aligned_traces_ds{ci, di};
            if isempty(tr) || all(isnan(tr)), continue; end

            % Replace NaN with 0 for xcorr computation
            tr_clean = tr(:);
            tr_clean(isnan(tr_clean)) = 0;

            if all(tr_clean == 0), continue; end

            % Cross-correlate cell trace with group mean template
            [r, lags] = xcorr(tr_clean, template_clean, XCORR_MAX_LAG_DS);

            % Find lag at maximum correlation
            [max_r, max_r_idx] = max(r);
            best_lag = lags(max_r_idx);  % positive = cell leads template

            % Apply correction: shift trace by -best_lag (circshift in sample space)
            corrected = NaN(display_len_ds, 1);
            for out_idx = 1:display_len_ds
                src_idx = out_idx + best_lag;
                if src_idx >= 1 && src_idx <= numel(tr)
                    corrected(out_idx) = tr(src_idx);
                end
            end
            xcorr_aligned_traces_ds{ci, di} = corrected;

            % Log diagnostic info
            % Convert correction from DS samples to original samples
            correction_samples = best_lag * DS_FACTOR;
            % Original peak shift (from Step 4)
            peak_shift_orig = time_to_max(ci, di) - REF_SAMPLE;
            total_shift = peak_shift_orig + correction_samples;

            diag_rows{end+1} = sprintf('%s\t%s\t%d\t%.1f\t%d\t%d\t%d\t%.4f', ...
                all_cells(ci).date_str, groups(g).name, di, ...
                pd_aligned_angles(di), round(peak_shift_orig), ...
                correction_samples, round(total_shift), max_r); %#ok<AGROW>
        end
    end
end

% Write diagnostic table
diag_file = fullfile(preview_dir, 'diag_xcorr_corrections.txt');
fid = fopen(diag_file, 'w');
fprintf(fid, 'cell_date\tgroup\tdir_idx\tdir_angle\tpeak_shift_samples\txcorr_correction_samples\ttotal_shift_samples\txcorr_max_r\n');
for r = 1:numel(diag_rows)
    fprintf(fid, '%s\n', diag_rows{r});
end
fclose(fid);
fprintf('  Diagnostic table saved: %s\n', diag_file);

% Flag large corrections (> ±200ms = ±2000 samples)
n_large = 0;
for r = 1:numel(diag_rows)
    parts = strsplit(diag_rows{r}, '\t');
    corr_val = str2double(parts{6});
    if abs(corr_val) > 2000
        n_large = n_large + 1;
    end
end
fprintf('  Corrections > ±200ms: %d / %d entries\n', n_large, numel(diag_rows));

% Compute xcorr-refined population statistics (same structure as groups)
fprintf('  Computing xcorr-refined population statistics...\n');
groups_xcorr = groups;  % copy structure, overwrite traces

for g = 1:4
    idx = groups_xcorr(g).indices;
    n_g = numel(idx);

    % Re-initialize trace storage
    groups_xcorr(g).mean_traces   = cell(16, 1);
    groups_xcorr(g).sem_traces    = cell(16, 1);

    % Polar stats are unchanged (not affected by temporal alignment)
    % groups_xcorr(g).polar_* remain the same as groups(g).polar_*

    for di = 1:16
        traces_mat = NaN(display_len_ds, n_g);
        for k = 1:n_g
            tr = xcorr_aligned_traces_ds{idx(k), di};
            len = min(numel(tr), display_len_ds);
            traces_mat(1:len, k) = tr(1:len);
        end

        % Mean ± SEM
        groups_xcorr(g).mean_traces{di}   = mean(traces_mat, 2, 'omitnan');
        groups_xcorr(g).sem_traces{di}    = std(traces_mat, 0, 2, 'omitnan') / sqrt(n_g);
    end
end

fprintf('  xcorr refinement complete.\n');

%% Step 8: Generate population ring-of-traces figures (xcorr-refined, mean±SEM)
fprintf('\n--- Generating population ring-of-traces ---\n');

% Shared layout parameters (matching plot_slow_bar_sweep_polar conventions)
RING_LAYOUT = struct( ...
    'fig_size', [50 50 1100 1100], ...
    'centerX', 0.50, 'centerY', 0.50, ...
    'radius', 0.32, ...
    'subW', 0.11, 'subH', 0.13, ...
    'polar_scale', 0.50, ...  % central polar size = radius * polar_scale * 2
    'font_label', 7, ...
    'font_title', 11);

% --- Population figures: xcorr-refined, mean±SEM only (2 PNGs: ON + OFF) ---
for on_off = ["ON", "OFF"]
    if on_off == "ON"
        gx_ctrl = groups_xcorr(1);
        gx_ttl  = groups_xcorr(2);
    else
        gx_ctrl = groups_xcorr(3);
        gx_ttl  = groups_xcorr(4);
    end

    fig = plot_ring_of_traces_population(gx_ctrl, gx_ttl, on_off, ...
        pd_aligned_angles, all_cells, ...
        display_len_ds, ref_sample_ds, DS_FACTOR, ...
        xcorr_aligned_traces_ds, RING_LAYOUT);

    fname = sprintf('pop_ring_%s_ctrl_vs_ttl.png', lower(char(on_off)));
    exportgraphics(fig, fullfile(preview_dir, fname), 'Resolution', 200);
    fprintf('  Saved: %s\n', fname);
    close(fig);
end

%% Step 8b: Generate single-cell representative comparison figures
fprintf('\n--- Generating single-cell representative comparisons ---\n');

% Representative cells (selected from candidate review)
rep_cells = struct( ...
    'ON_ctrl',  '2025_11_13_11_03', ...   % rank 1 candidate
    'ON_ttl',   '2025_08_27_12_32', ...   % rank 3 candidate
    'OFF_ctrl', '2025_07_17_15_05', ...   % rank 1 candidate
    'OFF_ttl',  '2025_10_30_10_47');      % rank 3 candidate

for on_off = ["ON", "OFF"]
    if on_off == "ON"
        ctrl_date = rep_cells.ON_ctrl;
        ttl_date  = rep_cells.ON_ttl;
    else
        ctrl_date = rep_cells.OFF_ctrl;
        ttl_date  = rep_cells.OFF_ttl;
    end

    % Find these cells in all_cells
    ctrl_ci = find_cell_by_date(all_cells, ctrl_date);
    ttl_ci  = find_cell_by_date(all_cells, ttl_date);

    if isempty(ctrl_ci) || isempty(ttl_ci)
        fprintf('  WARNING: Could not find representative cells for %s\n', on_off);
        continue;
    end

    fig = plot_ring_of_traces_single_pair(all_cells, ctrl_ci, ttl_ci, on_off, ...
        pd_aligned_angles, xcorr_aligned_traces_ds, ...
        display_len_ds, ref_sample_ds, DS_FACTOR, RING_LAYOUT);

    fname = sprintf('rep_ring_%s_ctrl_vs_ttl.png', lower(char(on_off)));
    exportgraphics(fig, fullfile(preview_dir, fname), 'Resolution', 200);
    fprintf('  Saved: %s\n', fname);
    close(fig);
end

%% Step 8c: Generate representative candidate figures (5 per group)
fprintf('\n--- Generating representative candidate figures ---\n');

cand_base_dir = fullfile(preview_dir, 'rep_candidates');
if ~isfolder(cand_base_dir), mkdir(cand_base_dir); end

group_tags = {'ON_ctrl', 'ON_ttl', 'OFF_ctrl', 'OFF_ttl'};
group_indices = {on_ctrl_idx, on_ttl_idx, off_ctrl_idx, off_ttl_idx};
N_CANDIDATES = 5;

for g = 1:4
    idx = group_indices{g};
    n_g = numel(idx);
    if n_g == 0, continue; end

    % Create subdirectory
    subdir = fullfile(cand_base_dir, group_tags{g});
    if ~isfolder(subdir), mkdir(subdir); end

    % Find top N candidates closest to group mean DSI
    cand_idx = find_top_n_by_dsi(idx, all_cells, min(N_CANDIDATES, n_g));

    for rank = 1:numel(cand_idx)
        ci = cand_idx(rank);
        fig = plot_ring_of_traces_single_cell(all_cells, ci, ...
            pd_aligned_angles, xcorr_aligned_traces_ds, ...
            display_len_ds, ref_sample_ds, DS_FACTOR, RING_LAYOUT);

        fname = sprintf('rep_cand_%s_%d_%s.png', ...
            group_tags{g}, rank, all_cells(ci).date_str);
        exportgraphics(fig, fullfile(subdir, fname), 'Resolution', 200);
        fprintf('  Saved: %s/%s\n', group_tags{g}, fname);
        close(fig);
    end
end

%% Step 9: Alignment diagnostic figure
%   4 rows (ON ctrl, ON TTL, OFF ctrl, OFF TTL) × 16 columns (PD-aligned dirs)
%   Each panel: overlay individual cell traces (thin lines), color-coded per cell

fprintf('\n--- Generating alignment diagnostic ---\n');

% Time axis for diagnostic (ms relative to aligned peak)
t_ms = ((1:display_len_ds) - ref_sample_ds) * DS_FACTOR * 0.1;

fig_diag = figure('Position', [50 50 2000 800], 'Color', 'w');
t_diag = tiledlayout(4, 16, 'TileSpacing', 'compact', 'Padding', 'compact');

group_order = {on_ctrl_idx, on_ttl_idx, off_ctrl_idx, off_ttl_idx};
group_names = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};

for g = 1:4
    idx = group_order{g};
    n_g = numel(idx);

    % Generate distinct colors for each cell in this group
    if n_g > 0
        cell_colors = lines(n_g);
    else
        cell_colors = [];
    end

    for di = 1:16
        ax = nexttile(t_diag);
        hold(ax, 'on');

        if n_g == 0
            set(ax, 'XTick', [], 'YTick', []);
            if di == 1
                ylabel(ax, group_names{g}, 'FontSize', 7);
            end
            continue;
        end

        for k = 1:n_g
            tr = aligned_traces_ds{idx(k), di};
            if isempty(tr) || all(isnan(tr)), continue; end
            plot(ax, t_ms(1:min(end,numel(tr))), tr(1:min(end,numel(t_ms))), ...
                'Color', [cell_colors(k,:), 0.6], 'LineWidth', 0.5);
        end

        % Vertical line at peak alignment point
        xline(ax, 0, 'Color', [0.3 0.3 0.3], 'LineWidth', 1);

        xlim(ax, [t_ms(1), t_ms(end)]);
        set(ax, 'XTick', [], 'YTick', [], 'FontSize', 5);

        if di == 1
            ylabel(ax, sprintf('%s (n=%d)', group_names{g}, n_g), 'FontSize', 7);
        end
        if g == 1
            if di == 5
                title(ax, 'PD', 'FontSize', 7, 'FontWeight', 'bold');
            elseif di == 13
                title(ax, 'ND', 'FontSize', 7);
            else
                title(ax, sprintf('%.0f°', pd_aligned_angles(di)), 'FontSize', 6);
            end
        end
    end
end

title(t_diag, 'Alignment Diagnostic — Per-Cell Traces (peak-aligned, ±1400ms)', ...
    'FontSize', 12);

fname_diag = 'diag_alignment_quality.png';
exportgraphics(fig_diag, fullfile(preview_dir, fname_diag), 'Resolution', 150);
fprintf('  Saved: %s\n', fname_diag);
close(fig_diag);

%% Step 10: Save time_to_max data
save_data.time_to_max   = time_to_max;
save_data.align_valid   = align_valid;
save_data.cell_info     = rmfield(all_cells, 'traces_aligned');
save_data.pd_aligned_angles = pd_aligned_angles;
save_data.params = struct( ...
    'smooth_win', SMOOTH_WIN, ...
    'bl_start', BL_START, 'bl_end', BL_END, ...
    'stim_start', STIM_START, 'stim_trim_end', STIM_TRIM_END, ...
    'display_half', DISPLAY_HALF, 'ds_factor', DS_FACTOR);

save_path = fullfile(data_root, 'population_results', 'ring_of_traces_alignment.mat');
save(save_path, '-struct', 'save_data');
fprintf('\nAlignment data saved to: %s\n', save_path);

%% Step 11: SymNN vs VecNN — Compute SymNN PDs and per-cell metrics
%   For each cell, compute SymNN PD via compute_pd_four_methods.
%   Compare to current VecNN PD. Identify cells that differ.
%   Compute per-cell DSI (PD-ND) and Aspect Ratio (PD+ND method) for both methods.

fprintf('\n--- Step 11: SymNN vs VecNN comparison ---\n');

n_cells = numel(all_cells);
symnn_shift   = zeros(n_cells, 1);   % direction-slot shift needed
dsi_vecnn_all = NaN(n_cells, 1);
dsi_symnn_all = NaN(n_cells, 1);
ar_vecnn_all  = NaN(n_cells, 1);     % PD+ND method
ar_symnn_all  = NaN(n_cells, 1);     % PD+ND method
symnn_polar   = cell(n_cells, 1);    % 16x2 SymNN-aligned [angle, peak]
vecnn_pd_deg  = NaN(n_cells, 1);
symnn_pd_deg  = NaN(n_cells, 1);

for ci = 1:n_cells
    pa = all_cells(ci).peak_amps_995;  % 16x2 [angle_rad, peak_amp]
    if isempty(pa) || size(pa, 2) ~= 2
        continue;
    end

    methods = compute_pd_four_methods(pa(:, 1), pa(:, 2));

    % Per-cell DSI from both methods
    dsi_vecnn_all(ci) = methods.vec_nn.dsi_pdnd;
    dsi_symnn_all(ci) = methods.sym_nn.dsi_pdnd;

    % Per-cell Aspect Ratio (PD+ND method) from both aligned curves
    ar_vecnn_all(ci) = compute_ar_pdnd(methods.vec_nn.aligned);
    ar_symnn_all(ci) = compute_ar_pdnd(methods.sym_nn.aligned);

    % Store SymNN-aligned polar data
    symnn_polar{ci} = methods.sym_nn.aligned;

    vecnn_pd_deg(ci) = methods.vec_nn.pd_deg;
    symnn_pd_deg(ci) = methods.sym_nn.pd_deg;

    % Determine circular shift: how many direction slots to shift
    %   VecNN PD is already at position 5 (90°). SymNN PD may be different.
    %   The input peak_amps_995 is already VecNN-aligned, so VecNN PD = 90°.
    %   SymNN PD is expressed in that same aligned coordinate frame.
    %   We need to find which position slot SymNN PD falls in, then shift so
    %   it goes to position 5.
    circ_diff = mod(pd_aligned_angles - methods.sym_nn.pd_deg + 180, 360) - 180;
    [~, sym_pos] = min(abs(circ_diff));
    symnn_shift(ci) = 5 - sym_pos;  % shift needed to put SymNN PD at slot 5
end

% Report which cells differ
n_differ = sum(symnn_shift ~= 0);
fprintf('  Cells with different SymNN vs VecNN PD: %d / %d (%.1f%%)\n', ...
    n_differ, n_cells, 100 * n_differ / n_cells);
for ci = find(symnn_shift ~= 0)'
    fprintf('    %s: VecNN=%.1f° SymNN=%.1f° (shift=%+d slots)\n', ...
        all_cells(ci).date_str, vecnn_pd_deg(ci), symnn_pd_deg(ci), symnn_shift(ci));
end


%% Step 12: Circshift traces for SymNN, re-run xcorr, compute population stats
%   For affected cells, circshift the 16-direction trace arrays.
%   Then re-run xcorr refinement and compute population statistics.

fprintf('\n--- Step 12: SymNN xcorr refinement ---\n');

% Circshift aligned traces for SymNN
symnn_traces_ds = aligned_traces_ds;  % start from peak-aligned VecNN traces
for ci = find(symnn_shift ~= 0)'
    symnn_traces_ds(ci, :) = circshift(aligned_traces_ds(ci, :), [0, symnn_shift(ci)]);
end
fprintf('  Circshifted %d cells for SymNN alignment.\n', n_differ);

% xcorr refinement for SymNN (same algorithm as Step 7b)
xcorr_symnn_traces_ds = symnn_traces_ds;

for g = 1:4
    idx = groups(g).indices;
    n_g = numel(idx);
    if n_g < 2, continue; end

    for di = 1:16
        % Compute group mean template from SymNN peak-aligned traces
        traces_mat = NaN(display_len_ds, n_g);
        for k = 1:n_g
            tr = symnn_traces_ds{idx(k), di};
            len = min(numel(tr), display_len_ds);
            traces_mat(1:len, k) = tr(1:len);
        end
        template = mean(traces_mat, 2, 'omitnan');
        template_clean = template;
        template_clean(isnan(template_clean)) = 0;
        if all(template_clean == 0), continue; end

        for k = 1:n_g
            ci = idx(k);
            tr = symnn_traces_ds{ci, di};
            if isempty(tr) || all(isnan(tr)), continue; end
            tr_clean = tr(:);
            tr_clean(isnan(tr_clean)) = 0;
            if all(tr_clean == 0), continue; end

            [r, lags] = xcorr(tr_clean, template_clean, XCORR_MAX_LAG_DS);
            [~, max_r_idx] = max(r);
            best_lag = lags(max_r_idx);

            corrected = NaN(display_len_ds, 1);
            for out_idx = 1:display_len_ds
                src_idx = out_idx + best_lag;
                if src_idx >= 1 && src_idx <= numel(tr)
                    corrected(out_idx) = tr(src_idx);
                end
            end
            xcorr_symnn_traces_ds{ci, di} = corrected;
        end
    end
end
fprintf('  xcorr refinement for SymNN complete.\n');

% Compute SymNN population statistics (same structure as groups_xcorr)
groups_symnn = groups;  % copy structure
for g = 1:4
    idx = groups_symnn(g).indices;
    n_g = numel(idx);

    groups_symnn(g).mean_traces = cell(16, 1);
    groups_symnn(g).sem_traces  = cell(16, 1);

    % Polar stats: use SymNN-aligned peak amplitudes
    polar_amps = NaN(16, n_g);
    for k = 1:n_g
        sp = symnn_polar{idx(k)};
        if ~isempty(sp) && size(sp, 2) == 2
            polar_amps(:, k) = sp(:, 2);
        end
    end
    groups_symnn(g).polar_mean = mean(polar_amps, 2, 'omitnan');
    groups_symnn(g).polar_sem  = std(polar_amps, 0, 2, 'omitnan') / sqrt(n_g);

    for di = 1:16
        traces_mat = NaN(display_len_ds, n_g);
        for k = 1:n_g
            tr = xcorr_symnn_traces_ds{idx(k), di};
            len = min(numel(tr), display_len_ds);
            traces_mat(1:len, k) = tr(1:len);
        end
        groups_symnn(g).mean_traces{di} = mean(traces_mat, 2, 'omitnan');
        groups_symnn(g).sem_traces{di}  = std(traces_mat, 0, 2, 'omitnan') / sqrt(n_g);
    end
end


%% Step 13: Generate SymNN population ring-of-traces (2 PNGs)
fprintf('\n--- Step 13: SymNN population ring-of-traces ---\n');

% ON: ctrl vs TTL (SymNN)
fig_on_sym = plot_ring_of_traces_population( ...
    groups_symnn(1), groups_symnn(2), 'ON', ...
    pd_aligned_angles, all_cells, ...
    display_len_ds, ref_sample_ds, DS_FACTOR, xcorr_symnn_traces_ds, RING_LAYOUT);
sgtitle(fig_on_sym, 'ON cells — SymNN alignment — ctrl vs TTL (xcorr, mean±SEM)', ...
    'FontSize', RING_LAYOUT.font_title);
exportgraphics(fig_on_sym, fullfile(preview_dir, 'pop_ring_on_ctrl_vs_ttl_symnn.png'), ...
    'Resolution', 200);
close(fig_on_sym);
fprintf('  Saved: pop_ring_on_ctrl_vs_ttl_symnn.png\n');

% OFF: ctrl vs TTL (SymNN)
fig_off_sym = plot_ring_of_traces_population( ...
    groups_symnn(3), groups_symnn(4), 'OFF', ...
    pd_aligned_angles, all_cells, ...
    display_len_ds, ref_sample_ds, DS_FACTOR, xcorr_symnn_traces_ds, RING_LAYOUT);
sgtitle(fig_off_sym, 'OFF cells — SymNN alignment — ctrl vs TTL (xcorr, mean±SEM)', ...
    'FontSize', RING_LAYOUT.font_title);
exportgraphics(fig_off_sym, fullfile(preview_dir, 'pop_ring_off_ctrl_vs_ttl_symnn.png'), ...
    'Resolution', 200);
close(fig_off_sym);
fprintf('  Saved: pop_ring_off_ctrl_vs_ttl_symnn.png\n');


%% Step 14: DSI and AR (PD+ND) comparison — VecNN vs SymNN
fprintf('\n--- Step 14: DSI and AR comparison figures ---\n');

% --- Figure A: DSI comparison (4 groups × 2 methods) ---
fig_dsi = figure('Position', [50 50 900 500], 'Color', 'w');

group_names = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
group_idx_list = {on_ctrl_idx, on_ttl_idx, off_ctrl_idx, off_ttl_idx};
group_colors = {[0.3 0.3 0.3], [0.9 0.1 0.1], [0.3 0.3 0.3], [0.9 0.1 0.1]};

ax_dsi = axes(fig_dsi);
hold(ax_dsi, 'on');

x_positions = [];
x_labels_all = {};
x_tick_pos = [];

for g = 1:4
    idx = group_idx_list{g};
    x_base = (g - 1) * 3 + 1;  % 1, 4, 7, 10

    % VecNN
    scatter_with_mean(ax_dsi, x_base, dsi_vecnn_all(idx), group_colors{g}, 'o');
    % SymNN
    scatter_with_mean(ax_dsi, x_base + 1, dsi_symnn_all(idx), group_colors{g}, 's');

    x_tick_pos = [x_tick_pos, x_base, x_base + 1]; %#ok<AGROW>
    x_labels_all = [x_labels_all, {'Vec', 'Sym'}]; %#ok<AGROW>

    % Add group name above
    text(ax_dsi, x_base + 0.5, max(ax_dsi.YLim) * 0.99, group_names{g}, ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end

ax_dsi.XTick = x_tick_pos;
ax_dsi.XTickLabel = x_labels_all;
ax_dsi.XLim = [0, 13];
ylabel(ax_dsi, 'DSI (PD-ND)');
title(ax_dsi, 'Direction Selectivity Index: VecNN vs SymNN');
box(ax_dsi, 'on');
grid(ax_dsi, 'on');

% Add group-level summary statistics text
summary_lines = {};
for g = 1:4
    idx = group_idx_list{g};
    v_vals = dsi_vecnn_all(idx);
    s_vals = dsi_symnn_all(idx);
    v_vals = v_vals(~isnan(v_vals));
    s_vals = s_vals(~isnan(s_vals));
    summary_lines{end+1} = sprintf('%s: Vec %.3f±%.3f, Sym %.3f±%.3f (n=%d)', ...
        group_names{g}, mean(v_vals), std(v_vals)/sqrt(numel(v_vals)), ...
        mean(s_vals), std(s_vals)/sqrt(numel(s_vals)), numel(v_vals)); %#ok<AGROW>
end
annotation(fig_dsi, 'textbox', [0.15 0.01 0.75 0.12], ...
    'String', strjoin(summary_lines, '\n'), 'FontSize', 7, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FitBoxToText', 'on');

exportgraphics(fig_dsi, fullfile(preview_dir, 'compare_dsi_vecnn_vs_symnn.png'), ...
    'Resolution', 200);
close(fig_dsi);
fprintf('  Saved: compare_dsi_vecnn_vs_symnn.png\n');

% Print summary
fprintf('\n  DSI summary (mean±SEM):\n');
for g = 1:4
    idx = group_idx_list{g};
    v = dsi_vecnn_all(idx); v = v(~isnan(v));
    s = dsi_symnn_all(idx); s = s(~isnan(s));
    fprintf('    %s (n=%d): VecNN %.3f±%.3f, SymNN %.3f±%.3f\n', ...
        group_names{g}, numel(v), mean(v), std(v)/sqrt(numel(v)), ...
        mean(s), std(s)/sqrt(numel(s)));
end


% --- Figure B: Aspect Ratio (PD+ND) comparison (4 groups × 2 methods) ---
fig_ar = figure('Position', [50 50 900 500], 'Color', 'w');
ax_ar = axes(fig_ar);
hold(ax_ar, 'on');

x_tick_pos_ar = [];
x_labels_ar = {};

for g = 1:4
    idx = group_idx_list{g};
    x_base = (g - 1) * 3 + 1;

    % VecNN
    scatter_with_mean(ax_ar, x_base, ar_vecnn_all(idx), group_colors{g}, 'o');
    % SymNN
    scatter_with_mean(ax_ar, x_base + 1, ar_symnn_all(idx), group_colors{g}, 's');

    x_tick_pos_ar = [x_tick_pos_ar, x_base, x_base + 1]; %#ok<AGROW>
    x_labels_ar = [x_labels_ar, {'Vec', 'Sym'}]; %#ok<AGROW>

    text(ax_ar, x_base + 0.5, max(ax_ar.YLim) * 0.99, group_names{g}, ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end

ax_ar.XTick = x_tick_pos_ar;
ax_ar.XTickLabel = x_labels_ar;
ax_ar.XLim = [0, 13];
ylabel(ax_ar, 'Aspect Ratio (PD+ND)/(ortho CW + ortho CCW)');
title(ax_ar, 'Aspect Ratio (PD+ND method): VecNN vs SymNN');
box(ax_ar, 'on');
grid(ax_ar, 'on');

% Summary text
summary_lines_ar = {};
for g = 1:4
    idx = group_idx_list{g};
    v_vals = ar_vecnn_all(idx); v_vals = v_vals(~isnan(v_vals));
    s_vals = ar_symnn_all(idx); s_vals = s_vals(~isnan(s_vals));
    summary_lines_ar{end+1} = sprintf('%s: Vec %.2f±%.2f, Sym %.2f±%.2f (n=%d)', ...
        group_names{g}, mean(v_vals), std(v_vals)/sqrt(numel(v_vals)), ...
        mean(s_vals), std(s_vals)/sqrt(numel(s_vals)), numel(v_vals)); %#ok<AGROW>
end
annotation(fig_ar, 'textbox', [0.15 0.01 0.75 0.12], ...
    'String', strjoin(summary_lines_ar, '\n'), 'FontSize', 7, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FitBoxToText', 'on');

exportgraphics(fig_ar, fullfile(preview_dir, 'compare_ar_pdnd_vecnn_vs_symnn.png'), ...
    'Resolution', 200);
close(fig_ar);
fprintf('  Saved: compare_ar_pdnd_vecnn_vs_symnn.png\n');

% Print summary
fprintf('\n  Aspect Ratio (PD+ND) summary (mean±SEM):\n');
for g = 1:4
    idx = group_idx_list{g};
    v = ar_vecnn_all(idx); v = v(~isnan(v));
    s = ar_symnn_all(idx); s = s(~isnan(s));
    fprintf('    %s (n=%d): VecNN %.2f±%.2f, SymNN %.2f±%.2f\n', ...
        group_names{g}, numel(v), mean(v), std(v)/sqrt(numel(v)), ...
        mean(s), std(s)/sqrt(numel(s)));
end


fprintf('\n=== Phase 0A Complete ===\n');
fprintf('Inspect PNGs in: %s\n', preview_dir);
fprintf('New SymNN comparison PNGs:\n');
fprintf('  pop_ring_on_ctrl_vs_ttl_symnn.png\n');
fprintf('  pop_ring_off_ctrl_vs_ttl_symnn.png\n');
fprintf('  compare_dsi_vecnn_vs_symnn.png\n');
fprintf('  compare_ar_pdnd_vecnn_vs_symnn.png\n');


%% ========================= Helper Functions ============================

function fig = plot_ring_of_traces_population(g_ctrl, g_ttl, on_off, ...
    pd_aligned_angles, all_cells, ...
    display_len_ds, ref_sample_ds, DS_FACTOR, aligned_traces_ds, layout)
% PLOT_RING_OF_TRACES_POPULATION  Population ring-of-traces with PD pointing up.
%   16 radial timeseries subplots (ctrl=black, TTL=red, mean±SEM) + central
%   polar with filled patch shading (99.5th percentile peaks).

    fig = figure('Position', layout.fig_size, 'Color', 'w');

    % Time axis (ms relative to aligned peak)
    t_ms = ((1:display_len_ds) - ref_sample_ds) * DS_FACTOR * 0.1;

    % Auto y-limits from all mean traces
    all_center = [];
    for di = 1:16
        all_center = [all_center; g_ctrl.mean_traces{di}; ...
            g_ttl.mean_traces{di}]; %#ok<AGROW>
    end
    y_range = [min(all_center, [], 'omitnan'), max(all_center, [], 'omitnan')];
    y_pad = 0.15 * diff(y_range);
    y_lim = [y_range(1) - y_pad, y_range(2) + y_pad];

    % --- 16 radial timeseries subplots ---
    ax_null = [];  % save null-direction axes for scale bar
    for di = 1:16
        angle_deg = pd_aligned_angles(di);
        angle_rad = deg2rad(180 - angle_deg);

        x_pos = layout.centerX + layout.radius * cos(angle_rad);
        y_pos = layout.centerY + layout.radius * sin(angle_rad);
        ax = axes('Position', [x_pos - layout.subW/2, y_pos - layout.subH/2, ...
            layout.subW, layout.subH]); %#ok<LAXES>
        hold(ax, 'on');

        % Ctrl: black mean ± SEM
        plot_shaded(ax, t_ms, g_ctrl.mean_traces{di}, ...
            g_ctrl.sem_traces{di}, [0 0 0], 0.15);
        % TTL: red mean ± SEM
        plot_shaded(ax, t_ms, g_ttl.mean_traces{di}, ...
            g_ttl.sem_traces{di}, [1 0 0], 0.15);

        ylim(ax, y_lim);
        xlim(ax, [t_ms(1), t_ms(end)]);
        axis(ax, 'off');

        if di == 13, ax_null = ax; end  % null direction (270°, bottom)
    end

    % --- Central polar plot (99.5th percentile, filled patch shading) ---
    cs = layout.radius * layout.polar_scale;
    polar_pos = [layout.centerX - cs, layout.centerY - cs, 2*cs, 2*cs];

    theta = all_cells(g_ctrl.indices(1)).peak_amps_995(:, 1);

    polar_opts = struct( ...
        'n_ctrl', g_ctrl.n, 'n_ttl', g_ttl.n, ...
        'ctrl_label', 'control', 'ttl_label', 'TTL');
    plot_polar_with_patch(polar_pos, theta, ...
        g_ctrl.polar_mean, g_ctrl.polar_sem, ...
        g_ttl.polar_mean, g_ttl.polar_sem, polar_opts);

    % --- Scale bars on null-direction subplot ---
    add_scale_bar_on_axes(ax_null, t_ms, y_lim);

    % --- Title ---
    sgtitle(sprintf('%s Cells — Population mean\\pmSEM (ctrl n=%d, TTL n=%d)', ...
        on_off, g_ctrl.n, g_ttl.n), 'FontSize', layout.font_title);
end


function fig = plot_ring_of_traces_single_pair(all_cells, ctrl_ci, ttl_ci, on_off, ...
    pd_aligned_angles, aligned_traces_ds, ...
    display_len_ds, ref_sample_ds, DS_FACTOR, layout)
% PLOT_RING_OF_TRACES_SINGLE_PAIR  Single ctrl vs TTL cell ring-of-traces.
%   Same layout as population figure but showing individual cell traces.
%   Ctrl = black, TTL = red.

    fig = figure('Position', layout.fig_size, 'Color', 'w');

    % Time axis
    t_ms = ((1:display_len_ds) - ref_sample_ds) * DS_FACTOR * 0.1;

    % Get traces for both cells
    ctrl_traces = cell(16, 1);
    ttl_traces  = cell(16, 1);
    all_mean = [];
    for di = 1:16
        ctrl_traces{di} = aligned_traces_ds{ctrl_ci, di};
        ttl_traces{di}  = aligned_traces_ds{ttl_ci, di};
        all_mean = [all_mean; ctrl_traces{di}(:); ttl_traces{di}(:)]; %#ok<AGROW>
    end
    y_range = [min(all_mean, [], 'omitnan'), max(all_mean, [], 'omitnan')];
    y_pad = 0.15 * diff(y_range);
    y_lim = [y_range(1) - y_pad, y_range(2) + y_pad];

    % --- 16 radial timeseries subplots ---
    ax_null = [];  % save null-direction axes for scale bar
    for di = 1:16
        angle_deg = pd_aligned_angles(di);
        angle_rad = deg2rad(180 - angle_deg);

        x_pos = layout.centerX + layout.radius * cos(angle_rad);
        y_pos = layout.centerY + layout.radius * sin(angle_rad);
        ax = axes('Position', [x_pos - layout.subW/2, y_pos - layout.subH/2, ...
            layout.subW, layout.subH]); %#ok<LAXES>
        hold(ax, 'on');

        % Ctrl trace (black)
        tr_c = ctrl_traces{di};
        len_c = min(numel(tr_c), numel(t_ms));
        if ~isempty(tr_c) && ~all(isnan(tr_c))
            plot(ax, t_ms(1:len_c), tr_c(1:len_c), 'Color', [0 0 0], 'LineWidth', 1.0);
        end

        % TTL trace (red)
        tr_t = ttl_traces{di};
        len_t = min(numel(tr_t), numel(t_ms));
        if ~isempty(tr_t) && ~all(isnan(tr_t))
            plot(ax, t_ms(1:len_t), tr_t(1:len_t), 'Color', [1 0 0], 'LineWidth', 1.0);
        end

        ylim(ax, y_lim);
        xlim(ax, [t_ms(1), t_ms(end)]);
        axis(ax, 'off');

        if di == 13, ax_null = ax; end  % null direction (270°, bottom)
    end

    % --- Central polar plot ---
    cs = layout.radius * layout.polar_scale;
    ax_polar = polaraxes('Position', ...
        [layout.centerX - cs, layout.centerY - cs, 2*cs, 2*cs]);
    hold(ax_polar, 'on');

    % Individual cell polar tuning
    theta_c = all_cells(ctrl_ci).max_v_aligned(:, 1);
    rho_c   = all_cells(ctrl_ci).max_v_aligned(:, 2);
    theta_t = all_cells(ttl_ci).max_v_aligned(:, 1);
    rho_t   = all_cells(ttl_ci).max_v_aligned(:, 2);

    polarplot(ax_polar, [theta_c; theta_c(1)], [rho_c; rho_c(1)], ...
        'k-', 'LineWidth', 1.5);
    polarplot(ax_polar, [theta_t; theta_t(1)], [rho_t; rho_t(1)], ...
        'Color', [1 0 0], 'LineWidth', 1.5);

    ax_polar.ThetaZeroLocation = 'right';
    ax_polar.ThetaDir = 'counterclockwise';
    ax_polar.ThetaTick = 0:22.5:337.5;   % 8 diameters at 22.5° spacing
    ax_polar.ThetaTickLabel = {};
    ax_polar.FontSize = 6;
    ax_polar.LineWidth = 1;

    % --- Scale bars on null-direction subplot ---
    add_scale_bar_on_axes(ax_null, t_ms, y_lim);

    % --- Title ---
    ctrl_dsi = all_cells(ctrl_ci).dsi;
    ttl_dsi  = all_cells(ttl_ci).dsi;
    sgtitle(sprintf('%s Cells — Representative (ctrl: %s, DSI=%.2f | TTL: %s, DSI=%.2f)', ...
        on_off, all_cells(ctrl_ci).date_str, ctrl_dsi, ...
        all_cells(ttl_ci).date_str, ttl_dsi), 'FontSize', layout.font_title - 1);
end


function add_scale_bar_on_axes(ax, t_ms, y_lim)
% ADD_SCALE_BAR_ON_AXES  Draw time + voltage scale bars on the given subplot.
%   Uses data coordinates (ms, mV) so the bar lengths are visually correct.
%   Placed in the lower-right corner of the axes.

    if isempty(ax) || ~isvalid(ax), return; end

    % Scale bar sizes in data units
    bar_t = 500;  % 500 ms time bar
    v_range = diff(y_lim);
    bar_v = max(2, round(v_range * 0.3));  % ~30% of y-range, min 2 mV

    % Position: lower-right corner of the axes
    x_right = t_ms(end) - 0.02 * (t_ms(end) - t_ms(1));  % 2% inset from right
    y_bottom = y_lim(1) + 0.05 * v_range;                  % 5% above bottom

    x_left = x_right - bar_t;
    y_top = y_bottom + bar_v;

    % Draw L-shaped scale bar
    plot(ax, [x_left, x_right], [y_bottom, y_bottom], 'k-', 'LineWidth', 1.5, ...
        'Clipping', 'off');
    plot(ax, [x_left, x_left], [y_bottom, y_top], 'k-', 'LineWidth', 1.5, ...
        'Clipping', 'off');

    % Labels
    text(ax, (x_left + x_right)/2, y_bottom - 0.06*v_range, '500 ms', ...
        'HorizontalAlignment', 'center', 'FontSize', 7, 'Clipping', 'off');
    text(ax, x_left - 0.02*(t_ms(end)-t_ms(1)), (y_bottom + y_top)/2, ...
        sprintf('%d mV', bar_v), ...
        'HorizontalAlignment', 'right', 'FontSize', 7, 'Rotation', 90, ...
        'Clipping', 'off');
end


function ci = find_cell_by_date(all_cells, date_str)
% FIND_CELL_BY_DATE  Find cell index matching a date string prefix.
    ci = [];
    for i = 1:numel(all_cells)
        if startsWith(all_cells(i).date_str, date_str)
            ci = i;
            return;
        end
    end
end


function [traces_16, pd_shift, lut_dirs_ordered] = extract_and_align_traces( ...
    bar_data, plot_order, lut_directions, pd_direction)
% EXTRACT_AND_ALIGN_TRACES  Extract 16 mean traces and PD-align by circular shift.
%   Uses externally-provided pd_direction (from batch results) for alignment.
%   Returns traces in PD-aligned order: position 5 = PD (90° in aligned frame).

    n_reps = size(bar_data, 2) - 1;  % last column is mean

    % Extract mean traces in subplot order
    traces_subplot_order = cell(16, 1);
    for si = 1:16
        data_row = plot_order(si);
        traces_subplot_order{si} = bar_data{data_row, n_reps + 1};
    end

    % Sort by direction angle (ascending)
    lut_dirs_subplot = lut_directions(plot_order);
    [sorted_angles, sort_idx] = sort(lut_dirs_subplot(:));
    traces_sorted = traces_subplot_order(sort_idx);

    % Find which sorted position is closest to the validated PD direction
    angle_diffs = abs(mod(sorted_angles - pd_direction + 180, 360) - 180);
    [~, pd_sorted_idx] = min(angle_diffs);

    % Circular shift so PD goes to position 5 (= 90° in aligned frame)
    pd_shift = 5 - pd_sorted_idx;
    traces_16 = circshift(traces_sorted, pd_shift);

    lut_dirs_ordered = circshift(sorted_angles, pd_shift);
end


function shifted = shift_and_window(tr, peak_time, ref_sample, display_len)
% SHIFT_AND_WINDOW  Shift trace so peak_time → ref_sample, window to display_len.
%   Pads with NaN at edges.

    shifted = NaN(display_len, 1);
    offset = round(peak_time) - ref_sample;

    for out_idx = 1:display_len
        src_idx = out_idx + offset;
        if src_idx >= 1 && src_idx <= numel(tr)
            shifted(out_idx) = tr(src_idx);
        end
    end
end


function shift_time = borrow_nearest_shift(time_to_max_row, valid_row, target_di)
% BORROW_NEAREST_SHIFT  Find nearest valid direction's alignment time (circular).

    shift_time = NaN;
    for dist = 1:8
        % Check both CW and CCW neighbors
        cw  = mod(target_di - 1 + dist, 16) + 1;
        ccw = mod(target_di - 1 - dist, 16) + 1;

        if valid_row(cw)
            shift_time = time_to_max_row(cw);
            return;
        end
        if valid_row(ccw)
            shift_time = time_to_max_row(ccw);
            return;
        end
    end
end


function plot_shaded(ax, x, m, s, color, alpha)
% PLOT_SHADED  Plot mean line with SEM shading.
    m = m(:)'; s = s(:)'; x = x(:)';

    % Trim to common length
    len = min([numel(x), numel(m), numel(s)]);
    x = x(1:len); m = m(1:len); s = s(1:len);

    % Remove NaN for fill
    valid = ~isnan(m) & ~isnan(s);
    if sum(valid) < 2, return; end

    xv = x(valid); mv = m(valid); sv = s(valid);

    fill(ax, [xv, fliplr(xv)], [mv+sv, fliplr(mv-sv)], ...
        color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    plot(ax, xv, mv, 'Color', color, 'LineWidth', 1.2);
end


function group = classify_group(is_on, is_ttl)
    if is_on && ~is_ttl
        group = 'on_control';
    elseif is_on && is_ttl
        group = 'on_ttl';
    elseif ~is_on && ~is_ttl
        group = 'off_control';
    else
        group = 'off_ttl';
    end
end


function exp_list = discover_early_experiments(data_root)
% DISCOVER_EARLY_EXPERIMENTS  Walk {control,ttl}/{ON,OFF}/ tree.

    exp_list = struct('folder', {}, 'date_str', {}, ...
        'treatment', {}, 'cell_type', {});

    for treatment = ["control", "ttl"]
        for cell_type = ["ON", "OFF"]
            subdir = fullfile(data_root, treatment, cell_type);
            if ~isfolder(subdir), continue; end

            dd = dir(subdir);
            dd = dd([dd.isdir]);
            dd = dd(~startsWith({dd.name}, '.'));

            for k = 1:numel(dd)
                exp_folder = fullfile(subdir, dd(k).name);
                if ~isfile(fullfile(exp_folder, 'currentExp.mat'))
                    continue;
                end

                ei.folder    = exp_folder;
                ei.date_str  = dd(k).name;
                ei.treatment = char(treatment);
                ei.cell_type = char(cell_type);
                exp_list(end + 1) = ei; %#ok<AGROW>
            end
        end
    end
end


function entry = find_batch_entry(batch_results, folder_name_or_path)
% FIND_BATCH_ENTRY  Look up a batch results entry by folder name or full path.
%   Matches by checking if the batch entry's folder ends with the query,
%   or if the query ends with the batch folder's basename.

    entry = [];
    [~, query_base] = fileparts(folder_name_or_path);
    if isempty(query_base)
        query_base = folder_name_or_path;
    end

    for i = 1:numel(batch_results)
        batch_folder = batch_results(i).folder;
        [~, batch_base] = fileparts(batch_folder);

        if strcmp(batch_base, query_base) || ...
           endsWith(batch_folder, folder_name_or_path) || ...
           strcmp(batch_folder, folder_name_or_path)
            entry = batch_results(i);
            return;
        end
    end
end


function [axPolar, axFill] = plot_polar_with_patch(ax_position, ...
    theta, center_ctrl, spread_ctrl, center_ttl, spread_ttl, opts)
% PLOT_POLAR_WITH_PATCH  Central polar plot with filled patch SEM/MAD shading.
%   Uses dual-axis Cartesian overlay technique (from plot_polar_population.m).
%   Polar axes provide grid/ticks; transparent Cartesian overlay draws patch
%   objects and center lines, guaranteeing perfect alignment.
%
%   INPUTS:
%     ax_position   - [x y w h] normalized figure position for the polar axes
%     theta         - 16×1 angles in radians (PD-aligned)
%     center_ctrl   - 16×1 center values for control group (mean or median)
%     spread_ctrl   - 16×1 spread values for control group (SEM or MAD)
%     center_ttl    - 16×1 center values for TTL group
%     spread_ttl    - 16×1 spread values for TTL group
%     opts          - Structure with optional fields:
%                       .ctrl_line_color  - [r g b] (default: [0 0 0])
%                       .ctrl_fill_color  - [r g b] (default: [0.80 0.80 0.80])
%                       .ttl_line_color   - [r g b] (default: [1 0 0])
%                       .ttl_fill_color   - [r g b] (default: [1 0.70 0.70])
%                       .alpha            - patch alpha (default: 0.40)
%                       .ctrl_label       - string (default: 'control')
%                       .ttl_label        - string (default: 'TTL')
%                       .n_ctrl           - count for legend
%                       .n_ttl            - count for legend

    if nargin < 7, opts = struct(); end
    if ~isfield(opts, 'ctrl_line_color'), opts.ctrl_line_color = [0 0 0]; end
    if ~isfield(opts, 'ctrl_fill_color'), opts.ctrl_fill_color = [0.80 0.80 0.80]; end
    if ~isfield(opts, 'ttl_line_color'),  opts.ttl_line_color  = [1 0 0]; end
    if ~isfield(opts, 'ttl_fill_color'),  opts.ttl_fill_color  = [1 0.70 0.70]; end
    if ~isfield(opts, 'alpha'),           opts.alpha           = 0.40; end
    if ~isfield(opts, 'ctrl_label'),      opts.ctrl_label      = 'control'; end
    if ~isfield(opts, 'ttl_label'),       opts.ttl_label       = 'TTL'; end

    % Create polar axes for grid and ticks
    axPolar = polaraxes('Position', ax_position);
    hold(axPolar, 'on');

    % Estimate rmax from both groups
    all_upper = [];
    if ~isempty(center_ctrl)
        all_upper = [all_upper; center_ctrl(:) + spread_ctrl(:)];
    end
    if ~isempty(center_ttl)
        all_upper = [all_upper; center_ttl(:) + spread_ttl(:)];
    end
    if ~isempty(all_upper)
        rmax = max(all_upper, [], 'omitnan');
        if isfinite(rmax) && rmax > 0
            rlim(axPolar, [0, rmax * 1.1]);
        end
    end

    % Create transparent Cartesian overlay for patch + line drawing
    axFill = axes('Position', axPolar.Position, 'Color', 'none', ...
                  'XColor', 'none', 'YColor', 'none', 'HitTest', 'off');
    axis(axFill, 'equal');
    hold(axFill, 'on');

    % Sync Cartesian limits to polar r-limits
    rL = rlim(axPolar);
    set(axFill, 'XLim', [-rL(2) rL(2)], 'YLim', [-rL(2) rL(2)]);

    % Plot control group
    hCtrl = gobjects(1, 1);
    if ~isempty(center_ctrl) && any(~isnan(center_ctrl))
        hCtrl = draw_polar_patch(axFill, theta, center_ctrl, spread_ctrl, ...
            opts.ctrl_line_color, opts.ctrl_fill_color, opts.alpha, 2);
    end

    % Plot TTL group
    hTtl = gobjects(1, 1);
    if ~isempty(center_ttl) && any(~isnan(center_ttl))
        hTtl = draw_polar_patch(axFill, theta, center_ttl, spread_ttl, ...
            opts.ttl_line_color, opts.ttl_fill_color, opts.alpha, 2);
    end

    % Re-sync Cartesian overlay after plotting
    rL = rlim(axPolar);
    set(axFill, 'XLim', [-rL(2) rL(2)], 'YLim', [-rL(2) rL(2)]);
    set(axFill, 'Position', axPolar.Position);

    % Keep polar axes underneath for grid/ticks
    axPolar.Color = 'none';
    uistack(axFill, 'top');

    % Legend on Cartesian overlay
    hs = []; labs = {};
    if ~isempty(center_ctrl) && any(~isnan(center_ctrl))
        hs(end+1)  = hCtrl;
        n_str = '';
        if isfield(opts, 'n_ctrl'), n_str = sprintf(' (n=%d)', opts.n_ctrl); end
        labs{end+1} = [opts.ctrl_label, n_str];
    end
    if ~isempty(center_ttl) && any(~isnan(center_ttl))
        hs(end+1)  = hTtl;
        n_str = '';
        if isfield(opts, 'n_ttl'), n_str = sprintf(' (n=%d)', opts.n_ttl); end
        labs{end+1} = [opts.ttl_label, n_str];
    end
    if ~isempty(hs)
        legend(axFill, hs, labs, 'Location', 'best', 'FontSize', 6);
    end

    % Orient polar so PD (π/2 = 90°) points up
    axPolar.ThetaZeroLocation = 'right';
    axPolar.ThetaDir = 'counterclockwise';
    axPolar.ThetaTick = 0:22.5:337.5;   % 16 radial spokes = 8 diameters at 22.5° spacing
    axPolar.ThetaTickLabel = {};
    axPolar.FontSize = 6;
    axPolar.LineWidth = 1;
end


function hLine = draw_polar_patch(axFill, theta, centerVals, bandVals, ...
    lineColor, fillColor, alphaVal, lw)
% DRAW_POLAR_PATCH  Draw center line + spread patch on Cartesian overlay.
%   Reference: plot_polar_population.m plot_with_shade() function.

    th = theta(:);
    m  = centerVals(:);
    b  = bandVals(:);

    % Ensure non-negative radii
    upper = max(m + b, 0);
    lower = max(m - b, 0);

    % Close the loop
    th_closed    = [th; th(1) + 2*pi];
    upper_closed = [upper; upper(1)];
    lower_closed = [lower; lower(1)];
    m_closed     = [m; m(1)];

    % Build polygon: outer ring forward, inner ring backward
    th_poly = [th_closed; flipud(th_closed)];
    r_poly  = [upper_closed; flipud(lower_closed)];
    [x_poly, y_poly] = pol2cart(th_poly, r_poly);

    % Shaded spread patch
    patch('XData', x_poly, 'YData', y_poly, ...
          'FaceColor', fillColor, 'FaceAlpha', alphaVal, ...
          'EdgeColor', 'none', 'Parent', axFill, 'HitTest', 'off');

    % Center line
    [x_line, y_line] = pol2cart(th_closed, m_closed);
    hLine = plot(axFill, x_line, y_line, '-', 'Color', lineColor, 'LineWidth', lw);
end


function cand_idx = find_top_n_by_dsi(group_idx, all_cells, n)
% FIND_TOP_N_BY_DSI  Select N cells closest to group mean DSI.
%   Returns cell indices from all_cells, sorted by |dsi - mean_dsi|.

    dsi_vals = arrayfun(@(i) all_cells(i).dsi, group_idx);
    mean_dsi = mean(dsi_vals, 'omitnan');
    dsi_dist = abs(dsi_vals - mean_dsi);
    [~, sort_order] = sort(dsi_dist);
    n = min(n, numel(group_idx));
    cand_idx = group_idx(sort_order(1:n));
end


function fig = plot_ring_of_traces_single_cell(all_cells, ci, ...
    pd_aligned_angles, aligned_traces_ds, ...
    display_len_ds, ref_sample_ds, DS_FACTOR, layout)
% PLOT_RING_OF_TRACES_SINGLE_CELL  Single-cell ring-of-traces figure.
%   Same radial layout as population figures. Shows individual cell trace
%   per direction (no SEM shading). Central polar shows tuning curve.
%   Title includes cell date, group, and DSI.

    fig = figure('Position', layout.fig_size, 'Color', 'w');

    % Time axis (ms relative to aligned peak)
    t_ms = ((1:display_len_ds) - ref_sample_ds) * DS_FACTOR * 0.1;

    % Group color
    if all_cells(ci).is_ttl
        trace_color = [1 0 0];      % red for TTL
    else
        trace_color = [0 0 0];      % black for control
    end

    % Auto y-limits from this cell's traces
    all_vals = [];
    for di = 1:16
        tr = aligned_traces_ds{ci, di};
        if ~isempty(tr) && ~all(isnan(tr))
            all_vals = [all_vals; tr(:)]; %#ok<AGROW>
        end
    end
    if isempty(all_vals)
        close(fig);
        return;
    end
    y_range = [min(all_vals, [], 'omitnan'), max(all_vals, [], 'omitnan')];
    y_pad = 0.15 * diff(y_range);
    y_lim = [y_range(1) - y_pad, y_range(2) + y_pad];

    % --- 16 radial timeseries subplots ---
    ax_null = [];  % save null-direction axes for scale bar
    for di = 1:16
        angle_deg = pd_aligned_angles(di);
        angle_rad = deg2rad(180 - angle_deg);

        x_pos = layout.centerX + layout.radius * cos(angle_rad);
        y_pos = layout.centerY + layout.radius * sin(angle_rad);
        ax = axes('Position', [x_pos - layout.subW/2, y_pos - layout.subH/2, ...
            layout.subW, layout.subH]); %#ok<LAXES>
        hold(ax, 'on');

        tr = aligned_traces_ds{ci, di};
        if ~isempty(tr) && ~all(isnan(tr))
            len = min(numel(tr), numel(t_ms));
            plot(ax, t_ms(1:len), tr(1:len), 'Color', trace_color, 'LineWidth', 1.0);
        end

        ylim(ax, y_lim);
        xlim(ax, [t_ms(1), t_ms(end)]);
        axis(ax, 'off');

        if di == 13, ax_null = ax; end  % null direction (270°, bottom)
    end

    % --- Central polar plot (individual cell tuning) ---
    cs = layout.radius * layout.polar_scale;
    ax_polar = polaraxes('Position', ...
        [layout.centerX - cs, layout.centerY - cs, 2*cs, 2*cs]);
    hold(ax_polar, 'on');

    theta_c = all_cells(ci).max_v_aligned(:, 1);
    rho_c   = all_cells(ci).max_v_aligned(:, 2);

    polarplot(ax_polar, [theta_c; theta_c(1)], [rho_c; rho_c(1)], ...
        'Color', trace_color, 'LineWidth', 1.5);

    ax_polar.ThetaZeroLocation = 'right';
    ax_polar.ThetaDir = 'counterclockwise';
    ax_polar.ThetaTick = 0:22.5:337.5;   % 8 diameters at 22.5° spacing
    ax_polar.ThetaTickLabel = {};
    ax_polar.FontSize = 6;
    ax_polar.LineWidth = 1;

    % --- Scale bars on null-direction subplot ---
    add_scale_bar_on_axes(ax_null, t_ms, y_lim);

    % --- Title ---
    group_str = classify_group(all_cells(ci).is_on, all_cells(ci).is_ttl);
    sgtitle(sprintf('%s — %s — DSI=%.3f', ...
        all_cells(ci).date_str, upper(strrep(group_str, '_', ' ')), ...
        all_cells(ci).dsi), 'FontSize', layout.font_title - 1);
end


function ar = compute_ar_pdnd(aligned)
% COMPUTE_AR_PDND  Aspect ratio using PD+ND method.
%   AR = (R_PD + R_ND) / (R_ortho_CW + R_ortho_CCW)
%   PD at pi/2, ND at 3*pi/2, ortho at 0 and pi.
%   Uses circular interpolation for grids that may not have exact values.

    angles    = aligned(:, 1);
    responses = aligned(:, 2);

    % Build circular interpolation (both-end wrap)
    angles_wrap    = [angles(end) - 2*pi; angles; angles(1) + 2*pi];
    responses_wrap = [responses(end);     responses; responses(1)];

    r_pd = interp1(angles_wrap, responses_wrap, pi/2, 'linear');
    r_nd = interp1(angles_wrap, responses_wrap, 3*pi/2, 'linear');
    r_o1 = interp1(angles_wrap, responses_wrap, 0, 'linear');
    r_o2 = interp1(angles_wrap, responses_wrap, pi, 'linear');

    denom = r_o1 + r_o2;
    if denom > 0
        ar = (r_pd + r_nd) / denom;
    else
        ar = NaN;
    end
end


function scatter_with_mean(ax, x_center, values, color, marker)
% SCATTER_WITH_MEAN  Jittered scatter + mean±SEM error bar.
%   Plots individual data points with horizontal jitter, plus a
%   filled marker and vertical error bar for the group mean±SEM.

    values = values(~isnan(values));
    if isempty(values), return; end

    n = numel(values);
    jitter = 0.15 * (rand(n, 1) - 0.5);  % ±0.075 horizontal spread

    % Individual points
    scatter(ax, x_center + jitter, values, 20, color, marker, ...
        'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0.6);

    % Mean ± SEM bar
    m = mean(values);
    s = std(values) / sqrt(n);
    plot(ax, [x_center, x_center], [m - s, m + s], '-', ...
        'Color', color, 'LineWidth', 2);
    plot(ax, x_center, m, marker, 'MarkerSize', 10, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', 'LineWidth', 1);
end

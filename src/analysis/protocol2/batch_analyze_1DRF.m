function results = batch_analyze_1DRF(data_root, opts)
% BATCH_ANALYZE_1DRF  Process all 1DRF experiments and generate population plots.
%
%   RESULTS = BATCH_ANALYZE_1DRF(DATA_ROOT) processes every experiment
%   folder in DATA_ROOT, classifies cells as ON/OFF and control/ttl, then
%   generates population plots with median +/- MAD (default) or mean +/-
%   SEM (configurable via opts.stat_method).
%
%   RESULTS = BATCH_ANALYZE_1DRF(DATA_ROOT, OPTS) uses the options
%   structure to override default parameters.
%
%   INPUTS:
%     data_root - Path to the 1DRF data directory containing experiment
%                 folders (e.g. '/path/to/protocol2/data/1DRF')
%     opts      - (Optional) structure. Inherits all fields from
%                 analyze_single_experiment opts, plus:
%                   .on_threshold  - Frame value threshold for ON vs OFF
%                                    classification (default: 129)
%                   .save_figs     - Save population figures as PDF
%                                    (default: true)
%                   .save_dir      - Output directory for population PDFs
%                                    (default: <data_root>/population_results)
%                   .flash_baseline - Sample indices for bar flash baseline
%                                     (default: 1:5000)
%                   .stat_method   - 'median_mad' (default) or 'mean_sem'
%
%   OUTPUT:
%     results - Structure array (one element per cell) with fields:
%       .folder            - Experiment folder name
%       .date_str          - Experiment date
%       .strain            - Strain name from metadata
%       .frame             - Frame number from metadata
%       .is_on             - Logical, true if ON cell (Frame > on_threshold)
%       .is_ttl            - Logical, true if Strain contains 'ttl'
%       .group             - String: 'on_control', 'on_ttl', 'off_control',
%                            or 'off_ttl'
%       .max_v_aligned     - 16x2 PD-aligned [angles, responses] from
%                            find_PD_and_order_idx
%       .pd_direction      - Preferred direction in degrees
%       .pd_orientation    - Bar orientation at PD in degrees
%       .pd_flash_bl       - 11xN baseline-subtracted mean flash traces
%                            along PD-ND axis
%       .ortho_flash_bl    - 11xN baseline-subtracted mean flash traces
%                            along orthogonal axis
%       .peak_pos          - Position of peak response (1-11, integer)
%       .peak_amplitudes   - 1x11 response amplitudes (mV)
%       .bump_width        - RF FWHM width in positions
%       .rest_voltage_median  - Median voltage during grey screen (mV)
%       .rest_voltage_mean    - Mean voltage during grey screen (mV)
%       .stim_voltage_median  - Median voltage during stimuli (mV)
%       .stim_voltage_mean    - Mean voltage during stimuli (mV)
%       .pd_flash_peak_aligned   - 11xN peak-aligned PD-ND traces
%       .ortho_flash_peak_aligned - 11xN peak-aligned orthogonal traces
%       .temporal_metrics  - Struct with rise_start, rise_time, decay_time
%                            and _delta variants (1x11 each)
%
%   FIGURES GENERATED:
%     For each of ON and OFF cells (6 figures total):
%       1. PD-aligned polar tuning curve (control=black, ttl=red, +/- spread)
%       2. 1x11 PD-ND bar flash (control=black, ttl=red, +/- spread)
%       3. 1x11 orthogonal bar flash (control=black, ttl=red, +/- spread)
%     Spread is MAD (default) or SEM, set via opts.stat_method.
%
%   ANALYSIS PIPELINE:
%     1. Discover experiment folders in data_root
%     2. For each experiment:
%        a. Load data and parse bar sweeps/flashes
%        b. Compute bar sweep responses, find PD, align to pi/2
%        c. Extract baseline-subtracted bar flash traces
%        d. Classify cell as ON/OFF and control/ttl
%     3. Group cells and generate population-averaged plots
%
%   EXAMPLE:
%     results = batch_analyze_1DRF('/path/to/protocol2/data/1DRF');
%
%   See also ANALYZE_SINGLE_EXPERIMENT, FIND_PD_AND_ORDER_IDX,
%            PLOT_POLAR_POPULATION, PLOT_FLASH_1X11_POPULATION,
%            VERIFY_LUT_DIRECTIONS, COMPUTE_BAR_SWEEP_RESPONSES,
%            FIND_PD_FROM_LUT
% ________________________________________________________________________

    %% Set defaults
    if nargin < 2, opts = struct(); end
    opts = set_batch_defaults(opts, data_root);

    %% Step 1: Load LUT (shared across experiments)
    S_lut = load(opts.lut_path, 'Tbl');
    Tbl = S_lut.Tbl;

    %% Step 2: Discover experiment folders
    d = dir(data_root);
    d = d([d.isdir]);
    d = d(~startsWith({d.name}, '.'));

    % Filter to folders that contain currentExp.mat
    valid = false(numel(d), 1);
    for i = 1:numel(d)
        valid(i) = isfile(fullfile(data_root, d(i).name, 'currentExp.mat'));
    end
    d = d(valid);
    n_exp = numel(d);

    fprintf('Found %d experiment folders in %s\n', n_exp, data_root);

    %% Step 3: Process each experiment
    results = struct([]);

    for exp_idx = 1:n_exp
        folder = d(exp_idx).name;
        exp_folder = fullfile(data_root, folder);
        fprintf('\n[%d/%d] Processing %s...\n', exp_idx, n_exp, folder);

        try
            r = process_single_cell(exp_folder, Tbl, opts);
            r.folder = folder;

            if isempty(results)
                results = r;
            else
                results(end + 1) = r; %#ok<AGROW>
            end

            fprintf('  -> %s | PD: %.0f° | %s\n', ...
                r.group, r.pd_direction, r.strain);

        catch ME
            fprintf('  ERROR: %s\n', ME.message);
            continue;
        end
    end

    fprintf('\n=== Processing Complete ===\n');
    fprintf('Total cells: %d\n', numel(results));

    % Count per group
    groups = {results.group};
    for g = ["on_control", "on_ttl", "off_control", "off_ttl"]
        fprintf('  %s: %d\n', g, sum(strcmp(groups, g)));
    end

    %% Step 4: Save results
    if opts.save_figs
        if ~isfolder(opts.save_dir)
            mkdir(opts.save_dir);
        end
        save_path = fullfile(opts.save_dir, 'batch_results.mat');
        save(save_path, 'results');
        fprintf('\nResults saved to: %s\n', save_path);
    end

    %% Step 5: Generate population plots
    fprintf('\nGenerating population plots...\n');

    for on_off_label = ["ON", "OFF"]
        if on_off_label == "ON"
            mask = [results.is_on];
        else
            mask = ~[results.is_on];
        end

        ctrl_mask = mask & ~[results.is_ttl];
        ttl_mask  = mask & [results.is_ttl];

        % Skip if no cells in this ON/OFF group
        if ~any(mask)
            fprintf('  No %s cells found — skipping plots.\n', on_off_label);
            continue;
        end

        % --- Polar tuning curves ---
        aligned_ctrl = {results(ctrl_mask).max_v_aligned};
        aligned_ttl  = {results(ttl_mask).max_v_aligned};

        polar_opts.stat_method = opts.stat_method;
        polar_title = sprintf('%s Cells — PD-Aligned Polar Tuning', on_off_label);
        fig_polar = plot_polar_population(aligned_ctrl, aligned_ttl, ...
            polar_title, polar_opts);

        % --- PD-ND bar flash traces ---
        pd_traces_ctrl = {results(ctrl_mask).pd_flash_bl};
        pd_traces_ttl  = {results(ttl_mask).pd_flash_bl};

        flash_opts.y_limits    = opts.flash_ylim;
        flash_opts.plot_type   = 'pd_nd';
        flash_opts.stat_method = opts.stat_method;
        pd_flash_title = sprintf('%s Cells — PD-ND Bar Flash', on_off_label);
        fig_pd = plot_flash_1x11_population(...
            pd_traces_ctrl, pd_traces_ttl, pd_flash_title, flash_opts);

        % --- Orthogonal bar flash traces ---
        ortho_traces_ctrl = {results(ctrl_mask).ortho_flash_bl};
        ortho_traces_ttl  = {results(ttl_mask).ortho_flash_bl};

        flash_opts.fig_position = [50 100 1800 300];
        flash_opts.plot_type    = 'orthogonal';
        ortho_flash_title = sprintf('%s Cells — Orthogonal Bar Flash', on_off_label);
        fig_ortho = plot_flash_1x11_population(...
            ortho_traces_ctrl, ortho_traces_ttl, ortho_flash_title, flash_opts);

        % --- Common aligned plot options ---
        aligned_flash_opts.y_limits          = opts.flash_ylim;
        aligned_flash_opts.stat_method       = opts.stat_method;
        aligned_flash_opts.show_n_per_pos    = true;
        aligned_flash_opts.require_both_n2   = true;
        aligned_flash_opts.stim_onset_sample  = 5001;
        aligned_flash_opts.stim_offset_sample = 5801;
        aligned_flash_opts.resp_end_sample    = 6551;
        aligned_flash_opts.show_ordinal_ranks = true;
        aligned_flash_opts.col_labels        = arrayfun(@(x) sprintf('%+d', x), ...
            -5:5, 'UniformOutput', false);

        % --- M2-aligned PD-ND bar flash traces ---
        pd_aligned_ctrl = {results(ctrl_mask).pd_flash_peak_aligned};
        pd_aligned_ttl  = {results(ttl_mask).pd_flash_peak_aligned};

        m2_flash_opts = aligned_flash_opts;
        m2_flash_opts.plot_type       = 'pd_nd';
        m2_flash_opts.col_labels{6}   = '0 (M2 peak)';

        pd_aligned_title = sprintf('%s Cells — M2-Aligned PD-ND Bar Flash', on_off_label);
        fig_pd_aligned = plot_flash_1x11_population(...
            pd_aligned_ctrl, pd_aligned_ttl, pd_aligned_title, m2_flash_opts);

        % --- M2-aligned orthogonal bar flash traces ---
        ortho_aligned_ctrl = {results(ctrl_mask).ortho_flash_peak_aligned};
        ortho_aligned_ttl  = {results(ttl_mask).ortho_flash_peak_aligned};

        m2_flash_opts.plot_type    = 'orthogonal';
        m2_flash_opts.fig_position = [50 100 1800 300];
        ortho_aligned_title = sprintf('%s Cells — M2-Aligned Orthogonal Bar Flash', on_off_label);
        fig_ortho_aligned = plot_flash_1x11_population(...
            ortho_aligned_ctrl, ortho_aligned_ttl, ortho_aligned_title, m2_flash_opts);

        % --- M5-aligned PD-ND bar flash traces ---
        pd_m5_ctrl = {results(ctrl_mask).pd_flash_m5_aligned};
        pd_m5_ttl  = {results(ttl_mask).pd_flash_m5_aligned};

        m5_flash_opts = aligned_flash_opts;
        m5_flash_opts.plot_type       = 'pd_nd';
        m5_flash_opts.col_labels{6}   = '0 (M5 centroid)';

        pd_m5_title = sprintf('%s Cells — M5-Aligned PD-ND Bar Flash', on_off_label);
        fig_pd_m5 = plot_flash_1x11_population(...
            pd_m5_ctrl, pd_m5_ttl, pd_m5_title, m5_flash_opts);

        % --- M5-aligned orthogonal bar flash traces ---
        ortho_m5_ctrl = {results(ctrl_mask).ortho_flash_m5_aligned};
        ortho_m5_ttl  = {results(ttl_mask).ortho_flash_m5_aligned};

        m5_flash_opts.plot_type    = 'orthogonal';
        m5_flash_opts.fig_position = [50 100 1800 300];
        ortho_m5_title = sprintf('%s Cells — M5-Aligned Orthogonal Bar Flash', on_off_label);
        fig_ortho_m5 = plot_flash_1x11_population(...
            ortho_m5_ctrl, ortho_m5_ttl, ortho_m5_title, m5_flash_opts);

        % Save figures
        if opts.save_figs
            save_population_figures(opts.save_dir, on_off_label, ...
                opts.stat_method, fig_polar, fig_pd, fig_ortho);
            save_aligned_figures(opts.save_dir, on_off_label, ...
                opts.stat_method, fig_pd_aligned, fig_ortho_aligned, ...
                fig_pd_m5, fig_ortho_m5);
        end
    end

    fprintf('\nDone.\n');

end


%% ========================= Core Processing ============================

function r = process_single_cell(exp_folder, Tbl, opts)
% PROCESS_SINGLE_CELL  Run the full per-cell pipeline and return results.

    % Save/restore working directory (load_protocol2_data uses cd)
    orig_dir = pwd;
    cleanup = onCleanup(@() cd(orig_dir));

    % Load data
    [date_str, ~, Log, ~, ~] = load_protocol2_data(exp_folder);

    f_data   = Log.ADC.Volts(1, :);
    v_data   = Log.ADC.Volts(2, :) * 10;
    median_v = median(v_data);

    ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
        'pattern_order', 'func_order', 'metadata');
    metadata = ce.metadata;

    % Classify cell
    r.date_str = date_str;
    r.strain   = metadata.Strain;
    r.frame    = metadata.Frame;
    r.is_on    = metadata.Frame > opts.on_threshold;
    r.is_ttl   = contains(metadata.Strain, 'ttl');

    if r.is_on && ~r.is_ttl
        r.group = 'on_control';
    elseif r.is_on && r.is_ttl
        r.group = 'on_ttl';
    elseif ~r.is_on && ~r.is_ttl
        r.group = 'off_control';
    else
        r.group = 'off_ttl';
    end

    % Parse bar sweep data and verify LUT
    bar_data = parse_bar_data(f_data, v_data);

    [lut_directions, lut_orientations, lut_patterns, lut_functions] = ...
        verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, opts.plot_order);

    % Compute bar sweep responses
    sweep_opts.baseline_range = opts.baseline_range;
    sweep_opts.stim_trim_end  = opts.stim_trim_end;
    sweep_opts.percentile     = opts.percentile;
    max_v = compute_bar_sweep_responses(bar_data, opts.plot_order, sweep_opts);

    % PD-align the tuning curve for population averaging
    % find_PD_and_order_idx expects data sorted by angle (0, 22.5, ..., 337.5)
    lut_dirs_ordered = lut_directions(opts.plot_order);
    [~, sort_idx] = sort(lut_dirs_ordered);
    max_v_sorted = max_v(sort_idx);
    max_v_polar = [max_v_sorted; max_v_sorted(1)];  % 17x1 for circular

    % median_voltage = 0 because max_v is already baseline-relative
    [d_aligned, ~, ~, ~, dir_fwhm, dir_cv, ~, dir_kappa] = ...
        find_PD_and_order_idx(max_v_polar, 0);
    r.max_v_aligned = d_aligned;

    % --- Direction selectivity metrics ---
    [r.sym_ratio, r.dsi_vector, r.dsi_pdnd, ~] = ...
        compute_bar_response_metrics(d_aligned);
    r.dir_tuning_fwhm = dir_fwhm;    % direction tuning width (degrees)
    r.dir_tuning_cv   = dir_cv;       % circular variance (0=sharp, 1=broad)
    r.dir_tuning_kappa = dir_kappa;   % von Mises concentration

    % --- Tuning aspect ratio: PD response / mean(ortho_cw, ortho_ccw) ---
    r.tuning_aspect_ratio = compute_tuning_aspect_ratio(d_aligned);

    % Find PD and map to bar flash columns
    pd_info = find_pd_from_lut(max_v, lut_directions, lut_orientations, ...
        lut_patterns, lut_functions, opts.plot_order, Tbl, opts.pattern_offset);

    r.pd_direction  = pd_info.pd_direction;
    r.pd_orientation = pd_info.pd_orientation;

    % Parse bar flash data (both speeds)
    prop_int = 0.5;  % matches old hardcoded gap_between_flashes = 5000
    [data_slow_bf, ~, mean_slow_bf, mean_fast_bf] = parse_bar_flash_data(f_data, v_data, prop_int);

    % Extract baseline-subtracted mean flash traces (11 x N_timepoints)
    bl_samples = opts.flash_baseline;

    r.pd_flash_bl    = extract_flash_traces(mean_slow_bf, ...
        pd_info.bar_flash_col, pd_info.pos_order, bl_samples);
    r.ortho_flash_bl = extract_flash_traces(mean_slow_bf, ...
        pd_info.ortho_flash_col, pd_info.ortho_pos_order, bl_samples);

    % --- RF metrics: peak position, amplitudes, bump width ---
    rf = compute_peak_metrics(mean_slow_bf, pd_info, bl_samples);
    r.peak_pos        = rf.peak_pos;
    r.peak_amplitudes = rf.peak_amplitudes;
    r.bump_width      = rf.bump_width;

    % --- M5 centroid (FWHM bump centroid) ---
    A_for_m5 = max(rf.peak_amplitudes, 0);
    m5 = compute_m5_centroid(A_for_m5);
    r.centroid_m5         = m5.centroid;       % fractional (e.g. 5.55)
    r.centroid_m5_rounded = m5.centroid_int;   % integer for alignment
    r.m5_bump_range       = m5.bump_range;     % [left, right]

    % --- Orthogonal RF bump width ---
    rf_ortho = compute_peak_metrics_for_col(mean_slow_bf, ...
        pd_info.ortho_flash_col, pd_info.ortho_pos_order, bl_samples);
    r.ortho_peak_pos        = rf_ortho.peak_pos;
    r.ortho_bump_width      = rf_ortho.bump_width;
    r.ortho_peak_amplitudes = rf_ortho.peak_amplitudes;

    % --- Voltage metrics ---
    grey_mask = f_data == 0;
    stim_mask = f_data ~= 0;
    r.rest_voltage_median = median(v_data(grey_mask));
    r.rest_voltage_mean   = mean(v_data(grey_mask));
    r.stim_voltage_median = median(v_data(stim_mask));
    r.stim_voltage_mean   = mean(v_data(stim_mask));

    % --- Peak-aligned flash traces (M2: shift so peak_pos -> row 6) ---
    r.pd_flash_peak_aligned    = reindex_to_peak(r.pd_flash_bl, r.peak_pos);
    r.ortho_flash_peak_aligned = reindex_to_peak(r.ortho_flash_bl, r.ortho_peak_pos);

    % --- M5-aligned flash traces (shift so centroid_m5_rounded -> row 6) ---
    r.pd_flash_m5_aligned    = reindex_to_peak(r.pd_flash_bl, r.centroid_m5_rounded);
    % Ortho M5: compute ortho-specific centroid for alignment
    A_ortho_m5 = max(rf_ortho.peak_amplitudes, 0);
    m5_ortho = compute_m5_centroid(A_ortho_m5);
    r.ortho_centroid_m5         = m5_ortho.centroid;
    r.ortho_centroid_m5_rounded = m5_ortho.centroid_int;
    r.ortho_flash_m5_aligned = reindex_to_peak(r.ortho_flash_bl, r.ortho_centroid_m5_rounded);

    % --- M6-aligned flash traces (68%-area centroid -> row 6) ---
    A_pd_m6 = max(r.peak_amplitudes, 0);
    m6_pd = compute_m6_centroid(A_pd_m6, 0.68);
    r.centroid_m6         = m6_pd.centroid;
    r.centroid_m6_rounded = m6_pd.centroid_int;
    r.m6_bump_range       = m6_pd.bump_range;
    r.m6_bump_width       = m6_pd.bump_width;
    r.pd_flash_m6_aligned = reindex_to_peak(r.pd_flash_bl, r.centroid_m6_rounded);

    A_ortho_m6 = max(rf_ortho.peak_amplitudes, 0);
    m6_ortho = compute_m6_centroid(A_ortho_m6, 0.68);
    r.ortho_centroid_m6         = m6_ortho.centroid;
    r.ortho_centroid_m6_rounded = m6_ortho.centroid_int;
    r.ortho_m6_bump_range       = m6_ortho.bump_range;
    r.ortho_m6_bump_width       = m6_ortho.bump_width;
    r.ortho_flash_m6_aligned = reindex_to_peak(r.ortho_flash_bl, r.ortho_centroid_m6_rounded);

    % --- Temporal metrics (Gruntman-style) ---
    r.temporal_metrics    = extract_temporal_metrics(r.pd_flash_bl, r.peak_pos, 10000);
    r.temporal_metrics_m5 = extract_temporal_metrics(r.pd_flash_bl, r.centroid_m5_rounded, 10000);

    % ===================== Fast flash data (14ms) ============================
    bl_samples_fast = 1:2500;  % prop_int=0.5 → 2500 sample baseline

    % Baseline-subtracted fast flash traces (11 x N)
    r.pd_flash_fast_bl    = extract_flash_traces(mean_fast_bf, ...
        pd_info.bar_flash_col, pd_info.pos_order, bl_samples_fast);
    r.ortho_flash_fast_bl = extract_flash_traces(mean_fast_bf, ...
        pd_info.ortho_flash_col, pd_info.ortho_pos_order, bl_samples_fast);

    % Fast RF metrics (peak amplitudes with fast timing)
    rf_fast = compute_peak_metrics_fast(mean_fast_bf, pd_info, bl_samples_fast);
    r.fast_peak_pos        = rf_fast.peak_pos;
    r.fast_peak_amplitudes = rf_fast.peak_amplitudes;
    r.fast_bump_width      = rf_fast.bump_width;

    rf_fast_ortho = compute_peak_metrics_fast_col(mean_fast_bf, ...
        pd_info.ortho_flash_col, pd_info.ortho_pos_order, bl_samples_fast);
    r.fast_ortho_peak_pos        = rf_fast_ortho.peak_pos;
    r.fast_ortho_peak_amplitudes = rf_fast_ortho.peak_amplitudes;
    r.fast_ortho_bump_width      = rf_fast_ortho.bump_width;

    % --- Fast-M6: M6 computed from fast flash amplitudes ---
    A_fast_pd = max(rf_fast.peak_amplitudes, 0);
    m6_fast_pd = compute_m6_centroid(A_fast_pd, 0.68);
    r.fast_centroid_m6         = m6_fast_pd.centroid;
    r.fast_centroid_m6_rounded = m6_fast_pd.centroid_int;
    r.fast_m6_bump_range       = m6_fast_pd.bump_range;
    r.fast_m6_bump_width       = m6_fast_pd.bump_width;

    A_fast_ortho = max(rf_fast_ortho.peak_amplitudes, 0);
    m6_fast_ortho = compute_m6_centroid(A_fast_ortho, 0.68);
    r.fast_ortho_centroid_m6         = m6_fast_ortho.centroid;
    r.fast_ortho_centroid_m6_rounded = m6_fast_ortho.centroid_int;

    % --- Combined-M6: average of normalized slow + fast profiles ---
    A_slow_norm = A_pd_m6 / max(max(A_pd_m6), 1e-6);
    A_fast_norm = A_fast_pd / max(max(A_fast_pd), 1e-6);
    A_comb = (A_slow_norm + A_fast_norm) / 2;
    m6_comb_pd = compute_m6_centroid(A_comb, 0.68);
    r.comb_centroid_m6         = m6_comb_pd.centroid;
    r.comb_centroid_m6_rounded = m6_comb_pd.centroid_int;

    A_slow_ortho_norm = max(rf_ortho.peak_amplitudes, 0);
    A_slow_ortho_norm = A_slow_ortho_norm / max(max(A_slow_ortho_norm), 1e-6);
    A_fast_ortho_norm = A_fast_ortho / max(max(A_fast_ortho), 1e-6);
    A_comb_ortho = (A_slow_ortho_norm + A_fast_ortho_norm) / 2;
    m6_comb_ortho = compute_m6_centroid(A_comb_ortho, 0.68);
    r.comb_ortho_centroid_m6         = m6_comb_ortho.centroid;
    r.comb_ortho_centroid_m6_rounded = m6_comb_ortho.centroid_int;

    % --- Fast-M6 aligned traces ---
    r.pd_flash_fast_m6fast_aligned    = reindex_to_peak(r.pd_flash_fast_bl, r.fast_centroid_m6_rounded);
    r.ortho_flash_fast_m6fast_aligned = reindex_to_peak(r.ortho_flash_fast_bl, r.fast_ortho_centroid_m6_rounded);
    % Also align slow traces with fast-M6 center (for comparison)
    r.pd_flash_slow_m6fast_aligned    = reindex_to_peak(r.pd_flash_bl, r.fast_centroid_m6_rounded);
    r.ortho_flash_slow_m6fast_aligned = reindex_to_peak(r.ortho_flash_bl, r.fast_ortho_centroid_m6_rounded);

    % --- Combined-M6 aligned traces ---
    r.pd_flash_fast_m6comb_aligned    = reindex_to_peak(r.pd_flash_fast_bl, r.comb_centroid_m6_rounded);
    r.ortho_flash_fast_m6comb_aligned = reindex_to_peak(r.ortho_flash_fast_bl, r.comb_ortho_centroid_m6_rounded);
    r.pd_flash_slow_m6comb_aligned    = reindex_to_peak(r.pd_flash_bl, r.comb_centroid_m6_rounded);
    r.ortho_flash_slow_m6comb_aligned = reindex_to_peak(r.ortho_flash_bl, r.comb_ortho_centroid_m6_rounded);

end


function traces = extract_flash_traces(mean_slow_bf, flash_col, pos_order, bl_samples)
% EXTRACT_FLASH_TRACES  Get baseline-subtracted mean traces for one orientation.
%
%   Returns an 11 x N matrix (positions x timepoints), ordered ND to PD.
%   All traces are truncated to the minimum length across positions to
%   handle minor length variations between flash stimuli.

    n_pos = 11;

    % First pass: find minimum trace length across all positions
    min_pts = Inf;
    for i = 1:n_pos
        ts = mean_slow_bf{i, flash_col};
        if ~isempty(ts)
            min_pts = min(min_pts, numel(ts));
        end
    end

    if isinf(min_pts)
        traces = [];
        return;
    end

    traces = NaN(n_pos, min_pts);

    for pos_idx = 1:n_pos
        flash_pos = pos_order(pos_idx);
        ts = mean_slow_bf{flash_pos, flash_col};
        if ~isempty(ts)
            ts_trunc = ts(1:min_pts);
            bl = mean(ts_trunc(bl_samples(bl_samples <= min_pts)));
            traces(pos_idx, :) = ts_trunc(:)' - bl;
        end
    end

end


%% ========================= Default Options ============================

function opts = set_batch_defaults(opts, data_root)
% SET_BATCH_DEFAULTS  Fill in default values for batch analysis options.

    % LUT path: same directory as this script
    if ~isfield(opts, 'lut_path')
        script_dir = fileparts(mfilename('fullpath'));
        opts.lut_path = fullfile(script_dir, 'bar_lut.mat');
    end
    if ~isfield(opts, 'plot_order')
        opts.plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
    end
    if ~isfield(opts, 'baseline_range')
        opts.baseline_range = [1000 9000];
    end
    if ~isfield(opts, 'stim_trim_end')
        opts.stim_trim_end = 7000;
    end
    if ~isfield(opts, 'percentile')
        opts.percentile = 98;
    end
    if ~isfield(opts, 'flash_baseline')
        opts.flash_baseline = 1:5000;
    end
    if ~isfield(opts, 'flash_ylim')
        opts.flash_ylim = [-15 35];
    end
    if ~isfield(opts, 'pattern_offset')
        opts.pattern_offset = 2;
    end
    if ~isfield(opts, 'on_threshold')
        opts.on_threshold = 129;
    end
    if ~isfield(opts, 'save_figs')
        opts.save_figs = true;
    end
    if ~isfield(opts, 'save_dir')
        opts.save_dir = fullfile(data_root, 'population_results');
    end
    if ~isfield(opts, 'stat_method')
        opts.stat_method = 'mean_sem';
    end

end


%% ========================= RF Metrics ================================

function rf = compute_peak_metrics(mean_slow_bf, pd_info, bl_samples)
% COMPUTE_PEAK_METRICS  Extract peak position, amplitudes, and bump width.
%   Uses the same amplitude extraction as compute_rf_centroid in
%   analyze_single_experiment_mr.m: 99.5th percentile in the response
%   window (samples 5001:6551), baseline-subtracted, negatives -> 0.

    n_pos = 11;
    col = pd_info.bar_flash_col;
    pos_order = pd_info.pos_order;
    peak_amplitudes = zeros(1, n_pos);

    resp_start = 5001;
    resp_end   = 5801 + 750;  % = 6551

    for pos_idx = 1:n_pos
        flash_pos = pos_order(pos_idx);
        ts_mean = mean_slow_bf{flash_pos, col};
        if ~isempty(ts_mean)
            bl_mean = mean(ts_mean(bl_samples(bl_samples <= numel(ts_mean))));
            win_end = min(resp_end, numel(ts_mean));
            resp_win = resp_start : win_end;
            if ~isempty(resp_win)
                peak_amplitudes(pos_idx) = prctile(ts_mean(resp_win) - bl_mean, 99.5);
            end
        end
    end

    A = max(peak_amplitudes, 0);
    [~, rf.peak_pos] = max(A);
    rf.peak_amplitudes = peak_amplitudes;
    rf.bump_width = compute_bump_width(A);

end


function rf = compute_peak_metrics_for_col(mean_slow_bf, flash_col, pos_order, bl_samples)
% COMPUTE_PEAK_METRICS_FOR_COL  Peak amplitudes and bump width for any flash column.
%   Same algorithm as compute_peak_metrics but takes column directly.

    n_pos = 11;
    peak_amplitudes = zeros(1, n_pos);
    resp_start = 5001;
    resp_end   = 5801 + 750;  % = 6551

    for pos_idx = 1:n_pos
        flash_pos = pos_order(pos_idx);
        ts_mean = mean_slow_bf{flash_pos, flash_col};
        if ~isempty(ts_mean)
            bl_mean = mean(ts_mean(bl_samples(bl_samples <= numel(ts_mean))));
            win_end = min(resp_end, numel(ts_mean));
            resp_win = resp_start : win_end;
            if ~isempty(resp_win)
                peak_amplitudes(pos_idx) = prctile(ts_mean(resp_win) - bl_mean, 99.5);
            end
        end
    end

    A = max(peak_amplitudes, 0);
    [~, rf.peak_pos] = max(A);
    rf.peak_amplitudes = peak_amplitudes;
    rf.bump_width = compute_bump_width(A);

end


function rf = compute_peak_metrics_fast(mean_fast_bf, pd_info, bl_samples)
% COMPUTE_PEAK_METRICS_FAST  Peak amplitudes for fast (14ms) flash data.
%   Same as compute_peak_metrics but with fast flash timing:
%   onset at sample 2501, response window 2501:4051.

    n_pos = 11;
    col = pd_info.bar_flash_col;
    pos_order = pd_info.pos_order;
    peak_amplitudes = zeros(1, n_pos);
    resp_start = 2501;
    resp_end   = 4051;  % onset + 1550

    for pos_idx = 1:n_pos
        flash_pos = pos_order(pos_idx);
        ts_mean = mean_fast_bf{flash_pos, col};
        if ~isempty(ts_mean)
            bl_mean = mean(ts_mean(bl_samples(bl_samples <= numel(ts_mean))));
            win_end = min(resp_end, numel(ts_mean));
            resp_win = resp_start : win_end;
            if ~isempty(resp_win)
                peak_amplitudes(pos_idx) = prctile(ts_mean(resp_win) - bl_mean, 99.5);
            end
        end
    end

    A = max(peak_amplitudes, 0);
    [~, rf.peak_pos] = max(A);
    rf.peak_amplitudes = peak_amplitudes;
    rf.bump_width = compute_bump_width(A);

end


function rf = compute_peak_metrics_fast_col(mean_fast_bf, flash_col, pos_order, bl_samples)
% COMPUTE_PEAK_METRICS_FAST_COL  Fast flash peak metrics for any column.
%   Same as compute_peak_metrics_fast but takes column directly.

    n_pos = 11;
    peak_amplitudes = zeros(1, n_pos);
    resp_start = 2501;
    resp_end   = 4051;

    for pos_idx = 1:n_pos
        flash_pos = pos_order(pos_idx);
        ts_mean = mean_fast_bf{flash_pos, flash_col};
        if ~isempty(ts_mean)
            bl_mean = mean(ts_mean(bl_samples(bl_samples <= numel(ts_mean))));
            win_end = min(resp_end, numel(ts_mean));
            resp_win = resp_start : win_end;
            if ~isempty(resp_win)
                peak_amplitudes(pos_idx) = prctile(ts_mean(resp_win) - bl_mean, 99.5);
            end
        end
    end

    A = max(peak_amplitudes, 0);
    [~, rf.peak_pos] = max(A);
    rf.peak_amplitudes = peak_amplitudes;
    rf.bump_width = compute_bump_width(A);

end


function bw = compute_bump_width(A)
% COMPUTE_BUMP_WIDTH  FWHM of the RF amplitude profile.
%   Contiguous region around peak where A >= peak/2.

    [peak_val, peak_pos] = max(A);
    if peak_val <= 0
        bw = NaN;
        return;
    end
    threshold = peak_val / 2;
    left = peak_pos;
    while left > 1 && A(left - 1) >= threshold
        left = left - 1;
    end
    right = peak_pos;
    while right < numel(A) && A(right + 1) >= threshold
        right = right + 1;
    end
    bw = right - left + 1;

end


function m5 = compute_m5_centroid(A)
% COMPUTE_M5_CENTROID  FWHM bump centroid (M5 method).
%   Returns fractional centroid, rounded integer for alignment, and bump range.
%   A is a 1x11 amplitude vector with negatives already thresholded to 0.
%   Algorithm matches compute_rf_centroid in analyze_single_experiment_mr.m.

    n_pos = numel(A);
    [peak_val, peak_pos] = max(A);

    if peak_val <= 0
        m5.centroid     = 6;
        m5.centroid_int = 6;
        m5.bump_range   = [1, n_pos];
        m5.bump_width   = n_pos;
        return;
    end

    threshold = peak_val / 2;

    left = peak_pos;
    while left > 1 && A(left - 1) >= threshold
        left = left - 1;
    end
    right = peak_pos;
    while right < n_pos && A(right + 1) >= threshold
        right = right + 1;
    end

    bump = left:right;
    A_bump = A(bump);
    m5.centroid     = sum(A_bump .* bump) / sum(A_bump);
    m5.centroid_int = round(m5.centroid);
    m5.centroid_int = max(1, min(n_pos, m5.centroid_int));  % clamp to [1, 11]
    m5.bump_range   = [left, right];
    m5.bump_width   = right - left + 1;

end


function m6 = compute_m6_centroid(A, area_fraction)
% COMPUTE_M6_CENTROID  68%-area bump centroid (M6 method).
%   Greedy expansion from peak: at each step, add the adjacent position
%   (left-1 or right+1) with larger amplitude, until contiguous region
%   contains >= area_fraction of total positive area.
%   Returns fractional centroid, rounded integer, bump range, and width.

    n_pos = numel(A);
    [peak_val, peak_pos] = max(A);
    total_area = sum(A);

    if peak_val <= 0 || total_area <= 0
        m6.centroid     = 6;
        m6.centroid_int = 6;
        m6.bump_range   = [1, n_pos];
        m6.bump_width   = n_pos;
        return;
    end

    target = area_fraction * total_area;

    left  = peak_pos;
    right = peak_pos;
    cum   = A(peak_pos);

    while cum < target && (left > 1 || right < n_pos)
        can_left  = (left > 1);
        can_right = (right < n_pos);

        if can_left && can_right
            if A(left - 1) >= A(right + 1)
                left = left - 1;
                cum  = cum + A(left);
            else
                right = right + 1;
                cum   = cum + A(right);
            end
        elseif can_left
            left = left - 1;
            cum  = cum + A(left);
        else
            right = right + 1;
            cum   = cum + A(right);
        end
    end

    bump   = left:right;
    A_bump = A(bump);
    m6.centroid     = sum(A_bump .* bump) / sum(A_bump);
    m6.centroid_int = round(m6.centroid);
    m6.centroid_int = max(1, min(n_pos, m6.centroid_int));
    m6.bump_range   = [left, right];
    m6.bump_width   = right - left + 1;

end


function aligned = reindex_to_peak(traces_11xN, peak_pos)
% REINDEX_TO_PEAK  Shift traces so peak_pos maps to row 6 (center of 11).
%   Positions that fall outside 1-11 become NaN rows.

    n_pos = 11;
    center = 6;
    shift = center - peak_pos;
    N = size(traces_11xN, 2);
    aligned = NaN(n_pos, N);

    for i = 1:n_pos
        src = i - shift;
        if src >= 1 && src <= n_pos
            aligned(i, :) = traces_11xN(src, :);
        end
    end

end


%% ========================= Figure Saving ==============================

function save_population_figures(save_dir, on_off_label, stat_method, ...
    fig_polar, fig_pd, fig_ortho)
% SAVE_POPULATION_FIGURES  Export population figures as 300 dpi PDFs.

    if ~isfolder(save_dir)
        mkdir(save_dir);
    end

    export_opts = {'ContentType', 'image', 'Resolution', 300};
    prefix = lower(char(on_off_label));

    exportgraphics(fig_polar, ...
        fullfile(save_dir, sprintf('%s_polar_population_%s.pdf', prefix, stat_method)), ...
        export_opts{:});
    exportgraphics(fig_pd, ...
        fullfile(save_dir, sprintf('%s_pd_nd_flash_population_%s.pdf', prefix, stat_method)), ...
        export_opts{:});
    exportgraphics(fig_ortho, ...
        fullfile(save_dir, sprintf('%s_ortho_flash_population_%s.pdf', prefix, stat_method)), ...
        export_opts{:});

    fprintf('  %s figures saved to: %s\n', on_off_label, save_dir);

end


function save_aligned_figures(save_dir, on_off_label, stat_method, ...
    fig_pd_aligned, fig_ortho_aligned, fig_pd_m5, fig_ortho_m5)
% SAVE_ALIGNED_FIGURES  Export M2 and M5 aligned population figures as 300 dpi PDFs.

    if ~isfolder(save_dir)
        mkdir(save_dir);
    end

    export_opts = {'ContentType', 'image', 'Resolution', 300};
    prefix = lower(char(on_off_label));

    % M2-aligned
    exportgraphics(fig_pd_aligned, ...
        fullfile(save_dir, sprintf('%s_pd_nd_M2_aligned_%s.pdf', prefix, stat_method)), ...
        export_opts{:});
    exportgraphics(fig_ortho_aligned, ...
        fullfile(save_dir, sprintf('%s_ortho_M2_aligned_%s.pdf', prefix, stat_method)), ...
        export_opts{:});

    % M5-aligned
    if nargin >= 6 && ~isempty(fig_pd_m5)
        exportgraphics(fig_pd_m5, ...
            fullfile(save_dir, sprintf('%s_pd_nd_M5_aligned_%s.pdf', prefix, stat_method)), ...
            export_opts{:});
    end
    if nargin >= 7 && ~isempty(fig_ortho_m5)
        exportgraphics(fig_ortho_m5, ...
            fullfile(save_dir, sprintf('%s_ortho_M5_aligned_%s.pdf', prefix, stat_method)), ...
            export_opts{:});
    end

    fprintf('  %s M2+M5 aligned figures saved to: %s\n', on_off_label, save_dir);

end


function ar = compute_tuning_aspect_ratio(d_aligned)
% COMPUTE_TUNING_ASPECT_RATIO  PD response / mean of two orthogonal responses.
%   d_aligned is 16x2 [angle_rad, response] with PD aligned to pi/2.
%   Orthogonal directions are at 0 and pi (90 deg CW and CCW from PD).

    angles = d_aligned(:, 1);
    resps  = d_aligned(:, 2);

    % Find indices closest to PD (pi/2), ortho_cw (0), ortho_ccw (pi)
    [~, pd_idx]     = min(abs(angles - pi/2));
    [~, ortho1_idx] = min(abs(angles - 0));
    [~, ortho2_idx] = min(abs(angles - pi));

    pd_resp    = resps(pd_idx);
    ortho_mean = mean([resps(ortho1_idx), resps(ortho2_idx)]);

    if ortho_mean > 0
        ar = pd_resp / ortho_mean;
    else
        ar = NaN;
    end

end

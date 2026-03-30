function cell_info = analyze_single_experiment_mr(exp_folder, opts)
% ANALYZE_SINGLE_EXPERIMENT_MR  Enhanced bar sweep and bar flash analysis.
%
%   ANALYZE_SINGLE_EXPERIMENT_MR(EXP_FOLDER) runs the full analysis pipeline
%   for a single Protocol 2 experiment, identical to ANALYZE_SINGLE_EXPERIMENT.
%
%   CELL_INFO = ANALYZE_SINGLE_EXPERIMENT_MR(EXP_FOLDER, OPTS) returns a
%   structure with per-cell metadata and RF centering metrics.
%
%   In addition to the options accepted by ANALYZE_SINGLE_EXPERIMENT, this
%   variant supports:
%
%     opts.visualize_everything  (default: false)
%       When true, adds diagnostic enhancements:
%         1. A voltage histogram figure showing the full-recording and
%            stimulus-only voltage distributions.
%         2. Vertical lines on all bar flash figures marking flash ON/OFF.
%         3. A thicker border and annotation on the PD row of the heatmap.
%         4. A red diametric line on the polar plot showing PD-ND axis.
%         5. Relabels PD-ND endpoints to "Leading" / "Trailing".
%         6. RF centroid marker and degree-labeled axes on PD-ND figure.
%
%   INPUTS:
%     exp_folder - Full path to an experiment directory
%                  e.g. '/path/to/data/1DRF/2025_11_10_10_17'
%     opts       - (Optional) structure with fields:
%                    .lut_path              - Path to bar_lut.mat
%                    .plot_order            - 1x16 data-row-to-subplot mapping
%                    .baseline_range        - [start end] samples for bar sweep
%                    .stim_trim_end         - Samples to trim from end of sweep
%                    .percentile            - Percentile for peak detection
%                    .flash_baseline        - Sample indices for flash baseline
%                    .flash_ylim            - [ymin ymax] for 1x11 plots
%                    .pattern_offset        - Experiment-to-fullfield offset
%                    .save_figs             - Boolean, save figures as PDF
%                    .save_dir              - Output directory for PDFs
%                    .visualize_everything  - Boolean, enable all enhancements
%
%   OUTPUTS:
%     cell_info - (Optional) Structure with fields:
%       .exp_folder      - Path to experiment
%       .date_str        - Experiment date string
%       .strain          - Fly strain
%       .frame           - Recording frame
%       .on_off          - 'ON' or 'OFF' classification
%       .pd_direction    - Preferred direction in degrees
%       .pd_orientation  - Bar orientation at PD in degrees
%       .centroid_idx    - RF centroid as fractional position (1-11)
%       .centroid_deg    - RF centroid offset from center in degrees
%       .peak_amplitudes - 1x11 peak amplitude profile
%
%     Figures:
%       Figure 1: Slow bar sweep polar timeseries with vector sum arrow
%       Figure 2: Full 8x11 bar flash heatmap (all orientations)
%       Figure 3: 1x11 bar flash subplots along PD-ND axis
%       Figure 4: 1x11 bar flash subplots along orthogonal axis
%       Figure 5: (visualize_everything only) Voltage histograms
%
%   EXAMPLE:
%     opts.visualize_everything = true;
%     opts.save_figs = true;
%     opts.save_dir  = fullfile(exp_folder, 'analysis_output_mr');
%     cell_info = analyze_single_experiment_mr(exp_folder, opts);
%
%   See also ANALYZE_SINGLE_EXPERIMENT, PLOT_VOLTAGE_HISTOGRAMS,
%            ADD_STIM_TIMING_LINES, LOAD_PROTOCOL2_DATA, PARSE_BAR_DATA,
%            PARSE_BAR_FLASH_DATA, PLOT_SLOW_BAR_SWEEP_POLAR,
%            PLOT_BAR_FLASH_HEATMAP, PLOT_BAR_FLASH_1X11
% ________________________________________________________________________

    %% Set defaults
    if nargin < 2, opts = struct(); end
    opts = set_default_opts(opts, exp_folder);

    %% Step 1: Load data
    S_lut = load(opts.lut_path, 'Tbl');
    Tbl = S_lut.Tbl;

    % Restore working directory after load_protocol2_data (which uses cd)
    orig_dir = pwd;
    cleanup = onCleanup(@() cd(orig_dir));

    [date_str, ~, Log, params, ~] = load_protocol2_data(exp_folder);

    f_data   = Log.ADC.Volts(1, :);       % frame position data
    v_data   = Log.ADC.Volts(2, :) * 10;  % voltage data (scaled to mV)
    median_v = median(v_data);

    ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
        'pattern_order', 'func_order', 'metadata');
    metadata = ce.metadata;

    %% Enhancement A: Voltage histogram figure
    fig_hist = [];
    if opts.visualize_everything
        hist_title = sprintf('Voltage Distribution — %s — %s — %s', ...
            strrep(date_str, '_', '-'), ...
            strrep(metadata.Strain, '_', ' '), ...
            params.on_off);
        fig_hist = plot_voltage_histograms(v_data, f_data, median_v, hist_title);
    end

    %% Step 2: Parse bar sweep data and verify LUT directions
    bar_data = parse_bar_data(f_data, v_data);

    [lut_directions, lut_orientations, lut_patterns, lut_functions] = ...
        verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, ...
        opts.plot_order);

    %% Step 3: Compute bar sweep responses and create polar figure
    sweep_opts.baseline_range = opts.baseline_range;
    sweep_opts.stim_trim_end  = opts.stim_trim_end;
    sweep_opts.percentile     = opts.percentile;
    max_v = compute_bar_sweep_responses(bar_data, opts.plot_order, sweep_opts);

    polar_title = sprintf('28 dps — %s — %s — %s', ...
        strrep(date_str, '_', '-'), ...
        strrep(metadata.Strain, '_', ' '), ...
        params.on_off);
    fig_polar = plot_slow_bar_sweep_polar(bar_data, max_v, lut_directions, ...
        opts.plot_order, median_v, polar_title);

    %% Step 4: Find preferred direction and map to bar flash columns
    %  (computed here so Enhancement D can use pd_info immediately after)
    pd_info = find_pd_from_lut(max_v, lut_directions, lut_orientations, ...
        lut_patterns, lut_functions, opts.plot_order, Tbl, opts.pattern_offset);

    % Compute direction selectivity metrics (DSI, tuning width)
    lut_dirs_ordered = lut_directions(opts.plot_order);
    [~, sort_idx] = sort(lut_dirs_ordered);
    max_v_sorted = max_v(sort_idx);
    max_v_polar_17 = [max_v_sorted; max_v_sorted(1)];  % 17x1 for circular
    [d_aligned, ~, ~, ~, dir_fwhm, ~, ~, ~] = ...
        find_PD_and_order_idx(max_v_polar_17, 0);
    [~, dsi_vector, dsi_pdnd, ~] = compute_bar_response_metrics(d_aligned);

    %% Enhancement D: PD and orthogonal direction lines on polar plot
    if opts.visualize_everything
        add_polar_direction_lines(fig_polar, ...
            pd_info.pd_direction, pd_info.ortho_orientation, max(max_v));
    end

    %% Step 5: Parse bar flash data and generate figures
    % prop_int = 0.5 matches old hardcoded gap_between_flashes = 5000 (slow)
    % and 2500 (fast). Laura's updated parse_bar_flash_data now requires this
    % argument: gap = 10000 * prop_int for slow, 5000 * prop_int for fast.
    prop_int = 0.5;
    [data_slow_bf, ~, mean_slow_bf, ~] = parse_bar_flash_data(f_data, v_data, prop_int);

    % Build orientation labels from LUT
    orient_labels = build_orient_labels(Tbl, opts.pattern_offset);

    % 8x11 heatmap
    heatmap_title = sprintf(...
        'Bar Flashes 80ms — PD: Dir %.0f° Orient %.0f° (row %d) — %s — %s', ...
        pd_info.pd_direction, pd_info.pd_orientation, pd_info.bar_flash_col, ...
        strrep(date_str, '_', '-'), strrep(metadata.Strain, '_', ' '));
    fig_heatmap = plot_bar_flash_heatmap(data_slow_bf, mean_slow_bf, ...
        median_v, pd_info.bar_flash_col, pd_info.pos_order, ...
        orient_labels, heatmap_title);

    % 1x11 PD-ND axis
    flash_opts.baseline_samples = opts.flash_baseline;
    flash_opts.y_limits         = opts.flash_ylim;

    pd_title = sprintf('Bar Flash PD-ND — Dir:%.0f° Orient:%.0f° — %s — %s', ...
        pd_info.pd_direction, pd_info.pd_orientation, ...
        strrep(date_str, '_', '-'), strrep(metadata.Strain, '_', ' '));
    fig_pd = plot_bar_flash_1x11(...
        data_slow_bf(:, pd_info.bar_flash_col, :), ...
        mean_slow_bf(:, pd_info.bar_flash_col), ...
        pd_info.pos_order, pd_title, flash_opts);

    %% Enhancement E: Relabel PD-ND endpoints to Leading / Trailing
    if opts.visualize_everything
        pd_axes = findobj(fig_pd, 'Type', 'axes');
        for k = 1:numel(pd_axes)
            t = get(pd_axes(k), 'Title');
            if strcmp(t.String, 'ND'),  t.String = 'Leading';  end
            if strcmp(t.String, 'PD'),  t.String = 'Trailing'; end
        end
    end

    %% Enhancement G: RF centroid and degree-labeled axes
    %  Compute response-weighted centroid from PD-ND bar flash amplitudes.
    %  This runs unconditionally (needed for the output struct), but the
    %  figure annotations are gated by opts.visualize_everything.
    bl_samples = opts.flash_baseline;
    [centroid_idx, centroid_deg, rf_metrics] = compute_rf_centroid( ...
        mean_slow_bf, pd_info, bl_samples);

    fprintf('\n=== RF Centroid (PD-ND axis) ===\n');
    fprintf('  Amplitudes (mV):      %s\n', ...
        mat2str(round(rf_metrics.peak_amplitudes, 1)));
    fprintf('  Peak position:        %d  (%.1f mV)\n', ...
        rf_metrics.peak_pos, rf_metrics.peak_val);
    fprintf('  Bump (FWHM):          positions %d-%d  (%d wide)\n', ...
        rf_metrics.bump_range(1), rf_metrics.bump_range(2), rf_metrics.bump_width);
    fprintf('  Bump centroid:        %.2f  (of 1-11, center=6)\n', centroid_idx);
    fprintf('  Centroid offset:      %+.1f°\n', centroid_deg);
    fprintf('  Centroid-peak delta:  %.2f positions\n', rf_metrics.centroid_peak_delta);

    if opts.visualize_everything
        annotate_centroid(fig_pd, centroid_idx, centroid_deg, rf_metrics);
    end

    % 1x11 Orthogonal axis
    flash_opts.fig_position = [50 100 1800 300];
    ortho_title = sprintf('Bar Flash Orthogonal — Orient:%.0f° — %s — %s', ...
        pd_info.ortho_orientation, ...
        strrep(date_str, '_', '-'), strrep(metadata.Strain, '_', ' '));
    fig_ortho = plot_bar_flash_1x11(...
        data_slow_bf(:, pd_info.ortho_flash_col, :), ...
        mean_slow_bf(:, pd_info.ortho_flash_col), ...
        pd_info.pos_order, ortho_title, flash_opts);

    % Compute orthogonal RF metrics (bump width)
    ortho_rf = compute_ortho_rf_metrics(mean_slow_bf, pd_info, opts.flash_baseline);
    ortho_bump_width = ortho_rf.bump_width;
    if rf_metrics.bump_width > 0 && ~isnan(ortho_bump_width) && ortho_bump_width > 0
        aspect_ratio = rf_metrics.bump_width / ortho_bump_width;
    else
        aspect_ratio = NaN;
    end

    fprintf('\n=== Direction & RF Shape Metrics ===\n');
    fprintf('  DSI (vector sum):    %.3f\n', dsi_vector);
    fprintf('  DSI (PD-ND):         %.3f\n', dsi_pdnd);
    fprintf('  Dir tuning FWHM:     %.0f°\n', dir_fwhm);
    fprintf('  PD bump width:       %d positions\n', rf_metrics.bump_width);
    fprintf('  Ortho bump width:    %d positions\n', ortho_bump_width);
    fprintf('  Aspect ratio (PD/O): %.2f\n', aspect_ratio);

    %% Annotate polar plot with DSI, ortho width, aspect ratio
    if opts.visualize_everything
        annotate_polar_metrics(fig_polar, dsi_pdnd, ortho_bump_width, aspect_ratio);
    end

    %% Enhancement F: Time scale bars
    %  (Must run BEFORE Enhancement B because querying axes positions in a
    %   tiledlayout triggers a layout recalculation that resets the Children
    %   stack ordering. By drawing scale bars first, the timing-line z-order
    %   applied in Enhancement B is preserved.)
    if opts.visualize_everything
        % Bar flash figures: trace is ~10801 samples at 10 kHz.
        % 0.5s = 5000 samples is a good scale.
        add_timebar_to_figure(fig_heatmap, 5000, '0.5 s', 'last_row_last_col');
        add_timebar_to_figure(fig_pd, 5000, '0.5 s', 'last_tile');
        add_timebar_to_figure(fig_ortho, 5000, '0.5 s', 'last_tile');

        % Polar figure: bar sweep traces are ~27000 samples at 10 kHz.
        % 1s = 10000 samples is a good scale.
        add_timebar_to_figure(fig_polar, 10000, '1 s', 'bottom_right_axes');
    end

    %% Enhancement B: Stimulus timing lines on bar flash figures
    if opts.visualize_everything
        % Compute stimulus timing from the extracted trace length
        % parse_bar_flash_data clips: v_data(st - gap : nd + gap)
        % So stimulus onset is at sample gap+1, offset at trace_len - gap
        sample_trace = mean_slow_bf{1, pd_info.bar_flash_col};
        if isempty(sample_trace)
            % Fallback: find any non-empty trace
            for ii = 1:size(mean_slow_bf, 1)
                sample_trace = mean_slow_bf{ii, pd_info.bar_flash_col};
                if ~isempty(sample_trace), break; end
            end
        end
        trace_len = numel(sample_trace);
        gap = 5000;  % gap_between_flashes for slow (80ms) bar flashes
        stim_onset  = gap + 1;          % sample 5001
        stim_offset = trace_len - gap;  % last stimulus sample

        fprintf('\n=== Stimulus timing (bar flash) ===\n');
        fprintf('  stim_onset  = sample %d (%.1f ms)\n', stim_onset, stim_onset / 10);
        fprintf('  stim_offset = sample %d (%.1f ms)\n', stim_offset, stim_offset / 10);
        fprintf('  flash duration = %d samples (%.1f ms)\n', ...
            stim_offset - stim_onset + 1, (stim_offset - stim_onset + 1) / 10);

        add_stim_timing_lines(fig_heatmap, stim_onset, stim_offset);
        add_stim_timing_lines(fig_pd, stim_onset, stim_offset);
        add_stim_timing_lines(fig_ortho, stim_onset, stim_offset);
    end

    %% Enhancement C: Highlight PD row on heatmap
    if opts.visualize_everything
        highlight_heatmap_pd_row(fig_heatmap, pd_info.bar_flash_col);
    end

    %% Step 6: Save figures
    if opts.save_figs
        save_analysis_figures(opts.save_dir, opts.visualize_everything, ...
            fig_polar, fig_heatmap, fig_pd, fig_ortho, fig_hist);
    end

    %% Output struct assembly
    if nargout > 0
        cell_info.exp_folder      = exp_folder;
        cell_info.date_str        = date_str;
        cell_info.strain          = metadata.Strain;
        cell_info.frame           = metadata.Frame;
        cell_info.on_off          = params.on_off;
        cell_info.pd_direction    = pd_info.pd_direction;
        cell_info.pd_orientation  = pd_info.pd_orientation;
        cell_info.centroid_idx    = centroid_idx;
        cell_info.centroid_deg    = centroid_deg;
        cell_info.peak_amplitudes = rf_metrics.peak_amplitudes;
        cell_info.peak_pos        = rf_metrics.peak_pos;
        cell_info.peak_val        = rf_metrics.peak_val;
        cell_info.bump_range      = rf_metrics.bump_range;
        cell_info.bump_width      = rf_metrics.bump_width;
        cell_info.centroid_peak_delta = rf_metrics.centroid_peak_delta;
        cell_info.dsi_vector      = dsi_vector;
        cell_info.dsi_pdnd        = dsi_pdnd;
        cell_info.dir_tuning_fwhm = dir_fwhm;
        cell_info.ortho_bump_width = ortho_bump_width;
        cell_info.aspect_ratio    = aspect_ratio;
    end

end


%% ========================= Local Functions ============================

function opts = set_default_opts(opts, exp_folder)
% SET_DEFAULT_OPTS  Fill in default values for any unset options.

    % Default LUT path: bar_lut.mat in the same directory as this script
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
    if ~isfield(opts, 'save_figs')
        opts.save_figs = true;
    end
    if ~isfield(opts, 'save_dir')
        opts.save_dir = fullfile(exp_folder, 'analysis_output');
    end
    if ~isfield(opts, 'visualize_everything')
        opts.visualize_everything = false;
    end

end


function orient_labels = build_orient_labels(Tbl, pattern_offset)
% BUILD_ORIENT_LABELS  Create orientation label strings from the LUT.

    orient_labels = cell(8, 1);
    for k = 1:8
        exp_pat = k + pattern_offset;
        mask = Tbl.pattern == exp_pat & Tbl.function == 3;
        orient_labels{k} = sprintf('Orient: %.0f°', Tbl.orientation(mask));
    end

end


function highlight_heatmap_pd_row(fig, ~)
% HIGHLIGHT_HEATMAP_PD_ROW  Thicken PD row borders and add arrow annotation.
%
%   HIGHLIGHT_HEATMAP_PD_ROW(FIG, ~) identifies the PD row in the 8x11
%   heatmap by looking for axes with XColor == [0.8 0 0] (set by
%   plot_bar_flash_heatmap for the PD row). Increases border width and
%   adds a left-side arrow annotation pointing to the PD row.

    all_axes = findobj(fig, 'Type', 'axes');

    % Collect positions of PD-row axes to place annotation
    pd_positions = [];

    for k = 1:numel(all_axes)
        ax = all_axes(k);
        xc = get(ax, 'XColor');

        % PD row axes have XColor = [0.8 0 0]
        if isnumeric(xc) && numel(xc) == 3 && all(abs(xc - [0.8 0 0]) < 0.01)
            set(ax, 'LineWidth', 2.5);
            pd_positions = [pd_positions; ax.Position]; %#ok<AGROW>
        end
    end

    % Add a bracket/arrow annotation to the left of the PD row
    if ~isempty(pd_positions)
        % Find leftmost and vertical extent of PD row axes
        left_edge = min(pd_positions(:, 1));
        bottom    = min(pd_positions(:, 2));
        top       = max(pd_positions(:, 2) + pd_positions(:, 4));
        mid_y     = (bottom + top) / 2;

        % Draw a small arrow pointing right at the PD row
        arrow_x = left_edge - 0.02;
        if arrow_x > 0
            annotation(fig, 'textarrow', ...
                [arrow_x, left_edge - 0.005], [mid_y, mid_y], ...
                'String', 'PD', 'FontSize', 11, 'FontWeight', 'bold', ...
                'Color', [0.8 0 0], 'HeadStyle', 'vback2', ...
                'HeadWidth', 8, 'HeadLength', 6);
        end
    end

end


function save_analysis_figures(save_dir, visualize_everything, ...
    fig_polar, fig_heatmap, fig_pd, fig_ortho, fig_hist)
% SAVE_ANALYSIS_FIGURES  Export all analysis figures as 300 dpi PDFs.

    if ~isfolder(save_dir)
        mkdir(save_dir);
    end

    export_opts = {'ContentType', 'image', 'Resolution', 300};

    exportgraphics(fig_polar, ...
        fullfile(save_dir, 'slow_bar_sweep_polar.pdf'), export_opts{:});
    exportgraphics(fig_heatmap, ...
        fullfile(save_dir, 'bar_flash_all_orientations.pdf'), export_opts{:});
    exportgraphics(fig_pd, ...
        fullfile(save_dir, 'bar_flash_PD_ND_axis.pdf'), export_opts{:});
    exportgraphics(fig_ortho, ...
        fullfile(save_dir, 'bar_flash_orthogonal_axis.pdf'), export_opts{:});

    if visualize_everything && ~isempty(fig_hist)
        exportgraphics(fig_hist, ...
            fullfile(save_dir, 'voltage_histograms.pdf'), export_opts{:});
    end

    fprintf('\nFigures saved to: %s\n', save_dir);

end


function add_polar_direction_lines(fig, pd_direction, ortho_orientation, max_rho)
% ADD_POLAR_DIRECTION_LINES  Draw PD and orthogonal lines on the polar plot.
%
%   Draws two diametric lines on the central polar plot:
%     - Red line at the PD motion direction (e.g. 225°/45°)
%     - Black line at the orthogonal motion direction (derived from
%       ortho_orientation + 90°, e.g. 135°/315°)
%   Both are drawn behind the tuning curve.
%
%   The polar plot axes represent motion direction, not bar orientation.
%   Motion is perpendicular to bar orientation, so we add 90° to convert
%   from bar orientation to one of the two possible motion directions.

    pax = findobj(fig, 'Type', 'PolarAxes');
    if isempty(pax), return; end

    r_max = max_rho * 1.05;

    hold(pax, 'on');

    % Orthogonal motion direction (black) — bar orientation + 90°
    % This is one of the two motion directions for the orthogonal bar.
    ortho_motion_deg = ortho_orientation + 90;
    ortho_rad = deg2rad(ortho_motion_deg);
    h_ortho = polarplot(pax, [ortho_rad, ortho_rad + pi], [r_max, r_max], ...
        '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

    % PD motion direction line (red)
    pd_rad = deg2rad(pd_direction);
    h_pd = polarplot(pax, [pd_rad, pd_rad + pi], [r_max, r_max], ...
        '-', 'Color', [0.8 0 0], 'LineWidth', 2);

    % Push both lines behind the tuning curve and arrow
    ch = get(pax, 'Children');
    is_guide = ismember(ch, [h_pd; h_ortho]);
    set(pax, 'Children', [ch(~is_guide); ch(is_guide)]);

end


function add_timebar_to_figure(fig, n_samples, label_str, mode)
% ADD_TIMEBAR_TO_FIGURE  Add a time scale bar to one panel of a figure.
%
%   Draws a vertical line and text label in a single representative panel
%   to indicate the time scale. No new axes are created.
%
%   INPUTS:
%     fig        - Figure handle
%     n_samples  - Length of the scale bar in samples (e.g. 5000 = 0.5s)
%     label_str  - Text label (e.g. '0.5 s')
%     mode       - Which panel to annotate:
%                    'last_tile'           - last (rightmost) tile in a 1xN layout
%                    'last_row_last_col'   - bottom-right tile in an MxN layout
%                    'bottom_right_axes'   - axes with lowest y-position, rightmost x

    all_axes = findobj(fig, 'Type', 'axes');

    % Filter out polar axes
    keep = true(size(all_axes));
    for k = 1:numel(all_axes)
        if isa(all_axes(k), 'matlab.graphics.axis.PolarAxes')
            keep(k) = false;
        end
    end
    all_axes = all_axes(keep);
    if isempty(all_axes), return; end

    % Select the target axes based on mode
    switch mode
        case 'last_tile'
            % In a tiled layout, findobj returns axes in reverse order,
            % so the first one is the last tile. Pick the one with the
            % rightmost position.
            positions = cell2mat(get(all_axes, 'Position'));
            [~, idx] = max(positions(:, 1));  % rightmost
            target_ax = all_axes(idx);

        case 'last_row_last_col'
            % Bottom-right: lowest Y and highest X
            positions = cell2mat(get(all_axes, 'Position'));
            % Find the bottom row (lowest y-position within tolerance)
            min_y = min(positions(:, 2));
            bottom_row = abs(positions(:, 2) - min_y) < 0.02;
            bottom_axes = all_axes(bottom_row);
            bottom_pos  = positions(bottom_row, :);
            [~, idx] = max(bottom_pos(:, 1));
            target_ax = bottom_axes(idx);

        case 'bottom_right_axes'
            % For the polar figure: axes positioned at various locations.
            % Pick the one in the lower-right quadrant.
            positions = cell2mat(get(all_axes, 'Position'));
            % Score: high x, low y
            scores = positions(:, 1) - positions(:, 2);
            [~, idx] = max(scores);
            target_ax = all_axes(idx);

        otherwise
            return;
    end

    % Draw the scale bar in the target axes
    hold(target_ax, 'on');
    xl = get(target_ax, 'XLim');
    yl = get(target_ax, 'YLim');

    % Position: right side of the axes, near the bottom.
    % The bar must fit entirely within the existing x-limits.
    x_range = xl(2) - xl(1);
    margin  = x_range * 0.02;               % 2% padding from right edge
    bar_x   = xl(2) - n_samples - margin;    % bar ends at xl(2) - margin
    bar_y   = yl(1) + (yl(2) - yl(1)) * 0.08;

    % Horizontal bar showing time span
    plot(target_ax, [bar_x, bar_x + n_samples], [bar_y, bar_y], ...
        '-k', 'LineWidth', 2, 'Clipping', 'on');

    % Small vertical end-caps
    cap_h = (yl(2) - yl(1)) * 0.04;
    plot(target_ax, [bar_x, bar_x], [bar_y - cap_h, bar_y + cap_h], ...
        '-k', 'LineWidth', 1.5, 'Clipping', 'on');
    plot(target_ax, [bar_x + n_samples, bar_x + n_samples], ...
        [bar_y - cap_h, bar_y + cap_h], '-k', 'LineWidth', 1.5, ...
        'Clipping', 'on');

    % Label sitting directly on the bar (white background to avoid clutter)
    text(target_ax, bar_x + n_samples / 2, bar_y, ...
        label_str, 'FontSize', 8, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'BackgroundColor', 'w', 'Margin', 1);

    % Ensure axes limits are unchanged (don't auto-expand for the bar)
    set(target_ax, 'XLim', xl, 'YLim', yl);

end


function [centroid_idx, centroid_deg, rf_metrics] = compute_rf_centroid( ...
    mean_slow_bf, pd_info, bl_samples)
% COMPUTE_RF_CENTROID  Bump-based centroid of the PD-ND bar flash profile.
%
%   [CENTROID_IDX, CENTROID_DEG, RF_METRICS] = COMPUTE_RF_CENTROID(
%       MEAN_SLOW_BF, PD_INFO, BL_SAMPLES)
%   measures the response amplitude at each of the 11 spatial positions
%   along the PD-ND axis, finds the contiguous depolarising bump (FWHM
%   region around the peak), and computes the response-weighted centroid
%   over the bump positions only.
%
%   Amplitude is the 99.5th percentile of the baseline-subtracted trace
%   in the response window (flash onset through offset + 75 ms), following
%   Gruntman et al. (eLife 2019). Negative amplitudes are thresholded to
%   zero so that only depolarising responses contribute.
%
%   INPUTS:
%     mean_slow_bf - cell array from parse_bar_flash_data (mean traces)
%     pd_info      - struct from find_pd_from_lut (needs .bar_flash_col, .pos_order)
%     bl_samples   - sample indices for baseline (e.g. 1:5000)
%
%   OUTPUTS:
%     centroid_idx - Fractional position (1-11), 6 = geometric center
%     centroid_deg - Offset from center in degrees (2.5°/position)
%     rf_metrics   - Structure with fields:
%       .peak_amplitudes    - 1x11 response amplitudes (mV)
%       .peak_pos           - Position of maximum amplitude (1-11)
%       .peak_val           - Maximum amplitude (mV)
%       .bump_range         - [left, right] FWHM boundary positions
%       .bump_width         - Number of positions in the bump
%       .centroid_peak_delta - |centroid_idx - peak_pos|
%
%   AMPLITUDE METHOD:
%     For each position, response amplitude is the 99.5th percentile of
%     the baseline-subtracted trace during the response window (flash
%     onset at sample 5001 through flash offset + 75 ms at sample 6551).
%     This is a robust estimate of the peak depolarisation, less sensitive
%     to single-sample noise spikes than max().
%
%   BUMP FINDING:
%     The RF bump is defined as the contiguous region around the peak
%     position where amplitude >= peak/2 (full width at half maximum).
%     The centroid is computed over bump positions only, excluding
%     low-amplitude flanking positions that would bias the estimate.
%     If the bump is genuinely broad (e.g. all 11 positions above
%     half-max), the centroid uses all positions — bump width itself
%     is a measured phenotype.
%
%   See also FIND_PD_FROM_LUT, ANALYZE_SINGLE_EXPERIMENT_MR

    n_pos = 11;
    col = pd_info.bar_flash_col;
    pos_order = pd_info.pos_order;
    peak_amplitudes = zeros(1, n_pos);

    % Response window: flash onset (5001) through flash offset + 75 ms
    % Flash offset = 5801 (80.1 ms at 10 kHz), tail = 750 samples (75 ms)
    resp_start = 5001;
    resp_end   = 5801 + 750;  % = 6551

    for pos_idx = 1:n_pos
        flash_pos = pos_order(pos_idx);
        ts_mean = mean_slow_bf{flash_pos, col};
        if ~isempty(ts_mean)
            bl_mean = mean(ts_mean(bl_samples));
            win_end = min(resp_end, numel(ts_mean));
            resp_win = resp_start : win_end;
            peak_amplitudes(pos_idx) = prctile(ts_mean(resp_win) - bl_mean, 99.5);
        end
    end

    % Threshold negative amplitudes to zero — only depolarizations count
    A = max(peak_amplitudes, 0);

    % --- Find the RF bump (contiguous FWHM around peak) ---
    [peak_val, peak_pos] = max(A);

    if peak_val > 0
        threshold = peak_val / 2;

        % Grow outward from peak to find contiguous half-max region
        left = peak_pos;
        while left > 1 && A(left - 1) >= threshold
            left = left - 1;
        end
        right = peak_pos;
        while right < n_pos && A(right + 1) >= threshold
            right = right + 1;
        end
        bump_positions = left:right;

        % Response-weighted centroid over bump positions only
        A_bump = A(bump_positions);
        centroid_idx = sum(A_bump .* bump_positions) / sum(A_bump);
    else
        % Fallback: no depolarising response at any position
        peak_pos = 6;
        left = 1;
        right = n_pos;
        bump_positions = left:right;
        centroid_idx = 6;
    end

    % Convert to degrees: position 6 = 0°, spacing = 2.5°
    centroid_deg = (centroid_idx - 6) * 2.5;

    % Validation: centroid should be near the peak
    centroid_peak_delta = abs(centroid_idx - peak_pos);
    if centroid_peak_delta > 1.0
        warning('compute_rf_centroid:largeDelta', ...
            'Centroid (%.2f) is >1 position from peak (%d) — check RF profile.', ...
            centroid_idx, peak_pos);
    end

    % Pack metrics
    rf_metrics.peak_amplitudes     = peak_amplitudes;
    rf_metrics.peak_pos            = peak_pos;
    rf_metrics.peak_val            = peak_val;
    rf_metrics.bump_range          = [left, right];
    rf_metrics.bump_width          = right - left + 1;
    rf_metrics.centroid_peak_delta = centroid_peak_delta;

end


function annotate_centroid(fig_pd, centroid_idx, centroid_deg, rf_metrics)
% ANNOTATE_CENTROID  Add centroid marker, bump range and degree labels.
%
%   ANNOTATE_CENTROID(FIG_PD, CENTROID_IDX, CENTROID_DEG, RF_METRICS)
%   1. Replaces tile titles with degree labels (−12.5° to +12.5°)
%   2. Highlights bump range tiles with a blue border
%   3. Highlights the tile nearest the centroid with a red border
%   4. Adds a text annotation showing centroid, peak, and bump width
%
%   The figure must be a 1x11 tiled layout from plot_bar_flash_1x11.

    n_pos = 11;

    % Degree labels for positions 1-11 (position 6 = center = 0°)
    deg_labels = cell(1, n_pos);
    for p = 1:n_pos
        deg_val = (p - 6) * 2.5;
        if p == 1
            deg_labels{p} = sprintf('%.1f° (Lead)', deg_val);
        elseif p == 6
            deg_labels{p} = '0° (Center)';
        elseif p == 11
            deg_labels{p} = sprintf('+%.1f° (Trail)', deg_val);
        elseif deg_val > 0
            deg_labels{p} = sprintf('+%.1f°', deg_val);
        else
            deg_labels{p} = sprintf('%.1f°', deg_val);
        end
    end

    % Find all Cartesian axes (sorted left-to-right by position)
    all_axes = findobj(fig_pd, 'Type', 'axes');
    % Filter out any polar axes
    keep = true(size(all_axes));
    for k = 1:numel(all_axes)
        if isa(all_axes(k), 'matlab.graphics.axis.PolarAxes')
            keep(k) = false;
        end
    end
    all_axes = all_axes(keep);

    % Sort by x-position (left to right) to match tile order
    positions = zeros(numel(all_axes), 4);
    for k = 1:numel(all_axes)
        positions(k, :) = all_axes(k).Position;
    end
    [~, sort_idx] = sort(positions(:, 1));
    all_axes = all_axes(sort_idx);

    % Nearest tile to centroid
    nearest_tile = round(centroid_idx);
    nearest_tile = max(1, min(n_pos, nearest_tile));

    % Bump range from rf_metrics
    bump_left  = rf_metrics.bump_range(1);
    bump_right = rf_metrics.bump_range(2);

    % Update tile titles and highlight centroid + bump tiles
    for k = 1:min(numel(all_axes), n_pos)
        ax = all_axes(k);
        title(ax, deg_labels{k}, 'FontSize', 9);

        if k == nearest_tile
            % Centroid tile: red border (highest priority)
            set(ax, 'XColor', [0.8 0 0], 'YColor', [0.8 0 0], 'LineWidth', 2.5);
        elseif k >= bump_left && k <= bump_right
            % Bump tiles: blue border
            set(ax, 'XColor', [0.2 0.4 0.8], 'YColor', [0.2 0.4 0.8], 'LineWidth', 1.5);
        end
    end

    % Add centroid + bump info text annotation at top of figure
    info_str = sprintf('Centroid: %+.1f° (pos %.1f)  |  Peak: pos %d  |  Bump: %d-%d (%d wide)', ...
        centroid_deg, centroid_idx, rf_metrics.peak_pos, ...
        bump_left, bump_right, rf_metrics.bump_width);
    annotation(fig_pd, 'textbox', [0.01 0.92 0.6 0.06], ...
        'String', info_str, ...
        'FitBoxToText', 'on', 'FontSize', 10, 'FontWeight', 'bold', ...
        'Color', [0.8 0 0], 'EdgeColor', 'none', 'BackgroundColor', 'w');

end


function ortho_rf = compute_ortho_rf_metrics(mean_slow_bf, pd_info, bl_samples)
% COMPUTE_ORTHO_RF_METRICS  Peak amplitudes and bump width for orthogonal axis.

    n_pos = 11;
    col = pd_info.ortho_flash_col;
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
    ortho_rf.peak_amplitudes = peak_amplitudes;

    % FWHM bump width
    [peak_val, ~] = max(A);
    if peak_val <= 0
        ortho_rf.bump_width = NaN;
        return;
    end
    threshold = peak_val / 2;
    [~, pk] = max(A);
    left = pk;
    while left > 1 && A(left - 1) >= threshold
        left = left - 1;
    end
    right = pk;
    while right < numel(A) && A(right + 1) >= threshold
        right = right + 1;
    end
    ortho_rf.bump_width = right - left + 1;

end


function annotate_polar_metrics(fig_polar, dsi_pdnd, ortho_width, aspect_ratio)
% ANNOTATE_POLAR_METRICS  Add DSI, ortho width, and aspect ratio to polar figure.

    figure(fig_polar);

    % Build annotation string
    lines = {};
    lines{end+1} = sprintf('DSI = %.2f', dsi_pdnd);
    if ~isnan(ortho_width)
        lines{end+1} = sprintf('Ortho width = %d pos', ortho_width);
    end
    if ~isnan(aspect_ratio)
        lines{end+1} = sprintf('Aspect ratio = %.2f', aspect_ratio);
    end

    annotation(fig_polar, 'textbox', [0.01 0.01 0.25 0.12], ...
        'String', strjoin(lines, '\n'), ...
        'FitBoxToText', 'on', 'FontSize', 9, 'FontWeight', 'bold', ...
        'Color', [0 0.3 0.7], 'EdgeColor', [0.7 0.7 0.7], ...
        'BackgroundColor', 'w', 'Interpreter', 'none');

end

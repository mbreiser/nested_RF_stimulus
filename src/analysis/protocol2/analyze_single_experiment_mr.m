function analyze_single_experiment_mr(exp_folder, opts)
% ANALYZE_SINGLE_EXPERIMENT_MR  Enhanced bar sweep and bar flash analysis.
%
%   ANALYZE_SINGLE_EXPERIMENT_MR(EXP_FOLDER) runs the full analysis pipeline
%   for a single Protocol 2 experiment, identical to ANALYZE_SINGLE_EXPERIMENT.
%
%   ANALYZE_SINGLE_EXPERIMENT_MR(EXP_FOLDER, OPTS) uses the options structure
%   to override default parameters. In addition to the options accepted by
%   ANALYZE_SINGLE_EXPERIMENT, this variant supports:
%
%     opts.visualize_everything  (default: false)
%       When true, adds three diagnostic enhancements:
%         1. A voltage histogram figure showing the full-recording and
%            stimulus-only voltage distributions, so you can see what the
%            baseline calculations are based on.
%         2. Vertical lines on all bar flash figures (8x11 heatmap, 1x11
%            PD-ND, 1x11 orthogonal) marking the flash ON and OFF times.
%         3. A thicker border and annotation on the PD row of the 8x11
%            heatmap, visually linking it to the 1x11 PD-ND figure.
%         4. A red diametric line on the polar plot showing the PD-ND
%            bar orientation axis.
%         5. Relabels PD-ND endpoints to "Leading" / "Trailing".
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
%     Figure 1: Slow bar sweep polar timeseries with vector sum arrow
%     Figure 2: Full 8x11 bar flash heatmap (all orientations)
%     Figure 3: 1x11 bar flash subplots along PD-ND axis
%     Figure 4: 1x11 bar flash subplots along orthogonal axis
%     Figure 5: (visualize_everything only) Voltage histograms
%
%   EXAMPLE:
%     opts.visualize_everything = true;
%     opts.save_figs = true;
%     opts.save_dir  = fullfile(exp_folder, 'analysis_output_mr');
%     analyze_single_experiment_mr(exp_folder, opts);
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

    % 1x11 Orthogonal axis
    flash_opts.fig_position = [50 100 1800 300];
    ortho_title = sprintf('Bar Flash Orthogonal — Orient:%.0f° — %s — %s', ...
        pd_info.ortho_orientation, ...
        strrep(date_str, '_', '-'), strrep(metadata.Strain, '_', ' '));
    fig_ortho = plot_bar_flash_1x11(...
        data_slow_bf(:, pd_info.ortho_flash_col, :), ...
        mean_slow_bf(:, pd_info.ortho_flash_col), ...
        pd_info.pos_order, ortho_title, flash_opts);

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

function analyze_single_experiment(exp_folder, opts)
% ANALYZE_SINGLE_EXPERIMENT  Bar sweep and bar flash analysis for one experiment.
%
%   ANALYZE_SINGLE_EXPERIMENT(EXP_FOLDER) runs the full analysis pipeline
%   for a single Protocol 2 experiment: loads data, verifies LUT directions,
%   computes the preferred direction from slow bar sweeps via vector sum,
%   and generates bar flash response figures along the PD-ND and orthogonal
%   axes.
%
%   ANALYZE_SINGLE_EXPERIMENT(EXP_FOLDER, OPTS) uses the options structure
%   to override default parameters.
%
%   INPUTS:
%     exp_folder - Full path to an experiment directory
%                  e.g. '/path/to/protocol2/data/1DRF/2025_11_10_10_17'
%     opts       - (Optional) structure with fields:
%                    .lut_path         - Path to bar_lut.mat
%                                        (default: bar_lut.mat in the same
%                                        directory as this script)
%                    .plot_order       - 1x16 data-row-to-subplot mapping
%                                        (default: [1,3,5,7,9,11,13,15,
%                                                    2,4,6,8,10,12,14,16])
%                    .baseline_range   - [start end] samples for bar sweep
%                                        pre-stimulus baseline
%                                        (default: [1000 9000])
%                    .stim_trim_end    - Samples to trim from end of the
%                                        bar sweep stimulus period
%                                        (default: 7000)
%                    .percentile       - Percentile for peak detection
%                                        (default: 98)
%                    .flash_baseline   - Sample indices for bar flash
%                                        baseline subtraction
%                                        (default: 1:5000, 500ms at 10kHz)
%                    .flash_ylim       - [ymin ymax] for 1x11 bar flash plots
%                                        (default: [-15 25])
%                    .pattern_offset   - Offset from experiment pattern number
%                                        to full-field pattern index
%                                        (default: 2)
%                    .save_figs        - Boolean, save figures as PDF
%                                        (default: true)
%                    .save_dir         - Output directory for PDFs
%                                        (default: <exp_folder>/analysis_output)
%
%   OUTPUTS:
%     Figure 1: Slow bar sweep polar timeseries with vector sum arrow
%     Figure 2: Full 8x11 bar flash heatmap (all orientations)
%     Figure 3: 1x11 bar flash subplots along PD-ND axis
%     Figure 4: 1x11 bar flash subplots along orthogonal axis
%
%   ANALYSIS PIPELINE:
%     1. Load data (LUT, TDMS recording, currentExp.mat)
%     2. Parse bar sweep data and verify LUT directions
%     3. Compute slow-speed (28 dps) bar sweep responses and create polar
%        timeseries figure
%     4. Find preferred direction via vector sum, map to bar flash columns
%     5. Parse bar flash data and generate 8x11 heatmap and 1x11 response
%        figures for PD-ND and orthogonal axes
%     6. Save all figures as PDFs
%
%   EXAMPLE:
%     % Basic usage with defaults:
%     analyze_single_experiment('/path/to/data/1DRF/2025_11_10_10_17')
%
%     % Override y-limits and disable saving:
%     opts.flash_ylim = [-20 30];
%     opts.save_figs  = false;
%     analyze_single_experiment('/path/to/data/1DRF/2025_11_10_10_17', opts)
%
%   See also LOAD_PROTOCOL2_DATA, PARSE_BAR_DATA, PARSE_BAR_FLASH_DATA,
%            VERIFY_LUT_DIRECTIONS, COMPUTE_BAR_SWEEP_RESPONSES,
%            FIND_PD_FROM_LUT, PLOT_SLOW_BAR_SWEEP_POLAR,
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
    v_data   = Log.ADC.Volts(2, :) * 10;  % voltage data (scaled)
    median_v = median(v_data);

    ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
        'pattern_order', 'func_order', 'metadata');
    metadata = ce.metadata;

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
    pd_info = find_pd_from_lut(max_v, lut_directions, lut_orientations, ...
        lut_patterns, lut_functions, opts.plot_order, Tbl, opts.pattern_offset);

    %% Step 5: Parse bar flash data and generate figures
    [data_slow_bf, ~, mean_slow_bf, ~] = parse_bar_flash_data(f_data, v_data);

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

    % 1x11 Orthogonal axis
    flash_opts.fig_position = [50 100 1800 300];
    ortho_title = sprintf('Bar Flash Orthogonal — Orient:%.0f° — %s — %s', ...
        pd_info.ortho_orientation, ...
        strrep(date_str, '_', '-'), strrep(metadata.Strain, '_', ' '));
    fig_ortho = plot_bar_flash_1x11(...
        data_slow_bf(:, pd_info.ortho_flash_col, :), ...
        mean_slow_bf(:, pd_info.ortho_flash_col), ...
        pd_info.ortho_pos_order, ortho_title, flash_opts);

    %% Step 6: Save figures
    if opts.save_figs
        save_analysis_figures(opts.save_dir, ...
            fig_polar, fig_heatmap, fig_pd, fig_ortho);
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


function save_analysis_figures(save_dir, fig_polar, fig_heatmap, fig_pd, fig_ortho)
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

    fprintf('\nFigures saved to: %s\n', save_dir);

end

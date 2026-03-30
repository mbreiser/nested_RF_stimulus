function fig = plot_flash_1x11_population(traces_ctrl, traces_ttl, title_str, opts)
% PLOT_FLASH_1X11_POPULATION  Population 1x11 bar flash plot with shaded spread.
%
%   FIG = PLOT_FLASH_1X11_POPULATION(TRACES_CTRL, TRACES_TTL, TITLE_STR)
%   creates a 1x11 tiled layout showing the population center bar flash
%   response at each spatial position, with a shaded spread band. Control
%   and ttl groups are overlaid on the same axes.
%
%   FIG = PLOT_FLASH_1X11_POPULATION(..., OPTS) uses the options structure
%   to override default plotting parameters.
%
%   INPUTS:
%     traces_ctrl - Cell array (one per control cell). Each element is an
%                   11xN numeric matrix of baseline-subtracted mean flash
%                   traces (11 positions x N timepoints). Positions are
%                   ordered ND (row 1) to PD (row 11).
%     traces_ttl  - Cell array (one per ttl cell), same format.
%     title_str   - Figure title string
%     opts        - (Optional) structure with fields:
%                     .y_limits     - [ymin ymax] in mV (default: [-15 35])
%                     .fig_position - [x y w h] in pixels
%                                     (default: [50 400 1800 300])
%                     .plot_type    - 'pd_nd' or 'orthogonal'
%                                     (default: 'pd_nd'). Controls whether
%                                     tiles 1/11 are labelled PD/ND.
%                     .stat_method  - 'median_mad' (default), 'mean_sem',
%                                     or 'mean_sd'
%                     .col_labels   - Cell array of column header strings
%                                     (default: auto-generated ND/Center/PD)
%                     .show_n_per_pos - Boolean, show per-position n-counts
%                                       below each tile (default: false)
%                     .stim_onset_sample  - Sample for stim onset green line
%                     .stim_offset_sample - Sample for stim offset green line
%                     .resp_end_sample - Sample index for response window end
%                                        line (thin dark). When set (e.g.
%                                        6551), marks the peak detection
%                                        window boundary.
%                     .show_ordinal_ranks - Boolean, show ordinal amplitude
%                                           rank (#1, #2, ...) for each
%                                           position (default: false)
%
%   OUTPUT:
%     fig - Figure handle
%
%   FIGURE LAYOUT:
%     1x11 tiled layout. Each tile shows the center trace for control
%     (black) and ttl (red) with shaded spread (gray and light red
%     respectively). The shaded band represents MAD (default) or SEM
%     depending on stat_method. For 'pd_nd' plots, tile 1 is labelled
%     'ND', tile 6 'Center', and tile 11 'PD'. For 'orthogonal' plots,
%     tile 1 is '1', tile 6 'Center', and tile 11 '11'.
%
%   See also PLOT_BAR_FLASH_1X11, BATCH_ANALYZE_1DRF, ANALYZE_SINGLE_EXPERIMENT
% ________________________________________________________________________

    % Set defaults
    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'y_limits'),     opts.y_limits     = [-15 35]; end
    if ~isfield(opts, 'fig_position'), opts.fig_position = [50 400 1800 300]; end
    if ~isfield(opts, 'plot_type'),    opts.plot_type    = 'pd_nd'; end
    if ~isfield(opts, 'stat_method'),  opts.stat_method  = 'median_mad'; end
    if ~isfield(opts, 'col_labels'),  opts.col_labels   = {}; end
    if ~isfield(opts, 'show_n_per_pos'), opts.show_n_per_pos = false; end
    if ~isfield(opts, 'require_both_n2'), opts.require_both_n2 = false; end
    if ~isfield(opts, 'resp_end_sample'), opts.resp_end_sample = []; end
    if ~isfield(opts, 'stim_onset_sample'), opts.stim_onset_sample = []; end
    if ~isfield(opts, 'stim_offset_sample'), opts.stim_offset_sample = []; end
    if ~isfield(opts, 'show_ordinal_ranks'), opts.show_ordinal_ranks = false; end

    % Derive n_pos from data (generalise beyond 11)
    n_pos = 0;
    all_inputs = [traces_ctrl(:); traces_ttl(:)];
    for kk = 1:numel(all_inputs)
        if ~isempty(all_inputs{kk})
            n_pos = size(all_inputs{kk}, 1);
            break;
        end
    end
    if n_pos == 0
        fig = figure('Name', title_str);
        return;
    end

    % Colours
    col_ctrl_line = [0 0 0];            % black
    col_ctrl_fill = [0.80 0.80 0.80];   % light gray
    col_ttl_line  = [1 0 0];            % red
    col_ttl_fill  = [1 0.70 0.70];      % light red
    alpha_val     = 0.40;

    % --- Pre-compute stats for all positions (needed for ordinal ranks) ---
    ctr_ctrl_all = cell(1, n_pos);
    spr_ctrl_all = cell(1, n_pos);
    n_ctrl_all   = zeros(1, n_pos);
    ctr_ttl_all  = cell(1, n_pos);
    spr_ttl_all  = cell(1, n_pos);
    n_ttl_all    = zeros(1, n_pos);

    for pos_idx = 1:n_pos
        [ctr_ctrl_all{pos_idx}, spr_ctrl_all{pos_idx}, n_ctrl_all(pos_idx)] = ...
            compute_trace_stats(traces_ctrl, pos_idx, opts.stat_method, n_pos);
        [ctr_ttl_all{pos_idx}, spr_ttl_all{pos_idx}, n_ttl_all(pos_idx)] = ...
            compute_trace_stats(traces_ttl, pos_idx, opts.stat_method, n_pos);
    end

    % --- Compute ordinal amplitude ranks if requested ---
    ranks_ctrl = NaN(1, n_pos);
    ranks_ttl  = NaN(1, n_pos);
    if opts.show_ordinal_ranks
        peaks_ctrl = NaN(1, n_pos);
        peaks_ttl  = NaN(1, n_pos);
        for pos_idx = 1:n_pos
            if ~isempty(ctr_ctrl_all{pos_idx})
                peaks_ctrl(pos_idx) = max(ctr_ctrl_all{pos_idx});
            end
            if ~isempty(ctr_ttl_all{pos_idx})
                peaks_ttl(pos_idx) = max(ctr_ttl_all{pos_idx});
            end
        end
        ranks_ctrl = compute_ordinal_ranks(peaks_ctrl);
        ranks_ttl  = compute_ordinal_ranks(peaks_ttl);
    end

    % --- Plot ---
    fig = figure('Name', title_str);
    tiledlayout(1, n_pos, 'TileSpacing', 'compact', 'Padding', 'compact');

    for pos_idx = 1:n_pos
        nexttile;
        hold on;

        ctr_ctrl = ctr_ctrl_all{pos_idx};
        spr_ctrl = spr_ctrl_all{pos_idx};
        n_ctrl   = n_ctrl_all(pos_idx);
        ctr_ttl  = ctr_ttl_all{pos_idx};
        spr_ttl  = spr_ttl_all{pos_idx};
        n_ttl    = n_ttl_all(pos_idx);

        % Joint n>=2 filter: require BOTH groups have n>=2 at this position
        if opts.require_both_n2 && ~(n_ctrl >= 2 && n_ttl >= 2)
            ylim(opts.y_limits);
            set(gca, 'XTick', []);
            if pos_idx == 1, ylabel('\DeltamV'); else, set(gca, 'YTickLabel', []); end
            box off;
            % Still label the position
            if ~isempty(opts.col_labels) && pos_idx <= numel(opts.col_labels)
                title(opts.col_labels{pos_idx}, 'FontSize', 9);
            end
            continue;
        end

        % Stimulus timing lines (green)
        if ~isempty(opts.stim_onset_sample)
            xline(opts.stim_onset_sample, '-', 'Color', [0.2 0.7 0.2], ...
                'LineWidth', 0.8, 'Alpha', 0.7);
        end
        if ~isempty(opts.stim_offset_sample)
            xline(opts.stim_offset_sample, '-', 'Color', [0.2 0.7 0.2], ...
                'LineWidth', 0.8, 'Alpha', 0.7);
        end

        % Response window end line (thin dark vertical line)
        if ~isempty(opts.resp_end_sample)
            xline(opts.resp_end_sample, '-', 'Color', [0.3 0.3 0.3], ...
                'LineWidth', 0.5, 'Alpha', 0.6);
        end

        % Plot control: shaded spread then center line (require n >= 2)
        if ~isempty(ctr_ctrl) && n_ctrl >= 2
            x = 1:numel(ctr_ctrl);
            plot_shaded_trace(x, ctr_ctrl, spr_ctrl, ...
                col_ctrl_line, col_ctrl_fill, alpha_val, 1.5);
        end

        % Plot ttl: shaded spread then center line (require n >= 2)
        if ~isempty(ctr_ttl) && n_ttl >= 2
            x = 1:numel(ctr_ttl);
            plot_shaded_trace(x, ctr_ttl, spr_ttl, ...
                col_ttl_line, col_ttl_fill, alpha_val, 1.5);
        end

        ylim(opts.y_limits);
        set(gca, 'XTick', []);

        % Y-axis label on first tile only
        if pos_idx == 1
            ylabel('\DeltamV');
        else
            set(gca, 'YTickLabel', []);
        end

        % Position labels — use custom labels if provided
        if ~isempty(opts.col_labels) && pos_idx <= numel(opts.col_labels)
            title(opts.col_labels{pos_idx}, 'FontSize', 9);
        elseif strcmpi(opts.plot_type, 'pd_nd')
            if pos_idx == 1
                title('ND');
            elseif pos_idx == 6
                title('Center');
            elseif pos_idx == n_pos
                title('PD');
            else
                title(sprintf('%d', pos_idx));
            end
        else
            % Orthogonal plot: no PD/ND labels
            if pos_idx == ceil(n_pos / 2)
                title('Center');
            else
                title(sprintf('%d', pos_idx));
            end
        end

        % Ordinal rank annotations (upper-right corner of each tile)
        if opts.show_ordinal_ranks
            yl = ylim;
            xl = xlim;
            x_rank = xl(2) - 0.03 * diff(xl);
            y_top  = yl(2) - 0.04 * diff(yl);
            y_step = 0.10 * diff(yl);
            if ~isnan(ranks_ctrl(pos_idx)) && n_ctrl >= 2
                text(x_rank, y_top, sprintf('#%d', ranks_ctrl(pos_idx)), ...
                    'FontSize', 7, 'FontWeight', 'bold', 'Color', col_ctrl_line, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
            end
            if ~isnan(ranks_ttl(pos_idx)) && n_ttl >= 2
                text(x_rank, y_top - y_step, sprintf('#%d', ranks_ttl(pos_idx)), ...
                    'FontSize', 7, 'FontWeight', 'bold', 'Color', col_ttl_line, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
            end
        end

        % Per-position n-count labels (stacked vertically at bottom)
        if opts.show_n_per_pos
            yl = ylim;
            xl = xlim;
            x_mid = mean(xl);
            y_base = yl(1) + (yl(2) - yl(1)) * 0.02;
            y_step = (yl(2) - yl(1)) * 0.07;
            if n_ctrl > 0
                text(x_mid, y_base, sprintf('n=%d', n_ctrl), ...
                    'FontSize', 6, 'Color', col_ctrl_line, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            end
            if n_ttl > 0
                text(x_mid, y_base + y_step, sprintf('n=%d', n_ttl), ...
                    'FontSize', 6, 'Color', col_ttl_line, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            end
        end

        box off;
        ax = gca;
        ax.LineWidth = 1.2;
        ax.TickDir = 'out';
        ax.TickLength = [0.015 0.015];
        ax.FontSize = 14;

        % Add legend to first tile only
        if pos_idx == 1
            h = [];
            labs = {};
            if ~isempty(ctr_ctrl)
                h(end+1) = plot(NaN, NaN, '-', 'Color', col_ctrl_line, 'LineWidth', 1.5);
                labs{end+1} = sprintf('control (n=%d)', n_ctrl);
            end
            if ~isempty(ctr_ttl)
                h(end+1) = plot(NaN, NaN, '-', 'Color', col_ttl_line, 'LineWidth', 1.5);
                labs{end+1} = sprintf('ttl (n=%d)', n_ttl);
            end
            if ~isempty(h)
                legend(h, labs, 'Location', 'northwest', 'FontSize', 8);
            end
        end
    end

    sgtitle(title_str, 'FontSize', 16);
    set(fig, 'Position', opts.fig_position);

end


%% ========================= Local Functions ============================

function [center_trace, spread_trace, n] = compute_trace_stats(traces_cell, pos_idx, stat_method, n_pos)
% COMPUTE_TRACE_STATS  Stack traces at one position across cells, compute center/spread.
%   stat_method: 'median_mad' (default), 'mean_sem', or 'mean_sd'.
%   n_pos: total number of positions (for bounds checking).
%   n: number of cells with non-NaN data at this position.

    center_trace = [];
    spread_trace = [];
    n            = 0;

    if isempty(traces_cell)
        return;
    end
    if nargin < 4, n_pos = 11; end

    % First pass: find minimum trace length across cells
    min_len = Inf;
    for k = 1:numel(traces_cell)
        mat = traces_cell{k};
        if ~isempty(mat) && pos_idx <= size(mat, 1)
            min_len = min(min_len, size(mat, 2));
        end
    end
    if isinf(min_len)
        return;
    end

    % Second pass: collect traces truncated to common length
    all_traces = [];
    for k = 1:numel(traces_cell)
        mat = traces_cell{k};
        if ~isempty(mat) && pos_idx <= size(mat, 1)
            row = mat(pos_idx, 1:min_len);
            % Skip rows that are entirely NaN (from canvas padding)
            if ~all(isnan(row))
                all_traces = [all_traces; row]; %#ok<AGROW>
            end
        end
    end

    if isempty(all_traces)
        return;
    end

    n = size(all_traces, 1);

    if strcmpi(stat_method, 'mean_sd')
        center_trace = mean(all_traces, 1, 'omitnan');
        spread_trace = std(all_traces, 0, 1, 'omitnan');
    elseif strcmpi(stat_method, 'mean_sem')
        center_trace = mean(all_traces, 1, 'omitnan');
        sd           = std(all_traces, 0, 1, 'omitnan');
        spread_trace = sd ./ sqrt(n);
    else  % 'median_mad'
        center_trace = median(all_traces, 1, 'omitnan');
        spread_trace = mad(all_traces, 1, 1);  % median absolute deviation (flag=1)
    end

end


function ranks = compute_ordinal_ranks(peaks)
% COMPUTE_ORDINAL_RANKS  Rank peaks from 1 (highest) to n (lowest).
%   NaN values in peaks get NaN ranks.

    ranks = NaN(size(peaks));
    valid = ~isnan(peaks);
    if ~any(valid), return; end

    vals = peaks(valid);
    [~, sort_idx] = sort(vals, 'descend');
    r = zeros(size(vals));
    r(sort_idx) = 1:numel(vals);
    ranks(valid) = r;

end


function plot_shaded_trace(x, center_vals, spread_vals, line_col, fill_col, alpha_val, lw)
% PLOT_SHADED_TRACE  Plot center line with shaded spread band on current axes.

    x = x(:)';
    m = center_vals(:)';
    s = spread_vals(:)';

    upper = m + s;
    lower = m - s;

    % Shaded spread band
    fill([x, fliplr(x)], [upper, fliplr(lower)], fill_col, ...
        'FaceAlpha', alpha_val, 'EdgeColor', 'none');

    % Center line on top
    plot(x, m, '-', 'Color', line_col, 'LineWidth', lw);

end

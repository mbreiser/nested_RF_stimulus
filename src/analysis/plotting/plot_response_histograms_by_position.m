function fig = plot_response_histograms_by_position(results, opts)
% PLOT_RESPONSE_HISTOGRAMS_BY_POSITION  Response distributions per position.
%
%   FIG = PLOT_RESPONSE_HISTOGRAMS_BY_POSITION(RESULTS) creates a 2x6 grid
%   of histogram panels (11 positions + 1 empty), each showing overlaid
%   distributions of per-cell mean response (stim onset to end of trace)
%   for control (black) vs TTL (red).
%
%   FIG = PLOT_RESPONSE_HISTOGRAMS_BY_POSITION(RESULTS, OPTS) uses options:
%     opts.cell_type    - 'ON', 'OFF', or 'both' (default: 'ON')
%     opts.trace_field  - field in results containing 11xN traces
%                         (default: 'pd_flash_peak_aligned')
%     opts.stim_onset   - sample index for stim onset (default: 5001)
%     opts.n_bins       - number of histogram bins (default: 20)
%     opts.fig_position - [x y w h] in pixels (default: [50 50 1400 500])
%     opts.reject_threshold - minimum |mean response| to include (mV, default: 0)
%
%   Each histogram shows the distribution of per-cell mean trace values
%   from stim onset to end of trace. This captures both depolarization and
%   hyperpolarization naturally for each cell at each position.
%
%   Wilcoxon rank-sum p-values are annotated when both groups have n >= 2.
%
%   INPUTS:
%     results - Struct array from batch_analyze_1DRF
%     opts    - (Optional) structure
%
%   OUTPUT:
%     fig - Figure handle
%
%   See also PLOT_AMPLITUDE_BY_POSITION, BATCH_ANALYZE_1DRF

    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'cell_type'),        opts.cell_type        = 'ON'; end
    if ~isfield(opts, 'trace_field'),      opts.trace_field      = 'pd_flash_peak_aligned'; end
    if ~isfield(opts, 'stim_onset'),       opts.stim_onset       = 5001; end
    if ~isfield(opts, 'n_bins'),           opts.n_bins           = 20; end
    if ~isfield(opts, 'fig_position'),     opts.fig_position     = [50 50 1400 500]; end
    if ~isfield(opts, 'reject_threshold'), opts.reject_threshold = 0; end

    % Select cells by type
    if strcmpi(opts.cell_type, 'ON')
        type_mask = [results.is_on];
    elseif strcmpi(opts.cell_type, 'OFF')
        type_mask = ~[results.is_on];
    else
        type_mask = true(size(results));
    end

    ctrl_mask = type_mask & ~[results.is_ttl];
    ttl_mask  = type_mask &  [results.is_ttl];

    col_ctrl = [0 0 0];
    col_ttl  = [1 0 0];

    % Extract per-cell mean response at each position
    ctrl_means = extract_mean_responses(results(ctrl_mask), opts);
    ttl_means  = extract_mean_responses(results(ttl_mask), opts);

    n_ctrl = size(ctrl_means, 1);
    n_ttl  = size(ttl_means, 1);

    % Determine common bin edges across all positions
    all_vals = [ctrl_means(:); ttl_means(:)];
    all_vals = all_vals(~isnan(all_vals));
    if isempty(all_vals)
        fig = figure('Name', 'No data');
        return;
    end
    bin_edges = linspace(min(all_vals), max(all_vals), opts.n_bins + 1);

    positions = -5:5;
    fig = figure('Name', sprintf('Response Histograms — %s Cells', opts.cell_type), ...
        'Position', opts.fig_position);
    tiledlayout(2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

    for pos_idx = 1:11
        % Map 11 positions into 2x6 grid (row-major: pos1=tile1, ..., pos11=tile11)
        nexttile; hold on;

        ctrl_vals = ctrl_means(:, pos_idx);
        ctrl_vals = ctrl_vals(~isnan(ctrl_vals));
        ttl_vals  = ttl_means(:, pos_idx);
        ttl_vals  = ttl_vals(~isnan(ttl_vals));

        % Histogram for control
        if ~isempty(ctrl_vals)
            histogram(ctrl_vals, bin_edges, 'FaceColor', col_ctrl, ...
                'FaceAlpha', 0.4, 'EdgeColor', col_ctrl, 'EdgeAlpha', 0.6);
        end

        % Histogram for TTL
        if ~isempty(ttl_vals)
            histogram(ttl_vals, bin_edges, 'FaceColor', col_ttl, ...
                'FaceAlpha', 0.4, 'EdgeColor', col_ttl, 'EdgeAlpha', 0.6);
        end

        % Vertical line at 0
        xline(0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

        % Position label
        pos_label = sprintf('%+d', positions(pos_idx));
        if pos_idx == 6
            pos_label = '0 (center)';
        end
        title(pos_label, 'FontSize', 9);

        % n-counts
        text(0.02, 0.95, sprintf('c:%d t:%d', numel(ctrl_vals), numel(ttl_vals)), ...
            'Units', 'normalized', 'FontSize', 7, 'VerticalAlignment', 'top');

        % Rank-sum test annotation
        if numel(ctrl_vals) >= 2 && numel(ttl_vals) >= 2
            p = ranksum(ctrl_vals, ttl_vals);
            if p < 0.05
                if p < 0.001
                    p_str = 'p<0.001';
                elseif p < 0.01
                    p_str = sprintf('p=%.3f', p);
                else
                    p_str = sprintf('p=%.2f', p);
                end
                text(0.98, 0.95, p_str, 'Units', 'normalized', ...
                    'FontSize', 7, 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'top', 'Color', [0.2 0.2 0.8]);
            end
        end

        set(gca, 'FontSize', 8, 'TickDir', 'out');
        if pos_idx == 1 || pos_idx == 7
            ylabel('Count');
        end
        if pos_idx > 5
            xlabel('mV');
        end
        box off;
    end

    % 12th tile: legend
    nexttile; axis off;
    hold on;
    h1 = patch(NaN, NaN, col_ctrl, 'FaceAlpha', 0.4, 'EdgeColor', col_ctrl);
    h2 = patch(NaN, NaN, col_ttl, 'FaceAlpha', 0.4, 'EdgeColor', col_ttl);
    legend([h1 h2], {sprintf('control (n=%d)', n_ctrl), ...
                      sprintf('TTL (n=%d)', n_ttl)}, ...
        'Location', 'north', 'FontSize', 10);

    sgtitle(sprintf('%s Cells — Response Distributions by Position (%s)', ...
        opts.cell_type, strrep(opts.trace_field, '_', ' ')), 'FontSize', 13);

end


function mean_matrix = extract_mean_responses(results_subset, opts)
% EXTRACT_MEAN_RESPONSES  Per-cell mean response at each position.
%   Returns n_cells x 11 matrix. Mean is from stim_onset to end of trace.

    n = numel(results_subset);
    mean_matrix = NaN(n, 11);

    if n == 0, return; end

    for k = 1:n
        if ~isfield(results_subset(k), opts.trace_field)
            continue;
        end
        traces = results_subset(k).(opts.trace_field);  % 11 x N
        if isempty(traces), continue; end

        n_samples = size(traces, 2);
        stim_start = opts.stim_onset;
        if stim_start > n_samples, continue; end

        for pos = 1:11
            trace = traces(pos, :);
            if all(isnan(trace)), continue; end

            resp_mean = mean(trace(stim_start:end), 'omitnan');

            % Rejection threshold
            if abs(resp_mean) < opts.reject_threshold
                continue;
            end

            mean_matrix(k, pos) = resp_mean;
        end
    end

end

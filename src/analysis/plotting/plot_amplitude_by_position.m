function fig = plot_amplitude_by_position(results, opts)
% PLOT_AMPLITUDE_BY_POSITION  Robust peak metrics vs RF position.
%
%   FIG = PLOT_AMPLITUDE_BY_POSITION(RESULTS) creates a 2x1 figure showing
%   depolarization (99.5th percentile) and hyperpolarization (0.5th
%   percentile) amplitudes at each aligned RF position.
%
%   FIG = PLOT_AMPLITUDE_BY_POSITION(RESULTS, OPTS) uses options:
%     opts.cell_type       - 'ON', 'OFF', or 'both' (default: 'ON')
%     opts.trace_field     - field in results containing 11xN traces
%                            (default: 'pd_flash_peak_aligned')
%     opts.dep_window      - [start, end] sample range for depol. pctile
%                            (default: [5001, 6551])
%     opts.hyp_window      - [start, end] sample range for hyperpol. pctile
%                            (default: [5001, Inf])  — Inf = end of trace
%     opts.dep_pctile      - percentile for depolarization (default: 99.5)
%     opts.hyp_pctile      - percentile for hyperpolarization (default: 0.5)
%     opts.reject_threshold - minimum |mean response| to include a cell at
%                             a given position (mV, default: 0.5)
%     opts.fig_position    - [x y w h] in pixels (default: [100 100 700 600])
%     opts.alignment_label - center position label (default: '0')
%
%   Each panel shows jittered dots for individual cells (ctrl=black,
%   TTL=red) with mean +/- SEM overlaid. Wilcoxon rank-sum asterisks are
%   shown at positions where ctrl vs TTL p < 0.05 (both n >= 2).
%
%   INPUTS:
%     results - Struct array from batch_analyze_1DRF with fields:
%               .group, .is_on, .is_ttl, and the trace field
%     opts    - (Optional) structure
%
%   OUTPUT:
%     fig - Figure handle
%
%   See also BATCH_ANALYZE_1DRF, PLOT_TIMING_BY_POSITION

    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'cell_type'),        opts.cell_type        = 'ON'; end
    if ~isfield(opts, 'trace_field'),      opts.trace_field      = 'pd_flash_peak_aligned'; end
    if ~isfield(opts, 'dep_window'),       opts.dep_window       = [5001, 6551]; end
    if ~isfield(opts, 'hyp_window'),       opts.hyp_window       = [5001, Inf]; end
    if ~isfield(opts, 'dep_pctile'),       opts.dep_pctile       = 99.5; end
    if ~isfield(opts, 'hyp_pctile'),       opts.hyp_pctile       = 0.5; end
    if ~isfield(opts, 'reject_threshold'), opts.reject_threshold = 0.5; end
    if ~isfield(opts, 'fig_position'),     opts.fig_position     = [100 100 700 600]; end
    if ~isfield(opts, 'alignment_label'),  opts.alignment_label  = '0'; end

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

    positions = -5:5;  % 11 positions relative to aligned center
    col_ctrl = [0 0 0];
    col_ttl  = [1 0 0];

    % Extract per-cell robust metrics at each position
    [ctrl_dep, ctrl_hyp] = extract_robust_amplitudes(results(ctrl_mask), opts);
    [ttl_dep, ttl_hyp]   = extract_robust_amplitudes(results(ttl_mask), opts);

    fig = figure('Name', sprintf('Amplitude — %s Cells', opts.cell_type), ...
        'Position', opts.fig_position);
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    % --- Panel 1: Depolarization ---
    nexttile; hold on;
    plot_dots_and_stats(positions, ctrl_dep, ttl_dep, col_ctrl, col_ttl);
    xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
    ylabel('Depolarization (mV)');
    xlabel('Position relative to center');
    xlim([-5.5 5.5]);
    box off; set(gca, 'TickDir', 'out', 'FontSize', 11);
    title(sprintf('Depolarization (%.1f%% pctile)', opts.dep_pctile));

    % Legend
    n_ctrl = size(ctrl_dep, 1);
    n_ttl  = size(ttl_dep, 1);
    h = [];
    labs = {};
    if n_ctrl > 0
        h(end+1) = plot(NaN, NaN, 'o', 'Color', col_ctrl, 'MarkerFaceColor', col_ctrl, 'MarkerSize', 5);
        labs{end+1} = sprintf('control (n=%d)', n_ctrl);
    end
    if n_ttl > 0
        h(end+1) = plot(NaN, NaN, 'o', 'Color', col_ttl, 'MarkerFaceColor', col_ttl, 'MarkerSize', 5);
        labs{end+1} = sprintf('TTL (n=%d)', n_ttl);
    end
    if ~isempty(h)
        legend(h, labs, 'Location', 'best', 'FontSize', 9);
    end

    % --- Panel 2: Hyperpolarization ---
    nexttile; hold on;
    plot_dots_and_stats(positions, ctrl_hyp, ttl_hyp, col_ctrl, col_ttl);
    xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
    ylabel('Hyperpolarization (mV)');
    xlabel('Position relative to center');
    xlim([-5.5 5.5]);
    box off; set(gca, 'TickDir', 'out', 'FontSize', 11);
    title(sprintf('Hyperpolarization (%.1f%% pctile)', opts.hyp_pctile));

    sgtitle(sprintf('%s Cells — Per-Position Amplitudes (%s)', ...
        opts.cell_type, strrep(opts.trace_field, '_', ' ')), 'FontSize', 14);

end


function [dep_matrix, hyp_matrix] = extract_robust_amplitudes(results_subset, opts)
% EXTRACT_ROBUST_AMPLITUDES  Compute robust depol/hyperpol per position.
%   Returns n_cells x 11 matrices. Positions failing rejection get NaN.

    n = numel(results_subset);
    dep_matrix = NaN(n, 11);
    hyp_matrix = NaN(n, 11);

    if n == 0, return; end

    for k = 1:n
        if ~isfield(results_subset(k), opts.trace_field)
            continue;
        end
        traces = results_subset(k).(opts.trace_field);  % 11 x N
        if isempty(traces), continue; end

        n_samples = size(traces, 2);

        for pos = 1:11
            trace = traces(pos, :);
            if all(isnan(trace)), continue; end

            % Rejection: skip if mean response is too close to baseline
            stim_start = 5001;
            if stim_start > n_samples, continue; end
            mean_resp = mean(trace(stim_start:end), 'omitnan');
            if abs(mean_resp) < opts.reject_threshold
                continue;
            end

            % Depolarization: percentile in dep_window
            dep_s = opts.dep_window(1);
            dep_e = opts.dep_window(2);
            if isinf(dep_e), dep_e = n_samples; end
            dep_e = min(dep_e, n_samples);
            if dep_s <= n_samples
                dep_win = trace(dep_s:dep_e);
                dep_win = dep_win(~isnan(dep_win));
                if ~isempty(dep_win)
                    dep_matrix(k, pos) = prctile(dep_win, opts.dep_pctile);
                end
            end

            % Hyperpolarization: percentile in hyp_window
            hyp_s = opts.hyp_window(1);
            hyp_e = opts.hyp_window(2);
            if isinf(hyp_e), hyp_e = n_samples; end
            hyp_e = min(hyp_e, n_samples);
            if hyp_s <= n_samples
                hyp_win = trace(hyp_s:hyp_e);
                hyp_win = hyp_win(~isnan(hyp_win));
                if ~isempty(hyp_win)
                    hyp_matrix(k, pos) = prctile(hyp_win, opts.hyp_pctile);
                end
            end
        end
    end

end


function plot_dots_and_stats(positions, data_ctrl, data_ttl, col_ctrl, col_ttl)
% PLOT_DOTS_AND_STATS  Jittered dots + mean/SEM + significance asterisks.

    for pos_idx = 1:numel(positions)
        x = positions(pos_idx);

        ctrl_vals = data_ctrl(:, pos_idx);
        ctrl_vals = ctrl_vals(~isnan(ctrl_vals));
        ttl_vals  = data_ttl(:, pos_idx);
        ttl_vals  = ttl_vals(~isnan(ttl_vals));

        % Jittered dots
        if ~isempty(ctrl_vals)
            jx = x - 0.15 + 0.1 * (rand(size(ctrl_vals)) - 0.5);
            scatter(jx, ctrl_vals, 30, col_ctrl, 'filled', ...
                'MarkerFaceAlpha', 0.5);
        end
        if ~isempty(ttl_vals)
            jx = x + 0.15 + 0.1 * (rand(size(ttl_vals)) - 0.5);
            scatter(jx, ttl_vals, 30, col_ttl, 'filled', ...
                'MarkerFaceAlpha', 0.5);
        end

        % Mean +/- SEM
        if numel(ctrl_vals) >= 2
            m = mean(ctrl_vals);
            s = std(ctrl_vals) / sqrt(numel(ctrl_vals));
            errorbar(x - 0.15, m, s, 'o', 'Color', col_ctrl, ...
                'MarkerFaceColor', col_ctrl, 'MarkerSize', 6, ...
                'LineWidth', 1.5, 'CapSize', 4);
        end
        if numel(ttl_vals) >= 2
            m = mean(ttl_vals);
            s = std(ttl_vals) / sqrt(numel(ttl_vals));
            errorbar(x + 0.15, m, s, 'o', 'Color', col_ttl, ...
                'MarkerFaceColor', col_ttl, 'MarkerSize', 6, ...
                'LineWidth', 1.5, 'CapSize', 4);
        end

        % Wilcoxon rank-sum test
        if numel(ctrl_vals) >= 2 && numel(ttl_vals) >= 2
            p = ranksum(ctrl_vals, ttl_vals);
            if p < 0.05
                yl = ylim;
                y_star = max([ctrl_vals; ttl_vals]) + 0.08 * diff(yl);
                if p < 0.001
                    star_str = '***';
                elseif p < 0.01
                    star_str = '**';
                else
                    star_str = '*';
                end
                text(x, y_star, star_str, 'FontSize', 12, ...
                    'HorizontalAlignment', 'center', 'Color', [0.2 0.2 0.8]);
            end
        end
    end

end

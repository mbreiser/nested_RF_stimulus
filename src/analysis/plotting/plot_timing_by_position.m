function fig = plot_timing_by_position(results, opts)
% PLOT_TIMING_BY_POSITION  Temporal metrics vs RF position.
%
%   FIG = PLOT_TIMING_BY_POSITION(RESULTS) creates a 4x1 figure showing
%   rise start, rise time, decay time, and time-to-90% as a function of
%   spatial position, reported as raw (absolute) values in ms.
%
%   FIG = PLOT_TIMING_BY_POSITION(RESULTS, OPTS) uses options:
%     opts.cell_type    - 'ON', 'OFF', or 'both' (default: 'ON')
%     opts.fig_position - [x y w h] in pixels (default: [100 100 600 900])
%
%   Each panel shows mean +/- SEM for control (black) and TTL (red).
%   X-axis: position relative to peak (-5 to +5).
%   Y-axis: time in ms (absolute values).
%
%   INPUTS:
%     results - Struct array from batch_analyze_1DRF with fields:
%               .group, .temporal_metrics, .is_on, .is_ttl
%     opts    - (Optional) structure
%
%   OUTPUT:
%     fig - Figure handle
%
%   See also EXTRACT_TEMPORAL_METRICS, BATCH_ANALYZE_1DRF

    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'cell_type'),      opts.cell_type      = 'ON'; end
    if ~isfield(opts, 'fig_position'),   opts.fig_position   = [100 100 600 900]; end
    if ~isfield(opts, 'temporal_field'), opts.temporal_field  = 'temporal_metrics'; end

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

    positions = -5:5;  % 11 positions relative to peak
    col_ctrl = [0 0 0];
    col_ttl  = [1 0 0];

    metrics_fields = {'rise_start', 'rise_time', 'decay_time', 'time_to_90'};
    y_labels = {'Rise start (ms)', 'Rise time (ms)', ...
                'Decay time (ms)', 'Time to 90% (ms)'};

    fig = figure('Name', sprintf('Timing — %s Cells', opts.cell_type), ...
        'Position', opts.fig_position);
    tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    for m = 1:4
        nexttile; hold on;

        field = metrics_fields{m};

        % Collect data matrices: n_cells x 11
        ctrl_data = collect_metric_matrix(results(ctrl_mask), field, opts.temporal_field);
        ttl_data  = collect_metric_matrix(results(ttl_mask), field, opts.temporal_field);

        % Plot control
        if ~isempty(ctrl_data)
            plot_mean_sem(positions, ctrl_data, col_ctrl);
        end

        % Plot TTL
        if ~isempty(ttl_data)
            plot_mean_sem(positions, ttl_data, col_ttl);
        end

        xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        xlabel('Position relative to peak');
        ylabel(y_labels{m});
        xlim([-5.5 5.5]);
        box off;
        set(gca, 'TickDir', 'out', 'FontSize', 11);

        % Legend on first panel
        if m == 1
            n_ctrl = size(ctrl_data, 1);
            n_ttl  = size(ttl_data, 1);
            h = [];
            labs = {};
            if n_ctrl > 0
                h(end+1) = plot(NaN, NaN, '-o', 'Color', col_ctrl, 'MarkerFaceColor', col_ctrl);
                labs{end+1} = sprintf('control (n=%d)', n_ctrl);
            end
            if n_ttl > 0
                h(end+1) = plot(NaN, NaN, '-o', 'Color', col_ttl, 'MarkerFaceColor', col_ttl);
                labs{end+1} = sprintf('TTL (n=%d)', n_ttl);
            end
            if ~isempty(h)
                legend(h, labs, 'Location', 'best', 'FontSize', 9);
            end
        end
    end

    sgtitle(sprintf('%s Cells — Temporal Metrics by Position', opts.cell_type), ...
        'FontSize', 14);

end


function data_matrix = collect_metric_matrix(results_subset, field, temporal_field)
% COLLECT_METRIC_MATRIX  Stack a 1x11 metric from each cell into n_cells x 11.
%   temporal_field: name of the struct field in results containing temporal
%                   metrics (default: 'temporal_metrics').

    if nargin < 3, temporal_field = 'temporal_metrics'; end

    n = numel(results_subset);
    if n == 0
        data_matrix = [];
        return;
    end

    data_matrix = NaN(n, 11);
    for k = 1:n
        if isfield(results_subset(k), temporal_field)
            tm = results_subset(k).(temporal_field);
        else
            continue;
        end
        if isstruct(tm) && isfield(tm, field)
            vals = tm.(field);
            if numel(vals) == 11
                data_matrix(k, :) = vals;
            end
        end
    end

end


function plot_mean_sem(x, data_matrix, col)
% PLOT_MEAN_SEM  Plot mean +/- SEM as line with error bars.

    m = mean(data_matrix, 1, 'omitnan');
    n_valid = sum(~isnan(data_matrix), 1);
    s = std(data_matrix, 0, 1, 'omitnan') ./ sqrt(max(n_valid, 1));

    % Only plot positions with at least 2 cells
    valid = n_valid >= 2;
    x_v = x(valid);
    m_v = m(valid);
    s_v = s(valid);

    errorbar(x_v, m_v, s_v, '-o', 'Color', col, 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', col, 'CapSize', 4);

end

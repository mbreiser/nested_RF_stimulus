function fig = plot_three_epoch_distributions(pooled_hist, cell_stats, bin_edges, title_str, opts)
% PLOT_THREE_EPOCH_DISTRIBUTIONS  Histograms + box-dot plots for 3 voltage epochs.
%
%   FIG = PLOT_THREE_EPOCH_DISTRIBUTIONS(POOLED_HIST, CELL_STATS, BIN_EDGES,
%         TITLE_STR) creates a 2x3 figure:
%     Row 1: Pooled voltage histograms (normalized to probability)
%     Row 2: Per-cell box-whisker + dot plots with Wilcoxon rank-sum brackets
%     Columns: Pre-stimulus | During stimulus | Post-stimulus
%
%   INPUTS:
%     pooled_hist - 4x3 cell array of histogram count vectors.
%                   Rows: groups (T4 ctrl, T4 TTL, T5 ctrl, T5 TTL).
%                   Cols: epochs (pre, during, post).
%                   Each cell is a 1xB count vector matching bin_edges.
%     cell_stats  - Struct array (one per cell) with fields:
%                     .group       - string: 'on_control','on_ttl','off_control','off_ttl'
%                     .epoch_median - 1x3 [pre, during, post] median voltage (mV)
%     bin_edges   - 1x(B+1) histogram bin edge vector
%     title_str   - Suptitle string (e.g., 'Bar Sweeps (n=48)')
%     opts        - (Optional) structure with fields:
%                     .fig_position  - [x y w h] in pixels (default: [50 50 1400 700])
%                     .epoch_names   - 1x3 cell of strings
%
%   OUTPUT:
%     fig - Figure handle
%
%   See also PLOT_BASELINE_VOLTAGE_COMPARISON, PLOT_VOLTAGE_HISTOGRAMS

    if nargin < 5, opts = struct(); end
    if ~isfield(opts, 'fig_position'), opts.fig_position = [50 50 1400 700]; end
    if ~isfield(opts, 'epoch_names')
        opts.epoch_names = {'Pre-stimulus', 'During stimulus', 'Post-stimulus'};
    end

    %% Group definitions
    groups      = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
    group_names = {'T4 ctrl', 'T4 TTL^{-}', 'T5 ctrl', 'T5 TTL^{-}'};
    colors = [0 0 0; 1 0 0; 0.4 0.4 0.4; 0.8 0.2 0.2];

    bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

    %% Create figure
    fig = figure('Name', 'Three Epoch Voltage', 'Position', opts.fig_position, 'Color', 'w');
    tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    %% Row 1: Pooled histograms
    for ep = 1:3
        ax = nexttile; hold(ax, 'on');

        for g = 1:4
            counts = pooled_hist{g, ep};
            if isempty(counts) || sum(counts) == 0, continue; end
            prob = counts / sum(counts);

            % Plot as filled stairs
            x_stair = [bin_edges(1), repelem(bin_edges(2:end), 1, 1)];
            y_stair = [repelem(prob, 1, 1), 0];
            % Use stairs for clean look
            stairs(ax, bin_edges(1:end-1), prob, '-', ...
                'Color', colors(g,:), 'LineWidth', 1.2);
            % Filled area
            patch(ax, [bin_centers, fliplr(bin_centers)], ...
                [prob, zeros(size(prob))], colors(g,:), ...
                'FaceAlpha', 0.15, 'EdgeColor', 'none');
        end

        xlabel(ax, 'Vm (mV)');
        if ep == 1, ylabel(ax, 'Probability'); end
        title(ax, opts.epoch_names{ep}, 'FontWeight', 'bold');
        box(ax, 'off');
        set(ax, 'TickDir', 'out', 'FontSize', 10);
    end

    % Build legend from Row 1 (attach to first axis)
    ax1 = tl.Children(end);  % first tile
    % Add invisible lines for legend
    hold(ax1, 'on');
    h_leg = gobjects(4, 1);
    for g = 1:4
        h_leg(g) = plot(ax1, NaN, NaN, '-', 'Color', colors(g,:), 'LineWidth', 2);
    end
    legend(ax1, h_leg, group_names, 'Location', 'northeast', 'FontSize', 8, ...
        'Interpreter', 'tex', 'Box', 'off');

    %% Row 2: Box-whisker + dot plots
    for ep = 1:3
        ax = nexttile; hold(ax, 'on');
        plot_group_boxes(ax, cell_stats, groups, group_names, colors, ep);
        if ep == 1, ylabel(ax, 'Median Vm (mV)'); end
        title(ax, opts.epoch_names{ep}, 'FontWeight', 'bold');
    end

    sgtitle(title_str, 'FontSize', 14, 'FontWeight', 'bold');

end


%% ========================= Local Functions ===============================

function plot_group_boxes(ax, cell_stats, groups, group_names, colors, epoch_idx)
% PLOT_GROUP_BOXES  Box+dot plot with Wilcoxon rank-sum brackets.

    all_vals = [];
    all_grp_idx = [];
    group_data = cell(1, numel(groups));

    for g = 1:numel(groups)
        mask = strcmp({cell_stats.group}, groups{g});
        medians = [cell_stats(mask).epoch_median];
        % epoch_median is 1x3, so when concatenated: 3*n vector
        % Reshape: each cell contributes a 1x3, pick epoch_idx
        n_cells = sum(mask);
        if n_cells == 0
            group_data{g} = [];
            continue;
        end
        vals_mat = reshape(medians, 3, [])';  % n_cells x 3
        vals = vals_mat(:, epoch_idx);
        vals = vals(~isnan(vals));
        group_data{g} = vals;
        n = numel(vals);
        all_vals = [all_vals; vals(:)]; %#ok<AGROW>
        all_grp_idx = [all_grp_idx; repmat(g, n, 1)]; %#ok<AGROW>
    end

    if isempty(all_vals), return; end

    % Box chart
    boxchart(ax, all_grp_idx, all_vals, ...
        'BoxFaceColor', [0.9 0.9 0.9], 'MarkerStyle', 'none', ...
        'BoxWidth', 0.5);

    % Overlay individual dots with group color and jitter
    for g = 1:numel(groups)
        vals = group_data{g};
        if isempty(vals), continue; end
        x = g + 0.15 * (rand(size(vals)) - 0.5);
        scatter(ax, x, vals, 40, colors(g,:), 'filled', ...
            'MarkerFaceAlpha', 0.8);
    end

    set(ax, 'XTick', 1:numel(groups), 'XTickLabel', group_names, ...
        'TickLabelInterpreter', 'tex');
    xlim(ax, [0.5, numel(groups) + 0.5]);
    box(ax, 'off');
    set(ax, 'TickDir', 'out', 'FontSize', 10);

    % Wilcoxon rank-sum: ON ctrl vs TTL, OFF ctrl vs TTL
    yl = ylim(ax);
    y_stat = yl(2) + 0.05 * diff(yl);

    % ON: groups 1 vs 2
    if numel(group_data{1}) >= 2 && numel(group_data{2}) >= 2
        p_on = ranksum(group_data{1}, group_data{2});
        draw_bracket(ax, 1, 2, y_stat, p_on);
    end

    % OFF: groups 3 vs 4
    if numel(group_data{3}) >= 2 && numel(group_data{4}) >= 2
        p_off = ranksum(group_data{3}, group_data{4});
        draw_bracket(ax, 3, 4, y_stat + 0.08 * diff(yl), p_off);
    end

    % Expand ylim for brackets
    ylim(ax, [yl(1), y_stat + 0.20 * diff(yl)]);
end


function draw_bracket(ax, x1, x2, y, p)
% DRAW_BRACKET  Draw a bracket between two x positions with p-value.

    line(ax, [x1 x1 x2 x2], [y - 0.02*abs(y), y, y, y - 0.02*abs(y)], ...
        'Color', 'k', 'LineWidth', 1);
    if p < 0.001
        p_str = 'p < 0.001';
    elseif p < 0.01
        p_str = sprintf('p = %.3f', p);
    else
        p_str = sprintf('p = %.2f', p);
    end
    text(ax, mean([x1 x2]), y, p_str, 'FontSize', 9, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

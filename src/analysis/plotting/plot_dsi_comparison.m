function fig = plot_dsi_comparison(results, opts)
% PLOT_DSI_COMPARISON  Direction selectivity metrics by cell group.
%
%   FIG = PLOT_DSI_COMPARISON(RESULTS) creates a 1x3 figure showing:
%     Left:   DSI (vector sum method, 0-1)
%     Center: DSI (PD-ND method, -1 to 1)
%     Right:  Direction tuning FWHM (degrees)
%   Each panel has 4 groups with box+dot and Wilcoxon rank-sum brackets.
%
%   FIG = PLOT_DSI_COMPARISON(RESULTS, OPTS) uses options:
%     opts.fig_position - [x y w h] in pixels (default: [100 200 1100 400])
%
%   INPUTS:
%     results - Struct array from batch_analyze_1DRF with fields:
%               .group, .dsi_vector, .dsi_pdnd, .dir_tuning_fwhm
%     opts    - (Optional) structure
%
%   OUTPUT:
%     fig - Figure handle
%
%   See also BATCH_ANALYZE_1DRF, COMPUTE_BAR_RESPONSE_METRICS

    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'fig_position'), opts.fig_position = [100 200 1100 400]; end

    groups = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
    group_names = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
    colors = [0 0 0; 1 0 0; 0.4 0.4 0.4; 0.8 0.2 0.2];

    fig = figure('Name', 'Direction Selectivity Comparison', ...
        'Position', opts.fig_position);
    t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % --- Panel 1: DSI (vector sum) ---
    nexttile; hold on;
    plot_group_with_stats(results, groups, group_names, colors, 'dsi_vector');
    ylabel('DSI (vector sum)');
    title('Direction Selectivity');

    % --- Panel 2: DSI (PD-ND) ---
    nexttile; hold on;
    plot_group_with_stats(results, groups, group_names, colors, 'dsi_pdnd');
    ylabel('DSI (PD-ND)');
    title('PD vs ND Selectivity');

    % --- Panel 3: Direction tuning width ---
    nexttile; hold on;
    plot_group_with_stats(results, groups, group_names, colors, 'dir_tuning_fwhm');
    ylabel('FWHM (\circ)');
    title('Direction Tuning Width');

    title(t, 'Direction Selectivity by Group', 'FontSize', 14);

end


function plot_group_with_stats(results, groups, group_names, colors, field)
% PLOT_GROUP_WITH_STATS  Box+dot plot with Wilcoxon rank-sum brackets.

    all_vals = [];
    all_grp_idx = [];
    group_data = cell(1, numel(groups));

    for g = 1:numel(groups)
        mask = strcmp({results.group}, groups{g});
        vals = [results(mask).(field)];
        vals = vals(~isnan(vals));
        group_data{g} = vals;
        n = numel(vals);
        all_vals = [all_vals, vals]; %#ok<AGROW>
        all_grp_idx = [all_grp_idx, repmat(g, 1, n)]; %#ok<AGROW>
    end

    if isempty(all_vals)
        return;
    end

    % Box chart
    boxchart(all_grp_idx(:), all_vals(:), ...
        'BoxFaceColor', [0.9 0.9 0.9], 'MarkerStyle', 'none', ...
        'BoxWidth', 0.5);

    % Overlay individual dots with group color and jitter
    for g = 1:numel(groups)
        idx = all_grp_idx == g;
        vals = all_vals(idx);
        x = g + 0.15 * (rand(size(vals)) - 0.5);
        scatter(x, vals, 50, colors(g, :), 'filled', ...
            'MarkerFaceAlpha', 0.8);
    end

    set(gca, 'XTick', 1:numel(groups), 'XTickLabel', group_names);
    xlim([0.5, numel(groups) + 0.5]);
    box off;
    set(gca, 'TickDir', 'out', 'FontSize', 11);

    % Statistical annotation: Wilcoxon rank-sum ctrl vs TTL
    yl = ylim;
    y_stat = yl(2) + 0.05 * diff(yl);

    % ON: ctrl (1) vs TTL (2)
    if numel(group_data{1}) >= 2 && numel(group_data{2}) >= 2
        p_on = ranksum(group_data{1}, group_data{2});
        draw_bracket(1, 2, y_stat, p_on);
    end

    % OFF: ctrl (3) vs TTL (4)
    if numel(group_data{3}) >= 2 && numel(group_data{4}) >= 2
        p_off = ranksum(group_data{3}, group_data{4});
        draw_bracket(3, 4, y_stat + 0.08 * diff(yl), p_off);
    end

    % Expand ylim to accommodate brackets
    ylim([yl(1), y_stat + 0.18 * diff(yl)]);

end


function draw_bracket(x1, x2, y, p)
% DRAW_BRACKET  Draw a bracket between two x positions with p-value.

    line([x1 x1 x2 x2], [y - 0.02*abs(y), y, y, y - 0.02*abs(y)], ...
        'Color', 'k', 'LineWidth', 1);
    if p < 0.001
        p_str = 'p < 0.001';
    elseif p < 0.01
        p_str = sprintf('p = %.3f', p);
    else
        p_str = sprintf('p = %.2f', p);
    end
    text(mean([x1 x2]), y, p_str, 'FontSize', 9, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

end

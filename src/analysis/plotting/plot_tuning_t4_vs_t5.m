function fig = plot_tuning_t4_vs_t5(results, opts)
% PLOT_TUNING_T4_VS_T5  T4 vs T5 tuning width and DSI comparison.
%
%   FIG = PLOT_TUNING_T4_VS_T5(RESULTS) creates a 1x2 figure comparing
%   direction tuning FWHM and DSI between T4 (ON) and T5 (OFF) cells,
%   split by treatment (control vs TTL).
%
%   FIG = PLOT_TUNING_T4_VS_T5(RESULTS, OPTS) uses options:
%     opts.fig_position - [x y w h] in pixels (default: [100 200 900 400])
%
%   Plot layout:
%     Panel 1: Direction tuning FWHM — 4 boxes: ON-ctrl, OFF-ctrl, ON-TTL, OFF-TTL
%     Panel 2: DSI (vector sum) — same 4 boxes
%
%   Statistical brackets: Wilcoxon rank-sum comparing ON vs OFF within
%   each treatment (ctrl: ON vs OFF, TTL: ON vs OFF).
%
%   Colors: T4 (ON) = blue, T5 (OFF) = orange. Ctrl = filled, TTL = open.
%
%   INPUTS:
%     results - Struct array from batch_analyze_1DRF with fields:
%               .group, .is_on, .is_ttl, .dir_tuning_fwhm, .dsi_vector
%     opts    - (Optional) structure
%
%   OUTPUT:
%     fig - Figure handle
%
%   See also PLOT_DSI_COMPARISON, PLOT_POLAR_POPULATION, BATCH_ANALYZE_1DRF

    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'fig_position'), opts.fig_position = [100 200 900 400]; end

    % Group order: ON-ctrl, OFF-ctrl, ON-TTL, OFF-TTL
    % This puts ON vs OFF side by side within each treatment for the bracket
    group_masks = {
        [results.is_on] & ~[results.is_ttl], ...   % ON ctrl
        ~[results.is_on] & ~[results.is_ttl], ...  % OFF ctrl
        [results.is_on] &  [results.is_ttl], ...   % ON TTL
        ~[results.is_on] &  [results.is_ttl]       % OFF TTL
    };

    group_names = {'T4 ctrl', 'T5 ctrl', 'T4 TTL', 'T5 TTL'};

    % Colors: T4=blue, T5=orange. Ctrl=filled, TTL=open (lighter)
    col_t4 = [0.12 0.47 0.71];      % blue
    col_t5 = [0.85 0.37 0.01];      % orange
    col_t4_light = [0.55 0.73 0.90]; % light blue
    col_t5_light = [0.95 0.65 0.40]; % light orange
    dot_colors = [col_t4; col_t5; col_t4_light; col_t5_light];

    fig = figure('Name', 'T4 vs T5 Tuning Comparison', ...
        'Position', opts.fig_position);
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % --- Panel 1: Direction tuning FWHM ---
    nexttile; hold on;
    plot_t4t5_metric(results, group_masks, group_names, dot_colors, 'dir_tuning_fwhm');
    ylabel('FWHM (\circ)');
    title('Direction Tuning Width');

    % --- Panel 2: DSI (vector sum) ---
    nexttile; hold on;
    plot_t4t5_metric(results, group_masks, group_names, dot_colors, 'dsi_vector');
    ylabel('DSI (vector sum)');
    title('Direction Selectivity Index');

    title(t, 'T4 (ON) vs T5 (OFF) — Tuning Comparison', 'FontSize', 14);

end


function plot_t4t5_metric(results, group_masks, group_names, dot_colors, field)
% PLOT_T4T5_METRIC  Box+dot with ON-vs-OFF rank-sum within treatment.

    n_groups = numel(group_masks);
    all_vals = [];
    all_grp_idx = [];
    group_data = cell(1, n_groups);

    for g = 1:n_groups
        mask = group_masks{g};
        vals = [results(mask).(field)];
        vals = vals(~isnan(vals));
        group_data{g} = vals;
        n = numel(vals);
        all_vals = [all_vals, vals]; %#ok<AGROW>
        all_grp_idx = [all_grp_idx, repmat(g, 1, n)]; %#ok<AGROW>
    end

    if isempty(all_vals), return; end

    % Box chart
    boxchart(all_grp_idx(:), all_vals(:), ...
        'BoxFaceColor', [0.9 0.9 0.9], 'MarkerStyle', 'none', ...
        'BoxWidth', 0.5);

    % Overlay individual dots with group color and jitter
    for g = 1:n_groups
        idx = all_grp_idx == g;
        vals = all_vals(idx);
        x = g + 0.15 * (rand(size(vals)) - 0.5);
        scatter(x, vals, 50, dot_colors(g, :), 'filled', ...
            'MarkerFaceAlpha', 0.8);
    end

    set(gca, 'XTick', 1:n_groups, 'XTickLabel', group_names);
    xlim([0.5, n_groups + 0.5]);
    box off;
    set(gca, 'TickDir', 'out', 'FontSize', 11);

    % Statistical annotation: ON vs OFF within each treatment
    yl = ylim;
    y_stat = yl(2) + 0.05 * diff(yl);

    % Ctrl: ON (1) vs OFF (2)
    if numel(group_data{1}) >= 2 && numel(group_data{2}) >= 2
        p_ctrl = ranksum(group_data{1}, group_data{2});
        draw_bracket(1, 2, y_stat, p_ctrl);
    end

    % TTL: ON (3) vs OFF (4)
    if numel(group_data{3}) >= 2 && numel(group_data{4}) >= 2
        p_ttl = ranksum(group_data{3}, group_data{4});
        draw_bracket(3, 4, y_stat + 0.08 * diff(yl), p_ttl);
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

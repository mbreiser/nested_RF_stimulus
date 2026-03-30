function fig = plot_dsi_ar_summary(all_cells, method_name, opts)
% PLOT_DSI_AR_SUMMARY  Publication-quality DSI + aspect ratio box/dot plots.
%
%   FIG = PLOT_DSI_AR_SUMMARY(ALL_CELLS, METHOD_NAME) creates a 1x2 figure:
%     Left:  DSI = (PD_max - ND_max) / PD_max
%     Right: Aspect ratio = PD_max / mean(ortho)
%   Each panel shows 4 groups (ON ctrl, ON TTL, OFF ctrl, OFF TTL) with
%   box charts, jittered dots, and 3 Wilcoxon rank-sum brackets:
%     1) ON ctrl vs ON TTL
%     2) OFF ctrl vs OFF TTL
%     3) ON ctrl vs OFF ctrl
%
%   FIG = PLOT_DSI_AR_SUMMARY(ALL_CELLS, METHOD_NAME, OPTS) uses:
%     opts.fig_position - [x y w h] in pixels (default: [100 200 900 420])
%
%   INPUTS:
%     all_cells   - Struct array with .group and .methods.(method_name)
%     method_name - String: 'vec_nn', 'sym_nn', 'vec_interp', 'sym_interp'
%     opts        - (Optional) structure
%
%   OUTPUT:
%     fig - Figure handle
%
%   See also PLOT_POLAR_SUMMARY_SINGLE, PLOT_DSI_COMPARISON

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'fig_position'), opts.fig_position = [100 200 900 420]; end

    % --- Group definitions ---
    groups     = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
    grp_labels = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
    dot_colors = [0    0    0;       % black
                  0.85 0    0;       % red
                  0.15 0.35 0.70;    % blue
                  0.90 0.50 0.10];   % orange

    % --- Extract per-cell DSI and AR ---
    grp_dsi = cell(1, 4);
    grp_ar  = cell(1, 4);
    grp_n   = zeros(1, 4);

    for g = 1:4
        mask = strcmp({all_cells.group}, groups{g});
        cg   = all_cells(mask);
        ng   = numel(cg);
        grp_n(g) = ng;

        dsi_v = NaN(1, ng);
        ar_v  = NaN(1, ng);

        for k = 1:ng
            m = cg(k).methods.(method_name);
            ar_v(k) = m.aspect_ratio;

            % DSI = (PD_max - ND_max) / PD_max
            aln  = m.aligned;
            ang  = aln(:, 1);
            resp = aln(:, 2);
            aw   = [ang(end) - 2*pi; ang; ang(1) + 2*pi];
            rw   = [resp(end);        resp; resp(1)];
            rpd  = interp1(aw, rw, pi/2,   'linear');
            rnd  = interp1(aw, rw, 3*pi/2, 'linear');
            if rpd > 0
                dsi_v(k) = (rpd - rnd) / rpd;
            end
        end

        grp_dsi{g} = dsi_v(~isnan(dsi_v));
        grp_ar{g}  = ar_v(~isnan(ar_v));
    end

    % --- Create figure ---
    method_display = strrep(method_name, '_', ' ');
    method_display(1) = upper(method_display(1));

    fig = figure('Name', sprintf('DSI & AR Summary (%s)', method_display), ...
        'Position', opts.fig_position, 'Color', 'w');
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % --- Panel 1: DSI ---
    nexttile; hold on;
    plot_metric_panel(grp_dsi, grp_labels, grp_n, dot_colors);
    ylabel('DSI  (PD_{max} - ND_{max}) / PD_{max}', ...
        'FontSize', 11, 'Interpreter', 'tex');
    title('Direction Selectivity Index', 'FontSize', 13);

    % --- Panel 2: Aspect Ratio ---
    nexttile; hold on;
    plot_metric_panel(grp_ar, grp_labels, grp_n, dot_colors);
    ylabel('Aspect Ratio  PD_{max} / mean(ortho)', ...
        'FontSize', 11, 'Interpreter', 'tex');
    title('Tuning Aspect Ratio', 'FontSize', 13);

    % --- Super title ---
    title(t, sprintf('Direction Selectivity — %s alignment', method_display), ...
        'FontSize', 14, 'FontWeight', 'bold');

end


%% ========================= Local Helpers ============================

function plot_metric_panel(grp_data, grp_labels, grp_n, dot_colors)
% PLOT_METRIC_PANEL  Box+dot+3-bracket panel for one metric.

    n_groups = numel(grp_data);

    % Assemble data for boxchart
    all_vals    = [];
    all_grp_idx = [];
    for g = 1:n_groups
        vals = grp_data{g};
        n = numel(vals);
        all_vals    = [all_vals, vals]; %#ok<AGROW>
        all_grp_idx = [all_grp_idx, repmat(g, 1, n)]; %#ok<AGROW>
    end

    if isempty(all_vals), return; end

    % Box chart
    boxchart(all_grp_idx(:), all_vals(:), ...
        'BoxFaceColor', [0.9 0.9 0.9], 'MarkerStyle', 'none', ...
        'BoxWidth', 0.5);

    % Jittered dots
    for g = 1:n_groups
        idx  = all_grp_idx == g;
        vals = all_vals(idx);
        x = g + 0.15 * (rand(size(vals)) - 0.5);
        scatter(x, vals, 50, dot_colors(g, :), 'filled', ...
            'MarkerFaceAlpha', 0.8);
    end

    % X-axis labels with n-count (single-line to avoid MATLAB tick rendering issues)
    tick_labels = cell(1, n_groups);
    for g = 1:n_groups
        tick_labels{g} = sprintf('%s (n=%d)', grp_labels{g}, grp_n(g));
    end
    set(gca, 'XTick', 1:n_groups, 'XTickLabel', tick_labels);
    xlim([0.3, n_groups + 0.7]);
    box off;
    set(gca, 'TickDir', 'out', 'FontSize', 11);

    % --- 3 statistical brackets ---
    yl = ylim;
    gap = 0.06 * diff(yl);

    % Tier 1 (low): within-treatment comparisons
    y_tier1 = yl(2) + gap;

    % Bracket 1: ON ctrl (1) vs ON TTL (2)
    if numel(grp_data{1}) >= 2 && numel(grp_data{2}) >= 2
        p = ranksum(grp_data{1}(:), grp_data{2}(:));
        draw_bracket(1, 2, y_tier1, p);
    end

    % Bracket 2: OFF ctrl (3) vs OFF TTL (4)
    if numel(grp_data{3}) >= 2 && numel(grp_data{4}) >= 2
        p = ranksum(grp_data{3}(:), grp_data{4}(:));
        draw_bracket(3, 4, y_tier1, p);
    end

    % Tier 2 (high): cross-type comparison
    y_tier2 = y_tier1 + 2.5 * gap;

    % Bracket 3: ON ctrl (1) vs OFF ctrl (3)
    if numel(grp_data{1}) >= 2 && numel(grp_data{3}) >= 2
        p = ranksum(grp_data{1}(:), grp_data{3}(:));
        draw_bracket(1, 3, y_tier2, p);
    end

    % Expand y-axis to fit brackets
    ylim([yl(1), y_tier2 + 2.5 * gap]);

end


function draw_bracket(x1, x2, y, p)
% DRAW_BRACKET  Bracket between two x-positions with p-value annotation.

    tick_h = 0.015 * abs(y);
    if tick_h == 0, tick_h = 0.01; end

    line([x1 x1 x2 x2], [y - tick_h, y, y, y - tick_h], ...
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

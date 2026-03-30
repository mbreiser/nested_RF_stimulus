function fig = plot_aspect_ratio_comparison(results, opts)
% PLOT_ASPECT_RATIO_COMPARISON  Tuning aspect ratio: T4 vs T5, ctrl vs TTL.
%
%   FIG = PLOT_ASPECT_RATIO_COMPARISON(RESULTS) creates a figure comparing
%   the directional tuning aspect ratio across T4 (ON) and T5 (OFF) cells,
%   split by treatment (control vs TTL).
%
%   Aspect ratio = PD response / mean(ortho_cw, ortho_ccw)
%     PD is the preferred direction (aligned to pi/2)
%     ortho_cw is 90 deg clockwise from PD (0 rad)
%     ortho_ccw is 90 deg counter-clockwise from PD (pi rad)
%
%   Values > 1 indicate the tuning is elongated along PD-ND; ~1 means
%   roughly circular tuning.
%
%   FIG = PLOT_ASPECT_RATIO_COMPARISON(RESULTS, OPTS) uses options:
%     opts.fig_position - [x y w h] in pixels (default: [100 200 500 450])
%
%   Statistical brackets (Wilcoxon rank-sum):
%     T4 ctrl vs T5 ctrl  (cell type within control)
%     T4 TTL  vs T5 TTL   (cell type within TTL)
%     T4 ctrl vs T4 TTL   (treatment within ON)
%     T5 ctrl vs T5 TTL   (treatment within OFF)
%
%   Colors: T4 (ON) = blue, T5 (OFF) = orange.
%           Control = filled, TTL = lighter shade.
%
%   INPUTS:
%     results - Struct array from batch_analyze_1DRF with field:
%               .max_v_aligned (16x2 PD-aligned tuning data)
%               .is_on, .is_ttl
%     opts    - (Optional) structure
%
%   OUTPUT:
%     fig - Figure handle
%
%   See also PLOT_TUNING_T4_VS_T5, PLOT_DSI_COMPARISON, BATCH_ANALYZE_1DRF

    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'fig_position'), opts.fig_position = [100 200 500 450]; end

    % Group order: T4 ctrl, T5 ctrl, T4 TTL, T5 TTL
    group_masks = {
        [results.is_on] & ~[results.is_ttl], ...   % T4 ctrl
        ~[results.is_on] & ~[results.is_ttl], ...  % T5 ctrl
        [results.is_on] &  [results.is_ttl], ...   % T4 TTL
        ~[results.is_on] &  [results.is_ttl]       % T5 TTL
    };

    group_names = {'T4 ctrl', 'T5 ctrl', 'T4 TTL', 'T5 TTL'};

    % Colors: T4=blue, T5=orange. Ctrl=filled, TTL=lighter
    col_t4       = [0.12 0.47 0.71];
    col_t5       = [0.85 0.37 0.01];
    col_t4_light = [0.55 0.73 0.90];
    col_t5_light = [0.95 0.65 0.40];
    dot_colors = [col_t4; col_t5; col_t4_light; col_t5_light];

    % Compute aspect ratio for each cell from max_v_aligned
    n = numel(results);
    ar_vals = NaN(1, n);
    for k = 1:n
        d = results(k).max_v_aligned;
        if isnumeric(d) && size(d, 1) == 16 && size(d, 2) == 2
            ar_vals(k) = compute_ar(d);
        end
    end

    % Build data vectors per group
    n_groups = numel(group_masks);
    all_vals = [];
    all_grp_idx = [];
    group_data = cell(1, n_groups);

    for g = 1:n_groups
        mask = group_masks{g};
        vals = ar_vals(mask);
        vals = vals(~isnan(vals));
        group_data{g} = vals;
        ng = numel(vals);
        all_vals = [all_vals, vals]; %#ok<AGROW>
        all_grp_idx = [all_grp_idx, repmat(g, 1, ng)]; %#ok<AGROW>
    end

    if isempty(all_vals)
        fig = figure('Name', 'No data');
        return;
    end

    fig = figure('Name', 'Tuning Aspect Ratio Comparison', ...
        'Position', opts.fig_position);
    hold on;

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

    % Reference line at aspect ratio = 1 (circular tuning)
    yline(1, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

    set(gca, 'XTick', 1:n_groups, 'XTickLabel', group_names);
    xlim([0.5, n_groups + 0.5]);
    ylabel('Aspect Ratio (PD / mean ortho)');
    title('Directional Tuning Aspect Ratio');
    box off;
    set(gca, 'TickDir', 'out', 'FontSize', 11);

    % --- Statistical brackets (4 comparisons) ---
    yl = ylim;
    dy = 0.07 * diff(yl);
    y0 = yl(2) + 0.03 * diff(yl);

    % Row 1: within-treatment cell type comparisons (adjacent groups)
    % T4 ctrl (1) vs T5 ctrl (2)
    if numel(group_data{1}) >= 2 && numel(group_data{2}) >= 2
        p = ranksum(group_data{1}, group_data{2});
        draw_bracket(1, 2, y0, p);
    end

    % T4 TTL (3) vs T5 TTL (4)
    if numel(group_data{3}) >= 2 && numel(group_data{4}) >= 2
        p = ranksum(group_data{3}, group_data{4});
        draw_bracket(3, 4, y0, p);
    end

    % Row 2: within-cell-type treatment comparisons (spanning groups)
    % T4 ctrl (1) vs T4 TTL (3)
    if numel(group_data{1}) >= 2 && numel(group_data{3}) >= 2
        p = ranksum(group_data{1}, group_data{3});
        draw_bracket(1, 3, y0 + dy, p);
    end

    % T5 ctrl (2) vs T5 TTL (4)
    if numel(group_data{2}) >= 2 && numel(group_data{4}) >= 2
        p = ranksum(group_data{2}, group_data{4});
        draw_bracket(2, 4, y0 + 2*dy, p);
    end

    % Expand ylim to accommodate all brackets
    ylim([yl(1), y0 + 3.2*dy]);

end


function ar = compute_ar(d_aligned)
% COMPUTE_AR  Aspect ratio from 16x2 PD-aligned tuning data.
%   PD at pi/2; orthogonal directions at 0 and pi.

    angles = d_aligned(:, 1);
    resps  = d_aligned(:, 2);

    [~, pd_idx]     = min(abs(angles - pi/2));
    [~, ortho1_idx] = min(abs(angles - 0));
    [~, ortho2_idx] = min(abs(angles - pi));

    pd_resp    = resps(pd_idx);
    ortho_mean = mean([resps(ortho1_idx), resps(ortho2_idx)]);

    if ortho_mean > 0
        ar = pd_resp / ortho_mean;
    else
        ar = NaN;
    end

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

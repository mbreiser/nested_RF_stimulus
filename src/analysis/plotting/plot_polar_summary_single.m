function fig = plot_polar_summary_single(all_cells, method_name, opts)
% PLOT_POLAR_SUMMARY_SINGLE  Publication-quality 1x2 polar tuning summary.
%
%   FIG = PLOT_POLAR_SUMMARY_SINGLE(ALL_CELLS, METHOD_NAME) draws a 1x2
%   figure with polar population mean +/- SEM tuning curves:
%     Left panel:  ON cells  — ctrl (black) vs TTL (red)
%     Right panel: OFF cells — ctrl (blue)  vs TTL (orange)
%
%   FIG = PLOT_POLAR_SUMMARY_SINGLE(ALL_CELLS, METHOD_NAME, OPTS) uses:
%     opts.fig_position - [x y w h] in pixels (default: [100 100 1100 520])
%
%   INPUTS:
%     all_cells   - Struct array with .group and .methods.(method_name).aligned
%     method_name - String: 'vec_nn', 'sym_nn', 'vec_interp', 'sym_interp'
%     opts        - (Optional) structure
%
%   OUTPUT:
%     fig - Figure handle
%
%   See also PLOT_POLAR_4METHOD_SUMMARY, PLOT_DSI_AR_SUMMARY

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'fig_position'), opts.fig_position = [100 100 1100 520]; end

    method_display = strrep(method_name, '_', ' ');
    method_display(1) = upper(method_display(1));

    % --- Group definitions ---
    groups     = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
    grp_labels = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
    line_colors = [0    0    0;       % black  (ON ctrl)
                   0.85 0    0;       % red    (ON TTL)
                   0.15 0.35 0.70;    % blue   (OFF ctrl)
                   0.90 0.50 0.10];   % orange (OFF TTL)
    fill_colors = [0.60 0.60 0.60;    % light gray
                   1.00 0.70 0.70;    % light red
                   0.65 0.75 0.90;    % light blue
                   1.00 0.87 0.72];   % light orange

    % --- Collect per-group statistics ---
    grp_stats = cell(1, 4);
    grp_n     = zeros(1, 4);

    for g = 1:4
        mask  = strcmp({all_cells.group}, groups{g});
        cg    = all_cells(mask);
        ng    = numel(cg);
        grp_n(g) = ng;

        aln_list = cell(1, ng);
        for k = 1:ng
            aln_list{k} = cg(k).methods.(method_name).aligned;
        end
        grp_stats{g} = compute_mean_sem_local(aln_list);
    end

    % --- Determine shared rmax across all 4 groups ---
    rmax = 0;
    for g = 1:4
        s = grp_stats{g};
        if ~isempty(s.center)
            rmax = max(rmax, max(s.center + s.spread, [], 'omitnan'));
        end
    end
    rmax = max(rmax, 0.1) * 1.15;

    % --- Create figure ---
    fig = figure('Name', sprintf('Direction Tuning (%s)', method_display), ...
        'Position', opts.fig_position, 'Color', 'w');

    % Panel definitions: {panel_title, [ctrl_idx, ttl_idx]}
    panels = {
        'ON cells (T4)',  [1, 2];
        'OFF cells (T5)', [3, 4]
    };

    for pi = 1:2
        panel_grps = panels{pi, 2};   % [ctrl_idx, ttl_idx]
        g_ctrl = panel_grps(1);
        g_ttl  = panel_grps(2);

        % Manual axes positioning for 1x2 layout
        x0 = 0.04 + (pi - 1) * 0.50;
        ax = axes('Position', [x0, 0.06, 0.46, 0.84]); %#ok<LAXES>
        hold(ax, 'on');
        axis(ax, 'equal');
        set(ax, 'Clipping', 'off');

        % --- Draw polar grid ---
        draw_polar_grid_local(ax, rmax);

        % --- Legend handles (invisible) ---
        hleg = gobjects(2, 1);
        hleg(1) = plot(ax, NaN, NaN, '-', ...
            'Color', line_colors(g_ctrl, :), 'LineWidth', 2);
        hleg(2) = plot(ax, NaN, NaN, '-', ...
            'Color', line_colors(g_ttl, :), 'LineWidth', 2);

        % --- Draw shaded traces: ctrl behind, TTL in front ---
        s = grp_stats{g_ctrl};
        if ~isempty(s.center)
            draw_shaded_polar_local(ax, s.theta, s.center, s.spread, ...
                line_colors(g_ctrl, :), fill_colors(g_ctrl, :), 0.30, 2.0);
        end
        s = grp_stats{g_ttl};
        if ~isempty(s.center)
            draw_shaded_polar_local(ax, s.theta, s.center, s.spread, ...
                line_colors(g_ttl, :), fill_colors(g_ttl, :), 0.30, 2.0);
        end

        % --- Axis formatting ---
        lim = rmax * 1.12;
        set(ax, 'XLim', [-lim, lim], 'YLim', [-lim, lim]);
        set(ax, 'XColor', 'none', 'YColor', 'none');

        % --- Legend ---
        leg_labs = {
            sprintf('%s (n=%d)', grp_labels{g_ctrl}, grp_n(g_ctrl)), ...
            sprintf('%s (n=%d)', grp_labels{g_ttl},  grp_n(g_ttl))
        };
        legend(hleg, leg_labs, 'Location', 'northwest', 'FontSize', 10, ...
            'Box', 'off');

        % --- Panel title ---
        title(ax, panels{pi, 1}, 'FontSize', 13, 'FontWeight', 'bold');
    end

    % --- Super title ---
    sgtitle(sprintf('Direction Tuning — %s alignment', method_display), ...
        'FontSize', 14, 'FontWeight', 'bold');

end


%% ========================= Local Helpers ============================

function stats = compute_mean_sem_local(aligned_list)
% Stack aligned data and compute mean +/- SEM.

    stats.theta  = [];
    stats.center = [];
    stats.spread = [];
    stats.n      = 0;

    if isempty(aligned_list), return; end

    vals = [];
    for k = 1:numel(aligned_list)
        d = aligned_list{k};
        if isnumeric(d) && size(d, 2) == 2 && size(d, 1) >= 2
            if isempty(stats.theta)
                stats.theta = d(:, 1);
            end
            vals = [vals, d(:, 2)]; %#ok<AGROW>
        end
    end

    if isempty(vals), return; end

    stats.n      = size(vals, 2);
    stats.center = mean(vals, 2, 'omitnan');
    stats.spread = std(vals, 0, 2, 'omitnan') ./ sqrt(stats.n);

end


function draw_polar_grid_local(ax, rmax)
% Draw concentric circles and radial spokes on Cartesian axes.

    grid_color = [0.88 0.88 0.88];
    label_color = [0.50 0.50 0.50];

    % Concentric circles
    n_rings = 4;
    r_ticks = linspace(0, rmax, n_rings + 1);
    r_ticks = r_ticks(2:end);  % skip origin
    th_circle = linspace(0, 2*pi, 181);

    for ri = 1:numel(r_ticks)
        [xc, yc] = pol2cart(th_circle, r_ticks(ri));
        plot(ax, xc, yc, '-', 'Color', grid_color, 'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
    end

    % Radial spokes every 45 degrees
    spoke_ang = deg2rad(0:45:315);
    for si = 1:numel(spoke_ang)
        [xs, ys] = pol2cart(spoke_ang(si), rmax);
        plot(ax, [0 xs], [0 ys], '-', 'Color', grid_color, 'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
    end

    % Cardinal direction labels
    card_angles = [0, 90, 180, 270];
    card_labels = {'0', '90 (PD)', '180', '270 (ND)'};
    deg_sym = char(176);

    for li = 1:4
        a = deg2rad(card_angles(li));
        [xl, yl] = pol2cart(a, rmax * 1.08);

        ha = 'center';  va = 'middle';
        if card_angles(li) == 0,   ha = 'left';   end
        if card_angles(li) == 180, ha = 'right';  end
        if card_angles(li) == 90,  va = 'bottom'; end
        if card_angles(li) == 270, va = 'top';    end

        text(ax, xl, yl, [card_labels{li}, deg_sym], 'FontSize', 9, ...
            'HorizontalAlignment', ha, 'VerticalAlignment', va, ...
            'Color', label_color, 'HandleVisibility', 'off');
    end

    % R-tick labels along 0-degree spoke
    for ri = 1:numel(r_ticks)
        text(ax, r_ticks(ri), rmax * 0.06, sprintf('%.1f', r_ticks(ri)), ...
            'FontSize', 7, 'Color', [0.6 0.6 0.6], ...
            'HorizontalAlignment', 'center', 'HandleVisibility', 'off');
    end

end


function draw_shaded_polar_local(ax, theta, center, spread, ...
    line_col, fill_col, alpha_val, lw)
% Draw shaded polar band + center line on Cartesian axes.

    th = theta(:);
    m  = center(:);
    b  = spread(:);

    upper = max(m + b, 0);
    lower = max(m - b, 0);

    % Close loop
    th_c    = [th; th(1) + 2*pi];
    upper_c = [upper; upper(1)];
    lower_c = [lower; lower(1)];
    m_c     = [m; m(1)];

    % Polygon for SEM band
    th_poly = [th_c; flipud(th_c)];
    r_poly  = [upper_c; flipud(lower_c)];
    [xp, yp] = pol2cart(th_poly, r_poly);

    patch('XData', xp, 'YData', yp, ...
        'FaceColor', fill_col, 'FaceAlpha', alpha_val, ...
        'EdgeColor', 'none', 'Parent', ax, 'HandleVisibility', 'off');

    % Center line
    [xl, yl] = pol2cart(th_c, m_c);
    plot(ax, xl, yl, '-', 'Color', line_col, 'LineWidth', lw, ...
        'HandleVisibility', 'off');

end

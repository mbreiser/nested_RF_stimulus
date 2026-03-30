function fig = plot_polar_4method_summary(all_cells, method_names, method_labels)
% PLOT_POLAR_4METHOD_SUMMARY  2x2 polar summary with 4 groups per panel.
%
%   Creates a 2x2 figure with one panel per PD alignment method. Each panel
%   shows population mean+/-SEM polar tuning for 4 groups (ON ctrl, ON TTL,
%   OFF ctrl, OFF TTL). Annotated with DSI = (PDmax-NDmax)/PDmax and
%   aspect ratio per group, with Wilcoxon rank-sum significance for
%   ctrl vs TTL within ON and OFF.
%
%   INPUTS:
%     all_cells     - struct array with .group, .methods fields
%                     (from ds_alignment_exploration.m)
%     method_names  - cell array of method field names
%                     e.g. {'vec_nn', 'sym_nn', 'vec_interp', 'sym_interp'}
%     method_labels - cell array of display labels
%                     e.g. {'VecNN', 'SymNN', 'VecInterp', 'SymInterp'}
%
%   OUTPUT:
%     fig - Figure handle
%
%   METRICS:
%     DSI  = (PD_max - ND_max) / PD_max   (PD at pi/2, ND at 3*pi/2)
%     AR   = PD_max / mean(ortho)          (ortho at 0 and pi)
%
%   See also COMPUTE_PD_FOUR_METHODS, DS_ALIGNMENT_EXPLORATION

    groups    = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
    grp_short = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
    pm = char(177);  % plus-minus symbol

    % Color scheme: ON ctrl=black, ON TTL=red, OFF ctrl=blue, OFF TTL=orange
    line_colors = [0    0    0;        % ON ctrl
                   0.85 0    0;        % ON TTL
                   0.15 0.35 0.70;     % OFF ctrl
                   0.90 0.50 0.10];    % OFF TTL
    fill_colors = [0.78 0.78 0.78;     % gray
                   1.00 0.72 0.72;     % light red
                   0.72 0.82 0.95;     % light blue
                   1.00 0.87 0.72];    % light orange

    fig = figure('Name', 'DS Method Comparison — 2x2 Polar Summary', ...
        'Position', [50 30 1200 1400], 'Color', 'w');

    n_methods = numel(method_names);

    for mi = 1:n_methods
        % --- Manual subplot positioning: leave room below for annotation ---
        row = ceil(mi / 2);       % 1 or 2
        col = mod(mi - 1, 2) + 1; % 1 or 2
        x0 = 0.05 + (col - 1) * 0.50;
        y0 = 0.10 + (2 - row) * 0.44;
        ax = axes('Position', [x0, y0, 0.44, 0.42]); %#ok<LAXES>
        hold(ax, 'on');
        axis(ax, 'equal');
        set(ax, 'Clipping', 'off');

        % --- Collect per-group data for this method ---
        grp_stats = cell(1, 4);
        grp_dsi   = cell(1, 4);
        grp_ar    = cell(1, 4);
        grp_n     = zeros(1, 4);

        for g = 1:4
            mask  = strcmp({all_cells.group}, groups{g});
            cg    = all_cells(mask);
            ng    = numel(cg);
            grp_n(g) = ng;

            aln_list = cell(1, ng);
            dsi_v    = NaN(1, ng);
            ar_v     = NaN(1, ng);

            for k = 1:ng
                m = cg(k).methods.(method_names{mi});
                aln_list{k} = m.aligned;
                ar_v(k) = m.aspect_ratio;

                % User DSI: (PDmax - NDmax) / PDmax
                aln = m.aligned;
                ang = aln(:, 1);
                resp = aln(:, 2);
                aw = [ang(end) - 2*pi; ang; ang(1) + 2*pi];
                rw = [resp(end);        resp; resp(1)];
                rpd = interp1(aw, rw, pi/2,   'linear');
                rnd = interp1(aw, rw, 3*pi/2, 'linear');
                if rpd > 0
                    dsi_v(k) = (rpd - rnd) / rpd;
                end
            end

            grp_stats{g} = compute_mean_sem_local(aln_list);
            grp_dsi{g}   = dsi_v(~isnan(dsi_v));
            grp_ar{g}    = ar_v(~isnan(ar_v));
        end

        % --- Determine rmax ---
        rmax = 0;
        for g = 1:4
            s = grp_stats{g};
            if ~isempty(s.center)
                rmax = max(rmax, max(s.center + s.spread, [], 'omitnan'));
            end
        end
        rmax = max(rmax, 0.1) * 1.15;

        % --- Draw polar grid ---
        draw_polar_grid_local(ax, rmax);

        % --- Create invisible legend handles ---
        hleg = gobjects(4, 1);
        for g = 1:4
            hleg(g) = plot(ax, NaN, NaN, '-', ...
                'Color', line_colors(g, :), 'LineWidth', 2);
        end

        % --- Draw shaded traces: ctrl behind, TTL in front ---
        for g = [1, 3, 2, 4]
            s = grp_stats{g};
            if isempty(s.center), continue; end
            draw_shaded_polar_local(ax, s.theta, s.center, s.spread, ...
                line_colors(g, :), fill_colors(g, :), 0.30, 2.0);
        end

        % --- Axis formatting ---
        lim = rmax * 1.08;
        set(ax, 'XLim', [-lim, lim], 'YLim', [-lim, lim]);
        set(ax, 'XColor', 'none', 'YColor', 'none');
        title(ax, method_labels{mi}, 'FontSize', 14, 'FontWeight', 'bold');

        % --- Build annotation text ---
        ann = build_annotation_text(grp_short, grp_dsi, grp_ar, grp_n, pm);
        text(ax, -rmax * 1.05, -rmax * 1.08, ann, ...
            'FontName', 'FixedWidth', 'FontSize', 7, ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

        % --- Legend on first panel only ---
        if mi == 1
            leg_labs = cell(1, 4);
            for g = 1:4
                leg_labs{g} = sprintf('%s (n=%d)', grp_short{g}, grp_n(g));
            end
            lg = legend(ax, hleg, leg_labs, ...
                'Location', 'northwest', 'FontSize', 8, 'Box', 'on');
            try  % semi-transparent legend background
                lg.BoxFace.ColorType = 'truecoloralpha';
                lg.BoxFace.ColorData = uint8([255; 255; 255; 220]);
            catch
            end
        end
    end

    % Super-title
    sgtitle({'PD Alignment Methods — ctrl vs TTL (ON & OFF)', ...
             'DSI = (PD_{max} - ND_{max}) / PD_{max}    |    AR = PD_{max} / mean(ortho)'}, ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');

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

        text(ax, xl, yl, [card_labels{li}, deg_sym], 'FontSize', 8, ...
            'HorizontalAlignment', ha, 'VerticalAlignment', va, ...
            'Color', label_color, 'HandleVisibility', 'off');
    end

    % R-tick labels along 0-degree spoke
    for ri = 1:numel(r_ticks)
        text(ax, r_ticks(ri), rmax * 0.06, sprintf('%.1f', r_ticks(ri)), ...
            'FontSize', 6.5, 'Color', [0.6 0.6 0.6], ...
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


function ann = build_annotation_text(grp_short, grp_dsi, grp_ar, grp_n, pm)
% Build annotation text block with DSI, AR, and significance per group.

    L = {};

    % Header
    L{end+1} = sprintf('%-10s  %-11s  %-9s  %s', '', 'DSI', 'AR', 'n');

    % ON groups (1 = ON ctrl, 2 = ON TTL)
    for g = [1, 2]
        nd = numel(grp_dsi{g});
        na = numel(grp_ar{g});
        ds = metric_str(grp_dsi{g}, nd, pm, '%.2f');
        as = metric_str(grp_ar{g},  na, pm, '%.1f');
        L{end+1} = sprintf('%-10s  %-11s  %-9s  %d', ...
            grp_short{g}, ds, as, grp_n(g)); %#ok<AGROW>
    end

    % ON p-values
    [pd, pa] = ranksum_pair(grp_dsi{1}, grp_dsi{2}, grp_ar{1}, grp_ar{2});
    L{end+1} = sprintf('  ON p:     %-11s  %-9s', fmt_p(pd), fmt_p(pa));

    % OFF groups (3 = OFF ctrl, 4 = OFF TTL)
    for g = [3, 4]
        nd = numel(grp_dsi{g});
        na = numel(grp_ar{g});
        ds = metric_str(grp_dsi{g}, nd, pm, '%.2f');
        as = metric_str(grp_ar{g},  na, pm, '%.1f');
        L{end+1} = sprintf('%-10s  %-11s  %-9s  %d', ...
            grp_short{g}, ds, as, grp_n(g)); %#ok<AGROW>
    end

    % OFF p-values
    [pd, pa] = ranksum_pair(grp_dsi{3}, grp_dsi{4}, grp_ar{3}, grp_ar{4});
    L{end+1} = sprintf('  OFF p:    %-11s  %-9s', fmt_p(pd), fmt_p(pa));

    ann = strjoin(L, newline);

end


function s = metric_str(vals, n, pm, fmt)
% Format mean +/- SEM string.

    if n > 0
        m = mean(vals);
        se = std(vals) / sqrt(n);
        s = sprintf([fmt, '%s', fmt], m, pm, se);
    else
        s = '---';
    end

end


function [p_dsi, p_ar] = ranksum_pair(dsi1, dsi2, ar1, ar2)
% Wilcoxon rank-sum for DSI and AR between two groups.

    if numel(dsi1) >= 2 && numel(dsi2) >= 2
        p_dsi = ranksum(dsi1(:), dsi2(:));
    else
        p_dsi = NaN;
    end
    if numel(ar1) >= 2 && numel(ar2) >= 2
        p_ar = ranksum(ar1(:), ar2(:));
    else
        p_ar = NaN;
    end

end


function s = fmt_p(p)
% Format p-value with significance stars.

    if isnan(p)
        s = 'n/a';
    elseif p < 0.001
        s = '<.001***';
    elseif p < 0.01
        s = sprintf('.%03d**', round(p * 1000));
    elseif p < 0.05
        s = sprintf('.%02d*', round(p * 100));
    else
        s = sprintf('.%02d', round(p * 100));
    end

end

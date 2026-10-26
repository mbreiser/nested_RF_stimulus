function fig = plot_rapid_readout(P)
% PLOT_RAPID_READOUT  Placement figure of the rapid Protocol 2 read-out (MATLAB twin of the Python figure).
%
%   FIG = PLOT_RAPID_READOUT(P) draws, from the struct P that RAPID_READOUT
%   (or RAPID_READOUT_PY) assembles:
%     left   arena map (1-based columns and rows, row 1 at the top): the 88
%            bar centers, the PD-axis and orthogonal bar sets, the 4-px
%            lattice squares, the Protocol 4 core (blue) and line-scan (cyan)
%            centers, the PD-frame window s -12..+16 px / v +-8 px, the RF
%            center (M6 x M6), the PD arrow and the 11 Prot3b sites;
%     middle 80 ms bar-flash profile along the PD axis (s; + = trailing /
%            proximal) with the post-flash minimum, the window and the sites;
%     right  the same along the orthogonal axis (v).
%   The profile panels are drawn only when P.profiles is present (the Python
%   JSON does not carry the profiles).
%
%   Fields of P: cell, pd (the placement PD), pd_vector (28 deg/s vector), optional
%   pd_vector_56 and pd_snapped (for the title), center (1x2, col row, 1-based), center_xy
%   (1x2, 0-based x and y-up coordinates of the center), e_pd, e_or,
%   bar_cx, bar_cy (1x88, 1-based), g_pd, g_or (0-based group), lattice
%   (struct array: col_start, row_start), sites (struct array: region,
%   col_start, row_start, col_center, row_center, s_target, v_target, s_px,
%   v_px), p4_frames (struct array: item, col_center, row_center), window
%   (s_lead, s_trail, v_half), and optionally profiles.pd / profiles.or
%   (pos, peak, under, m6, optional set) and others (struct array: name, col, row).

    col = struct('center', [0.84 0.15 0.16], 'leading', [1.0 0.50 0.05], 'trailing', [0.58 0.40 0.74], ...
                 'orthogonal', [0.17 0.63 0.17], 'outside', [0.4 0.4 0.4]);
    has_prof = isfield(P, 'profiles') && ~isempty(P.profiles);
    if has_prof
        fig = figure('Color', 'w', 'Position', [40 80 1650 540], 'Name', ['rapid read-out ' P.cell]);
        tl = tiledlayout(fig, 1, 10, 'TileSpacing', 'compact', 'Padding', 'compact');
        a = nexttile(tl, [1 4]);
    else
        fig = figure('Color', 'w', 'Position', [40 80 760 640], 'Name', ['rapid read-out ' P.cell]);
        a = axes(fig);
    end
    hold(a, 'on');

    % lattice squares
    for k = 1:numel(P.lattice)
        rectangle(a, 'Position', [P.lattice(k).col_start - 0.5, P.lattice(k).row_start - 0.5, 4, 4], 'EdgeColor', [0.8 0.8 0.8], 'LineWidth', 0.4);
    end
    % Protocol 4 core and line centers
    if isfield(P, 'p4_frames') && ~isempty(P.p4_frames)
        it = {P.p4_frames.item};
        m = strcmp(it, 'core');
        plot(a, [P.p4_frames(m).col_center], [P.p4_frames(m).row_center], '.', 'Color', [0.12 0.47 0.71], 'MarkerSize', 5);
        m = strcmp(it, 'line');
        plot(a, [P.p4_frames(m).col_center], [P.p4_frames(m).row_center], '.', 'Color', [0.09 0.75 0.81], 'MarkerSize', 7);
    end
    h = gobjects(0); lab = {};
    h(end+1) = plot(a, P.bar_cx, P.bar_cy, '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7); lab{end+1} = '88 bar centers';
    ip = P.g_pd * 11 + (1:11); io = P.g_or * 11 + (1:11);
    h(end+1) = plot(a, P.bar_cx(ip), P.bar_cy(ip), 'o', 'Color', [0.12 0.47 0.71], 'MarkerSize', 5); lab{end+1} = 'perpendicular bar set (PD axis)';
    h(end+1) = plot(a, P.bar_cx(io), P.bar_cy(io), 'o', 'Color', [0.17 0.63 0.17], 'MarkerSize', 5); lab{end+1} = 'parallel bar set (orthogonal axis)';
    % window polygon
    W = P.window; corners = [W.s_lead -W.v_half; W.s_trail -W.v_half; W.s_trail W.v_half; W.s_lead W.v_half];
    poly = zeros(4, 2);
    for k = 1:4, poly(k, :) = to_arena(P, corners(k, 1), corners(k, 2)); end
    h(end+1) = plot(a, poly([1:4 1], 1), poly([1:4 1], 2), 'k-', 'LineWidth', 1.2); lab{end+1} = '4-px lattice window';
    cc = P.center(1); cr = P.center(2);
    h(end+1) = plot(a, cc, cr, 'k+', 'MarkerSize', 14, 'LineWidth', 2); lab{end+1} = 'RF center (M6 x M6)';
    tip = to_arena(P, 6, 0);
    quiver(a, cc, cr, tip(1) - cc, tip(2) - cr, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.8);
    text(a, tip(1), tip(2), ' PD (trailing)', 'FontSize', 8);
    % sites
    for k = 1:numel(P.sites)
        st = P.sites(k);
        rectangle(a, 'Position', [st.col_start - 0.5, st.row_start - 0.5, 4, 4], 'FaceColor', [col.(st.region) 0.35], 'EdgeColor', 'k', 'LineWidth', 0.8);
        switch st.region
            case 'outside', txt = 'out';
            case 'orthogonal', txt = sprintf('v%+d', st.v_target);
            case 'center', txt = sprintf('c%+d', st.v_target);
            otherwise, txt = sprintf('%+d', st.s_target);
        end
        text(a, st.col_center, st.row_center, txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 7);
    end
    if isfield(P, 'others')
        for k = 1:numel(P.others)
            if contains(P.others(k).name, 'protocol1'), mk = 'x'; else, mk = 'p'; end
            h(end+1) = plot(a, P.others(k).col, P.others(k).row, mk, 'Color', [0.84 0.15 0.16], 'MarkerSize', 10); lab{end+1} = strrep(P.others(k).name, '_', ' '); %#ok<AGROW>
        end
    end
    xlim(a, [cc - 26, cc + 26]); ylim(a, [cr - 26, cr + 26]); set(a, 'YDir', 'reverse'); axis(a, 'equal');
    xlabel(a, 'arena column (1-based)'); ylabel(a, 'arena row (1-based, top = 1)');
    if isfield(P, 'pd_snapped')
        ttl = sprintf('%s  placement PD %.1f%s (28 %s/s %.1f%s, 56 %s/s %.1f%s; snapped %g%s)  center col %.1f row %.1f  lattice %d', ...
            strrep(P.cell, '_', '\_'), P.pd, char(176), char(176), P.pd_vector, char(176), char(176), P.pd_vector_56, char(176), P.pd_snapped, char(176), cc, cr, numel(P.lattice));
    else
        ttl = sprintf('%s  PD %.1f%s (vector %.1f%s)  center col %.1f row %.1f  lattice %d', strrep(P.cell, '_', '\_'), P.pd, char(176), P.pd_vector, char(176), cc, cr, numel(P.lattice));
    end
    title(a, ttl, 'FontSize', 8.5);
    legend(a, h, lab, 'FontSize', 6.5, 'Location', 'northeast');

    if ~has_prof, return; end
    axn = {'pd', 'or'}; labs = {'PD axis (s; + = trailing / proximal)', 'orthogonal axis (v)'};
    for j = 1:2
        b = nexttile(tl, [1 3]); hold(b, 'on');
        pr = P.profiles.(axn{j}); [pos, o] = sort(pr.pos);
        hp = plot(b, pos, pr.peak(o), 'o-', 'Color', 'k');
        hu = plot(b, pos, pr.under(o), 's--', 'Color', [0.12 0.47 0.71], 'MarkerSize', 4);
        yline(b, 0, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); xline(b, 0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        if j == 1, span = [W.s_lead W.s_trail]; else, span = [-W.v_half W.v_half]; end
        yl = ylim(b);
        patch(b, span([1 2 2 1]), [yl(1) yl(1) yl(2) yl(2)], [0.92 0.92 0.92], 'EdgeColor', 'none', 'FaceAlpha', 1); uistack(findobj(b, 'Type', 'patch'), 'bottom');
        for k = 1:numel(P.sites)
            st = P.sites(k);
            if j == 1 && any(strcmp(st.region, {'center', 'leading', 'trailing'})) && st.v_target == 0
                xline(b, st.s_px, 'Color', [0.84 0.15 0.16], 'Alpha', 0.6, 'LineWidth', 0.8);
            elseif j == 2 && any(strcmp(st.region, {'orthogonal', 'outside'}))
                xline(b, st.v_px, 'Color', [0.84 0.15 0.16], 'Alpha', 0.6, 'LineWidth', 0.8);
            end
        end
        xlabel(b, [labs{j} ' (px)']); ylabel(b, 'mV');
        if isfield(pr, 'set'), ttl = sprintf('placement-frame set %d, M6 bar %d of 11; gray = window; red = Prot3b sites', pr.set, pr.m6);
        else, ttl = sprintf('M6 bar %d of 11; gray = window; red = Prot3b sites', pr.m6); end
        title(b, ttl, 'FontSize', 9);
        legend(b, [hp hu], {'80 ms peak', 'post-flash minimum (150-500 ms)'}, 'FontSize', 7, 'Location', 'best');
    end
end

function xy = to_arena(P, s, v)
    p = P.center_xy + s * P.e_pd + v * P.e_or;      % 0-based (x, y up)
    xy = [p(1) + 1, -p(2) + 1];
end

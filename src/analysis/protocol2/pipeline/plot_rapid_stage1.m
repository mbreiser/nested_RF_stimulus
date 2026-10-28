function fig = plot_rapid_stage1(P)
% PLOT_RAPID_STAGE1  Direction tuning of a Protocol 2 stage-1 session (the sweeps), as RAPID_READOUT computes it.
%
%   FIG = PLOT_RAPID_STAGE1(P) draws
%     left   the 28 and 56 deg/s sweep amplitudes per direction (polar, each normalized to its maximum, the
%            maximum in mV in the legend), the 28 and 56 deg/s vector PDs, the placement PD and the snapped PD
%     right  the repeat-mean Vm of the 28 deg/s sweep in every direction (baseline subtracted, motion shaded),
%            on one y scale, the panel nearest the placement PD outlined
%   P fields: cell, dirs (deg), amp28, amp56 (mV), pd_vector_28, pd_vector_56, pd_placement, pd_snapped,
%   vector_length_28, vector_length_56, pd_rule, pd_flag_gt5, usable_pd, dsi_28, dsi_56 (optional), traces28 (cell per direction:
%   struct t_s, vm_mV, motion_s), name (window name).
%   27 October 2026.

if ~isfield(P, 'name'), P.name = ['Protocol 2 stage 1 ' P.cell]; end
fig = findobj(groot, 'Type', 'figure', 'Name', P.name);
if isempty(fig), fig = figure('Name', P.name, 'NumberTitle', 'off', 'Color', 'w', 'Position', [80 80 1300 620]);
else, fig = fig(1); clf(fig); end
tl = tiledlayout(fig, 4, 8, 'TileSpacing', 'tight', 'Padding', 'compact');
flag = ''; if P.pd_flag_gt5, flag = '  [28/56 differ > 5 deg]'; end
title(tl, sprintf('%s  Protocol 2 stage 1: placement PD %.1f deg (rule %s), snapped %g deg%s', P.cell, P.pd_placement, P.pd_rule, P.pd_snapped, flag), ...
    'Interpreter', 'none', 'FontWeight', 'bold');
if ~isfield(P, 'dsi_28'), P.dsi_28 = NaN; end
if ~isfield(P, 'dsi_56'), P.dsi_56 = NaN; end
subtitle(tl, sprintf('28 deg/s vector %.1f deg (L %.2f, DSI %.2f)   56 deg/s vector %.1f deg (L %.2f, DSI %.2f)   usable PD: %s', ...
    P.pd_vector_28, P.vector_length_28, P.dsi_28, P.pd_vector_56, P.vector_length_56, P.dsi_56, tern(P.usable_pd, 'yes', 'NO (stage 2 plays all 88 bars)')), 'FontSize', 9);

% ---- polar tuning
pax = polaraxes(tl); pax.Layout.Tile = 1; pax.Layout.TileSpan = [4 4]; hold(pax, 'on');
th = deg2rad([P.dirs(:); P.dirs(1)]);
a28 = P.amp28(:); a56 = P.amp56(:);
m28 = max(a28, [], 'omitnan'); m56 = max(a56, [], 'omitnan');
polarplot(pax, th, [a28; a28(1)] / m28, '-o', 'Color', [0.1 0.3 0.75], 'LineWidth', 1.5, 'MarkerSize', 4, ...
    'DisplayName', sprintf('28 deg/s (max %.1f mV)', m28));
if any(isfinite(a56))
    polarplot(pax, th, [a56; a56(1)] / m56, '-s', 'Color', [0.85 0.4 0.1], 'LineWidth', 1.2, 'MarkerSize', 4, ...
        'DisplayName', sprintf('56 deg/s (max %.1f mV)', m56));
end
arrow = @(deg, len, c, w, s, nm) polarplot(pax, deg2rad([deg deg]), [0 len], s, 'Color', c, 'LineWidth', w, 'DisplayName', nm);
arrow(P.pd_vector_28, P.vector_length_28 / max(P.vector_length_28, P.vector_length_56) * 0.9, [0.1 0.3 0.75], 1.5, '-', '28 deg/s vector');
if isfinite(P.pd_vector_56)
    arrow(P.pd_vector_56, P.vector_length_56 / max(P.vector_length_28, P.vector_length_56) * 0.9, [0.85 0.4 0.1], 1.5, '-', '56 deg/s vector');
end
arrow(P.pd_placement, 1.1, [0 0 0], 3, '-', sprintf('placement PD %.1f', P.pd_placement));
arrow(P.pd_snapped, 1.1, [0.5 0.5 0.5], 1.5, '--', sprintf('snapped PD %g', P.pd_snapped));
pax.RLim = [0 1.15]; pax.RTickLabel = {}; pax.ThetaTick = 0:45:315;
legend(pax, 'Location', 'southoutside', 'NumColumns', 2, 'FontSize', 8);

% ---- 28 deg/s traces, 4 x 4 small multiples in direction order
n = numel(P.dirs); yl = [Inf -Inf];
for i = 1:n
    tr = P.traces28{i}; if isempty(tr), continue; end
    yl = [min(yl(1), min(tr.vm_mV)) max(yl(2), max(tr.vm_mV))];
end
if ~all(isfinite(yl)), yl = [-1 1]; end
yl = yl + [-0.05 0.05] * diff(yl);
[~, ipd] = min(abs(mod(P.dirs - P.pd_placement + 180, 360) - 180));
for i = 1:min(n, 16)
    r = ceil(i / 4); c = mod(i - 1, 4) + 1;
    ax = nexttile(tl, (r - 1) * 8 + 4 + c); hold(ax, 'on');
    tr = P.traces28{i};
    if ~isempty(tr)
        patch(ax, [0 tr.motion_s tr.motion_s 0], yl([1 1 2 2]), [0.93 0.93 0.93], 'EdgeColor', 'none');
        plot(ax, tr.t_s, tr.vm_mV, 'k', 'LineWidth', 0.6);
        xlim(ax, [tr.t_s(1) tr.t_s(end)]);
    end
    ylim(ax, yl); yline(ax, 0, ':', 'Color', [0.6 0.6 0.6]);
    text(ax, 0.03, 0.92, sprintf('%g%s', P.dirs(i), char(176)), 'Units', 'normalized', 'FontSize', 8, 'FontWeight', 'bold');
    ax.FontSize = 7; box(ax, 'on');
    if c > 1, ax.YTickLabel = {}; else, ylabel(ax, 'mV'); end
    if r < 4, ax.XTickLabel = {}; else, xlabel(ax, 's'); end
    if i == ipd, ax.LineWidth = 2; ax.XColor = [0 0 0]; ax.YColor = [0 0 0]; box(ax, 'on'); end
end
drawnow;
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end

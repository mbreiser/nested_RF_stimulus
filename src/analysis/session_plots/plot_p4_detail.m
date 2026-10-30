function fig = plot_p4_detail(A, varargin)
%PLOT_P4_DETAIL  Protocol 4 at the small scales and the hyperpolarization: one figure from ANALYZE_P4's output.
%
%   FIG = PLOT_P4_DETAIL(A) draws, from A = ANALYZE_P4(...) (A.squares, A.traces, A.t_ms):
%     row 1  position x time images of dV (mean over the kept passes): the 4-px lattice ON along the PD axis
%            (|v| <= 2 px), the 3-px core along the PD axis (|v| < 1), the 2-px line scans along s and along v;
%            the color scale is asymmetric (blue saturates at the panel's own minimum, red at its maximum), so
%            a -3 mV hyperpolarization is as visible as a +20 mV peak
%     row 2  maps with the same asymmetric scale: 4-px ON early (0-260 ms) with the late hyperpolarization
%            contoured (-1, -2, -4 mV, 150-500 ms), 4-px ON late, 3-px core early and late (150-400 ms:
%            the core and line squares repeat every 400 ms)
%     row 3  the three square sizes on common profiles along s (early; late) and along v (early), SE over
%            passes as bands; the 4-px lattice traces on the trailing side, zoomed
%   Positions are rounded to whole pixels and squares at the same rounded position averaged.
%   Options: 'Cell' (title), 'Png' ('' = none), 'Figure' (true: visible), 'Name' (window name).
%   29 October 2026.

opts = struct('Cell', '', 'Png', '', 'Figure', true, 'Name', 'Protocol 4 detail');
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
Q = A.squares; TR = A.traces; t = A.t_ms;
TR = movmean(TR, 9, 2, 'omitnan');                                        % 9-ms running mean, display only
on = Q.polarity == "ON";
sets = { ...                                                              % label, rows of Q, axis, color
    '4-px lattice', Q.item == "lattice" & on & abs(Q.v_px) <= 2, 's', [0.15 0.15 0.15]; ...
    '3-px core',    Q.item == "core" & abs(Q.v_px) < 1,           's', [0.85 0.4 0.1]; ...
    '2-px line',    Q.item == "line" & abs(Q.v_px) < 0.75,        's', [0.1 0.35 0.8]; ...
    '2-px line',    Q.item == "line" & abs(Q.s_px) < 0.75,        'v', [0.1 0.35 0.8]};
P = cellfun(@(m, ax) profile_of(Q, TR, m, ax), sets(:, 2), sets(:, 3));

fig = session_figure(opts.Name, opts.Figure, [70 50 1500 1000]);
tl = tiledlayout(fig, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s  Protocol 4 detail: small squares and the hyperpolarization (s + toward the PD)', opts.Cell), 'FontWeight', 'bold', 'Interpreter', 'none');
subtitle(tl, sprintf('passes kept %d/%d; color scales are asymmetric: blue saturates at each panel''s minimum, red at its maximum; white = 0 mV   (%s)', ...
    sum(A.passes.kept), height(A.passes), A.folder), 'FontSize', 8, 'Interpreter', 'none');

% ---- row 1: position x time
ttl = {'4-px ON lattice along the PD axis', '3-px ON core along the PD axis', '2-px ON line along s', '2-px ON line along v'};
for j = 1:4
    ax = nexttile(tl); hold(ax, 'on');
    Z = P(j).trace; per = 600 - 200 * (j > 1);                            % lattice 600 ms period, others 400 ms
    imagesc(ax, t, P(j).x, Z, 'AlphaData', ~isnan(Z));
    asym_scale(ax, Z(:, t >= 0 & t < per));
    xline(ax, [0 160], 'k:', 'HandleVisibility', 'off'); yline(ax, 0, '-', 'Color', [0.45 0.45 0.45]);
    patch(ax, [150 min(500, per) min(500, per) 150], P(j).x(end) + [0.3 0.3 0.6 0.6], [0.2 0.2 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    axis(ax, 'xy'); xlim(ax, [-100 per]); ylim(ax, [P(j).x(1) - 0.5, P(j).x(end) + 0.7]); box(ax, 'on');
    xlabel(ax, 'ms from flash onset'); ylabel(ax, sprintf('%s (px)', sets{j, 3}));
    title(ax, ttl{j}, 'FontSize', 9);
end

% ---- row 2: maps
L = Q(Q.item == "lattice" & on, :); C = Q(Q.item == "core", :);
specs = {L, 'early', '4-px ON, 0-260 ms (late contours)'; L, 'late', '4-px ON, 150-500 ms'; C, 'early', '3-px core, 0-260 ms'; C, 'late', '3-px core, 150-400 ms'};
for j = 1:4
    ax = nexttile(tl); hold(ax, 'on'); T = specs{j, 1}; val = T.(specs{j, 2});
    [vg, sg, Z] = grid_map(T.v_px, T.s_px, val);
    imagesc(ax, vg, sg, Z, 'AlphaData', ~isnan(Z)); asym_scale(ax, val);
    if j == 1 || j == 2
        [~, ~, Zl] = grid_map(T.v_px, T.s_px, T.late);
        [cm, hc] = contour(ax, vg, sg, Zl, [-4 -2 -1], 'LineColor', [0.05 0.15 0.6], 'LineWidth', 0.9);
        clabel(cm, hc, 'FontSize', 7, 'Color', [0.05 0.15 0.6]);
        if j == 1 && isfinite(A.fit.amp), th = linspace(0, 2 * pi, 100); plot(ax, A.fit.v0 + A.fit.sig_v * cos(th), A.fit.s0 + A.fit.sig_s * sin(th), 'k--', 'LineWidth', 0.8); end
        if A.trailing_extent_px > 0, yline(ax, A.trailing_extent_px, '--', 'extent', 'Color', [0.05 0.15 0.6], 'FontSize', 7, 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom'); end
    end
    scatter(ax, T.v_px, T.s_px, 4, [0.35 0.35 0.35], 'filled');
    plot(ax, 0, 0, 'k+', 'MarkerSize', 9, 'LineWidth', 1.3);
    axis(ax, 'xy'); axis(ax, 'equal'); axis(ax, 'tight'); xlabel(ax, 'v (px)'); ylabel(ax, 's (px)  PD \uparrow'); box(ax, 'on');
    title(ax, sprintf('%s  [%+.1f, %+.1f]', specs{j, 3}, min(val), max(val)), 'FontSize', 9);
end

% ---- row 3: profiles and trailing traces
pr = {[1 2 3], 'early', 's (px)', 'along s, 0-260 ms'; [1 2 3], 'late', 's (px)', 'along s, late (150-500 / 150-400 ms)'; [4], 'early', 'v (px)', 'along v, 0-260 ms'};
for j = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    if j == 3                                                             % lattice and core along v as well
        pv = [profile_of(Q, TR, Q.item == "lattice" & on & abs(Q.s_px) <= 2, 'v'), profile_of(Q, TR, Q.item == "core" & abs(Q.s_px) < 1, 'v'), P(4)];
        idx = 1:3; src = pv;
    else
        idx = pr{j, 1}; src = P(idx);
    end
    if strcmp(pr{j, 2}, 'late'), yregion = [-Inf 0]; else, yregion = []; end
    for i = 1:numel(src)
        k = tern(j == 3, i, idx(i)); c = sets{min(k, 3), 4};
        y = src(i).(pr{j, 2}); se = src(i).([pr{j, 2} '_se']); x = src(i).x;
        fill(ax, [x; flipud(x)], [y - se; flipud(y + se)], c, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax, x, y, '-o', 'Color', c, 'MarkerSize', 3, 'MarkerFaceColor', c, 'LineWidth', 1.1, 'DisplayName', sets{min(k, 3), 1});
    end
    yline(ax, 0, 'k:', 'HandleVisibility', 'off'); xline(ax, 0, 'k:', 'HandleVisibility', 'off');
    if ~isempty(yregion)
        yl = ylim(ax); patch(ax, [-30 30 30 -30], [yl(1) yl(1) 0 0], [0.88 0.92 1], 'EdgeColor', 'none', 'HandleVisibility', 'off');   % +-30: covers the +22 lattice
        to_back(ax); ylim(ax, yl);
        if A.trailing_extent_px > 0, xline(ax, A.trailing_extent_px, '--', 'Color', [0.05 0.15 0.6], 'HandleVisibility', 'off'); end
    end
    xlim(ax, [min(vertcat(src.x)) - 1, max(vertcat(src.x)) + 1]);
    xlabel(ax, pr{j, 3}); ylabel(ax, 'dV (mV), mean \pm SE over passes'); box(ax, 'on');
    if j == 1, legend(ax, 'Location', 'northwest', 'FontSize', 7); else, legend(ax, 'off'); end
    title(ax, pr{j, 4}, 'FontSize', 9);
end
ax = nexttile(tl); hold(ax, 'on');
p1 = P(1); want = [4 8 12 16 20]; cc = [0.85 0.2 0.15; 0.95 0.55 0.1; 0.6 0.4 0.8; 0.15 0.4 0.85; 0.05 0.55 0.45];   % s +20 only with the +22 lattice
for i = 1:numel(want)
    [d, r] = min(abs(p1.x - want(i))); if d > 1, continue; end
    plot(ax, t, p1.trace(r, :), 'Color', cc(i, :), 'LineWidth', 1.2, 'DisplayName', sprintf('s %+d (late %.1f mV)', p1.x(r), p1.late(r)));
end
r0 = find(p1.x == 0, 1); if ~isempty(r0), plot(ax, t, p1.trace(r0, :), 'Color', [0.6 0.6 0.6], 'DisplayName', 's 0 (clipped)'); end
patch(ax, [150 500 500 150], [-20 -20 20 20], [0.9 0.93 1], 'EdgeColor', 'none', 'HandleVisibility', 'off'); to_back(ax);
xline(ax, [0 160], 'k:', 'HandleVisibility', 'off'); yline(ax, 0, 'k:', 'HandleVisibility', 'off');
lo = min(p1.trace(p1.x >= 4, t >= 0 & t < 600), [], 'all'); hi = max(p1.trace(p1.x >= 4, t >= 0 & t < 600), [], 'all');
ylim(ax, [min(-3, 1.2 * lo), max(3, 1.2 * hi)]); xlim(ax, [-100 600]);
xlabel(ax, 'ms from flash onset'); ylabel(ax, 'dV (mV)'); box(ax, 'on'); legend(ax, 'Location', 'best', 'FontSize', 7);
title(ax, '4-px ON on the trailing side (|v| \leq 2), zoomed', 'FontSize', 9);
drawnow;
if ~isempty(opts.Png), exportgraphics(fig, opts.Png, 'Resolution', 120); end
end

% ------------------------------------------------------------------------------------------------
function p = profile_of(Q, TR, m, ax)
% squares in M, by whole-pixel position along AX: mean early/late, pooled SE, mean trace
x = round(Q.([ax '_px'])(m)); [g, xu] = findgroups(x);
e = Q.early(m); l = Q.late(m); es = Q.early_se(m); ls = Q.late_se(m); T = TR(m, :);
p.x = xu(:);
p.early = splitapply(@(v) mean(v, 'omitnan'), e, g); p.late = splitapply(@(v) mean(v, 'omitnan'), l, g);
p.early_se = splitapply(@(v) sqrt(mean(v.^2, 'omitnan') / numel(v)), es, g);
p.late_se = splitapply(@(v) sqrt(mean(v.^2, 'omitnan') / numel(v)), ls, g);
p.trace = cell2mat(arrayfun(@(k) mean(T(g == k, :), 1, 'omitnan'), (1:numel(xu)).', 'UniformOutput', false));
end

function [vg, sg, Z] = grid_map(v, s, val)
ok = isfinite(val); v = v(ok); s = s(ok); val = val(ok);
vg = min(v):0.25:max(v); sg = min(s):0.25:max(s); [V, S] = meshgrid(vg, sg);
if numel(val) < 4, Z = nan(size(V)); return; end
Z = scatteredInterpolant(v, s, val, 'natural', 'none'); Z = Z(V, S);
end

function asym_scale(ax, vals)
% color limits [min max] of VALS with white at 0: each side of zero uses its full blue / red range
vals = vals(isfinite(vals)); if isempty(vals), return; end
lo = min(min(vals), -0.5); hi = max(max(vals), 0.5);
n = 256; nb = max(2, round(n * -lo / (hi - lo))); nr = max(2, n - nb);
b = linspace(0, 1, nb).'; r = linspace(0, 1, nr).';
blue = [0.05 0.2 0.7] + b * ([1 1 1] - [0.05 0.2 0.7]);
red = [1 1 1] + r * ([0.7 0.05 0.05] - [1 1 1]);
colormap(ax, [blue; red(2:end, :)]); clim(ax, [lo hi]);
cb = colorbar(ax); cb.Label.String = 'mV'; cb.Ticks = unique([round(lo, 1), lo / 2, 0, hi / 2, round(hi, 1)]);
cb.TickLabels = compose('%+.1f', cb.Ticks);
end

function to_back(ax)
% the object drawn last goes behind everything else
ax.Children = ax.Children([2:end, 1]);
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end

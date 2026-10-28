function A = analyze_p4(folder, varargin)
%ANALYZE_P4  Protocol 4 read-outs and one compact figure: the combined RF map at three scales.
%
%   A = ANALYZE_P4(FOLDER) analyses a recorded Protocol 4 package (the sequencer's p4/ folder) as the
%   hand-over specifies (JY_PROTOCOL_HANDOVER_20261026.md, sections 5 and 7), with the windows of
%   datasets/20261022_jy_rnai_rf/multiscale_rf.py (square_maps):
%     flashes     onsets from the frame log (LOAD_STEP_EVENTS); dV = Vm minus the flash's own pre-flash
%                 mean (-100..0 ms), at 1 kHz
%     pass QC     per pass the median pre-flash Vm; a pass is kept when within 5 mV of pass 1
%     maps        per square (pattern frame), the mean over the kept passes of the mean dV in the early
%                 window 0-260 ms and the late window 150-500 ms (lattice; 440 ms gray keeps it intact),
%                 in the PD frame (s + toward the PD, v 90 deg clockwise) of the read-out center
%     fit         2-D Gaussian (amplitude, s0, v0, sigma_s, sigma_v; three starts sigma 3/5/8 px,
%                 fminsearch on the squared error) on the early lattice ON map
%     trailing    the trailing extent of the hyperpolarization: along |v| <= 2 px, in 0.5-px bins of s > 0,
%                 the far end of the first run of bins where the late lattice ON map is below -1 mV
%     core, lines the 3-px ON core (early window) and the 2-px ON line scans (s line at v = 0, v line at
%                 s = 0), profiles with SE over passes, the column period (3.6-4 px) marked
%     OFF series  the four 3-px OFF squares (s = 0, -4, -8, +4): mean dV traces
%   Options: 'Cell' (title), 'Png' ('' = none), 'Figure' (true: visible), 'Name' (window name).
%   A holds the per-square table (A.squares), A.passes, A.fit, A.trailing_extent_px, A.summary.
%   28 October 2026.

opts = struct('Cell', '', 'Png', '', 'Figure', true, 'Name', 'Protocol 4');
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
t_all = tic;
D = load_step_events(folder); E = D.E;
[X, t] = step_snips(D.vm, E.k_on, [-100 600]);
E.v_pre = mean(X(:, t < 0), 2); Xb = X - E.v_pre;
per = E.flash_duration_s + E.gray_after_s;                                 % s
for i = 1:height(E), Xb(i, t >= per(i) * 1000) = NaN; end                  % never past the next flash
E.early = mean(Xb(:, t >= 0 & t < 260), 2, 'omitnan');
E.late = mean(Xb(:, t >= 150 & t < 500), 2, 'omitnan');

% ---- pass QC
pl = unique(E.pass).'; base = arrayfun(@(p) median(E.v_pre(E.pass == p), 'omitnan'), pl);
Pt = table(pl(:), base(:), abs(base(:) - base(1)) <= 5, 'VariableNames', {'pass', 'base_mV', 'kept'});
keep = ismember(E.pass, Pt.pass(Pt.kept));

% ---- per square
[g, fr] = findgroups(E.pattern_frame);
Q = table(fr, 'VariableNames', {'frame'});
first = splitapply(@(i) i(1), (1:height(E)).', g);
for c = {'item', 'polarity', 'size_px', 's_px', 'v_px'}, Q.(c{1}) = E.(c{1})(first); end
Q.n = splitapply(@(k) sum(k), keep, g);
Q.early = splitapply(@(v, k) mean(v(k), 'omitnan'), E.early, keep, g);
Q.late = splitapply(@(v, k) mean(v(k), 'omitnan'), E.late, keep, g);
Q.early_se = splitapply(@(v, k) std(v(k), 'omitnan') / sqrt(max(sum(k), 1)), E.early, keep, g);

% ---- Gaussian fit on the early lattice ON map
L1 = Q(Q.item == "lattice" & Q.polarity == "ON", :);
F = fit_gauss2d(L1.s_px, L1.v_px, L1.early);
% ---- trailing extent of the hyperpolarization (late lattice ON along v = 0)
ax0 = L1(abs(L1.v_px) <= 2, :); trail = 0;
sx = unique(round(ax0.s_px(ax0.s_px > 0) * 2) / 2).';
for s = sx
    v = mean(ax0.late(abs(ax0.s_px - s) < 0.75), 'omitnan');
    if v < -1, trail = s; elseif trail > 0, break; end
end

A = struct('folder', folder, 'events', E, 'squares', Q, 'passes', Pt, 'fit', F, 'trailing_extent_px', trail, 't_ms', t);
A.summary = struct('passes_kept', sum(Pt.kept), 'passes', height(Pt), 'flashes_missing', D.n_missing, ...
    'lattice_ON_peak_mV', max(L1.early), 'gauss_amp_mV', F.amp, 'gauss_s0_px', F.s0, 'gauss_v0_px', F.v0, ...
    'gauss_sigma_s_px', F.sig_s, 'gauss_sigma_v_px', F.sig_v, 'trailing_extent_px', trail, ...
    'lattice_ON_late_min_mV', min(L1.late), 'vm_first_pass_mV', base(1), 'vm_last_pass_mV', base(end));

% =========================================================================================== figure
fig = session_figure(opts.Name, opts.Figure, [60 60 1500 900]);
tl = tiledlayout(fig, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s  Protocol 4: combined RF map (PD frame: s + toward the PD, up)', opts.Cell), 'FontWeight', 'bold', 'Interpreter', 'none');
subtitle(tl, sprintf('passes kept %d/%d   lattice ON peak %.1f mV   Gaussian center s %.1f v %.1f px, sigma %.1f x %.1f px   trailing hyperpolarization to s = %+.1f px   (%s)', ...
    sum(Pt.kept), height(Pt), max(L1.early), F.s0, F.v0, F.sig_s, F.sig_v, trail, folder), 'FontSize', 8, 'Interpreter', 'none');
specs = {"lattice", "ON", "early", '4-px ON, 0-260 ms'; "lattice", "ON", "late", '4-px ON, 150-500 ms'; ...
         "lattice", "OFF", "early", '4-px OFF, 0-260 ms'; "lattice", "OFF", "late", '4-px OFF, 150-500 ms'; "core", "ON", "early", '3-px ON core, 0-260 ms'};
for j = 1:size(specs, 1)
    ax = nexttile(tl);
    m = Q.item == specs{j, 1} & Q.polarity == specs{j, 2};
    map_panel(ax, Q.s_px(m), Q.v_px(m), Q.(specs{j, 3})(m), specs{j, 4});
    if j == 1, gauss_contour(ax, F); end
    if j == 2 && trail > 0, plot(ax, 0, trail, 'k^', 'MarkerFaceColor', 'y', 'MarkerSize', 7); end
end
% line scans
lines = Q(Q.item == "line", :); sl = lines(abs(lines.v_px) < 0.75, :); vl = lines(abs(lines.s_px) < 0.75, :);
core = Q(Q.item == "core", :);
for j = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    if j == 1, T = sortrows(sl, 's_px'); x = T.s_px; C = sortrows(core(abs(core.v_px) < 0.75, :), 's_px'); xc = C.s_px; xl = 's (px)';
    else, T = sortrows(vl, 'v_px'); x = T.v_px; C = sortrows(core(abs(core.s_px) < 0.75, :), 'v_px'); xc = C.v_px; xl = 'v (px)'; end
    errorbar(ax, x, T.early, T.early_se, '-o', 'Color', [0.1 0.3 0.75], 'MarkerSize', 3, 'MarkerFaceColor', [0.1 0.3 0.75], 'DisplayName', '2-px line');
    if ~isempty(C), plot(ax, xc, C.early, '-s', 'Color', [0.85 0.4 0.1], 'MarkerSize', 3, 'DisplayName', '3-px core'); end
    for p = -12:3.8:12, xline(ax, p, ':', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off'); end
    yline(ax, 0, ':', 'HandleVisibility', 'off'); xlabel(ax, xl); ylabel(ax, 'dV 0-260 ms (mV)'); box(ax, 'on');
    title(ax, sprintf('%s line scan (column period 3.8 px dotted)', tern(j == 1, 's', 'v')), 'FontSize', 9);
    legend(ax, 'Location', 'best', 'FontSize', 7);
end
% traces at chosen lattice ON positions
ax = nexttile(tl); hold(ax, 'on');
targets = [0 0; 4 0; 8 0; 12 0; -4 0]; cc = lines_colors(size(targets, 1));
for j = 1:size(targets, 1)
    [d, i] = min(hypot(L1.s_px - targets(j, 1), L1.v_px - targets(j, 2))); if d > 1.5, continue; end
    m = E.pattern_frame == L1.frame(i) & keep;
    plot(ax, t, mean(Xb(m, :), 1, 'omitnan'), 'Color', cc(j, :), 'DisplayName', sprintf('s %+.0f v %+.0f', L1.s_px(i), L1.v_px(i)));
end
xline(ax, [0 160], ':', 'HandleVisibility', 'off'); yline(ax, 0, ':', 'HandleVisibility', 'off'); xlim(ax, [-100 600]);
xlabel(ax, 'ms'); ylabel(ax, 'dV (mV)'); legend(ax, 'Location', 'best', 'FontSize', 7); box(ax, 'on'); title(ax, '4-px ON along the PD axis', 'FontSize', 9);
% OFF series
ax = nexttile(tl); hold(ax, 'on'); o3 = Q(Q.item == "off3", :); o3 = sortrows(o3, 's_px'); cc = lines_colors(height(o3));
for j = 1:height(o3)
    m = E.pattern_frame == o3.frame(j) & keep; mu = mean(Xb(m, :), 1, 'omitnan');
    plot(ax, t, mu, 'Color', cc(j, :), 'DisplayName', sprintf('s %+.0f (%.1f mV)', o3.s_px(j), o3.early(j)));
end
xline(ax, [0 160], ':', 'HandleVisibility', 'off'); yline(ax, 0, ':', 'HandleVisibility', 'off'); xlim(ax, [-100 400]);
xlabel(ax, 'ms'); ylabel(ax, 'dV (mV)'); legend(ax, 'Location', 'best', 'FontSize', 7); box(ax, 'on'); title(ax, '3-px OFF series (0-260 ms mean)', 'FontSize', 9);
% pass QC and Vm record
ax = nexttile(tl, [1 2]); hold(ax, 'on');
k1 = 1:100:numel(D.vm); plot(ax, (k1 - 1) / D.fs / 60, D.vm(k1), 'Color', [0.35 0.35 0.35], 'LineWidth', 0.4);
pk = arrayfun(@(p) median(E.k_on(E.pass == p), 'omitnan'), Pt.pass);
plot(ax, (pk - 1) / D.fs / 60, Pt.base_mV, 'o-', 'Color', [0.9 0.6 0], 'MarkerFaceColor', [1 0.8 0], 'LineWidth', 1.2);
for j = find(~Pt.kept).', plot(ax, (pk(j) - 1) / D.fs / 60, Pt.base_mV(j), 'rx', 'MarkerSize', 10, 'LineWidth', 1.5); end
xlabel(ax, 'min'); ylabel(ax, 'Vm (mV)'); box(ax, 'on');
title(ax, sprintf('Vm; per-pass median baseline (x: > 5 mV from pass 1); %d / %d passes kept', sum(Pt.kept), height(Pt)), 'FontSize', 9);
ax = nexttile(tl); axis(ax, 'off');
cnt = @(it, pol) sum(Q.item == it & Q.polarity == pol);
txt = {'\bfread-outs\rm', ...
    sprintf('lattice ON peak  %6.1f mV (0-260 ms)', max(L1.early)), ...
    sprintf('Gaussian amp     %6.1f mV, R^2 %.2f', F.amp, F.r2), ...
    sprintf('  center s %+5.1f  v %+5.1f px', F.s0, F.v0), ...
    sprintf('  sigma  s %5.1f  v %5.1f px', F.sig_s, F.sig_v), ...
    sprintf('late ON minimum  %6.1f mV (150-500 ms)', min(L1.late)), ...
    sprintf('trailing hyperpol. to s = %+.1f px', trail), '', '\bfsquares\rm', ...
    sprintf('lattice ON %d, OFF %d; core %d; line %d; OFF3 %d', cnt("lattice", "ON"), cnt("lattice", "OFF"), sum(Q.item == "core"), sum(Q.item == "line"), sum(Q.item == "off3")), ...
    sprintf('flashes %d, missing in the frame log %d', height(E), D.n_missing), ...
    sprintf('Vm pass 1 %.1f mV, last pass %.1f mV', base(1), base(end))};
text(ax, 0, 1, txt, 'VerticalAlignment', 'top', 'FontName', 'Consolas', 'FontSize', 8, 'Interpreter', 'tex');
drawnow;
if ~isempty(opts.Png), exportgraphics(fig, opts.Png, 'Resolution', 120); end
A.figure = fig; A.seconds = toc(t_all);
end

% ------------------------------------------------------------------------------------------------
function map_panel(ax, s, v, val, ttl)
% interpolated signed map in the PD frame (s up = toward the PD), symmetric diverging colors, squares marked
hold(ax, 'on');
ok = isfinite(val); s = s(ok); v = v(ok); val = val(ok);
if numel(val) >= 4
    [Vg, Sg] = meshgrid(min(v):0.25:max(v), min(s):0.25:max(s));
    Fi = scatteredInterpolant(v, s, val, 'natural', 'none'); Z = Fi(Vg, Sg);
    imagesc(ax, Vg(1, :), Sg(:, 1), Z, 'AlphaData', ~isnan(Z));
end
scatter(ax, v, s, 6, [0.3 0.3 0.3], 'filled');
plot(ax, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
m = max(abs(val)); if ~(m > 0), m = 1; end
clim(ax, [-m m]); colormap(ax, diverging_map()); cb = colorbar(ax); cb.Label.String = 'mV';
axis(ax, 'xy'); axis(ax, 'equal'); axis(ax, 'tight'); xlabel(ax, 'v (px)'); ylabel(ax, 's (px)  PD \uparrow');
title(ax, sprintf('%s  [%+.1f, %+.1f]', ttl, min(val), max(val)), 'FontSize', 9);
end

function gauss_contour(ax, F)
if ~isfinite(F.amp), return; end
th = linspace(0, 2 * pi, 100);
plot(ax, F.v0 + F.sig_v * cos(th), F.s0 + F.sig_s * sin(th), 'k-', 'LineWidth', 1);
plot(ax, F.v0, F.s0, 'kx', 'MarkerSize', 8, 'LineWidth', 1.5);
end

function F = fit_gauss2d(s, v, z)
% amp * exp(-((s-s0)^2/(2 sig_s^2) + (v-v0)^2/(2 sig_v^2))): three starts, bounded by transforms
F = struct('amp', NaN, 's0', NaN, 'v0', NaN, 'sig_s', NaN, 'sig_v', NaN, 'r2', NaN);
ok = isfinite(z); s = s(ok); v = v(ok); z = z(ok); if numel(z) < 6, return; end
[zm, i] = max(z); best = Inf;
model = @(p) p(1) * exp(-((s - p(2)).^2 / (2 * p(4)^2) + (v - p(3)).^2 / (2 * p(5)^2)));
unpack = @(q) [q(1), q(2), q(3), 0.5 + exp(q(4)), 0.5 + exp(q(5))];
for sg = [3 5 8]
    q0 = [zm, s(i), v(i), log(sg - 0.5), log(sg - 0.5)];
    [q, fv] = fminsearch(@(q) sum((model(unpack(q)) - z).^2), q0, optimset('MaxFunEvals', 4000, 'MaxIter', 4000, 'Display', 'off'));
    if fv < best, best = fv; p = unpack(q); end
end
F = struct('amp', p(1), 's0', p(2), 'v0', p(3), 'sig_s', p(4), 'sig_v', p(5), 'r2', 1 - best / sum((z - mean(z)).^2));
end

function c = diverging_map()
n = 128; a = linspace(0, 1, n).';
c = [[a, a, ones(n, 1)]; [ones(n, 1), flipud(a), flipud(a)]];      % blue - white - red
c = c([1:n, n+2:2*n], :);
end

function c = lines_colors(n)
c = [0 0 0; 0.85 0.2 0.15; 0.95 0.55 0.1; 0.6 0.4 0.8; 0.15 0.4 0.85; 0.1 0.6 0.3];
c = c(1 + mod(0:n - 1, size(c, 1)), :);
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end

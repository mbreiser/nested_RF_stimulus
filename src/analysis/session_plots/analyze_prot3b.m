function A = analyze_prot3b(folder, varargin)
%ANALYZE_PROT3B  Prot3b read-outs and one compact figure: apparent input resistance during site flashes.
%
%   A = ANALYZE_PROT3B(FOLDER) analyses a recorded Prot3b package (the sequencer's p3b/ folder) as the
%   hand-over specifies (JY_PROTOCOL_HANDOVER_20261026.md, sections 4 and 7), with the arithmetic of
%   datasets/20261025_jy_prot3_dr_vs_dv/dr_time.py adapted to the two-current pairs of Prot3b:
%     flashes     onsets from the frame log (LOAD_STEP_EVENTS); traces at 1 kHz from -100 to 500 ms
%     pairs       blocks 1-20: consecutive +0.5 / -0.5 pA blocks (ABBA: 1-2, 4-5, 6-7, ...), 8 pairs in 4
%                 quartets; QC per quartet on the block baselines (median pre-flash Vm, -100..0 ms): the two
%                 blocks of the same current within 2 mV, every block within 5 mV of the first block of its
%                 current (the hand-over's rule, read per current: opposite currents differ by R x 1 pA)
%     R(t)        per pair, site and polarity: (V+(t) - V-(t)) / (I+ - I-), nominal currents (mV/pA = GOhm);
%                 R_gray the same on the site's pre-flash window -100..0 ms; dR(t) = R(t) - R_gray, the
%                 mean over the pairs kept, SE across pairs; regions pool their sites
%     read-outs   ON trailing dR 300-500 ms, OFF center dR 100-200 ms, ON center dR 150-250 ms, and
%                 dR(t) in 50-ms bins to 500 ms per region and polarity
%     linearity   the 0 pA blocks: V0 - (V+ + V-)/2 per region, 50-300 ms (the mid-point residual)
%     Mi9 series  blocks 21-26 (3-px OFF squares at s = 0, -4, -8): baseline-subtracted response per
%                 current (-1, 0, +1 pA), its mean over 100-200 ms and 200-500 ms, and dR_OFF = slope of the
%                 response-window Vm on current minus the slope of the pre-flash Vm (100-200 ms)
%   The figure: the Vm record with the block currents and pair QC; the 0 pA responses per region (ON, OFF);
%   dR(t) per region (ON, OFF); the site map; the Mi9 traces and I-V; the a-priori read-outs.
%
%   Options: 'Cell' (title), 'Png' ('' = none), 'Figure' (true: visible), 'Name' (window name).
%   A holds the tables (A.events, A.blocks, A.pairs, A.dR_bins, A.readouts, A.mi9) and A.summary.
%   28 October 2026.

opts = struct('Cell', '', 'Png', '', 'Figure', true, 'Name', 'Prot3b');
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
t_all = tic;
D = load_step_events(folder); E = D.E; B = D.B;
pre = [-100 0];
[X, t] = step_snips(D.vm, E.k_on, [-100 500]);
E.v_pre = mean(X(:, t >= pre(1) & t < pre(2)), 2);
Xb = X - E.v_pre;                                                         % baseline-subtracted

% ---- blocks: baseline, measured current (settle window), pairs
nb = height(B); blk = table(B.block_index, B.phase, B.requested_current_pA, 'VariableNames', {'block', 'phase', 'I_pA'});
blk.base_mV = arrayfun(@(b) median(E.v_pre(E.block_index == b), 'omitnan'), blk.block);
blk.i_mon_pA = nan(nb, 1); blk.t0 = nan(nb, 1); blk.t1 = nan(nb, 1);
for b = 1:nb
    r = find(startsWith(string({D.rows.label}), sprintf('block %02d', blk.block(b))), 1);
    if isempty(r), continue; end
    k0 = D.rows(r).t0; blk.t0(b) = k0; blk.t1(b) = k0 + D.rows(r).n;
    if ~isempty(D.im), blk.i_mon_pA(b) = median(D.im(k0 + 5000:min(k0 + 30000, end))); end
end
rin = find(blk.phase == "rin");
pairs = zeros(0, 2); i = 1;
while i < numel(rin)
    a = rin(i); b = rin(i + 1);
    if b == a + 1 && blk.I_pA(a) ~= 0 && blk.I_pA(b) == -blk.I_pA(a), pairs(end+1, :) = [a b]; i = i + 2; else, i = i + 1; end %#ok<AGROW>
end
% QC per ABBA quartet (pairs 2q-1, 2q: + - | - +): the two blocks of the same current within 2 mV of each
% other, and every block within 5 mV of the first block of its current. (Blocks of opposite current differ by
% R x 1 pA, about 13 mV at 13 GOhm, by design, so "the two blocks of an ABBA pair" can only be the A-A and
% B-B blocks; and "the run's first block" is read per current.)
P = table(pairs(:, 1), pairs(:, 2), 'VariableNames', {'block_a', 'block_b'});
first_of = @(c) blk.base_mV(rin(find(blk.I_pA(rin) == c, 1)));
P.quartet = ceil((1:height(P)).' / 2);
P.d_same_mV = nan(height(P), 1); P.max_from_first_mV = nan(height(P), 1);
for q = unique(P.quartet).'
    j = find(P.quartet == q); bl = unique(pairs(j, :));
    d = 0; for c = unique(blk.I_pA(bl)).', bc = bl(blk.I_pA(bl) == c); if numel(bc) == 2, d = max(d, abs(diff(blk.base_mV(bc)))); end, end
    f = max(arrayfun(@(b) abs(blk.base_mV(b) - first_of(blk.I_pA(b))), bl));
    P.d_same_mV(j) = d; P.max_from_first_mV(j) = f;
end
P.kept = P.d_same_mV <= 2 & P.max_from_first_mV <= 5;

% ---- R(t) per pair x site x polarity (Rin blocks)
regions = ["center", "leading", "trailing", "orthogonal", "outside"];
isr = E.phase == "rin";
sites = unique(E(isr, {'s_px', 'v_px', 'region'}), 'rows');
pols = ["ON", "OFF"];
np = height(P); ns = height(sites);
R = nan(np, ns, 2, numel(t)); Rg = nan(np, ns, 2);
for p = 1:np
    ba = P.block_a(p); bb = P.block_b(p); dI = blk.I_pA(ba) - blk.I_pA(bb);
    for s = 1:ns
        for q = 1:2
            ma = isr & E.block_index == blk.block(ba) & E.s_px == sites.s_px(s) & E.v_px == sites.v_px(s) & E.polarity == pols(q);
            mb = isr & E.block_index == blk.block(bb) & E.s_px == sites.s_px(s) & E.v_px == sites.v_px(s) & E.polarity == pols(q);
            if ~any(ma) || ~any(mb), continue; end
            R(p, s, q, :) = (mean(X(ma, :), 1, 'omitnan') - mean(X(mb, :), 1, 'omitnan')) / dI;
            Rg(p, s, q) = (mean(E.v_pre(ma), 'omitnan') - mean(E.v_pre(mb), 'omitnan')) / dI;
        end
    end
end
dR = R - Rg;                                                              % pair x site x pol x time
kept = P.kept; qc_all = ~any(kept);
if qc_all, kept = true(np, 1); end                                        % nothing kept: show all, flagged
reg_of = sites.region;
bins = 0:50:450;
dRb = table();
for q = 1:2
    for g = regions
        sm = reg_of == g; if ~any(sm), continue; end
        cv = squeeze(mean(dR(kept, sm, q, :), 2, 'omitnan'));             % pairs x time
        if isvector(cv), cv = cv(:).'; end
        for b0 = bins
            w = t >= b0 & t < b0 + 50; v = mean(cv(:, w), 2, 'omitnan');
            dRb = [dRb; table(pols(q), g, b0, mean(v, 'omitnan'), std(v, 'omitnan') / sqrt(sum(~isnan(v))), sum(~isnan(v)), ...
                'VariableNames', {'polarity', 'region', 'bin_start_ms', 'dR_GOhm', 'se', 'n_pairs'})]; %#ok<AGROW>
        end
    end
end
win = {"ON", "trailing", [300 500]; "OFF", "center", [100 200]; "ON", "center", [150 250]};
RO = table();
for j = 1:size(win, 1)
    q = find(pols == win{j, 1}); sm = reg_of == win{j, 2}; w = t >= win{j, 3}(1) & t < win{j, 3}(2);
    v = squeeze(mean(mean(dR(kept, sm, q, w), 4, 'omitnan'), 2, 'omitnan'));
    RO = [RO; table(string(win{j, 1}), string(win{j, 2}), win{j, 3}, mean(v, 'omitnan'), std(v, 'omitnan') / sqrt(sum(~isnan(v))), sum(~isnan(v)), ...
        'VariableNames', {'polarity', 'region', 'window_ms', 'dR_GOhm', 'se', 'n_pairs'})]; %#ok<AGROW>
end
rg_all = Rg(kept, :, :); R_gray = mean(rg_all(:), 'omitnan');

% ---- linearity: 0 pA mid-point residual per region (50-300 ms, absolute Vm)
w = t >= 50 & t < 300; lin = table();
for q = 1:2
    for g = regions
        m0 = isr & E.region == g & E.polarity == pols(q);
        v0 = mean(mean(X(m0 & E.requested_current_pA == 0, w), 2, 'omitnan'), 'omitnan');
        vp = mean(mean(X(m0 & E.requested_current_pA > 0, w), 2, 'omitnan'), 'omitnan');
        vn = mean(mean(X(m0 & E.requested_current_pA < 0, w), 2, 'omitnan'), 'omitnan');
        lin = [lin; table(pols(q), g, v0 - (vp + vn) / 2, 'VariableNames', {'polarity', 'region', 'residual_mV'})]; %#ok<AGROW>
    end
end

% ---- Mi9 holding series
im9 = E.phase ~= "rin";
[Xm, tm] = step_snips(D.vm, E.k_on, [-100 1000]);
vpre_m = mean(Xm(:, tm < 0), 2); Xmb = Xm - vpre_m;
sq = unique(round(E.s_px(im9))).'; cur = unique(E.requested_current_pA(im9)).';
mi9 = table();
for s = sq
    for c = cur
        m = im9 & round(E.s_px) == s & E.requested_current_pA == c;
        mi9 = [mi9; table(s, c, sum(m), mean(vpre_m(m), 'omitnan'), mean(mean(Xm(m, tm >= 100 & tm < 200), 2), 'omitnan'), ...
            mean(mean(Xmb(m, tm >= 100 & tm < 200), 2), 'omitnan'), mean(mean(Xmb(m, tm >= 200 & tm < 500), 2), 'omitnan'), ...
            'VariableNames', {'s_px', 'I_pA', 'n', 'v_pre_mV', 'v_resp_100_200_mV', 'resp_100_200_mV', 'resp_200_500_mV'})]; %#ok<AGROW>
    end
end
dR_off = nan(size(sq));
for j = 1:numel(sq)
    m = mi9.s_px == sq(j); if sum(m) < 2, continue; end
    pr = polyfit(mi9.I_pA(m), mi9.v_resp_100_200_mV(m), 1); pg = polyfit(mi9.I_pA(m), mi9.v_pre_mV(m), 1);
    dR_off(j) = pr(1) - pg(1);
end

A = struct('folder', folder, 'events', E, 'blocks', blk, 'pairs', P, 'dR_bins', dRb, 'readouts', RO, 'linearity', lin, 'mi9', mi9, ...
    'mi9_s_px', sq, 'mi9_dR_off_GOhm', dR_off, 'R_gray_GOhm', R_gray, 't_ms', t);
A.summary = struct('pairs_kept', sum(P.kept), 'pairs', np, 'qc_none_kept_all_shown', qc_all, 'R_gray_GOhm', R_gray, 'flashes_missing', D.n_missing, ...
    'ON_trailing_dR_300_500', RO.dR_GOhm(1), 'OFF_center_dR_100_200', RO.dR_GOhm(2), 'ON_center_dR_150_250', RO.dR_GOhm(3), ...
    'mi9_dR_off_s0', dR_off(sq == 0), 'monitor_slope', NaN);
ok = ~isnan(blk.i_mon_pA) & blk.phase == "rin";
if sum(ok) > 2, pf = polyfit(blk.I_pA(ok), blk.i_mon_pA(ok), 1); A.summary.monitor_slope = pf(1); end

% =========================================================================================== figure
cols = struct('center', [0 0 0], 'leading', [0.15 0.4 0.85], 'trailing', [0.85 0.2 0.15], 'orthogonal', [0.1 0.6 0.3], 'outside', [0.6 0.6 0.6]);
fig = session_figure(opts.Name, opts.Figure, [60 60 1500 900]);
tl = tiledlayout(fig, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s  Prot3b: input resistance during site flashes', opts.Cell), 'FontWeight', 'bold', 'Interpreter', 'none');
subtitle(tl, sprintf('R_{gray} %.2f G\\Omega   pairs kept %d/%d%s   monitor slope %.2f   flashes missing %d   (%s)', R_gray, sum(P.kept), np, ...
    tern(qc_all, ' (NONE: all pairs shown, read with care)', ''), ...
    A.summary.monitor_slope, D.n_missing, strrep(strrep(folder, '\', '\\'), '_', '\_')), 'FontSize', 8, 'Interpreter', 'tex');

% 1. Vm record with currents
ax = nexttile(tl, [1 3]); hold(ax, 'on');
k1 = 1:100:numel(D.vm); tt = (k1 - 1) / D.fs / 60;
yl = [min(D.vm(k1)) max(D.vm(k1))] + [-2 2];
for b = 1:nb
    if isnan(blk.t0(b)), continue; end
    c = blk.I_pA(b); col = [1 1 1] - min(abs(c), 1) * 0.25 * [0 1 1] * (c > 0) - min(abs(c), 1) * 0.25 * [1 1 0] * (c < 0);
    x = ([blk.t0(b) blk.t1(b)] - 1) / D.fs / 60;
    patch(ax, x([1 2 2 1]), yl([1 1 2 2]), col, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
end
plot(ax, tt, D.vm(k1), 'Color', [0.2 0.2 0.2], 'LineWidth', 0.4);
bx = (blk.t0 + 3 * D.fs) / D.fs / 60;
plot(ax, bx, blk.base_mV, 'o', 'MarkerFaceColor', [1 0.8 0], 'MarkerEdgeColor', 'k', 'MarkerSize', 4);
for p = find(~P.kept).', plot(ax, bx(P{p, 1:2}), blk.base_mV(P{p, 1:2}), 'rx', 'MarkerSize', 9, 'LineWidth', 1.5); end
ylim(ax, yl); xlabel(ax, 'min'); ylabel(ax, 'Vm (mV)'); box(ax, 'on');
title(ax, 'Vm; blocks shaded red +pA, blue -pA; dots: block baselines; x: quartet failed QC', 'FontWeight', 'normal', 'FontSize', 9);

% 2. site map
ax = nexttile(tl); hold(ax, 'on');
for g = regions
    m = sites.region == g; if ~any(m), continue; end
    scatter(ax, sites.s_px(m), sites.v_px(m), 70, cols.(g), 'filled', 's', 'DisplayName', g);
end
sm9 = unique(E(im9, {'s_px', 'v_px'}), 'rows'); scatter(ax, sm9.s_px, sm9.v_px, 40, [0.5 0 0.6], 'd', 'LineWidth', 1.2, 'DisplayName', 'Mi9 series');
axis(ax, 'equal'); grid(ax, 'on'); xlabel(ax, 's (px, + toward PD)'); ylabel(ax, 'v (px)');
legend(ax, 'Location', 'southoutside', 'NumColumns', 3, 'FontSize', 7); title(ax, 'sites (PD frame)', 'FontWeight', 'normal', 'FontSize', 9);

% 3-4. 0 pA responses per region, ON / OFF
for q = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    for g = regions
        m = isr & E.region == g & E.polarity == pols(q) & E.requested_current_pA == 0;
        if ~any(m), continue; end
        plot(ax, t, mean(Xb(m, :), 1, 'omitnan'), 'Color', cols.(g), 'LineWidth', 1.2, 'DisplayName', g);
    end
    xline(ax, [0 160], ':', 'HandleVisibility', 'off'); yline(ax, 0, ':', 'HandleVisibility', 'off'); xlim(ax, [-100 500]); xlabel(ax, 'ms'); ylabel(ax, '\DeltaVm (mV)');
    title(ax, sprintf('%s flashes at 0 pA', pols(q)), 'FontWeight', 'normal', 'FontSize', 9); box(ax, 'on');
    if q == 1, legend(ax, 'Location', 'northeast', 'FontSize', 7); end
end
% 5-6. dR(t) per region
for q = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    for g = regions
        sm = reg_of == g; if ~any(sm), continue; end
        cv = squeeze(mean(dR(kept, sm, q, :), 2, 'omitnan')); if isvector(cv), cv = cv(:).'; end
        cv = movmean(cv, 20, 2); mu = mean(cv, 1, 'omitnan'); se = std(cv, 0, 1, 'omitnan') / sqrt(size(cv, 1));
        fill(ax, [t fliplr(t)], [mu + se fliplr(mu - se)], cols.(g), 'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HandleVisibility', 'off');
        plot(ax, t, mu, 'Color', cols.(g), 'LineWidth', 1.2);
    end
    xline(ax, [0 160], ':', 'HandleVisibility', 'off'); yline(ax, 0, ':', 'HandleVisibility', 'off'); xlim(ax, [-100 500]); xlabel(ax, 'ms'); ylabel(ax, '\DeltaR (G\Omega)');
    title(ax, sprintf('%s: \\DeltaR(t) = R(t) - R_{gray}, mean \\pm SE over pairs (20-ms smoothing)', pols(q)), 'FontWeight', 'normal', 'FontSize', 8); box(ax, 'on');
end

% 7. Mi9 traces at s = 0 per current; 8. Mi9 I-V
ax = nexttile(tl); hold(ax, 'on'); cc = [0.15 0.4 0.85; 0.3 0.3 0.3; 0.85 0.2 0.15];
for j = 1:numel(cur)
    m = im9 & round(E.s_px) == 0 & E.requested_current_pA == cur(j);
    plot(ax, tm, mean(Xmb(m, :), 1, 'omitnan'), 'Color', cc(min(j, 3), :), 'LineWidth', 1.2, 'DisplayName', sprintf('%+g pA', cur(j)));
end
xline(ax, [0 500], ':', 'HandleVisibility', 'off'); yline(ax, 0, ':', 'HandleVisibility', 'off'); xlim(ax, [-100 1000]); xlabel(ax, 'ms'); ylabel(ax, '\DeltaVm (mV)'); box(ax, 'on');
legend(ax, 'Location', 'best', 'FontSize', 7); title(ax, 'Mi9 series: 3-px OFF at s = 0', 'FontWeight', 'normal', 'FontSize', 9);
ax = nexttile(tl); hold(ax, 'on'); mk = 'osd';
for j = 1:numel(sq)
    m = mi9.s_px == sq(j);
    plot(ax, mi9.v_pre_mV(m), mi9.resp_100_200_mV(m), ['-' mk(min(j, 3))], 'LineWidth', 1.2, 'DisplayName', sprintf('s %+d: \\DeltaR_{OFF} %.2f G\\Omega', sq(j), dR_off(j)));
end
yline(ax, 0, ':', 'HandleVisibility', 'off'); xlabel(ax, 'pre-flash Vm (mV), at -1 / 0 / +1 pA'); ylabel(ax, 'OFF response 100-200 ms (mV)'); box(ax, 'on');
legend(ax, 'Location', 'best', 'FontSize', 7, 'Interpreter', 'tex'); title(ax, 'Mi9 series: response vs holding (sign reversal?)', 'FontWeight', 'normal', 'FontSize', 9);

% 9. a-priori read-outs
ax = nexttile(tl); hold(ax, 'on');
lab = compose("%s %s %d-%d", RO.polarity, RO.region, RO.window_ms(:, 1), RO.window_ms(:, 2));
bar(ax, 1:height(RO), RO.dR_GOhm, 0.6, 'FaceColor', [0.55 0.65 0.8]); errorbar(ax, 1:height(RO), RO.dR_GOhm, RO.se, 'k.', 'LineWidth', 1);
yline(ax, 0, '-', 'HandleVisibility', 'off'); xticks(ax, 1:height(RO)); xticklabels(ax, lab); xtickangle(ax, 20); ylabel(ax, '\DeltaR (G\Omega)'); box(ax, 'on');
title(ax, 'a-priori read-outs (mean \pm SE over pairs)', 'FontWeight', 'normal', 'FontSize', 9);

% 10. linearity and numbers
ax = nexttile(tl); axis(ax, 'off');
L = lin; txt = {'\bf0 pA mid-point residual (50-300 ms)\rm'};
for g = regions, txt{end+1} = sprintf('%-10s ON %+5.2f   OFF %+5.2f mV', g, L.residual_mV(L.region == g & L.polarity == "ON"), L.residual_mV(L.region == g & L.polarity == "OFF")); end %#ok<AGROW>
txt{end+1} = ''; txt{end+1} = '\bfblocks\rm';
for x = unique(blk.I_pA(blk.phase == "rin")).', m = blk.phase == "rin" & blk.I_pA == x; txt{end+1} = sprintf('%+4.1f pA: base %6.1f mV, monitor %+5.2f pA (n=%d)', x, median(blk.base_mV(m)), median(blk.i_mon_pA(m)), sum(m)); end %#ok<AGROW>
text(ax, 0, 1, txt, 'VerticalAlignment', 'top', 'FontName', 'Consolas', 'FontSize', 8, 'Interpreter', 'tex');
drawnow;
if ~isempty(opts.Png), exportgraphics(fig, opts.Png, 'Resolution', 120); end
A.figure = fig; A.seconds = toc(t_all);
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end

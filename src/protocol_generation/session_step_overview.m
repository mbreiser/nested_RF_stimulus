function O = session_step_overview(folder, varargin)
%SESSION_STEP_OVERVIEW  One-figure check of a recorded step: Vm, its drift, the frames shown and the current.
%
%   O = SESSION_STEP_OVERVIEW(FOLDER) reads the G4 TDMS log under FOLDER (a package folder, its 'Log Files'
%   or the time-stamped log folder; FAST_TDMS_FOLDER2STRUCT resolves it) and draws
%     1  Vm over the whole step (min/max envelope), the rows (Start-Display to the next) shaded by pattern ID
%     2  the median Vm of each row with its 5-95 % range: drift and stability across the step
%     3  the frame log (pattern frame shown over time): the stimuli reached the arena
%     4  the current monitor (pA) when it is recorded (ADC1 Vm / ADC2 current layout), else nothing
%   O.summary holds the numbers, O.line a one-line text for the session log, O.seconds the time taken.
%
%   Options (name/value): 'Title' (the folder name), 'Png' (file to write, '' = none), 'Figure' (true:
%   visible window; false: drawn off screen, only for the PNG), 'Name' (window name, reused between calls).
%   Channel layout as RAPID_READOUT: an integer-valued first ADC channel is the frame index (historic,
%   ADC0 frames / ADC1 Vm x 10), otherwise ADC1 is Vm (x 10 mV/V) and ADC2 the current (x 2000 pA/V).
%   27 October 2026.

t_all = tic;
opts = struct('Title', '', 'Png', '', 'Figure', true, 'Name', 'Session step overview');
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
folder = char(folder);
[~, fname] = fileparts(folder);
if isempty(opts.Title), opts.Title = fname; end

[Log, info] = fast_tdms_folder2struct(folder);
t_read = toc(t_all);

% ---- channels
av = Log.ADC.Volts; at = double(Log.ADC.Time);
v0 = av(1, :);
if mean(abs(v0 - round(v0)) < 1e-3) > 0.9
    vm = double(av(2, :)) * 10; tick = at(2, :); ipa = []; layout = 'ADC0 frames / ADC1 Vm';
else
    vm = double(av(1, :)) * 10; tick = at(1, :); layout = 'ADC1 Vm / ADC2 current';
    if size(av, 1) >= 2, ipa = double(av(2, :)) * 2000; else, ipa = []; end
end
ups = median(diff(tick)) * 1e4;                        % clock units per second (10 kHz sampling)
t0 = tick(1);
tmin = @(t) (double(t) - t0) / ups / 60;
dur_min = tmin(tick(end));

% ---- rows from the command log
nm = Log.Commands.Name; dt = Log.Commands.Data; tc = double(Log.Commands.Time(:).');
dec = @(d) [1 256] * hex2dec(strsplit(d, '/'));      % 'lo/hi' bytes in hex
is = find(strcmp(nm, 'Start-Display'));
rows = struct('t0', {}, 't1', {}, 'pid', {});
pid = NaN;
for j = 1:numel(nm)
    if strcmp(nm{j}, 'Set Pattern ID'), pid = dec(dt{j}); end
    if strcmp(nm{j}, 'Start-Display'), rows(end+1) = struct('t0', tc(j), 't1', NaN, 'pid', pid); end %#ok<AGROW>
end
for j = 1:numel(rows)
    if j < numel(rows), rows(j).t1 = rows(j + 1).t0; else, rows(j).t1 = rows(j).t0 + dec(dt{is(j)}) / 10 * ups; end
end
if isempty(rows), rows = struct('t0', t0, 't1', tick(end), 'pid', NaN); end

% ---- per-row Vm (1 kHz) and frame changes
v1k = vm(1:10:end); t1k = tick(1:10:end);
if isempty(Log.Frames.Time), fr_t = []; fr_p = []; else, fr_t = double(Log.Frames.Time(1, :)); fr_p = double(Log.Frames.Position(1, :)); end
nr = numel(rows); med = nan(1, nr); lo = med; hi = med; nchg = zeros(1, nr); rmid = med;
for j = 1:nr
    w = t1k >= rows(j).t0 & t1k < rows(j).t1;
    if any(w)
        x = sort(v1k(w)); n = numel(x);
        med(j) = x(ceil(n / 2)); lo(j) = x(max(1, round(0.05 * n))); hi(j) = x(max(1, round(0.95 * n)));
    end
    rmid(j) = tmin((rows(j).t0 + min(rows(j).t1, tick(end))) / 2);
    f = fr_p(fr_t >= rows(j).t0 & fr_t < rows(j).t1);
    nchg(j) = nnz(diff(f) ~= 0);
end

% ---- summary numbers
first = v1k(t1k < t0 + 60 * ups); last = v1k(t1k > tick(end) - 60 * ups);
sm = struct('folder', info.tdms_folder, 'layout', layout, 'duration_min', round(dur_min, 2), 'n_rows', numel(is), ...
    'vm_median_mV', round(median(v1k), 2), 'vm_first_min_mV', round(median(first), 2), 'vm_last_min_mV', round(median(last), 2), ...
    'vm_drift_mV', round(median(last) - median(first), 2), 'vm_min_mV', round(min(vm), 2), 'vm_max_mV', round(max(vm), 2), ...
    'noise_mV', round(median(abs(diff(vm))) / 0.6745 / sqrt(2), 3), ...
    'rows_with_frame_changes', nnz(nchg > 0), 'current_median_pA', [], 'current_first_last_pA', []);
if ~isempty(ipa)
    sm.current_median_pA = round(median(ipa(1:10:end)), 1);
    sm.current_first_last_pA = round([median(ipa(tick < t0 + 60 * ups)), median(ipa(tick > tick(end) - 60 * ups))], 1);
end
line = sprintf('%d rows, %.1f min; Vm median %.1f mV (first min %.1f, last min %.1f, drift %+.1f mV), range %.1f..%.1f, noise %.2f mV', ...
    sm.n_rows, sm.duration_min, sm.vm_median_mV, sm.vm_first_min_mV, sm.vm_last_min_mV, sm.vm_drift_mV, sm.vm_min_mV, sm.vm_max_mV, sm.noise_mV);
if sm.n_rows > 0, line = sprintf('%s; frames changed in %d/%d rows', line, sm.rows_with_frame_changes, sm.n_rows); end
if ~isempty(ipa), line = sprintf('%s; current %.0f pA (first/last min %.0f / %.0f)', line, sm.current_median_pA, sm.current_first_last_pA); end

% ---- figure
fig = findobj(groot, 'Type', 'figure', 'Name', opts.Name);
if isempty(fig)
    fig = figure('Name', opts.Name, 'NumberTitle', 'off', 'Color', 'w', 'Position', [620 80 1180 820], 'Visible', tern(opts.Figure, 'on', 'off'));
else
    fig = fig(1); clf(fig); fig.Visible = tern(opts.Figure, 'on', 'off');
end
tl = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s', opts.Title), 'Interpreter', 'none', 'FontWeight', 'bold');
subtitle(tl, line, 'FontSize', 9);

% 1: Vm envelope with the rows shaded by pattern ID
ax1 = nexttile(tl, [2 1]); hold(ax1, 'on');
[te, mn, mx] = envelope_of(tmin(tick), vm, 4000);
ylim_v = [min(mn) max(mx)] + [-1 1] * max(1, 0.05 * (max(mx) - min(mn)));
shade_rows(ax1, rows, tmin, ylim_v);
fill(ax1, [te fliplr(te)], [mn fliplr(mx)], [0.1 0.1 0.1], 'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 0.3);
ylim(ax1, ylim_v); xlim(ax1, [0 max(dur_min, eps)]); ylabel(ax1, 'Vm (mV)'); grid(ax1, 'on'); box(ax1, 'on');
ax1.XTickLabel = {};

% 2: per-row median Vm and 5-95 % range
ax2 = nexttile(tl); hold(ax2, 'on');
ok = isfinite(med);
if any(ok)
    errorbar(ax2, rmid(ok), med(ok), med(ok) - lo(ok), hi(ok) - med(ok), '.', 'Color', [0.7 0.7 0.7], 'CapSize', 0);
    plot(ax2, rmid(ok), med(ok), '.', 'Color', [0.1 0.3 0.7], 'MarkerSize', 9);
    yline(ax2, sm.vm_median_mV, ':', sprintf('%.1f mV', sm.vm_median_mV), 'Color', [0.1 0.3 0.7], 'LabelHorizontalAlignment', 'right');
    ylim(ax2, [min(lo(ok)) - 1, max(med(ok)) + max(3, 0.5 * (max(med(ok)) - min(lo(ok))))]);   % keep the medians readable
end
xlim(ax2, [0 max(dur_min, eps)]); ylabel(ax2, {'row median', 'Vm (mV)'}); grid(ax2, 'on'); box(ax2, 'on'); ax2.XTickLabel = {};

% 3: frames shown, and the current when recorded
ax3 = nexttile(tl); hold(ax3, 'on');
if ~isempty(fr_t)
    step = max(1, floor(numel(fr_t) / 20000));
    plot(ax3, tmin(fr_t(1:step:end)), fr_p(1:step:end), '.', 'Color', [0.2 0.55 0.3], 'MarkerSize', 3);
end
ylabel(ax3, 'frame'); grid(ax3, 'on'); box(ax3, 'on');
if ~isempty(ipa)
    yyaxis(ax3, 'right');
    [tb, ib] = bin_mean(tmin(tick), ipa, 10000);       % 1-s means: the monitor noise hides the steps otherwise
    plot(ax3, tb, ib, '-', 'Color', [0.8 0.35 0.1], 'LineWidth', 1.2);
    ylabel(ax3, 'current, 1-s mean (pA)'); ax3.YColor = [0.8 0.35 0.1];
    mi = median(ib); ylim(ax3, mi + max(10, max(abs(ib - mi))) * [-1.1 1.1]);   % at least +-10 pA: noise is not a step
    yyaxis(ax3, 'left'); ax3.YColor = [0.2 0.55 0.3];
end
xlim(ax3, [0 max(dur_min, eps)]); xlabel(ax3, 'time in step (min)');
linkaxes([ax1 ax2 ax3], 'x');

if ~isempty(opts.Png)
    d = fileparts(opts.Png); if ~isempty(d) && ~isfolder(d), mkdir(d); end
    exportgraphics(fig, opts.Png, 'Resolution', 90);
end
if opts.Figure, drawnow; end

O = struct('summary', sm, 'line', line, 'png', opts.Png, 'figure', fig, 'seconds', toc(t_all), 'read_seconds', t_read);
end

function [te, mn, mx] = envelope_of(t, x, nb)
% min/max of x in nb equal bins (fast reshape; the tail that does not fill a bin is dropped)
n = numel(x); k = max(1, floor(n / nb)); m = floor(n / k);
X = reshape(x(1:k * m), k, m); T = reshape(t(1:k * m), k, m);
mn = min(X, [], 1); mx = max(X, [], 1); te = T(1, :);
end

function [tb, xb] = bin_mean(t, x, k)
% means of x over consecutive bins of k samples
m = floor(numel(x) / k);
if m < 1, tb = t; xb = x; return; end
xb = mean(reshape(x(1:k * m), k, m), 1); tb = t(round(k / 2):k:k * m);
end

function shade_rows(ax, rows, tmin, yl)
% one patch: a band per row, colored by pattern ID (alternating lightness so neighbours of one ID separate)
if isempty(rows), return; end
ids = [rows.pid]; ids(isnan(ids)) = -1; u = unique(ids);
cmap = lines(max(numel(u), 1)); cmap = 0.25 * cmap + 0.75;   % pale
nr = numel(rows); V = zeros(4 * nr, 2); F = reshape(1:4 * nr, 4, nr).'; C = zeros(nr, 3);
for j = 1:nr
    a = tmin(rows(j).t0); b = tmin(rows(j).t1);
    V(4 * j - 3:4 * j, :) = [a yl(1); b yl(1); b yl(2); a yl(2)];
    C(j, :) = cmap(u == ids(j), :) * (1 - 0.06 * mod(j, 2));
end
patch(ax, 'Faces', F, 'Vertices', V, 'FaceVertexCData', C, 'FaceColor', 'flat', 'EdgeColor', 'none');
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end

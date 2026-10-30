function D = load_step_events(folder)
%LOAD_STEP_EVENTS  A recorded sequencer step (p2b, p3b, p4) as Vm, current and ledger flashes on the frame log.
%
%   D = LOAD_STEP_EVENTS(FOLDER) reads the G4 TDMS log of a step package written by RUN_SESSION
%   (FOLDER holds Functions/, Reports/event_ledger.csv, Reports/run_status.json and Log Files/) with
%   FAST_TDMS_FOLDER2STRUCT and returns
%     vm, im       Vm (mV, ADC1 x 10) and the current monitor (pA, ADC2 x 2000; relative), 10 kHz
%     fs           10000
%     rows         one per Start-Display: t0 (sample index), n (samples), pid, fid, label (run_status.json)
%     E            the event ledger, one row per flash, plus
%                    k_on      sample index of the flash onset in the FRAME LOG (NaN when not shown)
%                    dur_log_s its duration in the frame log
%                    row       the Start-Display row that played it
%     B            the block ledger (Reports/block_ledger.csv) when there is one, else []
%     n_missing    ledger flashes not found in the frame log (e.g. cut by a row end)
%   Each row's flashes are the runs of non-background frames of its position function (pfnparam.func,
%   500 Hz, 1-based; the frame log is 0-based) matched in order to the frame-log runs of that row, so the
%   onsets are the arena's, not the programmed ones. The ADC1 Vm / ADC2 current layout of October 2026.
%   28 October 2026.

[Log, ~] = fast_tdms_folder2struct(folder);
av = Log.ADC.Volts; at = double(Log.ADC.Time(1, :));
D.vm = double(av(1, :)) * 10; D.im = []; if size(av, 1) >= 2, D.im = double(av(2, :)) * 2000; end
D.fs = 10000; ups = median(diff(at)) * D.fs;
rs = jsondecode(fileread(fullfile(folder, 'Reports', 'run_status.json'))); if iscell(rs.rows), rs.rows = [rs.rows{:}]; end
nm = Log.Commands.Name; dt = Log.Commands.Data; tc = double(Log.Commands.Time(:).');
dec = @(d) [1 256] * hex2dec(strsplit(d, '/'));
is = find(strcmp(nm, 'Start-Display'));
nr = min(numel(is), numel(rs.rows));
tick2k = @(t) min(max(discretize(t, [at, Inf]), 1), numel(at));           % last ADC sample at or before t
rows = struct('t0', cell(1, nr), 'n', [], 'pid', [], 'fid', [], 'label', '');
for j = 1:nr
    q = find(strcmp(nm(1:is(j)), 'Set Pattern ID'), 1, 'last'); pid = NaN; if ~isempty(q), pid = dec(dt{q}); end
    q = find(strcmp(nm(1:is(j)), 'Set Pattern Function ID'), 1, 'last'); fid = NaN; if ~isempty(q), fid = dec(dt{q}); end
    rows(j) = struct('t0', NaN, 'n', round(dec(dt{is(j)}) / 10 * D.fs), 'pid', pid, 'fid', fid, 'label', rs.rows(j).label);
end
if nr > 0, t0 = num2cell(tick2k(tc(is(1:nr)))); [rows.t0] = t0{:}; end
D.rows = rows;

E = readtable(fullfile(folder, 'Reports', 'event_ledger.csv'), 'TextType', 'string', 'Delimiter', ',');
D.B = [];
if isfile(fullfile(folder, 'Reports', 'block_ledger.csv')), D.B = readtable(fullfile(folder, 'Reports', 'block_ledger.csv'), 'TextType', 'string', 'Delimiter', ','); end
if ~ismember('position_function', E.Properties.VariableNames)
    if ismember('block_index', E.Properties.VariableNames), E.position_function = D.B.position_function(E.block_index);
    else, E.position_function = repmat(D.B.position_function(1), height(E), 1); end
end
t_on = nan(height(E), 1); E.dur_log_s = nan(height(E), 1); E.row = nan(height(E), 1);

ft = double(Log.Frames.Time(1, :)); fp = double(Log.Frames.Position(1, :));
fb = dir(fullfile(folder, 'Functions', '*.mat')); fid_all = cellfun(@(s) str2double(s(1:4)), {fb.name});
taken = false(height(E), 1);
for r = 1:nr
    k = find(fid_all == rows(r).fid, 1); if isempty(k), continue; end
    stem = string(fb(k).name(1:end-4));
    li = find(E.position_function == stem & ~taken); if isempty(li), continue; end
    Fn = load(fullfile(folder, 'Functions', fb(k).name)); fv = double(Fn.pfnparam.func(:).');
    ef = fv([true, diff(fv) ~= 0]); ef = ef(ef ~= 1);
    li = li(1:min(numel(li), numel(ef)));
    a = tc(is(r)); e = a + rows(r).n / D.fs * ups; w = ft >= a & ft < e;
    v = fp(w); t = ft(w); st = [true, diff(v) ~= 0]; gv = v(st); gt = t(st); gd = diff([gt, e]) / ups;
    keep = gv ~= 0; gv = gv(keep); gt = gt(keep); gd = gd(keep);
    n = min(numel(li), numel(gv));
    ok = gv(1:n) + 1 == E.pattern_frame(li(1:n)).';                        % frame log 0-based
    if ~all(ok), n = find(~ok, 1) - 1; end                                  % stop at the first mismatch
    t_on(li(1:n)) = gt(1:n); E.dur_log_s(li(1:n)) = gd(1:n); E.row(li(1:n)) = r;
    taken(li) = true;
end
E.k_on = nan(height(E), 1); s = ~isnan(t_on); E.k_on(s) = tick2k(t_on(s));
D.E = E; D.n_missing = sum(isnan(E.k_on));
end

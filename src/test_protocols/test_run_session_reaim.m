% TEST_RUN_SESSION_REAIM  The arena re-aim after Protocol 1 (steps reaim, p1_confirm), mock controller (29 Oct 2026).
%
% With the bkg4 RHS2 12px_6px Protocol 1 package's patterns and a mock Protocol 1 read-out (pattern- and
% folder-dependent answers), dry run, scripted operator answers. Checks:
%   N  no edge (6-px peak at row 27.5): edge_check logged, nothing asked, p1_confirm skipped, the peak unchanged
%   T  top edge (row 9.5): suggested +18.5 deg (15 rows down), the confirmation package holds only the window's
%      squares and the Protocol 1 timing, the found peak replaces the Protocol 1 peak, decision confirm GOOD;
%      the operator typed a different angle -> that angle drives the prediction and is recorded
%   B  bottom edge (row 39.5) with reaim_sign -1: suggested +18.5 deg (the sign flips the rotation, not the
%      direction of the RF move: 15 rows up)
%   K  the operator skips ('s'): flagged, no confirmation, the peak unchanged
%   C  row 3.5: 21 rows wanted = 26.25 deg, clamped to +20 (and to +12 from a stage already at +8)
%   W  the RF moved the wrong way (mock read-out finds nothing in the window) -> stop; reaim_reset with the
%      other sign asks again and confirms
%   F  a full dry run of all steps (two-stage mode) with a top-edge peak: every step done or skipped as expected
% Writes only under tempdir.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(genpath(fullfile(repo, 'src')));
assert(startsWith(which('run_session'), repo), 'run_session is not this repository''s (%s)', which('run_session'));
proto = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\bkg4';
pkg_rhs = fullfile(proto, 'RHS2', 'protocol1_10kHz_4reps_12px_6px_RHS2_5sbkg_200msfl_150msint_81_180_07-10-25_15-47-98');
assert(isfolder(pkg_rhs), 'not found: %s', pkg_rhs);

scratch = fullfile(tempdir, ['run_session_reaim_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
mkdir(scratch); fprintf('Scratch folder: %s\n', scratch);
global NFAIL; NFAIL = 0;

pf = fullfile(scratch, 'p1_rhs'); mkdir(pf); copyfile(fullfile(pkg_rhs, 'Patterns'), fullfile(pf, 'Patterns'));
d6 = dir(fullfile(pf, 'Patterns', '*_6px_square*.mat')); q = load(fullfile(d6(1).folder, d6(1).name)); P6 = double(q.pattern.Pats);
on6 = on_frames(P6);
col = 131.5;                                            % a column of the RHS2 6-px grid (81..180)
f_at = @(row) frame_at(P6, on6, col, row);
fprintf('6-px ON frames at col %g: row 3.5 -> %d, 9.5 -> %d, 27.5 -> %d, 39.5 -> %d\n', col, f_at(3.5), f_at(9.5), f_at(27.5), f_at(39.5));

% mock read-out: Protocol 1 (folder pf) answers FIRST_ROW; the confirmation (folder p1_confirm) answers the
% square nearest the true RF after the move (TRUE_ROW_AFTER), or nothing when that is outside the window
mock = @(first_row, true_row_after) @(f, o, pat) mock_readout(f, pat, P6, on6, col, first_row, true_row_after);
base = {'dry_run', true, 'side', 'R', 'p1_folder', pf, 'overview_plots', false, 'steps', {'p1_peak', 'reaim', 'p1_confirm'}};
run1 = @(name, fn, varargin) run_session(fullfile(scratch, name), base{:}, 'cell', name, 'p1_readout_fn', fn, varargin{:});

%% N. no edge
S = run1('n_center', mock(27.5, NaN), 'operator_answers', {'should not be used'});
d = last_decision(S, 'edge_check');
check(strcmp(status(S, 'reaim'), 'done') && strcmp(status(S, 'p1_confirm'), 'skipped') && ~S.out.reaim.needed && ...
    contains(d.outcome, 'no re-aim') && ~d.fallback && isempty(find_decision(S, 'rotation')) && S.out.p1_peak.aim_xy(2) == 27.5, ...
    sprintf('N: peak at row 27.5 (%g / %g px room): no re-aim, nothing asked, p1_confirm skipped', S.out.reaim.room_px));

%% T. top edge, the operator moves 16 deg instead of the suggested 18.5
S = run1('t_top', mock(9.5, 9.5 + 16/1.25), 'operator_answers', {'', '16'});
A = S.out.reaim; r = last_decision(S, 'rotation'); c = last_decision(S, 'confirm');
check(A.needed && A.done && strcmp(A.edge, 'top') && A.suggested_deg == 18.5 && A.actual_deg == 16 && ~A.clamped && ...
    abs(A.predicted_xy(2) - (9.5 + 12.8)) < 1e-9 && contains(r.flags, 'moved +16 deg, suggested +18.5'), ...
    sprintf('T: top edge (row 9.5): suggested %+g deg, moved %+g deg recorded; predicted row %.1f', A.suggested_deg, A.actual_deg, A.predicted_xy(2)));
D = jsondecode(fileread(fullfile(S.exp_folder, 'decisions.json'))); if iscell(D), D = [D{:}]; end
rin = jsondecode(D(strcmp({D.decision}, 'rotation')).inputs);
check(rin.suggested_deg == 18.5 && rin.actual_deg == 16 && rin.sign == 1, 'T: decisions.json holds the suggested and the actual angle');
L = readtable(fullfile(S.exp_folder, 'p1_confirm', 'Reports', 'event_ledger.csv'), 'TextType', 'string', 'Delimiter', ',');
rows = jsondecode(fileread(fullfile(S.exp_folder, 'p1_confirm', 'Reports', 'rows.json'))).rows;
pfn = dir(fullfile(S.exp_folder, 'p1_confirm', 'Functions', '0002_*.mat')); fm = load(fullfile(pfn.folder, pfn.name)); fu = fm.pfnparam.func;
check(all(abs(L.square_center_y - A.predicted_xy(2)) <= 12) && all(abs(L.square_center_x - col) <= 12) && max(L.repeat) == 2 && ...
    sum(L.repeat == 1) == 2 * sum(L.repeat == 1 & L.polarity == "ON") && numel(rows) == 4 && all(fu(1:500) == 1) && sum(fu(501:600) == fu(501)) == 100 && ...
    isfile(fullfile(S.exp_folder, 'p1_confirm', 'Patterns', ['0001_' d6(1).name(6:end)])), ...
    sprintf('T: confirmation package: %d squares (ON and OFF) within 12 px of the prediction, 2 repeats, 200 ms flashes after 1 s gray, the P1 pattern as 0001, %.1f s', ...
    sum(L.repeat == 1 & L.polarity == "ON"), sum([rows.duration_s])));
pk = S.out.p1_peak;
check(strcmp(status(S, 'p1_confirm'), 'done') && contains(c.outcome, 'GOOD') && ~c.fallback && pk.aim_xy(2) == 21.5 && pk.before_reaim.aim_xy(2) == 9.5 && ...
    contains(pk.pattern, 'p1_confirm') && pk.reaim_deg == 16, ...
    sprintf('T: confirmation GOOD, found row %g (predicted %.1f); the Protocol 1 peak is now there, the old one kept', pk.aim_xy(2), A.predicted_xy(2)));

%% B. bottom edge, reaim_sign -1
S = run1('b_bottom', mock(39.5, 39.5 - 18.5/1.25), 'reaim_sign', -1);
A = S.out.reaim; c = last_decision(S, 'confirm');
check(strcmp(A.edge, 'bottom') && A.suggested_deg == 18.5 && A.actual_deg == 18.5 && abs(A.predicted_xy(2) - (39.5 - 14.8)) < 1e-9 && ...
    contains(c.outcome, 'GOOD') && S.out.p1_peak.aim_xy(2) == 27.5, ...
    sprintf('B: bottom edge (row 39.5), sign -1: suggested %+g deg, RF predicted %.1f px up at row %.1f, found row %g', A.suggested_deg, 14.8, A.predicted_xy(2), S.out.p1_peak.aim_xy(2)));
S = run1('b_bottom_plus', mock(39.5, NaN), 'steps', {'p1_peak', 'reaim'});
check(S.out.reaim.suggested_deg == -18.5, sprintf('B: same peak with sign +1: suggested %+g deg', S.out.reaim.suggested_deg));

%% K. the operator skips
S = run1('k_skip', mock(9.5, NaN), 'operator_answers', {'s'});
r = last_decision(S, 'rotation');
check(strcmp(status(S, 'reaim'), 'done') && strcmp(status(S, 'p1_confirm'), 'skipped') && ~S.out.reaim.done && S.out.reaim.skipped && ...
    r.fallback && strcmp(r.flags, 'operator_skip') && S.out.p1_peak.aim_xy(2) == 9.5 && ~isfield(S.out.p1_peak, 'before_reaim'), ...
    'K: skip at the prompt: flagged operator_skip, no confirmation, the peak unchanged');

%% C. clamped at the stage limit
S = run1('c_clamp', mock(3.5, 3.5 + 20/1.25));
A = S.out.reaim; r = last_decision(S, 'rotation');
check(A.clamped && A.suggested_deg == 20 && abs(A.wanted_deg - 26.25) < 1e-9 && contains(r.flags, 'clamped') && abs(A.predicted_xy(2) - 19.5) < 1e-9 && ...
    contains(last_decision(S, 'confirm').outcome, 'GOOD'), sprintf('C: row 3.5 wants %+.2f deg: clamped to %+g, predicted row %.1f', A.wanted_deg, A.suggested_deg, A.predicted_xy(2)));
S = run1('c_clamp_from8', mock(3.5, NaN), 'stage_angle_deg', 8, 'steps', {'p1_peak', 'reaim'});
A = S.out.reaim;
check(A.clamped && A.suggested_deg == 12 && A.stage_after_deg == 20, sprintf('C: stage already at +8 deg: suggested %+g deg (to +20)', A.suggested_deg));
S = run1('c_at_limit', mock(3.5, NaN), 'stage_angle_deg', 20, 'steps', {'p1_peak', 'reaim'}, 'operator_answers', {'should not be used'});
check(S.out.reaim.suggested_deg == 0 && strcmp(last_decision(S, 'rotation').flags, 'stage_at_limit') && ~S.out.reaim.done, 'C: stage at +20 deg already: no rotation possible, flagged, nothing asked');

%% W. wrong sign: the RF went toward the edge; stop, then reset with the other sign
name = 'w_wrong_sign';
S = run1(name, mock(9.5, 9.5 - 18.5/1.25));
st = S.steps(strcmp({S.steps.id}, 'p1_confirm'));
check(strcmp(st.status, 'interrupted') && contains(st.note, 'reaim_sign') && last_decision(S, 'confirm').fallback && S.out.p1_peak.aim_xy(2) == 9.5, ...
    'W: nothing in the window (the RF moved the other way): the session stops at p1_confirm, the peak unchanged');
S = run1(name, mock(9.5, 9.5 + 18.5/1.25), 'reaim_reset', true, 'reaim_sign', -1, 'operator_answers', {'', '-18.5'});
check(strcmp(status(S, 'p1_confirm'), 'done') && S.out.reaim.suggested_deg == -18.5 && contains(last_decision(S, 'confirm').outcome, 'GOOD') && ...
    any(strcmp({S.decisions.decision}, 'reset')) && S.out.p1_peak.aim_xy(2) == 21.5, ...
    'W: reaim_reset with reaim_sign -1: asked again (-18.5 deg), confirmed at row 21.5');

%% F. full dry run, two-stage mode, top-edge peak (Protocol 2 placeholders; read-out tables of JY1008)
analysis_repo = fullfile(fileparts(repo), 'ephys-standard-analysis');
tools_out = fullfile(analysis_repo, 'tools', 'out');
p1_package = fullfile(repo, 'protocols', 'LHS', 'protocol1_10kHz_4reps_12px_6px_LHS_2sbkg_200msfl_50msint_12-13-24_14-33-03');
if ~isfile(fullfile(tools_out, 'JY1008_prot3b_blocks.csv')) || ~isfolder(p1_package)
    fprintf('  SKIP  F: read-out tables or the LHS package not found\n');
else
    rd = fullfile(scratch, 'readout_src'); mkdir(rd);
    for f = dir(fullfile(tools_out, 'JY1008_*')).', copyfile(fullfile(f.folder, f.name), fullfile(rd, f.name)); end
    for nm = {'stage1.json', 'stage2_sequence.csv', 'stage2_bars.csv'}, copyfile(fullfile(tools_out, ['TUTL_20251023a_' nm{1}]), fullfile(rd, ['JY1008_' nm{1}])); end
    P1 = load(fullfile(p1_package, 'Patterns', dir(fullfile(p1_package, 'Patterns', '0002_*.mat')).name)); Pl = double(P1.pattern.Pats);
    fl = find(arrayfun(@(f) is_square(Pl(:, :, f), 9.5), 2:size(Pl, 3)), 1) + 1;    % a 6-px square at rows 7..12
    S = run_session(fullfile(scratch, 'f_full'), 'dry_run', true, 'cell', 'JY1008', 'side', 'L', 'p1_package', p1_package, 'peak_frame', fl, ...
        'readout_dir', rd, 'analysis_repo', analysis_repo, 'overview_plots', false);
    ids = {S.steps.id}; sts = {S.steps.status};
    check(isequal(ids, {'dark30', 'p1', 'p1_peak', 'reaim', 'p1_confirm', 'p2a', 'readout1', 'p2b', 'readout', 'build', 'p3b', 'p4'}) && all(strcmp(sts, 'done')) && ...
        S.out.p1_peak.before_reaim.aim_xy(2) == 9.5 && abs(S.out.p1_confirm.predicted_xy(2) - 24.3) < 1e-9 && ~isempty(S.finished), ...
        sprintf('F: full dry run with a top-edge peak (frame %d): 12 steps done, re-aimed +18.5 deg, the confirmation took row %g', fl, S.out.p1_peak.aim_xy(2)));
end
fprintf('\n%d check(s) failed.\n', NFAIL);

function Q = mock_readout(folder, pat, P6, on6, col, first_row, true_row_after)
[~, n] = fileparts(pat);
if contains(n, '12px'), Q = struct('frame', 2, 'xy', [], 'peak_mV', 0); return; end
if ~contains(folder, 'p1_confirm'), Q = struct('frame', frame_at(P6, on6, col, first_row), 'xy', [], 'peak_mV', 8); return; end
% the confirmation: the window's frames are those the package played (its ledger)
L = readtable(fullfile(folder, 'Reports', 'event_ledger.csv'), 'TextType', 'string', 'Delimiter', ',');
L = L(L.polarity == "ON", :);
[dmin, j] = min(hypot(L.square_center_x - col, L.square_center_y - true_row_after));
if isnan(true_row_after) || dmin > 4.5, Q = struct('frame', L.pattern_frame(1), 'xy', [], 'peak_mV', 0.4); return; end
Q = struct('frame', L.pattern_frame(j), 'xy', [], 'peak_mV', 7);
end

function on = on_frames(P)
on = false(1, size(P, 3));
for f = 2:size(P, 3), x = P(:, :, f); on(f) = max(x(:)) > mode(x(:)); end
end

function f = frame_at(P, on, col, row)
f = NaN;
for k = find(on)
    x = P(:, :, k); [r, c] = find(x ~= mode(x(:)));
    if median(c) == col && median(r) == row, f = k; return; end
end
end

function t = is_square(x, row)
[r, c] = find(x ~= mode(x(:))); t = ~isempty(r) && median(r) == row && max(x(:)) > mode(x(:)) && median(c) > 40 && median(c) < 90;
end

function s = status(S, id)
s = S.steps(strcmp({S.steps.id}, id)).status;
end

function d = find_decision(S, name)
d = S.decisions(strcmp({S.decisions.decision}, name));
end

function d = last_decision(S, name)
D = S.decisions; d = D(find(strcmp({D.decision}, name), 1, 'last'));
end

function check(cond, msg)
global NFAIL
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); NFAIL = NFAIL + 1; end
end

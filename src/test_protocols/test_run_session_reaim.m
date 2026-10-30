% TEST_RUN_SESSION_REAIM  The arena re-aim after Protocol 1 (steps reaim, p1_confirm), mock controller (29 Oct 2026).
%
% With the bkg4 RHS2 12px_6px Protocol 1 package's patterns and a mock Protocol 1 read-out (pattern- and
% folder-dependent answers), dry run, scripted operator answers. Checks:
%   N  no edge (6-px peak at row 27.5): edge_check logged, nothing asked, p1_confirm skipped, the peak unchanged
%   T  top edge (row 9.5): suggested +18.5 deg (15 rows down), the confirmation package holds every ON square of
%      the Protocol 1 pattern (no OFF) and the Protocol 1 timing, the found peak replaces the Protocol 1 peak,
%      decision confirm GOOD with its SNR and ratio; the operator typed a different angle -> that angle drives
%      the prediction and is recorded
%   V  reaim_confirm 'window': the 29 Oct window of ON and OFF squares around the prediction
%   Q  the QC (30 Oct 2026): noise only as JY1030a (above the 2 mV floor) -> NO RESPONSE; a small response
%      with little noise fails the amplitude test; a weak real RF passes; a peak on the scan's outer row ->
%      CHECK; a large offset from the prediction -> GOOD with a geometry flag; the NO RESPONSE choices p
%      (rotate back, the old peak), r (rotate back, redo the re-aim), c (repeat), an invalid answer, u (a clear
%      but weaker peak used anyway, flagged)
%   B  bottom edge (row 39.5) with reaim_sign -1: suggested +18.5 deg (the sign flips the rotation, not the
%      direction of the RF move: 15 rows up)
%   K  the operator skips ('s'): flagged, no confirmation, the peak unchanged
%   C  row 3.5: 21 rows wanted = 26.25 deg, clamped to +20 (and to +12 from a stage already at +8)
%   W  the RF moved the wrong way, off the screen (noise only in the scan) -> NO RESPONSE, s stops; reaim_reset
%      with the other sign asks again and confirms
%   G  the gate after the confirmation (figure saved; c repeats the check, a re-aims again, s stops and resumes)
%   F  a full dry run of all steps (two-stage mode) with a top-edge peak: every step done or skipped as expected
% Tests run_session_v2 (the re-aim version); run_session.m is the unchanged fallback.
% Writes only under tempdir.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(genpath(fullfile(repo, 'src')));
assert(startsWith(which('run_session_v2'), repo), 'run_session is not this repository''s (%s)', which('run_session_v2'));
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

% mock read-out: Protocol 1 (folder pf) answers FIRST_ROW (peak first_mV, default 8 mV); the confirmation
% (folder p1_confirm) simulates every square of its ledger, every repeat: amp mV (default 7) on the ON square
% nearest the true RF after the move (TRUE_ROW_AFTER; none when it is more than 4.5 px from a square of the
% scan) plus Gaussian noise of noise mV (default 1) per flash, then the map, the peak and snr_mean as
% rapid_p1_peak computes them
mock = @(first_row, true_row_after, varargin) @(f, o, pat) mock_readout(f, pat, P6, on6, col, first_row, true_row_after, varargin{:});
n6 = sum(on6);
base = {'dry_run', true, 'side', 'R', 'p1_folder', pf, 'overview_plots', false, 'steps', {'p1_peak', 'reaim', 'p1_confirm'}};
run1 = @(name, fn, varargin) run_session_v2(fullfile(scratch, name), base{:}, 'cell', name, 'p1_readout_fn', fn, varargin{:});

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
cfg = jsondecode(fileread(fullfile(S.exp_folder, 'p1_confirm', 'Reports', 'config.json')));
t_full = sum([rows.duration_s]);
check(all(L.polarity == "ON") && sum(L.repeat == 1) == n6 && sum(L.repeat == 2) == n6 && numel(unique(L.pattern_frame)) == n6 && max(L.repeat) == 2 && ...
    numel(rows) == 4 && all(fu(1:500) == 1) && sum(fu(501:600) == fu(501)) == 100 && strcmp(cfg.mode, 'full') && ...
    abs(S.steps(strcmp({S.steps.id}, 'p1_confirm')).est_min - t_full/60) < 1e-9 && ...
    isfile(fullfile(S.exp_folder, 'p1_confirm', 'Patterns', ['0001_' d6(1).name(6:end)])), ...
    sprintf('T: confirmation package: all %d ON squares of the Protocol 1 pattern, no OFF, 2 repeats, 200 ms flashes after 1 s gray, the P1 pattern as 0001: %.1f s (%.2f min)', ...
    n6, t_full, t_full/60));
pk = S.out.p1_peak; cin = jsondecode(c.inputs);
check(strcmp(status(S, 'p1_confirm'), 'done') && contains(c.outcome, 'GOOD') && ~c.fallback && pk.aim_xy(2) == 21.5 && pk.before_reaim.aim_xy(2) == 9.5 && ...
    contains(pk.pattern, 'p1_confirm') && pk.reaim_deg == 16 && cin.snr_mean >= 4 && cin.ratio_to_pre_move >= 0.5 && cin.qc.applied && ...
    contains(c.outcome, 'SNR') && strcmp(S.out.p1_confirm.verdict, 'GOOD'), ...
    sprintf('T: confirmation GOOD (%.1f mV, SNR %.1f, %.2f x the peak before), found row %g (predicted %.1f); the Protocol 1 peak is now there, the old one kept', ...
    cin.peak_mV, cin.snr_mean, cin.ratio_to_pre_move, pk.aim_xy(2), A.predicted_xy(2)));

%% V. the window option (the 29 Oct confirmation): ON and OFF squares within 12 px of the prediction
S = run1('v_window', mock(9.5, 9.5 + 18.5/1.25), 'reaim_confirm', 'window');
L = readtable(fullfile(S.exp_folder, 'p1_confirm', 'Reports', 'event_ledger.csv'), 'TextType', 'string', 'Delimiter', ',');
A = S.out.reaim; rows = jsondecode(fileread(fullfile(S.exp_folder, 'p1_confirm', 'Reports', 'rows.json'))).rows;
check(all(abs(L.square_center_y - A.predicted_xy(2)) <= 12) && all(abs(L.square_center_x - col) <= 12) && ...
    sum(L.repeat == 1) == 2 * sum(L.repeat == 1 & L.polarity == "ON") && contains(last_decision(S, 'confirm').outcome, 'GOOD') && S.out.p1_peak.aim_xy(2) == 21.5, ...
    sprintf('V: reaim_confirm window: %d squares (ON and OFF) within 12 px of the prediction, %.1f s; GOOD at row %g', ...
    sum(L.repeat == 1 & L.polarity == "ON"), sum([rows.duration_s]), S.out.p1_peak.aim_xy(2)));

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

%% Q. the QC of the scan's peak: (a) SNR, (b) amplitude vs the peak before the move, (c) the scan's edge
% Q1: noise only, as JY1030a on 29 Oct (single-flash noise 2 mV, the Protocol 1 peak 9.8 mV): the largest square
% of noise is above the 2 mV floor (the old rule's GOOD) but fails the SNR and the amplitude tests
S = run1('q_noise', mock(9.5, NaN, 'noise', 2, 'first_mV', 9.8), 'operator_answers', {'', '', 's'});
c = last_decision(S, 'confirm'); cin = jsondecode(c.inputs); st = S.steps(strcmp({S.steps.id}, 'p1_confirm'));
check(contains(c.outcome, 'NO RESPONSE') && c.fallback && cin.peak_mV >= 2 && cin.snr_mean < 4 && cin.ratio_to_pre_move < 0.5 && ~cin.qc.snr_ok && ~cin.qc.amplitude_ok && ...
    contains(c.flags, 'SNR') && strcmp(S.out.p1_confirm.verdict, 'NO RESPONSE') && S.out.p1_peak.aim_xy(2) == 9.5 && ~isfield(S.out.p1_peak, 'before_reaim') && ...
    strcmp(st.status, 'interrupted') && contains(st.note, 'reaim_reset') && strcmp(last_decision(S, 'gate').flags, 'operator_stop'), ...
    sprintf('Q1: noise only: largest square %.1f mV (above the 2 mV floor), SNR %.1f, %.2f x the peak before -> NO RESPONSE, Protocol 2 not aimed there; s stops', ...
    cin.peak_mV, cin.snr_mean, cin.ratio_to_pre_move));
% Q2: a small response with little noise: SNR passes, the amplitude (2.5 of 8 mV) does not
S = run1('q_small', mock(9.5, 9.5 + 18.5/1.25, 'amp', 2.5, 'noise', 0.2), 'operator_answers', {'', '', 's'});
c = last_decision(S, 'confirm'); cin = jsondecode(c.inputs);
check(contains(c.outcome, 'NO RESPONSE') && cin.qc.snr_ok && ~cin.qc.amplitude_ok && contains(c.flags, 'before the move') && S.out.p1_peak.aim_xy(2) == 9.5, ...
    sprintf('Q2: 2.5 mV at SNR %.1f but %.2f x the 8 mV before the move -> NO RESPONSE (amplitude)', cin.snr_mean, cin.ratio_to_pre_move));
% Q3: a weak but real RF (4 mV before, 3.5 mV after, single-flash noise 0.8 mV) passes
S = run1('q_weak', mock(9.5, 9.5 + 18.5/1.25, 'amp', 3.5, 'noise', 0.8, 'first_mV', 4));
c = last_decision(S, 'confirm'); cin = jsondecode(c.inputs);
check(contains(c.outcome, 'GOOD') && ~c.fallback && cin.snr_mean >= 4 && cin.ratio_to_pre_move >= 0.5 && S.out.p1_peak.aim_xy(2) == 21.5 && strcmp(status(S, 'p1_confirm'), 'done'), ...
    sprintf('Q3: weak real RF: %.1f mV, SNR %.1f, %.2f x the %g mV before -> GOOD, aimed at row %g', cin.peak_mV, cin.snr_mean, cin.ratio_to_pre_move, 4, S.out.p1_peak.aim_xy(2)));
% Q4: a real response on the scan's outer row (row 3.5) -> CHECK, used; the geometry flags; still near the edge
S = run1('q_edge', mock(9.5, 3.5), 'operator_answers', {'', '', ''});
c = last_decision(S, 'confirm');
check(contains(c.outcome, 'CHECK') && c.fallback && contains(c.flags, 'outer row') && contains(c.flags, 'OPPOSITE') && contains(c.flags, 'still near the edge') && ...
    S.out.p1_peak.aim_xy(2) == 3.5 && contains(last_decision(S, 'gate').outcome, 'continue'), ...
    'Q4: real peak on the outer row of the scan -> CHECK (flagged: outer row, moved the opposite way, still near the edge); Enter aims there');
% Q5: a real response 9 px from the prediction -> GOOD; the offset is a geometry flag only
S = run1('q_offset', mock(9.5, 33.5));
c = last_decision(S, 'confirm'); cin = jsondecode(c.inputs);
check(contains(c.outcome, 'GOOD') && ~c.fallback && contains(c.flags, 'geometry') && abs(cin.offset_px(2) - 9.2) < 1e-9 && S.out.p1_peak.aim_xy(2) == 33.5, ...
    sprintf('Q5: found %+.1f rows from the prediction -> GOOD with a geometry flag [%s]', cin.offset_px(2), c.flags));
% Q6: NO RESPONSE, p: the stage back to 0, go on with the peak from before the move
S = run1('q_premove', mock(9.5, NaN, 'noise', 2, 'first_mV', 9.8), 'operator_answers', {'', '', 'p', ''});
g = last_decision(S, 'gate');
check(strcmp(status(S, 'p1_confirm'), 'done') && strcmp(g.flags, 'operator_premove') && g.fallback && S.out.p1_peak.aim_xy(2) == 9.5 && ...
    S.out.stage_now_deg == 0 && S.out.reaim.reverted && ~contains(S.out.p1_peak.pattern, 'p1_confirm'), ...
    'Q6: NO RESPONSE, p: stage back at +0 deg, the session goes on with the Protocol 1 peak from before the move (row 9.5)');
% Q7: NO RESPONSE, r: the stage back to 0, the re-aim asked again (then s stops)
S = run1('q_redo', mock(9.5, NaN, 'noise', 2, 'first_mV', 9.8), 'operator_answers', {'', '', 'r', '', '', '', 's'});
D = S.decisions; rot = D(strcmp({D.decision}, 'rotation')); g = D(strcmp({D.decision}, 'gate'));
rin = jsondecode(rot(2).inputs);
check(numel(rot) == 2 && rin.stage_before_deg == 0 && rin.suggested_deg == 18.5 && strcmp(g(1).flags, 'operator_rotate_back') && strcmp(g(2).flags, 'operator_stop') && ...
    sum(strcmp({D.decision}, 'confirm')) == 2 && S.out.p1_peak.aim_xy(2) == 9.5, ...
    'Q7: NO RESPONSE, r: stage back at +0 deg, the re-aim asked again from the old peak (+18.5 deg), scanned again; s stops');
% Q8: the operator types something else at the NO RESPONSE choice: Enter is not a choice there
S = run1('q_enter', mock(9.5, NaN, 'noise', 2, 'first_mV', 9.8), 'operator_answers', {'', '', 'x', 'c', 's'});
D = S.decisions; g = D(strcmp({D.decision}, 'gate'));
check(numel(g) == 2 && strcmp(g(1).flags, 'operator_repeat') && strcmp(g(2).flags, 'operator_stop') && sum(strcmp({D.decision}, 'confirm')) == 2, ...
    'Q8: NO RESPONSE: an invalid answer is asked again; c repeats the scan; s stops');

% Q9: a clear but weaker peak (3.2 of 8 mV, SNR well over 4): NO RESPONSE by the amplitude test; u uses it anyway
S = run1('q_accept', mock(9.5, 9.5 + 18.5/1.25, 'amp', 3.2, 'noise', 0.3), 'operator_answers', {'', '', 'u'});
c = last_decision(S, 'confirm'); g = last_decision(S, 'gate'); pk = S.out.p1_peak;
check(contains(c.outcome, 'NO RESPONSE') && strcmp(g.flags, 'operator_accept') && g.fallback && strcmp(status(S, 'p1_confirm'), 'done') && ...
    pk.aim_xy(2) == 21.5 && pk.accepted_by_operator && contains(pk.source, 'accepted by the operator') && pk.before_reaim.aim_xy(2) == 9.5 && ...
    S.out.p1_confirm.accepted_by_operator && contains(g.outcome, 'use the peak anyway'), ...
    sprintf('Q9: NO RESPONSE (%.2f x the peak before), u: the scan''s peak used anyway at row %g, flagged operator_accept, the old peak kept', ...
    jsondecode(c.inputs).ratio_to_pre_move, pk.aim_xy(2)));

%% W. wrong sign: the RF went off the top edge; stop, then reset with the other sign
name = 'w_wrong_sign';
S = run1(name, mock(9.5, 9.5 - 18.5/1.25), 'operator_answers', {'', '', 's'});
st = S.steps(strcmp({S.steps.id}, 'p1_confirm'));
check(strcmp(st.status, 'interrupted') && contains(st.note, 'reaim_sign') && last_decision(S, 'confirm').fallback && ...
    contains(last_decision(S, 'confirm').outcome, 'NO RESPONSE') && S.out.p1_peak.aim_xy(2) == 9.5, ...
    'W: no response in the scan (the RF moved off the screen): NO RESPONSE, s stops at p1_confirm, the peak unchanged');
S = run1(name, mock(9.5, 9.5 + 18.5/1.25), 'reaim_reset', true, 'reaim_sign', -1, 'operator_answers', {'', '-18.5'});
check(strcmp(status(S, 'p1_confirm'), 'done') && S.out.reaim.suggested_deg == -18.5 && contains(last_decision(S, 'confirm').outcome, 'GOOD') && ...
    any(strcmp({S.decisions.decision}, 'reset')) && S.out.p1_peak.aim_xy(2) == 21.5, ...
    'W: reaim_reset with reaim_sign -1: asked again (-18.5 deg), confirmed at row 21.5');

%% G. the gate after the confirmation: the figure, 'c' repeats the check, 'a' re-aims again from the found peak, 's' stops
S = run1('g_repeat', mock(9.5, 9.5 + 18.5/1.25), 'operator_answers', {'', '', 'c', ''});
D = S.decisions; g = D(strcmp({D.decision}, 'gate'));
check(strcmp(status(S, 'p1_confirm'), 'done') && sum(strcmp({D.decision}, 'confirm')) == 2 && numel(g) == 2 && strcmp(g(1).flags, 'operator_repeat') && ...
    contains(g(2).outcome, 'continue') && S.out.p1_peak.before_reaim.aim_xy(2) == 9.5 && S.out.p1_peak.aim_xy(2) == 21.5 && ...
    isfile(S.out.p1_confirm.png) && numel(dir(fullfile(S.exp_folder, 'plots', 'p1_confirm_reaim_*.png'))) >= 1, ...
    'G: c at the gate repeats the check (2 confirmations, Protocol 1 peak kept as before), then Enter continues; the confirmation figure is saved');
S = run1('g_again', mock(3.5, 3.5 + 4.8), 'stage_angle_deg', 14, 'operator_answers', {'', '', 'a'});
D = S.decisions; rot = D(strcmp({D.decision}, 'rotation'));
check(strcmp(D(find(strcmp({D.decision}, 'gate'), 1)).flags, 'operator_reaim_again') && numel(rot) == 2 && contains(rot(1).outcome, 'moved +6 deg') && ...
    strcmp(rot(2).flags, 'stage_at_limit') && contains(rot(2).outcome, '+20') && S.out.p1_peak.aim_xy(2) == 9.5 && S.out.p1_peak.before_reaim.aim_xy(2) == 3.5 && ...
    strcmp(status(S, 'p1_confirm'), 'skipped'), ...
    'G: a at the gate (found row 9.5, still near the edge): re-aim again from there with the stage at +20 -> no rotation left, flagged; the found peak kept');
S = run1('g_stop', mock(9.5, 9.5 + 18.5/1.25), 'operator_answers', {'', '', 's'});
st = S.steps(strcmp({S.steps.id}, 'p1_confirm'));
check(strcmp(st.status, 'interrupted') && contains(st.note, 'stopped by the operator') && S.out.p1_peak.aim_xy(2) == 9.5 && ...
    strcmp(last_decision(S, 'gate').flags, 'operator_stop'), 's at the gate: the session stops at p1_confirm, the peak as before the check');
S = run1('g_stop', mock(9.5, 9.5 + 18.5/1.25));
check(strcmp(status(S, 'p1_confirm'), 'done') && S.out.p1_peak.aim_xy(2) == 21.5, 'G: the same command resumes: the check runs again, Enter continues');

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
    S = run_session_v2(fullfile(scratch, 'f_full'), 'dry_run', true, 'cell', 'JY1008', 'side', 'L', 'p1_package', p1_package, 'peak_frame', fl, ...
        'readout_dir', rd, 'analysis_repo', analysis_repo, 'overview_plots', false);
    ids = {S.steps.id}; sts = {S.steps.status};
    check(isequal(ids, {'dark30', 'p1', 'p1_peak', 'reaim', 'p1_confirm', 'p2a', 'readout1', 'p2b', 'readout', 'build', 'p3b', 'p4'}) && all(strcmp(sts, 'done')) && ...
        S.out.p1_peak.before_reaim.aim_xy(2) == 9.5 && abs(S.out.p1_confirm.predicted_xy(2) - 24.3) < 1e-9 && ~isempty(S.finished), ...
        sprintf('F: full dry run with a top-edge peak (frame %d): 12 steps done, re-aimed +18.5 deg, the confirmation took row %g', fl, S.out.p1_peak.aim_xy(2)));
end
fprintf('\n%d check(s) failed.\n', NFAIL);

function Q = mock_readout(folder, pat, P6, on6, col, first_row, true_row_after, varargin)
o = struct('amp', 7, 'noise', 1, 'first_mV', 8);
for k = 1:2:numel(varargin), o.(varargin{k}) = varargin{k+1}; end
[~, n] = fileparts(pat);
if contains(n, '12px'), Q = struct('frame', 2, 'xy', [], 'peak_mV', 0); return; end
if ~contains(folder, 'p1_confirm'), Q = struct('frame', frame_at(P6, on6, col, first_row), 'xy', [], 'peak_mV', o.first_mV); return; end
% the confirmation: the frames the package played (its ledger), every repeat, as rapid_p1_peak's map
L = readtable(fullfile(folder, 'Reports', 'event_ledger.csv'), 'TextType', 'string', 'Delimiter', ',');
fr = unique(L.pattern_frame); nrep = max(L.repeat);
on = L(L.polarity == "ON", :);
[dmin, j] = min(hypot(on.square_center_x - col, on.square_center_y - true_row_after));
rf_frame = NaN; if ~isnan(true_row_after) && dmin <= 4.5, rf_frame = on.pattern_frame(j); end
rng(11, 'twister');
rows = struct('frame_1based', {}, 'polarity', {}, 'x_center_1based', {}, 'y_center_1based', {}, 'n_repeats', {}, 'dV_mean_mV', {}, 'dV_sd_across_repeats_mV', {});
for i = 1:numel(fr)
    li = L(L.pattern_frame == fr(i), :);
    dv = o.amp * (fr(i) == rf_frame) + o.noise * randn(1, nrep);
    rows(i) = struct('frame_1based', fr(i), 'polarity', char(li.polarity(1)), 'x_center_1based', li.square_center_x(1), 'y_center_1based', li.square_center_y(1), ...
        'n_repeats', nrep, 'dV_mean_mV', mean(dv), 'dV_sd_across_repeats_mV', std(dv));
end
r = rows(strcmp({rows.polarity}, 'ON')); a = [r.dV_mean_mV]; sig = sqrt(mean([r.dV_sd_across_repeats_mV] .^ 2));
[pk, ip] = max(a);
Q = struct('frame', r(ip).frame_1based, 'xy', [r(ip).x_center_1based, r(ip).y_center_1based], 'peak_mV', pk, 'n_repeats', nrep, ...
    'ON', struct('snr_mean', pk / (sig / sqrt(nrep)), 'noise_single_mV', sig));
Q.rows = rows;
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

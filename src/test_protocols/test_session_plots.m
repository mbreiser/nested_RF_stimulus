% TEST_SESSION_PLOTS  The Prot3b / Protocol 4 analyses and the session comparison on recorded sessions (rig PC, 28 Oct 2026).
%
% Reads Jin Yong's cell of 28 Oct, recorded twice by Auto_T4_v1 (D:\nested_RF_sessions\JY1028a and JY1028b,
% the second after the cell had degraded); writes only under tempdir. Checks that every ledger flash is found
% in the frame log, the analyses' numbers against the independent checks of that day (session_check.m in
% runbook_20261026\jy1028: R_gray from the settle windows 13.6 / 5.5 GOhm; current monitor slope near 1;
% the Protocol 4 RF on the read-out center) and that the figures and metrics files are written.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(genpath(fullfile(repo, 'src')));
runs = {'D:\nested_RF_sessions\JY1028a', 'D:\nested_RF_sessions\JY1028b'};
for r = runs, assert(isfile(fullfile(r{1}, 'session_state.json')), 'session not found: %s', r{1}); end
scratch = fullfile(tempdir, ['session_plots_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]); mkdir(scratch);
fprintf('Scratch folder: %s\n', scratch);
global NFAIL; NFAIL = 0;
snap = @(d) dir(fullfile(d, '**', '*')); before = cellfun(snap, runs, 'UniformOutput', false);

A3 = analyze_prot3b(fullfile(runs{1}, 'p3b'), 'Cell', 'JY1028a', 'Figure', false, 'Png', fullfile(scratch, 'a_p3b.png'));
s = A3.summary;
check(s.flashes_missing == 0 && height(A3.events) == 548, sprintf('p3b: all %d ledger flashes found in the frame log', height(A3.events)));
check(abs(s.R_gray_GOhm - 13.6) < 1.5 && s.pairs == 8 && s.pairs_kept >= 4, sprintf('p3b: R_gray %.2f GOhm (settle-window estimate 13.6), %d/%d pairs kept', s.R_gray_GOhm, s.pairs_kept, s.pairs));
check(abs(s.monitor_slope - 1) < 0.35, sprintf('p3b: current monitor follows the requested levels (slope %.2f per pA)', s.monitor_slope));
check(s.ON_center_dR_150_250 < 0 && height(A3.mi9) == 9, sprintf('p3b: ON center dR %.2f GOhm (a conductance increase); Mi9 table 3 squares x 3 currents', s.ON_center_dR_150_250));
A4 = analyze_p4(fullfile(runs{1}, 'p4'), 'Cell', 'JY1028a', 'Figure', false, 'Png', fullfile(scratch, 'a_p4.png'));
s = A4.summary;
check(s.flashes_missing == 0 && height(A4.events) == 1710 && s.passes == 6, sprintf('p4: all %d flashes found, %d passes', height(A4.events), s.passes));
check(hypot(s.gauss_s0_px, s.gauss_v0_px) < 3 && s.lattice_ON_peak_mV > 10 && s.gauss_sigma_s_px > 1 && s.gauss_sigma_s_px < 8, ...
    sprintf('p4: Gaussian RF on the read-out center (s %+.1f v %+.1f px), sigma %.1f x %.1f px, ON peak %.1f mV', s.gauss_s0_px, s.gauss_v0_px, s.gauss_sigma_s_px, s.gauss_sigma_v_px, s.lattice_ON_peak_mV));
R = analyze_session(runs, 'Out', scratch, 'Labels', {'run1', 'run2'});
M = R.metrics;
check(height(M) == 2 && all(cellfun(@isfile, R.pngs)) && isfile(fullfile(scratch, 'run1_session_metrics.json')) && isfile(fullfile(scratch, 'compare_run1_vs_run2.png')), ...
    sprintf('session: metrics for both runs, %d figures listed and present, comparison written', numel(R.pngs)));
check(M.R_gray_GOhm(2) < M.R_gray_GOhm(1) && M.vm_p4_mV(2) > M.vm_p4_mV(1) + 20 && M.p4_lattice_ON_peak_mV(2) < M.p4_lattice_ON_peak_mV(1), ...
    'session: the second run shows the degradation (lower R_gray, Vm 20 mV higher in p4, smaller lattice peak)');
after = cellfun(snap, runs, 'UniformOutput', false);
check(isequal(cellfun(@numel, before), cellfun(@numel, after)), 'the session folders were only read (no files added)');
fprintf('\n%d check(s) failed.\n', NFAIL);

function check(cond, msg)
global NFAIL
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); NFAIL = NFAIL + 1; end
end

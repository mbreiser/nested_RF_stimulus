% TEST_RUN_SESSION_P1_AIM  Protocol 1 peak -> Protocol 2 aim with 12px_6px packages (rig PC, 28 Oct 2026).
%
% Two bugs fixed on 28 Oct: (1) with a package holding both grids, rapid_p1_peak analysed the 12-px
% pattern (the first *square*.mat) and generate_protocol2_package looked that frame up in the rig's fixed
% 6-px pattern, aiming Protocol 2 elsewhere; (2) the 12-px fallback could never fire (frame12/peak12_mV are
% [ON OFF] of the same grid). Checks:
%   R  rapid_p1_peak's default pattern on Jin Yong's 23 Oct recording (6-px grid, pattern ID 2) with a
%      Patterns folder holding BOTH grids: the 6-px one is chosen, peak frame 239 at [97.5 33.5] as before
%   D  run_session's p1_peak decision with a mock read-out on the bkg4 RHS2 12px_6px package: the 6-px
%      pattern is analysed first; 6-px above the floor -> 6-px aim; below -> the 12-px read-out, aimed at
%      that frame of the 12-px pattern; both below -> stop; 6pxONLY package below -> stop; a two-argument
%      read-out function still works
%   G  a stage-1 package generated from the 12-px fallback is centered on the 12-px frame's square (the
%      FLASHES bars' centroid), not where the old frame lookup would have put it
% Needs G4_Display_Tools on the path (G). Writes only under tempdir.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(genpath(fullfile(repo, 'src')));
for f = {'run_session', 'rapid_p1_peak', 'create_protocol2'}
    assert(startsWith(which(f{1}), repo), '%s is not this repository''s (%s)', f{1}, which(f{1}));
end
proto = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\bkg4';
pkg_lhs6 = fullfile(proto, 'LHS', 'protocol1_10kHz_4reps_6pxONLY_LHS_2sbkg_200msfl_150msint_09-18-26_12-59-40');
pkg_lhs12 = fullfile(proto, 'LHS', 'protocol1_10kHz_4reps_12px_6px_LHS_2sbkg_200msfl_150msint_07-10-25_15-45-60');
pkg_rhs = fullfile(proto, 'RHS2', 'protocol1_10kHz_4reps_12px_6px_RHS2_5sbkg_200msfl_150msint_81_180_07-10-25_15-47-98');
p1_src = fullfile(pkg_lhs6, '10_23_2026', 'SS02344XJFRC28 (T4) - GABAR RNAi-14_31_92');
for d = {pkg_lhs6, pkg_lhs12, pkg_rhs, p1_src}, assert(isfolder(d{1}), 'not found: %s', d{1}); end

scratch = fullfile(tempdir, ['run_session_p1_aim_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
mkdir(scratch); fprintf('Scratch folder: %s\n', scratch);
global NFAIL; NFAIL = 0;

%% R. rapid_p1_peak on a recording whose package folder holds both grids
root = fullfile(scratch, 'p1_both'); cellf = fullfile(root, '10_23_2026', 'cell'); mkdir(cellf);
copyfile(p1_src, cellf);
mkdir(fullfile(root, 'Patterns'));
copyfile(fullfile(pkg_lhs12, 'Patterns', '0001_*.mat'), fullfile(root, 'Patterns'));
copyfile(fullfile(pkg_lhs6, 'Patterns', '0002_*.mat'), fullfile(root, 'Patterns'));
R = rapid_p1_peak(cellf, 'Out', fullfile(scratch, 'p1_both_out'), 'Cell', 'JY1023', 'Figure', false, 'SavePng', false);
[~, pn] = fileparts(R.pattern);
check(contains(pn, '6px') && ~contains(pn, '12px') && R.frame == 239 && isequal(R.xy, [97.5 33.5]) && R.n_flashes > 0, ...
    sprintf('R: default pattern with both grids present is %s; ON peak frame %d at %s, %d flashes', pn, R.frame, mat2str(R.xy), R.n_flashes));

%% D. the p1_peak decision with a mock read-out (pattern-dependent answers)
P6 = load_pats(pkg_rhs, '6px'); P12 = load_pats(pkg_rhs, '12px');
% a 12-px ON frame whose square is well inside the Protocol 2 crop and far from the 6-px square of the same number
f12 = NaN;
for f = 2:size(P12, 3)
    [c12, on12] = med_xy(P12(:, :, f)); if ~on12 || f > size(P6, 3), continue; end
    c6 = med_xy(P6(:, :, f));
    if c12(2) >= 16 && c12(2) <= 33 && c12(1) >= 32 && c12(1) <= 177 && norm(c12 - c6) > 10, f12 = f; break; end
end
assert(isfinite(f12), 'no suitable 12-px frame found');
[c12, ~] = med_xy(P12(:, :, f12)); c6_same = med_xy(P6(:, :, f12));
f6 = find(arrayfun(@(f) on_frame(P6(:, :, f)), 2:size(P6, 3)), 1) + 1; c6 = med_xy(P6(:, :, f6));
fprintf('12-px ON frame %d at %s (the 6-px frame %d is at %s); 6-px ON frame %d at %s\n', f12, mat2str(c12), f12, mat2str(c6_same), f6, mat2str(c6));

mock = @(pk6, pk12) @(f, o, pat) mock_readout(pat, f6, pk6, f12, pk12);
pf = fullfile(scratch, 'p1_rhs'); mkdir(pf); copyfile(fullfile(pkg_rhs, 'Patterns'), fullfile(pf, 'Patterns'));
pf6 = fullfile(scratch, 'p1_lhs6'); mkdir(pf6); copyfile(fullfile(pkg_lhs6, 'Patterns'), fullfile(pf6, 'Patterns'));
run1 = @(name, fn, folder, side) run_session(fullfile(scratch, name), 'dry_run', true, 'cell', name, 'side', side, 'p1_folder', folder, ...
    'p1_readout_fn', fn, 'steps', {'p1_peak'}, 'overview_plots', false);

S = run1('d_strong6', mock(5, 9), pf, 'R'); pk = S.out.p1_peak; d = last_decision(S, 'aim');
check(strcmp(S.steps(strcmp({S.steps.id}, 'p1_peak')).status, 'done') && pk.grid_px == 6 && pk.frame == f6 && contains(pk.pattern, '6px') && ...
    isequal(pk.aim_xy, c6) && strcmpi(pk.aim_polarity, 'ON') && ~d.fallback, ...
    sprintf('D: 6-px peak above the floor -> 6-px frame %d, aimed at %s (%s)', pk.frame, mat2str(pk.aim_xy), pk.aim_polarity));

S = run1('d_weak6', mock(1, 9), pf, 'R'); pk = S.out.p1_peak; d = last_decision(S, 'aim');
check(strcmp(S.steps(strcmp({S.steps.id}, 'p1_peak')).status, 'done') && pk.grid_px == 12 && pk.frame == f12 && contains(pk.pattern, '12px') && ...
    isequal(pk.aim_xy, c12) && d.fallback && contains(d.flags, '12-px') && pk.peak_6px_mV == 1 && pk.peak_12px_mV == 9, ...
    sprintf('D: 6-px peak below the floor -> 12-px fallback fires: frame %d aimed at %s, not at %s (the 6-px frame of that number)', pk.frame, mat2str(pk.aim_xy), mat2str(c6_same)));
S12 = S;

S = run1('d_weak_both', mock(1, 1), pf, 'R');
check(strcmp(S.steps(strcmp({S.steps.id}, 'p1_peak')).status, 'interrupted') && isempty(S.finished), 'D: both grids below the floor -> the session stops at p1_peak');

S = run1('d_weak_6only', mock(1, 9), pf6, 'L'); st = S.steps(strcmp({S.steps.id}, 'p1_peak'));
check(strcmp(st.status, 'interrupted') && contains(st.note, 'no 12-px grid'), 'D: 6pxONLY package below the floor -> stops (no 12-px grid to fall back on)');

S = run1('d_two_arg', @(f, o) struct('frame', f6, 'xy', [], 'peak_mV', 5), pf, 'R'); pk = S.out.p1_peak;
check(pk.grid_px == 6 && isequal(pk.aim_xy, c6) && isequal(pk.xy, c6), 'D: a two-argument read-out function is still accepted (6-px grid, xy from the pattern)');

S = run_session(fullfile(scratch, 'd_override'), 'dry_run', true, 'cell', 'ov', 'side', 'R', 'p1_folder', pf, 'peak_frame', f6, 'steps', {'p1_peak'}, 'overview_plots', false);
pk = S.out.p1_peak;
check(pk.grid_px == 6 && pk.frame == f6 && isequal(pk.aim_xy, c6) && contains(pk.source, 'override'), 'D: peak_frame override is a 6-px frame, aimed through the package''s 6-px pattern');

%% G. stage-1 package generated from the 12-px fallback
if isempty(which('create_exp_dir_G4'))
    fprintf('  SKIP  G: G4_Display_Tools not on the path\n');
else
    S = run_session(fullfile(scratch, 'g_weak6'), 'dry_run', true, 'cell', 'g12', 'side', 'R', 'p1_folder', pf, 'p1_readout_fn', mock(1, 9), ...
        'p2_generate_in_dry_run', true, 'p2_package_root', fullfile(scratch, 'generated'), 'steps', {'p1_peak', 'p2a'}, 'overview_plots', false);
    fl = dir(fullfile(S.out.p2a.folder, 'Patterns', '*FLASHES*.mat')); q = load(fullfile(fl(1).folder, fl(1).name)); Pt = double(q.pattern.Pats);
    on = Pt(:, :, 2:end) > 4; [rr, cc] = ndgrid(1:size(Pt, 1), 1:size(Pt, 2)); n = squeeze(sum(on, [1 2]));
    cx = squeeze(sum(on .* cc, [1 2])) ./ n; cy = squeeze(sum(on .* rr, [1 2])) ./ n; ctr = [mean(cx), mean(cy)];
    check(norm(ctr - double(int16(c12))) < 1.5 && norm(ctr - c6_same) > 8, ...
        sprintf('G: stage-1 FLASHES bars centered at %s, on the 12-px square %s (the old lookup would aim at %s)', mat2str(round(ctr, 2)), mat2str(c12), mat2str(c6_same)));
end
fprintf('\n%d check(s) failed.\n', NFAIL);

function Q = mock_readout(pat, f6, pk6, f12, pk12)
[~, n] = fileparts(pat);
if contains(n, '12px'), Q = struct('frame', f12, 'xy', [], 'peak_mV', pk12);
else, Q = struct('frame', f6, 'xy', [], 'peak_mV', pk6); end
end

function Pt = load_pats(pkg, grid)
d = dir(fullfile(pkg, 'Patterns', ['*_' grid '_square*.mat'])); q = load(fullfile(d(1).folder, d(1).name)); Pt = double(q.pattern.Pats);
end

function [xy, on] = med_xy(f)
bkg = mode(f(:)); [r, c] = find(f ~= bkg); xy = [median(c), median(r)]; on = max(f(:)) > bkg;
end

function on = on_frame(f)
[xy, on] = med_xy(f); on = on && xy(2) >= 16 && xy(2) <= 33 && xy(1) >= 32 && xy(1) <= 177;
end

function d = last_decision(S, name)
D = S.decisions; d = D(find(strcmp({D.decision}, name), 1, 'last'));
end

function check(cond, msg)
global NFAIL
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); NFAIL = NFAIL + 1; end
end

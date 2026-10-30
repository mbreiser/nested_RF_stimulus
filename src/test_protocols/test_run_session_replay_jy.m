% TEST_RUN_SESSION_REPLAY_JY  run_session offline on a recorded Jin Yong session, rig PC (28 Oct 2026).
%
% The live paths that the empty-arena runs could not exercise, on Jin Yong's cell of 23 Oct 2026
% (Protocol 1 ...\10_23_2026\SS02344XJFRC28 (T4) - GABAR RNAi-14_31_92, Protocol 2
% nested_RF_protocol2\data\2026_10_23_14_40; earlier read-out: Protocol 1 ON peak frame 239 at
% col/row [97.5 33.5], center [97.8 33.35]). Mock controller, no arena, no Host.
%   A  two-stage, as Jin Yong's session will run: the Protocol 1 peak from run_session's default
%      rapid_p1_peak call; the stage-1 package GENERATED from that peak (p2_generate_in_dry_run, into
%      the scratch folder) and compared file by file with the package setup-pharma made on 23 Oct;
%      stage-1 rows from it; read-out calls 1 and 2 by rapid_readout's emulations on the recorded
%      log ('Stage1Test', 'Stage2Test': a real stage-1 / stage-2 log only comes from the arena);
%      stage 2, build, Prot3b, Protocol 4, allOff at the end.
%   B  full mode on the recorded Protocol 2 folder: run_session's default rapid_readout call
%      (call_readout, 'P1Peak') on the real log.
% Both sessions read COPIES of the recorded folders, made in the scratch folder; the originals are
% only read (by copyfile) and checked unchanged at the end. Needs G4_Display_Tools on the path.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(genpath(fullfile(repo, 'src')));
for f = {'run_session', 'create_protocol2', 'generate_bar_pos_fns', 'rapid_readout', 'rapid_p1_peak'}
    assert(startsWith(which(f{1}), repo), '%s is not this repository''s (%s)', f{1}, which(f{1}));
end
assert(~isempty(which('create_exp_dir_G4')) && ~isempty(which('G4_Function_Generator')), 'G4_Display_Tools not on the path');
if ~exist('p1_src', 'var')
    p1_root = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\bkg4\LHS\protocol1_10kHz_4reps_6pxONLY_LHS_2sbkg_200msfl_150msint_09-18-26_12-59-40';
    p1_src = fullfile(p1_root, '10_23_2026', 'SS02344XJFRC28 (T4) - GABAR RNAi-14_31_92');
end
if ~exist('p2_src', 'var'), p2_src = 'C:\matlabroot\G4_Protocols\nested_RF_protocol2\data\2026_10_23_14_40'; end
p1_root = fileparts(fileparts(p1_src));
assert(isfolder(fullfile(p1_src, 'Log Files')) || ~isempty(dir(fullfile(p1_src, '**', '*.tdms'))), 'Protocol 1 session not found: %s', p1_src);
assert(isfolder(fullfile(p2_src, 'Log Files')), 'Protocol 2 session not found: %s', p2_src);
tag = 'JY1023r';
ref = struct('frame', 239, 'xy', [97.5 33.5], 'peak_mV', 19.94, 'center', [97.8 33.35]);

scratch = fullfile(tempdir, ['run_session_replay_jy_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
mkdir(scratch);
fprintf('Scratch folder: %s\n', scratch);
global NFAIL; NFAIL = 0;
snap = @(d) dir(fullfile(d, '**', '*'));
before1 = snap(p1_src); before2 = snap(p2_src);

% copies: the Protocol 1 cell folder under <root>\<date>\ with the root's Patterns (rapid_p1_peak finds the
% square pattern two levels up), and the whole Protocol 2 session
t0 = tic;
p1 = fullfile(scratch, 'src_p1', '10_23_2026', 'p1_cell'); mkdir(p1);
copyfile(p1_src, p1);
copyfile(fullfile(p1_root, 'Patterns'), fullfile(scratch, 'src_p1', 'Patterns'));
p2 = fullfile(scratch, 'src_p2'); mkdir(p2);
copyfile(p2_src, p2);
fprintf('copies made in %.0f s\n', toc(t0));
gen_root = fullfile(scratch, 'generated');

%% A. Two-stage with the Protocol 2 package generated from the computed peak
t0 = tic;
expA = fullfile(scratch, 'cell_A_two_stage');
rfn = @(call, S, out_dir) replay_readout(call, p2, tag, out_dir, S);
SA = run_session(expA, 'dry_run', true, 'cell', tag, 'side', 'L', 'p1_folder', p1, ...
    'p2_generate_in_dry_run', true, 'p2_package_root', gen_root, 'readout_fn', rfn, ...
    'steps', {'p1_peak', 'p2a', 'readout1', 'p2b', 'readout', 'build', 'p3b', 'p4'});
tA = toc(t0);
stA = containers.Map({SA.steps.id}, {SA.steps.status});
sel = {'p1_peak', 'p2a', 'readout1', 'p2b', 'readout', 'build', 'p3b', 'p4'};
check(all(cellfun(@(k) strcmp(stA(k), 'done'), sel)), sprintf('A: all selected steps done (%s)', strjoin(cellfun(@(k) [k ' ' stA(k)], sel, 'UniformOutput', false), ', ')));
pk = SA.out.p1_peak;
check(strcmp(pk.source, 'rapid_p1_peak') && pk.frame == ref.frame && isequal(pk.xy(:).', ref.xy) && abs(pk.peak_mV - ref.peak_mV) < 0.05, ...
    sprintf('A: default rapid_p1_peak call: ON peak frame %d at col/row %s, %.2f mV (Jin Yong noted frame %d)', pk.frame, mat2str(pk.xy(:).'), pk.peak_mV, ref.frame));
gen = SA.out.p2a.folder;
check(startsWith(gen, fullfile(gen_root, 'data')) && isfile(fullfile(gen, 'currentExp.mat')), sprintf('A: package generated under the scratch root: %s', gen));
c23 = load(fullfile(p2_src, 'currentExp.mat')); cg = load(fullfile(gen, 'currentExp.mat'));
check(isfield(c23, 'metadata') && isfield(c23.metadata, 'Frame') && double(c23.metadata.Frame) == pk.frame, ...
    sprintf('A: the 23 Oct package was made for frame %s, the replay generated it for frame %d', mat2str(getf_(c23, 'Frame')), pk.frame));
check(isequal(double(cg.pattern_order), double(c23.pattern_order)) && isequal(double(cg.func_order), double(c23.func_order)) && ...
    max(abs(double(cg.trial_dur) - double(c23.trial_dur))) < 1e-9, ...
    sprintf('A: pattern_order, func_order, trial_dur as setup-pharma''s (%d conditions)', numel(cg.pattern_order)));
for d = {'Patterns', 'Functions'}
    [same, msg] = same_files(fullfile(gen, d{1}), fullfile(p2_src, d{1}), tern_(strcmp(d{1}, 'Patterns'), '*.pat', '*.pfn'));
    check(same, sprintf('A: %s: %s', d{1}, msg));
    [same, msg] = same_names(fullfile(gen, d{1}), fullfile(p2_src, d{1}));
    check(same, sprintf('A: %s file names: %s', d{1}, msg));
end
sa = jsondecode(fileread(fullfile(gen, 'Reports', 'run_status.json')));
pf = dir(fullfile(gen, 'Patterns', '*FLASHES*.mat')); fl = str2double(pf(1).name(1:4));
n_sweep = sum(cg.pattern_order ~= fl);
rr = sa.rows; if iscell(rr), rr = [rr{:}]; end
check(isscalar(pf) && sa.rows_total == 3 * n_sweep && all([rr.pattern_id] ~= fl) && strcmp(sa.state, 'completed'), ...
    sprintf('A: stage 1 = %d conditions without the FLASHES pattern %d x 3 = %d rows, %.1f min', n_sweep, fl, sa.rows_total, sum([rr.duration_s]) / 60));
J1 = jsondecode(fileread(fullfile(expA, 'readout', [tag '_stage1.json'])));
check(J1.stage == 1 && SA.out.readout1.usable_pd && strcmp(SA.out.readout1.stage2_mode, 'bars22'), 'A: read-out call 1: usable PD, stage 2 = 22 bars');
sb = jsondecode(fileread(fullfile(expA, 'p2b', 'Reports', 'run_status.json')));
dp = dir(fullfile(expA, 'p2b', 'Patterns', '0001_*.mat'));
check(SA.out.p2b.n_flashes == 132 && isscalar(dp) && contains(dp.name, 'FLASHES') && contains(dp.name, 'X98_Y34'), ...
    sprintf('A: stage-2 package: %s, %d rows, %d flashes', dp(1).name, sb.rows_total, SA.out.p2b.n_flashes));
rsA = SA.out.readout.summary;
check(norm(rsA.center_col_row(:).' - ref.center) < 0.1, sprintf('A: read-out call 2: center col/row %s (earlier read-out %s)', mat2str(rsA.center_col_row(:).', 4), mat2str(ref.center)));
dA = SA.decisions(strcmp({SA.decisions.decision}, 'center'));
check(~dA.fallback && isequal(SA.center_shift_px, [0 0]) && SA.out.build.p3b.n_blocks == 26, ...
    sprintf('A: center accepted; Prot3b %d blocks, Protocol 4 %d frames, %.1f min', SA.out.build.p3b.n_blocks, SA.out.build.p4.n_frames, SA.out.build.p4.programmed_s / 60));
LA = jsondecode(fileread(fullfile(expA, 'readout', [tag '_locations.json'])));
F4 = readtable(fullfile(expA, 'readout', [tag '_prot4_frames.csv']), 'TextType', 'string', 'Delimiter', ',');
lat = F4(F4.item == "lattice", :); W = LA.window_px;
check(W.s_lead == -14 && W.s_trail == 18 && W.v_half == 8 && all(lat.s_px > -15 & lat.s_px < 19 & abs(lat.v_px) < 9) && ...
    (min(lat.s_px) < -13 || LA.lattice.n_dropped_off_screen > 0) && (max(lat.s_px) > 17 || LA.lattice.n_dropped_off_screen > 0), ...
    sprintf('A: Protocol 4 lattice window s %d..%+d, v +-%d px: %d squares at s %.1f..%+.1f (%d off screen)', W.s_lead, W.s_trail, W.v_half, ...
    height(lat), min(lat.s_px), max(lat.s_px), LA.lattice.n_dropped_off_screen));
check(isfield(SA.out, 'arena_off') && SA.out.arena_off.ok && ~isempty(SA.finished), 'A: session finished, allOff at the end');
fprintf('A took %.0f s\n', tA);

%% B. Full mode: run_session's own rapid_readout call on the recorded Protocol 2 log
t0 = tic;
expB = fullfile(scratch, 'cell_B_full');
SB = run_session(expB, 'dry_run', true, 'cell', tag, 'side', 'L', 'p2_mode', 'full', 'p1_folder', p1, 'p2_folder', p2, ...
    'steps', {'p1_peak', 'p2', 'readout', 'build'});
tB = toc(t0);
stB = containers.Map({SB.steps.id}, {SB.steps.status});
check(all(cellfun(@(k) strcmp(stB(k), 'done'), {'p1_peak', 'p2', 'readout', 'build'})), 'B: all selected steps done');
rsB = SB.out.readout.summary; ck = SB.out.readout.checks;
check(strcmp(SB.out.readout.method, 'rapid_readout.m') && norm(rsB.center_col_row(:).' - ref.center) < 0.1, ...
    sprintf('B: default rapid_readout call: center col/row %s', mat2str(rsB.center_col_row(:).', 4)));
L = jsondecode(fileread(fullfile(expB, 'readout', [tag '_locations.json'])));
check(isequal(L.other_centers.protocol1_ON_peak(:).', ref.xy) && L.sequencer_inputs.usable_pd && ck.near_p1 && ck.inside_bar_set && ck.peak_above_floor, ...
    sprintf('B: P1Peak passed through; %.2f px from the Protocol 1 peak, inside the bar set, center peak %.1f mV', ck.dist_to_p1_px, rsB.center_peak_mV));
check(norm(rsA.center_col_row(:).' - rsB.center_col_row(:).') < 1e-6, 'A and B: the same center (two-stage emulation and one call)');
fprintf('B took %.0f s\n', tB);

%% The recorded folders were not written to
a1 = snap(p1_src); a2 = snap(p2_src);
check(isequal({a1.name}, {before1.name}) && isequal([a1.datenum], [before1.datenum]) && isequal({a2.name}, {before2.name}) && isequal([a2.datenum], [before2.datenum]), ...
    'the recorded Protocol 1 and Protocol 2 folders are unchanged');

fprintf('\n%d check(s) failed. A %.0f s, B %.0f s. Output in %s\n', NFAIL, tA, tB, scratch);

function replay_readout(call, p2, tag, out_dir, S)
    xy = double(S.out.p1_peak.xy(:).');
    if call == 1
        rapid_readout(p2, 'Cell', tag, 'Out', out_dir, 'Figure', false, 'Stage1Test', true, 'P1Peak', xy);
    else
        rapid_readout(p2, 'Cell', tag, 'Out', out_dir, 'Figure', false, 'Stage2Test', true, 'P1Peak', xy);
    end
end

function [ok, msg] = same_files(a, b, pat)
    fa = dir(fullfile(a, pat)); fb = dir(fullfile(b, pat));
    ok = isequal(sort({fa.name}), sort({fb.name}));
    if ~ok, msg = sprintf('%d generated vs %d recorded %s files, names differ', numel(fa), numel(fb), pat); return, end
    bad = {};
    for k = 1:numel(fa)
        if ~isequal(read_bytes(fullfile(a, fa(k).name)), read_bytes(fullfile(b, fa(k).name))), bad{end+1} = fa(k).name; end %#ok<AGROW>
    end
    ok = isempty(bad);
    if ok, msg = sprintf('all %d %s files byte-identical to the 23 Oct package', numel(fa), pat);
    else, msg = sprintf('%d of %d %s files differ: %s', numel(bad), numel(fa), pat, strjoin(bad, ', ')); end
end

function [ok, msg] = same_names(a, b)
    fa = dir(fullfile(a, '*.mat')); fb = dir(fullfile(b, '*.mat'));
    ok = isequal(sort({fa.name}), sort({fb.name}));
    msg = sprintf('%d .mat files, %s', numel(fa), tern_(ok, 'same names', ['differ: ' strjoin(setxor({fa.name}, {fb.name}), ', ')]));
end

function b = read_bytes(f)
    fid = fopen(f, 'r'); b = fread(fid, Inf, 'uint8=>uint8'); fclose(fid);
end

function v = getf_(c, f)
    v = []; if isfield(c, 'metadata') && isfield(c.metadata, f), v = c.metadata.(f); end
end

function v = tern_(c, a, b)
    if c, v = a; else, v = b; end
end

function check(cond, msg)
    global NFAIL
    if cond, fprintf('  PASS  %s\n', msg); else, fprintf(2, '  FAIL  %s\n', msg); NFAIL = NFAIL + 1; end
end

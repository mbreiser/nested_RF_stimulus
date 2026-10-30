% TEST_RUN_SESSION_V2_WIDEN  Auto_T4_v2's wider stage-2 bar flashes (n_flank 7) (rig PC, 30 Oct 2026).
%
% Why: JY1028a's PD-side hyperpolarization runs from s +5 to beyond +16 px (still -0.7 mV at +16); the 11 stage-2
% bars per set cover only +/-10 px. v2 passes p2_n_flank 7 (15 bars per set, +/-14 px, 120 frames); the shared
% files keep the old values as defaults, so Auto_T4_v1 is
% unchanged (test_run_session_replay_jy checks v1's generated package byte for byte against the 23 Oct one).
%   D  the defaults: generate_bar_flash_pos_fns with n_flank 5 draws randperm(88) with 89 for frame 1, as before
%   G  run_session_v2 generates the stage-1 package from the 23 Oct Protocol 1 peak (dry run): the FLASHES pattern
%      has 121 frames (8 x 15 bars + background), the six bar-flash functions flash each of the 120 bars once; the
%      sweep conditions and their durations are those of the 23 Oct package
%   S  read-out call 1 on the 23 Oct log with the generated 120-frame FLASHES pattern: 30 stage-2 bars, 15 per set,
%      one orientation each; the stage-2 package: 30 bars x 2 classes x 3 repeats = 180 flashes, its duration
% Writes only under tempdir. Needs G4_Display_Tools on the path.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(genpath(fullfile(repo, 'src')));
for f = {'run_session_v2', 'rapid_readout', 'generate_bar_flash_pos_fns', 'generate_bar_flash_stimulus_xy'}
    assert(startsWith(which(f{1}), repo), '%s is not this repository''s (%s)', f{1}, which(f{1}));
end
assert(~isempty(which('G4_Function_Generator')), 'G4_Display_Tools not on the path');
p1_root = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\bkg4\LHS\protocol1_10kHz_4reps_6pxONLY_LHS_2sbkg_200msfl_150msint_09-18-26_12-59-40';
p1_src = fullfile(p1_root, '10_23_2026', 'SS02344XJFRC28 (T4) - GABAR RNAi-14_31_92');
p2_src = 'C:\matlabroot\G4_Protocols\nested_RF_protocol2\data\2026_10_23_14_40';
jy_p2a = 'C:\matlabroot\G4_Protocols\nested_RF_protocol2\data\2026_10_28_14_14';
jy_p2b = 'D:\nested_RF_sessions\JY1028a\p2b';
tag = 'JY1023w';
scratch = fullfile(tempdir, ['run_session_v2_widen_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
mkdir(scratch); fprintf('Scratch folder: %s\n', scratch);
global NFAIL; NFAIL = 0;
here0 = pwd; restore_cd = onCleanup(@() cd(here0));          % generate_bar_flash_stimulus_xy cd's into the bar folder

p1 = fullfile(scratch, 'src_p1', '10_23_2026', 'p1_cell'); mkdir(p1);
copyfile(p1_src, p1); copyfile(fullfile(p1_root, 'Patterns'), fullfile(scratch, 'src_p1', 'Patterns'));
p2 = fullfile(scratch, 'src_p2'); mkdir(p2); copyfile(p2_src, p2);

%% D. the defaults (v1)
d5 = fullfile(scratch, 'fns_n5'); mkdir(d5);
rng(5, 'twister'); generate_bar_flash_pos_fns(d5, 5, 3);
rng(5, 'twister'); want = randperm(88); want(want == 1) = 89;
fm = dir(fullfile(d5, '*Rep1*.mat')); fm = fm(contains({fm.name}, '80ms'));
q = load(fullfile(fm(1).folder, fm(1).name)); v = q.pfnparam.val(2:2:end);
check(isequal(v, want) && q.pfnparam.frames == 89 && numel(dir(fullfile(d5, '*.mat'))) == 6, ...
    'D: n_flank 5: the 80 ms Rep1 function is randperm(88) with 89 for frame 1, 89 frames, as before 30 Oct');

%% G. run_session_v2 generates the stage-1 package with n_flank 7
gen_root = fullfile(scratch, 'generated');
p2w = fullfile(scratch, 'src_p2_120'); mkdir(p2w); copyfile(p2_src, p2w);     % the 23 Oct log with the generated FLASHES pattern
rfn = @(call, S, out_dir) widen_readout(call, p2w, tag, out_dir, S);
expG = fullfile(scratch, 'cell_G');
SG = run_session_v2(expG, 'dry_run', true, 'cell', tag, 'side', 'L', 'p1_folder', p1, 'p2_generate_in_dry_run', true, ...
    'p2_package_root', gen_root, 'readout_fn', rfn, 'overview_plots', false, 'steps', {'p1_peak', 'p2a', 'readout1', 'p2b'});
st = containers.Map({SG.steps.id}, {SG.steps.status});
check(all(cellfun(@(k) strcmp(st(k), 'done'), {'p1_peak', 'p2a', 'readout1', 'p2b'})), 'G: p1_peak, p2a, readout1, p2b done');
gen = SG.out.p2a.folder;
pf = dir(fullfile(gen, 'Patterns', '*FLASHES*.mat')); qp = load(fullfile(pf.folder, pf.name), 'pattern');
nfr = size(qp.pattern.Pats, 3);
fl = dir(fullfile(gen, 'Functions', '*bar_flashes*.mat')); okf = true;
for k = 1:numel(fl)
    q = load(fullfile(fl(k).folder, fl(k).name)); vals = q.pfnparam.val(2:2:end);
    okf = okf && isequal(sort(vals), 2:121) && q.pfnparam.frames == 121;
end
check(nfr == 121 && numel(fl) == 6 && okf && contains(fl(1).name, '120frames'), ...
    sprintf('G: FLASHES pattern %s: %d frames (8 x 15 bars + background); 6 bar-flash functions, each flashes bars 2..121 once', pf.name, nfr));
cg = load(fullfile(gen, 'currentExp.mat')); c23 = load(fullfile(p2_src, 'currentExp.mat'));
fid_g = str2double(pf.name(1:4)); f23 = dir(fullfile(p2_src, 'Patterns', '*FLASHES*.mat')); fid_23 = str2double(f23.name(1:4));
sw_g = cg.trial_dur(cg.pattern_order ~= fid_g); sw_23 = c23.trial_dur(c23.pattern_order ~= fid_23);
bf_g = cg.trial_dur(cg.pattern_order == fid_g); bf_23 = c23.trial_dur(c23.pattern_order == fid_23);
check(isequal(double(cg.pattern_order), double(c23.pattern_order)) && max(abs(double(sw_g) - double(sw_23))) < 1e-9 && fid_g == fid_23, ...
    sprintf('G: the sweep conditions and durations are the 23 Oct package''s (stage 1: %.1f s per repeat, unchanged)', sum(sw_g)));
t_full_g = 3 * sum(cg.trial_dur); t_full_23 = 3 * sum(c23.trial_dur);
fprintf('  INFO  bar-flash conditions per repeat: %s s (23 Oct, 88 bars: %s s); full-mode Protocol 2 x3: %.1f min (23 Oct: %.1f min, %+.1f min)\n', ...
    mat2str(round(bf_g, 2)), mat2str(round(bf_23, 2)), t_full_g / 60, t_full_23 / 60, (t_full_g - t_full_23) / 60);

%% S. read-out call 1 with the 120-frame pattern; the stage-2 package
B2 = readtable(fullfile(expG, 'readout', [tag '_stage2_bars.csv']), 'TextType', 'string', 'Delimiter', ',');
sets = unique(B2.set); okset = numel(sets) == 2;
for k = 1:numel(sets)
    fr = sort(B2.bar_pattern_frame_1based(B2.set == sets(k))).';
    okset = okset && numel(fr) == 15 && mod(fr(1) - 1, 15) == 0 && isequal(fr, fr(1):fr(1) + 14) && numel(unique(B2.bar_orientation_deg(B2.set == sets(k)))) == 1;
end
J1 = jsondecode(fileread(fullfile(expG, 'readout', [tag '_stage1.json'])));
check(height(B2) == 30 && okset && J1.stage2.n_bars == 30, sprintf('S: read-out call 1: %d stage-2 bars, 15 per set (one orientation, consecutive bars); %.1f s by the read-out''s count', height(B2), J1.stage2.total_s));
p2b = SG.out.p2b; cfg = jsondecode(fileread(fullfile(expG, 'p2b', 'Reports', 'config.json')));
rw = jsondecode(fileread(fullfile(expG, 'p2b', 'Reports', 'rows.json'))).rows;
t2 = sum([rw.duration_s]);
dp = dir(fullfile(expG, 'p2b', 'Patterns', '0001_*.mat')); qq = load(fullfile(dp.folder, dp.name), 'pattern');
check(p2b.n_bars == 30 && p2b.n_flashes == 180 && cfg.n_bars == 30 && size(qq.pattern.Pats, 3) == 121 && strcmp(SG.out.readout1.stage2_mode, 'bars22'), ...
    sprintf('S: stage-2 package: %d bars x 2 classes x %d repeats = %d flashes on the 121-frame pattern, %.1f s (%.2f min)', p2b.n_bars, p2b.n_reps, p2b.n_flashes, t2, t2 / 60));
t2_22 = 2 + 3 * (2 * 3 + 22 * (1.0 + 0.5)) + 2;                  % the same package with 22 bars (v1): settle 3 s per function
fprintf('  INFO  stage 2: %.1f s with 30 bars vs %.1f s with 22 (%+.1f s, %+.2f min)\n', t2, t2_22, t2 - t2_22, (t2 - t2_22) / 60);

close all force
fprintf('\n%d check(s) failed. Output in %s\n', NFAIL, scratch);

function widen_readout(call, p2w, tag, out_dir, S)
% the 23 Oct log, with the FLASHES pattern of the package that run_session_v2 generated (120 frames), for call 1;
% call 2 is not used (a 120-bar stage-2 log only comes from the arena)
    if call ~= 1, error('test:widen', 'call 2 is not emulated'); end
    old = dir(fullfile(p2w, 'Patterns', '*FLASHES*.mat')); new = dir(fullfile(S.out.p2a.folder, 'Patterns', '*FLASHES*.mat'));
    copyfile(fullfile(new.folder, new.name), fullfile(old.folder, old.name));      % the recorded name, the generated content
    rapid_readout(p2w, 'Cell', tag, 'Out', out_dir, 'Figure', false, 'Stage1Test', true, 'P1Peak', double(S.out.p1_peak.xy(:).'));
end

function check(cond, msg)
    global NFAIL
    if cond, fprintf('  PASS  %s\n', msg); else, fprintf(2, '  FAIL  %s\n', msg); NFAIL = NFAIL + 1; end
end

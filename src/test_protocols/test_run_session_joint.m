% TEST_RUN_SESSION_JOINT  run_session with the real MATLAB read-outs on a local TUTL session (26 Oct 2026).
%
% Joint dry run of the sequencer (mock controller) with rapid_readout.m and rapid_p1_peak.m
% (src/analysis/protocol2) on a complete Protocol 2 session kept on this Mac
% (~/Documents/ttl_1DRF/2025_10_23_10_31: historic layout, three sweep speeds, 88 bars). The
% session folder is read-only reference data: the sequencer's Protocol 2 package is a copy of its
% currentExp.mat, Patterns/ and Functions/ (no logs), and every read-out call writes into the scratch
% session folder. Call 1 and call 2 of the two-stage form use the read-out's own emulations
% ('Stage1Test': bar-flash trials dropped; 'Stage2Test': only the 22 placement-set bars kept), since a
% real stage-2 log can only come from the arena.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(genpath(fullfile(repo, 'src')));
analysis_repo = fullfile(fileparts(repo), 'ephys-standard-analysis');
session = fullfile(getenv('HOME'), 'Documents', 'ttl_1DRF', '2025_10_23_10_31');
assert(isfolder(fullfile(session, 'Log Files')), 'TUTL session not found: %s', session);
assert(exist('rapid_readout', 'file') == 2 && exist('rapid_p1_peak', 'file') == 2, 'rapid_readout / rapid_p1_peak not on the path');
tag = 'TUTL1023a';

scratch = fullfile(tempdir, ['run_session_joint_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
mkdir(scratch);
fprintf('Scratch folder: %s\n', scratch);
global NFAIL; NFAIL = 0;

% the Protocol 2 package without its logs (the sequencer's mock would otherwise move the session's Log Files aside)
p2pkg = fullfile(scratch, 'p2_pkg'); mkdir(p2pkg);
copyfile(fullfile(session, 'currentExp.mat'), fullfile(p2pkg, 'currentExp.mat'));
copyfile(fullfile(session, 'Patterns'), fullfile(p2pkg, 'Patterns'));
copyfile(fullfile(session, 'Functions'), fullfile(p2pkg, 'Functions'));
c = load(fullfile(p2pkg, 'currentExp.mat'));
bar_flash_id = max(c.pattern_order);              % this historic package has ONE bar-flash pattern (11), used at both speeds
p1_fn = @(f, o) rapid_p1_peak(session, 'Out', o, 'Cell', tag);

%% A. Full mode: Protocol 1 peak from rapid_p1_peak, read-out by rapid_readout (one call, P1Peak passed)
expA = fullfile(scratch, 'cell_A_full');
SA = run_session(expA, 'dry_run', true, 'cell', tag, 'side', 'L', 'p2_mode', 'full', 'p1_folder', session, 'p1_readout_fn', p1_fn, ...
    'p2_folder', session, 'steps', {'p1_peak', 'p2', 'readout', 'build', 'p3b', 'p4'}, 'analysis_repo', analysis_repo);
stA = containers.Map({SA.steps.id}, {SA.steps.status});
check(strcmp(stA('p1_peak'), 'done') && strcmp(stA('readout'), 'done') && strcmp(stA('build'), 'done') && strcmp(stA('p3b'), 'done') && strcmp(stA('p4'), 'done'), 'full mode with the real read-outs: all selected steps done');
pk = SA.out.p1_peak;
check(pk.frame == 303 && isequal(pk.xy(:).', [81.5 23.5]) && abs(pk.peak_mV - 14.77) < 0.05 && strcmp(pk.source, 'rapid_p1_peak'), ...
    sprintf('rapid_p1_peak: ON peak frame %d at col/row %s, %.2f mV', pk.frame, mat2str(pk.xy(:).'), pk.peak_mV));
rs = SA.out.readout.summary;
check(strcmp(SA.out.readout.method, 'rapid_readout.m') && abs(rs.center_col_row(1) - 79.1) < 0.1 && abs(rs.center_col_row(2) - 25.0) < 0.1, ...
    sprintf('rapid_readout one-call: center col/row %s', mat2str(rs.center_col_row(:).', 4)));
L = jsondecode(fileread(fullfile(expA, 'readout', [tag '_locations.json'])));
check(isfield(L, 'sequencer_inputs') && L.sequencer_inputs.usable_pd && isequal(L.other_centers.protocol1_ON_peak(:).', [81.5 23.5]), ...
    'locations.json: sequencer_inputs present, usable_pd true, the Protocol 1 peak passed through P1Peak');
ck = SA.out.readout.checks;
check(ck.near_p1 && abs(ck.dist_to_p1_px - norm([79.1 25.0] - [81.5 23.5])) < 0.15 && ck.inside_bar_set && ck.peak_above_floor && ischar(ck.unimodal), ...
    sprintf('center checks: %.2f px from the Protocol 1 peak, placement offset inside the bar set, center peak %.1f mV, unimodal n/a', ck.dist_to_p1_px, rs.center_peak_mV));
dA = SA.decisions(strcmp({SA.decisions.decision}, 'center'));
check(~dA.fallback && isequal(SA.center_shift_px, [0 0]), 'center accepted, no shift');
F4 = readtable(fullfile(expA, 'readout', [tag '_prot4_frames.csv']), 'TextType', 'string', 'Delimiter', ',');
% 399 frames and 13.3 min with the hand-over's lattice window (s -12..+16); the window is s -14..+18 since 29 Oct 2026
check(SA.out.build.p3b.n_frames == 26 && SA.out.build.p3b.n_blocks == 26 && SA.out.build.p4.n_frames == height(F4) && SA.out.build.p4.n_frames > 399 && ...
    SA.out.build.p4.programmed_s/60 > 13.3 && SA.out.build.p4.programmed_s/60 < 15.5, ...
    sprintf('packages from the real tables: Prot3b 26 frames / 26 blocks; Protocol 4 %d frames, %.1f min', SA.out.build.p4.n_frames, SA.out.build.p4.programmed_s/60));

%% B. Two-stage: stage 1 rows from the package, call 1 (Stage1Test), stage-2 package from the real sequence, call 2 (Stage2Test)
expB = fullfile(scratch, 'cell_B_two_stage');
p1xy = pk.xy(:).';
rfn = @(call, S, out_dir) joint_readout(call, session, tag, out_dir, p1xy);
SB = run_session(expB, 'dry_run', true, 'cell', tag, 'side', 'L', 'p1_folder', session, 'p1_readout_fn', p1_fn, 'p2_folder', p2pkg, ...
    'p2_bar_flash_pattern_id', bar_flash_id, 'readout_fn', rfn, 'steps', {'p1_peak', 'p2a', 'readout1', 'p2b', 'readout', 'build'}, 'analysis_repo', analysis_repo);
stB = containers.Map({SB.steps.id}, {SB.steps.status});
check(all(cellfun(@(k) strcmp(stB(k), 'done'), {'p1_peak', 'p2a', 'readout1', 'p2b', 'readout', 'build'})), 'two-stage with the real read-out: all selected steps done');
sa = jsondecode(fileread(fullfile(p2pkg, 'Reports', 'run_status.json')));
n_sweep_cond = sum(c.pattern_order ~= bar_flash_id);
check(sa.rows_total == 3 * n_sweep_cond && all([sa.rows.pattern_id] ~= bar_flash_id) && strcmp(sa.state, 'completed'), ...
    sprintf('stage 1: %d non-bar-flash conditions x 3 reps = %d rows (bar-flash pattern %d excluded), deciseconds as run_protocol2', n_sweep_cond, sa.rows_total, bar_flash_id));
J1 = jsondecode(fileread(fullfile(expB, 'readout', [tag '_stage1.json'])));
check(J1.stage == 1 && J1.sequencer_inputs.usable_pd && SB.out.readout1.usable_pd && strcmp(SB.out.readout1.stage2_mode, 'bars22'), ...
    'read-out call 1 (stage-1 emulation): stage1.json with usable_pd true -> stage 2 = 22 bars');
Sq = readtable(fullfile(expB, 'readout', [tag '_stage2_sequence.csv']), 'TextType', 'string', 'Delimiter', ',');
check(height(Sq) == 132 && numel(unique(Sq.bar_pattern_frame_1based)) == 22, 'stage-2 sequence from the read-out: 3 repeats x 2 classes x 22 bars');
sb = jsondecode(fileread(fullfile(expB, 'p2b', 'Reports', 'run_status.json')));
dp = dir(fullfile(expB, 'p2b', 'Patterns', '0001_*.mat'));
check(sb.rows_total == 8 && SB.out.p2b.n_flashes == 132 && isscalar(dp) && contains(dp.name, 'FLASHES'), ...
    sprintf('stage-2 package: pattern %d of the package copied as 0001 (%s), 8 rows, 132 flashes', bar_flash_id, dp.name));
rs2 = SB.out.readout.summary;
check(abs(rs2.center_col_row(1) - rs.center_col_row(1)) < 1e-6 && abs(rs2.center_col_row(2) - rs.center_col_row(2)) < 1e-6, ...
    'read-out call 2 (stage-2 emulation, 22 bars): the same center as the one-call form');
dB = SB.decisions(strcmp({SB.decisions.decision}, 'center'));
check(~dB.fallback && isequal(SB.center_shift_px, [0 0]) && SB.out.build.p3b.n_blocks == 26, 'two-stage: center accepted, Prot3b built');
check(~isfolder(fullfile(session, 'rapid_readout')) && isfolder(fullfile(session, 'Log Files')) && numel(dir(fullfile(session, 'Log Files'))) == 4, ...
    'the reference session folder was not written to');

fprintf('\n%d check(s) failed. Output in %s\n', NFAIL, scratch);

function joint_readout(call, session, tag, out_dir, p1xy)
    if call == 1
        rapid_readout(session, 'Cell', tag, 'Out', out_dir, 'Figure', false, 'Stage1Test', true, 'P1Peak', p1xy);
    else
        rapid_readout(session, 'Cell', tag, 'Out', out_dir, 'Figure', false, 'Stage2Test', true, 'P1Peak', p1xy);
    end
end

function check(cond, msg)
    global NFAIL
    if cond, fprintf('  PASS  %s\n', msg); else, fprintf(2, '  FAIL  %s\n', msg); NFAIL = NFAIL + 1; end
end

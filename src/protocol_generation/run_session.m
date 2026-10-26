function S = run_session(exp_folder, varargin)
%RUN_SESSION  One cell's protocol sequence on the G4 arena, automatic from one command (prototype, 26 Oct 2026).
%
%   S = RUN_SESSION(EXP_FOLDER)                       start, or RESUME, the session kept in EXP_FOLDER
%   S = RUN_SESSION(EXP_FOLDER, 'name', value, ...)   options below
%   S = RUN_SESSION(EXP_FOLDER, 'abandon', true)      the cell is lost: close the record (rest marked skipped)
%
%   Operator inputs: starting the session (break-in done, 0 pA), one key after the 30 s dark period
%   (the holding current is set once, by hand, and never changed by the sequencer) and Ctrl-C to
%   abort. Every other decision is rule-based with a defined fallback, is printed as a flag when it
%   fires and is logged with its inputs in EXP_FOLDER/decisions.json, so the analysis can re-make it
%   offline. There is no approval prompt. One override path exists and is never prompted for: the
%   option 'center' ([col row] or 'p1_peak') replaces the read-out center and is logged as an override.
%
%   The sequence (docs/briefs/JY_PROTOCOL_HANDOVER_20261026.md and M. Reiser's decisions of 26 Oct):
%     dark30   30 s in the dark at 0 pA right after break-in, logged (the resting Vm is on disk)
%     p1       Protocol 1 (the pre-made grids; rows driven from the .g4p, or 'p1_mode','conductor')
%     p1_peak  the Protocol 1 ON peak, computed by the Protocol 1 read-out (rapid_p1_peak.m), not by eye;
%              fallback: the 12-px grid; below the floor the session stops with a flag
%     p2a      Protocol 2 stage 1: the five sweep blocks, 16 directions, 3 repeats (about 7.3 min)
%     readout1 read-out call 1: placement PD (rule: mean of the 28 and 56 deg/s vectors within 10 deg,
%              otherwise the 28 deg/s vector with a warning; > 5 deg sets a flag), the two bar sets, the
%              stage-2 frame list and sequence; fallback when no usable PD: stage 2 plays all 88 bars
%              and the sites are placed from the Protocol 1 peak
%     p2b      Protocol 2 stage 2: the 22 bars at 80 ms and 14 ms, random order per repeat, 3 repeats
%              (about 1.7 min; 'n_reps_stage2', 5 as the option), a position function on the existing
%              88-frame bar pattern (no new pattern)
%     readout  read-out call 2: M6 center, window, 11 sites, Prot3b and Protocol 4 tables, the figure;
%              sanity checks (center within a few px of the Protocol 1 peak and inside the bar set, peak
%              response above a floor, unimodal profiles when reported); on failure the center falls back
%              to the Protocol 1 peak and the tables are shifted there, flagged
%     build    Prot3b and Protocol 4 packages from the tables
%     p3b      Prot3b: 26 blocks, five AO levels on AO2 (7.5 min)
%     p4       Protocol 4: one randomized sequence with mixed flash periods (about 14 min)
%   With 'p2_mode','full' Protocol 2 runs as today (all 88 bars in one run, run_protocol2) and the
%   read-out is called once at its end; the steps are then dark30, p1, p1_peak, p2, readout, build, p3b, p4.
%
%   State: EXP_FOLDER/session_state.json after every step and row, session_log.txt, decisions.json.
%   RUN_SESSION(EXP_FOLDER) after a crash or Ctrl-C continues from the first step that is not done; a step
%   left running, failed or interrupted is re-run (its old logs are moved aside). Each protocol runs in
%   its own package folder under EXP_FOLDER (p1/, p2a/ or p2/, p2b/, p3b/, p4/) with its own 'Log Files'
%   and Reports/ ledgers, so the analysis sees the same per-protocol layout as before.
%
%   OPTIONS (name/value or one struct):
%     dry_run        false   SessionMockController instead of PanelsController; no G4 Host, no arena, no
%                            Protocol 2 generation (needs the rig paths); files, ledgers and decisions are
%                            written and the controller calls are logged
%     cell, side, age, strain   metadata (cell tag such as 'JY1027'; side 'L'|'R')
%     p2_mode        'two_stage' | 'full'
%     steps          all     cellstr of step ids to run; others are marked 'skipped (not selected)'
%     p1_mode        'sequencer' | 'conductor' | 'skip';  p1_package  the pre-made Protocol 1 package folder
%     p1_folder      an existing Protocol 1 results folder (resume / conductor mode)
%     p1_readout_fn  @(p1_folder, out_dir) -> struct(frame, xy, peak_mV, polarity, frame12, peak12_mV, ...)
%                    default rapid_p1_peak when on the path
%     peak_frame     Protocol 1 ON peak frame (6-px grid) given by hand: an OVERRIDE, logged
%     p2_folder      an existing Protocol 2 package (stage 1 or full), p2b_folder an existing stage-2 package
%     readout_fn     @(call, S, out_dir) producing the read-out files (call 1: <cell>_stage1.json,
%                    <cell>_stage2_sequence.csv; call 2: <cell>_prot3b_frames.csv, _prot3b_blocks.csv,
%                    _prot4_frames.csv, _prot4_sequence.csv, _sites.csv, _locations.json, _locations.png);
%                    default rapid_readout(p2_folder, 'Cell', cell, 'Out', out_dir[, 'P1Peak', xy][, 'Stage2', p2b_folder]);
%                    the JSON object sequencer_inputs (usable_pd, centroid_offset_placement_px,
%                    center_expected_ON_80ms_peak_mV, profiles_unimodal) feeds the decisions when present
%     readout_dir    folder that already holds those files (dry run, resume)
%     center         [] | [col row] | 'p1_peak'   the override (never prompted for)
%     n_reps_stage2  3;  p2_bar_flash_pattern_id  [] = max(pattern_order)-1 (run_protocol2's slow bar-flash pattern)
%     thresholds: p1_peak_floor_mV 2, pd_vector_length_floor 0.1, center_max_dist_px 6, center_max_offset_px 4,
%                 peak_floor_mV 3
%     ao_channel 2, ai_channels [1 2], adc_rates_hz [0 10000 10000 0], host_config_file, pretrial_s 2,
%     ao_zero_tail_s 0.1, dark_s 30, controller_factory, sim_speed 0, show_figures (= ~dry_run), auto_accept (= dry_run)
%
%   Conventions: arena columns 17..192 and rows 1..48, 1-based, row 1 at the top; background 4, ON 15,
%   OFF 0; s is positive toward the PD (the trailing / proximal side), v is + 90 deg clockwise from the
%   PD; sides are leading (distal) / trailing (proximal).
%
%   Reuses the run_protocol2 call pattern and its generators (nested_RF_stimulus), the per-row controller
%   sequence and cleanup guard of run_protocol3 (TUTL p4_experiment_design, Oct 2026), the G4 file layout
%   of save_function_G4, and save_pattern_G4 / maDisplayTools for patterns. No existing runner is modified.
%
%   See also RUN_PROTOCOL2, GENERATE_PROTOCOL2, RAPID_READOUT, SESSIONMOCKCONTROLLER.

opts = parse_options(varargin{:});
exp_folder = char(exp_folder);
if ~isfolder(exp_folder), mkdir(exp_folder); end
state_file = fullfile(exp_folder, 'session_state.json');

% ---------------------------------------------------------------- state
if isfile(state_file)
    S = load_state(state_file);
    S.exp_folder = exp_folder;
    S = merge_metadata(S, opts);
    say(S, 'Resuming session %s (%s)', S.cell, exp_folder);
else
    S = new_state(exp_folder, opts);
    say(S, 'New session %s in %s%s  [p2_mode %s]', S.cell, exp_folder, tern(opts.dry_run, ' [DRY RUN]', ''), opts.p2_mode);
end
S.opts = opts_for_state(opts);
save_state(S, state_file);

if opts.abandon
    for k = 1:numel(S.steps)
        if any(strcmp(S.steps(k).status, {'pending', 'running', 'interrupted', 'failed'}))
            S.steps(k).status = 'skipped'; S.steps(k).note = 'cell lost (abandon)';
        end
    end
    S.finished = iso_now();
    save_state(S, state_file);
    say(S, 'Session abandoned by the operator.');
    print_checklist(S); print_decisions(S);
    return
end

print_checklist(S);

% ---------------------------------------------------------------- run the steps in order
for k = 1:numel(S.steps)
    st = S.steps(k);
    if ~ismember(st.id, opts.steps)
        if strcmp(st.status, 'pending')
            S.steps(k).status = 'skipped'; S.steps(k).note = 'not selected';
            save_state(S, state_file);
        end
        continue
    end
    if any(strcmp(st.status, {'done', 'skipped'})), continue; end
    if any(strcmp(st.status, {'running', 'interrupted', 'failed'}))
        say(S, 'FLAG step %s was left "%s" (%s); re-running it.', st.id, st.status, st.note);
        S = decide(S, st.id, 'resume', 'a step left running/failed/interrupted is re-run; its old logs are moved aside', ...
            struct('previous_status', st.status, 'note', st.note), 'rerun', false, 'rerun');
    end

    S.steps(k).status = 'running'; S.steps(k).t_start = iso_now(); S.steps(k).note = '';
    save_state(S, state_file);
    say(S, '---- %d/%d  %s  (%s, about %.1f min) ----', k, numel(S.steps), st.id, st.label, st.est_min);
    t0 = tic;
    try
        switch st.id
            case 'dark30',   S = step_dark30(S, opts);
            case 'p1',       S = step_p1(S, opts);
            case 'p1_peak',  S = step_p1_peak(S, opts);
            case 'p2',       S = step_p2_full(S, opts);
            case 'p2a',      S = step_p2_stage1(S, opts);
            case 'readout1', S = step_readout1(S, opts);
            case 'p2b',      S = step_p2_stage2(S, opts);
            case 'readout',  S = step_readout(S, opts);
            case 'build',    S = step_build(S, opts);
            case 'p3b',      S = step_run_package(S, opts, 'p3b');
            case 'p4',       S = step_run_package(S, opts, 'p4');
            otherwise, error('run_session:step', 'Unknown step %s', st.id);
        end
    catch err
        S.steps(k).elapsed_s = toc(t0);
        if strcmp(err.identifier, 'run_session:skip')
            S.steps(k).status = 'skipped'; S.steps(k).note = err.message;
            save_state(S, state_file);
            say(S, 'Step %s skipped: %s', st.id, err.message);
            continue
        elseif strcmp(err.identifier, 'run_session:stop')
            S.steps(k).status = 'interrupted'; S.steps(k).note = err.message;
            save_state(S, state_file);
            say(S, 'FLAG session stopped at %s: %s', st.id, err.message);
            print_checklist(S); print_decisions(S);
            return
        end
        S.steps(k).status = 'failed'; S.steps(k).note = sprintf('%s: %s', err.identifier, err.message);
        save_state(S, state_file);
        say(S, 'Step %s FAILED: %s', st.id, err.message);
        print_checklist(S);
        rethrow(err);
    end
    S.steps(k).elapsed_s = toc(t0);
    S.steps(k).t_end = iso_now();
    if strcmp(S.steps(k).status, 'running'), S.steps(k).status = 'done'; end
    save_state(S, state_file);
    say(S, 'Step %s done in %.1f min.', st.id, S.steps(k).elapsed_s/60);
end

S.finished = iso_now();
save_state(S, state_file);
print_checklist(S); print_decisions(S);
say(S, 'Session complete. State: %s', state_file);
end

% =========================================================================================
% options and state
% =========================================================================================
function o = parse_options(varargin)
d = struct('dry_run', false, 'abandon', false, 'cell', '', 'side', '', 'age', '', 'strain', '', ...
    'p2_mode', 'two_stage', 'steps', {{}}, 'p1_mode', 'sequencer', 'p1_package', '', 'p1_folder', '', ...
    'p1_readout_fn', [], 'peak_frame', [], 'p2_folder', '', 'p2b_folder', '', 'readout_fn', [], 'readout_dir', '', ...
    'analysis_repo', '', 'center', [], 'n_reps_stage2', 3, 'n_reps_p2', 3, 'p2_bar_flash_pattern_id', [], ...
    'p1_peak_floor_mV', 2, 'pd_vector_length_floor', 0.1, 'center_max_dist_px', 6, 'center_max_offset_px', 4, ...
    'peak_floor_mV', 3, 'auto_accept', [], 'show_figures', [], 'ao_channel', 2, 'ai_channels', [1 2], ...
    'adc_rates_hz', [0 10000 10000 0], 'host_config_file', '', 'pretrial_s', 2, 'ao_zero_tail_s', 0.1, 'dark_s', 30, ...
    'controller_factory', [], 'sim_speed', 0, 'verbose_mock', false);
if numel(varargin) == 1 && isstruct(varargin{1})
    given = varargin{1};
else
    if mod(numel(varargin), 2) ~= 0, error('run_session:options', 'Options come in name/value pairs'); end
    given = struct();
    for k = 1:2:numel(varargin), given.(char(varargin{k})) = varargin{k+1}; end
end
o = d;
for f = fieldnames(given).'
    if ~isfield(d, f{1}), error('run_session:options', 'Unknown option %s', f{1}); end
    o.(f{1}) = given.(f{1});
end
if ~ismember(o.p2_mode, {'two_stage', 'full'}), error('run_session:options', 'p2_mode must be two_stage or full'); end
if isempty(o.auto_accept), o.auto_accept = o.dry_run; end
if isempty(o.show_figures), o.show_figures = ~o.dry_run; end
if isempty(o.steps), o.steps = step_ids(o.p2_mode); end
o.steps = cellstr(o.steps);
if isstring(o.cell), o.cell = char(o.cell); end
if isstring(o.side), o.side = char(o.side); end
if isstring(o.center), o.center = char(o.center); end
if isempty(o.analysis_repo)
    o.analysis_repo = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))), 'ephys-standard-analysis');
end
if isempty(o.controller_factory)
    if o.dry_run
        mo = struct('sim_speed', o.sim_speed, 'verbose', o.verbose_mock);
        o.controller_factory = @() SessionMockController(mo);
    else
        o.controller_factory = @() PanelsController();
    end
end
end

function ids = step_ids(mode)
if nargin < 1, mode = 'two_stage'; end
if strcmp(mode, 'full')
    ids = {'dark30', 'p1', 'p1_peak', 'p2', 'readout', 'build', 'p3b', 'p4'};
else
    ids = {'dark30', 'p1', 'p1_peak', 'p2a', 'readout1', 'p2b', 'readout', 'build', 'p3b', 'p4'};
end
end

function [labels, est] = step_labels(ids)
L = struct('dark30', {{'30 s dark at 0 pA (break-in)', 0.5}}, 'p1', {{'Protocol 1 (coarse map)', 4.5}}, ...
    'p1_peak', {{'Protocol 1 ON peak (read-out)', 0.2}}, 'p2', {{'Protocol 2 (all 88 bars, as today)', 14.6}}, ...
    'p2a', {{'Protocol 2 stage 1 (sweeps)', 7.3}}, 'readout1', {{'Read-out 1 (placement PD, bar sets)', 0.3}}, ...
    'p2b', {{'Protocol 2 stage 2 (22 bars)', 1.7}}, 'readout', {{'Read-out 2 (center, sites, tables)', 0.5}}, ...
    'build', {{'Build Prot3b + Protocol 4', 0.5}}, 'p3b', {{'Prot3b (26 blocks, 5 AO levels)', 7.5}}, ...
    'p4', {{'Protocol 4 (combined map)', 14}});
labels = cellfun(@(id) L.(id){1}, ids, 'UniformOutput', false);
est = cellfun(@(id) L.(id){2}, ids);
end

function S = new_state(exp_folder, opts)
ids = step_ids(opts.p2_mode);
[labels, est] = step_labels(ids);
steps = struct('id', ids, 'label', labels, 'status', 'pending', 't_start', '', 't_end', '', ...
    'elapsed_s', num2cell(zeros(size(ids))), 'est_min', num2cell(est), 'note', '');
S = struct();
S.schema = 'run_session/2';
S.created = iso_now();
S.finished = '';
S.exp_folder = exp_folder;
S.dry_run = opts.dry_run;
S.p2_mode = opts.p2_mode;
S.cell = opts.cell; S.side = opts.side; S.age = opts.age; S.strain = opts.strain;
if isempty(S.cell), S.cell = 'cell'; end
if isempty(S.side), S.side = 'L'; end
S.steps = steps;
S.center_shift_px = [0 0];
S.decisions = empty_decisions();
S.out = struct();
for f = [{'setup'}, ids], S.out.(f{1}) = struct(); end
if ~isempty(opts.p1_folder), S.out.p1.folder = opts.p1_folder; end
if ~isempty(opts.p2_folder)
    if strcmp(opts.p2_mode, 'full'), S.out.p2.folder = opts.p2_folder; else, S.out.p2a.folder = opts.p2_folder; end
end
if ~isempty(opts.p2b_folder) && isfield(S.out, 'p2b'), S.out.p2b.folder = opts.p2b_folder; end
end

function S = merge_metadata(S, opts)
for f = {'cell', 'side', 'age', 'strain'}
    if ~isempty(opts.(f{1})), S.(f{1}) = opts.(f{1}); end
end
if ~isempty(opts.p1_folder), S.out.p1.folder = opts.p1_folder; end
if ~isempty(opts.p2_folder)
    if strcmp(S.p2_mode, 'full'), S.out.p2.folder = opts.p2_folder; else, S.out.p2a.folder = opts.p2_folder; end
end
if ~isempty(opts.p2b_folder) && isfield(S.out, 'p2b'), S.out.p2b.folder = opts.p2b_folder; end
end

function o = opts_for_state(opts)
o = rmfield(opts, {'controller_factory', 'readout_fn', 'p1_readout_fn'});
end

function S = load_state(f)
S = jsondecode(fileread(f));
if iscell(S.steps), S.steps = [S.steps{:}]; end
for k = 1:numel(S.steps)
    for g = {'note', 't_start', 't_end'}
        if isempty(S.steps(k).(g{1})), S.steps(k).(g{1}) = ''; end
    end
end
if ~isfield(S, 'p2_mode'), S.p2_mode = 'two_stage'; end
if ~isfield(S, 'out') || isempty(S.out), S.out = struct(); end
for f = [{'setup'}, step_ids(S.p2_mode)]
    if ~isfield(S.out, f{1}) || isempty(S.out.(f{1})), S.out.(f{1}) = struct(); end
end
if ~isfield(S, 'decisions') || isempty(S.decisions)
    S.decisions = empty_decisions();
elseif iscell(S.decisions)
    S.decisions = [S.decisions{:}];
end
if ~isfield(S, 'center_shift_px') || isempty(S.center_shift_px), S.center_shift_px = [0 0]; end
S.center_shift_px = S.center_shift_px(:).';
end

function save_state(S, f)
write_text(f, [jsonencode(S, 'PrettyPrint', true) newline]);
end

% =========================================================================================
% the decision log: every automatic decision with its inputs, rule, outcome and fallback flag
% =========================================================================================
function d = empty_decisions()
d = struct('time', {}, 'step', {}, 'decision', {}, 'rule', {}, 'inputs', {}, 'outcome', {}, 'fallback', {}, 'flags', {});
end

function S = decide(S, step, decision, rule, inputs, outcome, fallback, flags)
if isstruct(inputs), inputs = jsonencode(inputs); end
d = struct('time', iso_now(), 'step', step, 'decision', decision, 'rule', rule, 'inputs', char(inputs), ...
    'outcome', char(outcome), 'fallback', logical(fallback), 'flags', char(flags));
if isempty(S.decisions), S.decisions = d; else, S.decisions(end+1) = d; end
write_json(fullfile(S.exp_folder, 'decisions.json'), S.decisions);
say(S, '%sDECISION %s/%s: %s%s', tern(fallback, 'FLAG ', ''), step, decision, outcome, tern(isempty(flags), '', [' [' char(flags) ']']));
end

function print_decisions(S)
if isempty(S.decisions), return; end
fprintf('\n  Automatic decisions (%d; details in decisions.json):\n', numel(S.decisions));
for k = 1:numel(S.decisions)
    d = S.decisions(k);
    fprintf('  %s %-8s %-16s %s%s\n', tern(d.fallback, '!', ' '), d.step, d.decision, d.outcome, tern(isempty(d.flags), '', ['  [' d.flags ']']));
end
fprintf('\n');
end

% =========================================================================================
% steps
% =========================================================================================
function S = step_dark30(S, opts)
% Launching run_session is the first operator input: break-in done, current clamp, 0 pA. The 30 s of
% dark are logged. The second (and last) input follows: the holding current is set by hand once.
if ~opts.dry_run
    S = apply_host_config_once(S, opts);
end
folder = fullfile(S.exp_folder, 'dark30');
if ~isfolder(folder), mkdir(folder); end
say(S, 'Logging %g s of dark at 0 pA (break-in). Ctrl-C aborts.', opts.dark_s);
ctlr = opts.controller_factory();
ctx = containers.Map('KeyType', 'char', 'ValueType', 'any');
ctx('ctlr') = ctlr; ctx('opened') = false; ctx('log') = false; ctx('done') = false; ctx('ao') = [];
guard = onCleanup(@() controller_cleanup(ctx)); %#ok<NASGU>
ctlr.open(true); ctx('opened') = true;
checked('setRootDirectory', ctlr.setRootDirectory(folder));
ctx('log') = true;
checked('startLog', ctlr.startLog());
checked('allOff', ctlr.allOff());
t_wait = tern(opts.dry_run, opts.dark_s * opts.sim_speed, opts.dark_s);
for s = t_wait:-5:0.01
    fprintf('    dark: %4.0f s left\r', s); pause(min(5, s));
end
fprintf('\n');
ctx('done') = true;
controller_cleanup(ctx);
S.out.dark30 = struct('folder', folder, 'seconds', opts.dark_s, 'logged', true);
say(S, 'Dark period logged in %s.', folder);
% The one physical hand-off the sequencer cannot do: the holding current. Set once, logged, never changed.
hold = ask(opts, 'Set the holding current now (it stays fixed for the session). Type it in pA and press Enter: ', 'not recorded');
S.out.dark30.holding_pA = hold;
say(S, 'Holding current recorded as "%s" pA. Protocol 1 starts.', hold);
end

function S = step_p1(S, opts)
switch opts.p1_mode
    case 'skip'
        error('run_session:skip', 'p1_mode = skip');
    case 'conductor'
        say(S, 'Protocol 1 in the G4 Conductor (pre-made .g4p for side %s); the sequencer continues when its results folder exists.', S.side);
        f = ask(opts, 'Protocol 1 results folder (Enter when the Conductor has finished): ', '');
        S.out.p1 = struct('mode', 'conductor', 'folder', f, 'package_source', opts.p1_package);
    case 'sequencer'
        if isempty(opts.p1_package) || ~isfolder(opts.p1_package)
            error('run_session:p1', 'p1_mode sequencer needs p1_package (the pre-made Protocol 1 package folder)');
        end
        pkg = fullfile(S.exp_folder, 'p1');
        if ~isfolder(fullfile(pkg, 'Patterns')), copy_package(opts.p1_package, pkg); end
        rows = rows_from_g4p(pkg);
        S.steps(strcmp({S.steps.id}, 'p1')).est_min = sum([rows.duration_s]) / 60;
        say(S, 'Protocol 1: %d rows, %.1f min programmed, from %s', numel(rows), sum([rows.duration_s])/60, opts.p1_package);
        status = run_rows(opts, pkg, rows, [], S);
        S.out.p1 = struct('mode', 'sequencer', 'folder', pkg, 'package_source', opts.p1_package, ...
            'rows', numel(rows), 'state', status.state, 'log_folder', fullfile(pkg, 'Log Files'));
    otherwise
        error('run_session:p1', 'p1_mode must be sequencer, conductor or skip');
end
end

function S = step_p1_peak(S, opts)
% Decision (a): aim Protocol 2 from the Protocol 1 ON peak computed by the read-out.
% Rule: 6-px grid ON peak when its response is >= p1_peak_floor_mV; fallback the 12-px grid peak; below
% the floor on both grids the session stops (nothing can be aimed) unless 'peak_frame' overrides.
out_dir = fullfile(S.exp_folder, 'readout'); if ~isfolder(out_dir), mkdir(out_dir); end
pkg = '';
if isfield(S.out.p1, 'package_source') && ~isempty(S.out.p1.package_source), pkg = S.out.p1.package_source; end
if isfield(S.out.p1, 'folder') && ~isempty(S.out.p1.folder) && isfolder(fullfile(S.out.p1.folder, 'Patterns')), pkg = S.out.p1.folder; end
R = struct('frame', [], 'xy', [], 'peak_mV', [], 'polarity', '', 'frame12', [], 'peak12_mV', [], 'source', '');
if ~isempty(opts.peak_frame)
    R.frame = opts.peak_frame; R.source = 'override (peak_frame option)';
else
    fn = opts.p1_readout_fn;
    if isempty(fn) && exist('rapid_p1_peak', 'file') == 2, fn = @(f, o) rapid_p1_peak(f, 'Out', o, 'Cell', S.cell); end
    if isempty(fn)
        error('run_session:stop', 'no Protocol 1 read-out on the path (rapid_p1_peak.m) and no peak_frame override');
    end
    p1_folder = ''; if isfield(S.out.p1, 'folder'), p1_folder = S.out.p1.folder; end
    Q = fn(p1_folder, out_dir);
    for g = fieldnames(Q).', if isfield(R, g{1}), R.(g{1}) = Q.(g{1}); end, end
    % rapid_p1_peak returns frame12 and peak12_mV as [ON OFF]; the fallback uses the ON 12-px peak
    if numel(R.frame12) > 1, R.frame12 = R.frame12(1); end
    if numel(R.peak12_mV) > 1, R.peak12_mV = R.peak12_mV(1); end
    R.source = 'rapid_p1_peak';
end
flags = '';
fallback = false;
if ~isempty(R.peak_mV) && R.peak_mV < opts.p1_peak_floor_mV
    if ~isempty(R.frame12) && ~isempty(R.peak12_mV) && R.peak12_mV >= opts.p1_peak_floor_mV
        R.frame = R.frame12; R.xy = []; fallback = true; flags = sprintf('6-px peak %.1f mV < floor %.1f; using the 12-px grid peak', R.peak_mV, opts.p1_peak_floor_mV);
    else
        S = decide(S, 'p1_peak', 'aim', 'ON peak >= floor on the 6-px grid, else the 12-px grid, else stop', R, 'no usable Protocol 1 response', true, 'stop');
        error('run_session:stop', 'Protocol 1 peak %.1f mV is below the floor %.1f mV on both grids; nothing to aim at', R.peak_mV, opts.p1_peak_floor_mV);
    end
end
if isempty(R.frame) || ~isfinite(R.frame) || R.frame ~= round(R.frame) || R.frame < 2
    error('run_session:stop', 'no valid Protocol 1 peak frame (got %s)', mat2str(R.frame));
end
if isempty(R.xy) && ~isempty(pkg)
    [R.xy, R.polarity] = frame_to_xy(pkg, R.frame);
end
S.out.p1_peak = R;
S = decide(S, 'p1_peak', 'aim', 'ON peak >= floor on the 6-px grid, else the 12-px grid, else stop', R, ...
    sprintf('frame %d at col/row %s (%s)', R.frame, mat2str(R.xy, 4), R.source), fallback, flags);
end

function S = step_p2_full(S, opts)
% Protocol 2 as today: all 88 bars in one run (run_protocol2), one read-out at the end.
if isfield(S.out.p2, 'folder') && ~isempty(S.out.p2.folder) && isfolder(S.out.p2.folder)
    say(S, 'Protocol 2 folder given: %s', S.out.p2.folder); S.out.p2.mode = 'given'; return
end
if opts.dry_run
    S.out.p2 = dry_p2_placeholder(S, 'p2', 'full: generate_protocol2_package + run_protocol2 (all 88 bars)');
    return
end
[folder, pattern_order, func_order, trial_dur] = generate_protocol2_package(S.out.p1_peak.frame, S.side, p2_metadata(S));
S.out.p2 = struct('folder', folder, 'mode', 'run_protocol2', 'peak_frame', S.out.p1_peak.frame, 'n_reps', opts.n_reps_p2);
save_state(S, fullfile(S.exp_folder, 'session_state.json'));
say(S, 'Protocol 2 generated in %s; running (%d reps, all 88 bars).', folder, opts.n_reps_p2);
run_protocol2(folder, pattern_order, func_order, trial_dur, opts.n_reps_p2);   % as today; converts TDMS at the end
end

function S = step_p2_stage1(S, opts)
% Stage 1: the package is generated as today, but only the sweep conditions are run (every condition
% whose pattern is not one of run_protocol2's two bar-flash patterns), n_reps_p2 repeats.
if isfield(S.out.p2a, 'folder') && ~isempty(S.out.p2a.folder) && isfolder(S.out.p2a.folder) && ...
        isfield(S.out.p2a, 'state') && strcmp(S.out.p2a.state, 'completed')
    say(S, 'Stage 1 already completed in %s', S.out.p2a.folder); return
end
if isfield(S.out.p2a, 'folder') && ~isempty(S.out.p2a.folder) && isfolder(fullfile(S.out.p2a.folder, 'Patterns'))
    folder = S.out.p2a.folder; given = true;
elseif opts.dry_run
    S.out.p2a = dry_p2_placeholder(S, 'p2a', 'stage 1: generate_protocol2_package, then the sweep rows x n_reps');
    return
else
    folder = generate_protocol2_package(S.out.p1_peak.frame, S.side, p2_metadata(S)); given = false;
end
P = p2_plan(folder, opts);
rows = P.stage1_rows;
S.steps(strcmp({S.steps.id}, 'p2a')).est_min = sum([rows.duration_s]) / 60;
say(S, 'Protocol 2 stage 1: %d sweep rows x %d reps = %d rows, %.1f min (bar-flash patterns %s left for stage 2).', ...
    numel(P.stage1_idx), P.n_reps, numel(rows), sum([rows.duration_s])/60, mat2str(P.bar_flash_patterns));
status = run_rows(opts, folder, rows, [], S);
S.out.p2a = struct('folder', folder, 'mode', tern(given, 'given package, stage-1 rows', 'generated, stage-1 rows'), ...
    'peak_frame', S.out.p1_peak.frame, 'n_reps', P.n_reps, 'rows', numel(rows), 'state', status.state, ...
    'bar_flash_patterns', P.bar_flash_patterns, 'log_folder', fullfile(folder, 'Log Files'));
end

function S = step_readout1(S, opts)
% Read-out call 1 on the stage-1 session: placement PD, the two bar sets, the stage-2 frame list and
% sequence. Decision (b): the placement-PD rule lives in the read-out; the sequencer displays its
% outcome and flag. Decision (d): no usable PD -> stage 2 plays all 88 bars and the sites are placed
% from the Protocol 1 peak.
out_dir = fullfile(S.exp_folder, 'readout'); if ~isfolder(out_dir), mkdir(out_dir); end
tag = S.cell;
method = call_readout(1, S, opts, out_dir);
js = fullfile(out_dir, [tag '_stage1.json']); seq = fullfile(out_dir, [tag '_stage2_sequence.csv']);
if ~isfile(js), error('run_session:readout', 'Read-out call 1 did not write %s', js); end
J = jsondecode(fileread(js));
P = getf(J, 'placement', struct());
pdinfo = struct('pd_deg', getf(J, 'pd_deg', []), 'pd_vector_28', getf(J, 'pd_vector_28', []), 'pd_vector_56', getf(J, 'pd_vector_56', []), ...
    'pd_placement', getf(J, 'pd_placement', []), 'pd_rule', getf(J, 'pd_rule', ''), 'pd_flag_gt5', getf(J, 'pd_flag_gt5', false), ...
    'agree', getf(P, 'agree', []), 'diff_28_56_deg', getf(P, 'diff_28_56_deg', []), 'vector_length_28', getf(P, 'vector_length_28', []), ...
    'vector_length_56', getf(P, 'vector_length_56', []), 'margin_to_lut_boundary_deg', getf(P, 'margin_to_lut_boundary_deg', []));
flags = {};
if islogical(pdinfo.pd_flag_gt5) && pdinfo.pd_flag_gt5, flags{end+1} = '28 and 56 deg/s vectors differ by > 5 deg'; end
if islogical(pdinfo.agree) && ~pdinfo.agree, flags{end+1} = 'WARNING > 10 deg: placement uses the 28 deg/s vector alone'; end
S = decide(S, 'readout1', 'placement_pd', 'mean of the 28 and 56 deg/s vector PDs within 10 deg, else the 28 deg/s vector (warning); > 5 deg sets a flag; the TUTL snapped PD stays the reported metric', ...
    pdinfo, sprintf('placement PD %s deg (rule %s; snapped %s deg; 28: %s deg L %s; 56: %s deg L %s)', num2str(pdinfo.pd_placement), pdinfo.pd_rule, ...
    num2str(pdinfo.pd_deg), num2str(pdinfo.pd_vector_28), num2str(pdinfo.vector_length_28), num2str(pdinfo.pd_vector_56), num2str(pdinfo.vector_length_56)), ...
    false, strjoin(flags, '; '));
% usable PD?
si = getf(J, 'sequencer_inputs', struct());
usable = getf(si, 'usable_pd', getf(J, 'usable_pd', []));
if ischar(usable) || isstring(usable), usable = []; end
rule = 'sequencer_inputs.usable_pd from the read-out when given (all 16 directions measured and normalized 28 deg/s vector length >= 0.1), else vector length (28 deg/s) >= pd_vector_length_floor and a stage-2 sequence written';
if isempty(usable)
    L28 = pdinfo.vector_length_28;
    usable = ~isempty(L28) && L28 >= opts.pd_vector_length_floor && isfile(seq);
end
usable = logical(usable);
S.out.readout1 = struct('dir', out_dir, 'method', method, 'stage1_json', js, 'stage2_sequence', tern(isfile(seq), seq, ''), ...
    'placement', pdinfo, 'usable_pd', usable, 'stage2_mode', tern(usable, 'bars22', 'all88'));
S = decide(S, 'readout1', 'usable_pd', rule, struct('vector_length_28', pdinfo.vector_length_28, 'floor', opts.pd_vector_length_floor, ...
    'sequence_file', isfile(seq), 'usable_pd_from_readout', getf(si, 'usable_pd', 'n/a'), 'pd_inputs', getf(si, 'pd_inputs', 'n/a')), ...
    tern(usable, 'stage 2 plays the 22 bars of the two placement sets', 'NO usable PD: stage 2 plays all 88 bars; sites will be placed from the Protocol 1 peak'), ...
    ~usable, tern(usable, '', 'no_pd_fallback'));
end

function S = step_p2_stage2(S, opts)
% Stage 2 package p2b/: the existing bar-flash pattern (copied from the stage-1 package) and position
% functions written from the stage-2 sequence (22 bars at 80 ms + 920 ms and 14 ms + 486 ms, random order
% per repeat; one function per repeat and class, named like the Protocol 2 generator's so the read-out
% classifies the trials). Fallback: all 88 bars with the stage-1 package's own bar-flash rows, as run_protocol2 plays them.
if isfield(S.out.p2b, 'folder') && ~isempty(S.out.p2b.folder) && isfield(S.out.p2b, 'state') && strcmp(S.out.p2b.state, 'completed')
    say(S, 'Stage 2 already completed in %s', S.out.p2b.folder); return
end
p2a = S.out.p2a.folder;
out_dir = fullfile(S.exp_folder, 'p2b');
mode = S.out.readout1.stage2_mode;
if opts.dry_run && ~isfolder(fullfile(p2a, 'Patterns'))
    S.out.p2b = dry_p2_placeholder(S, 'p2b', sprintf('stage 2 (%s) needs the stage-1 package', mode)); S.out.p2b.mode = mode;
    return
end
P = p2_plan(p2a, opts);
if strcmp(mode, 'bars22')
    B = build_stage2_package(p2a, S.out.readout1.stage2_sequence, out_dir, P, opts, S.cell);
    rows = rows_from_package(out_dir);
    say(S, 'Stage 2: %d bars x %d repeats (%d flashes), %.1f min, pattern %d of the stage-1 package.', ...
        B.n_bars, B.n_reps, B.n_flashes, B.programmed_s/60, B.pattern_id);
else
    copy_package(p2a, out_dir);                 % all 88 bars: the package's own bar-flash rows and IDs
    rows = P.stage2_rows_all88;
    if ~isfolder(fullfile(out_dir, 'Reports')), mkdir(fullfile(out_dir, 'Reports')); end
    write_json(fullfile(out_dir, 'Reports', 'rows.json'), struct('rows', rows, 'cleanup', struct('pattern_id', rows(1).pattern_id, 'gray_pfn_id', rows(1).pfn_id, 'zero_afn_id', 0), 'programmed_s', sum([rows.duration_s])));
    B = struct('n_bars', 88, 'n_reps', P.n_reps, 'n_flashes', 88*2*P.n_reps, 'programmed_s', sum([rows.duration_s]), 'pattern_id', P.bar_flash_patterns(1));
    say(S, 'FLAG stage 2 fallback: all 88 bars as today (%d rows, %.1f min).', numel(rows), B.programmed_s/60);
end
S.steps(strcmp({S.steps.id}, 'p2b')).est_min = B.programmed_s / 60;
status = run_rows(opts, out_dir, rows, [], S);
S.out.p2b = struct('folder', out_dir, 'mode', mode, 'n_bars', B.n_bars, 'n_reps', B.n_reps, 'n_flashes', B.n_flashes, ...
    'rows', numel(rows), 'state', status.state, 'log_folder', fullfile(out_dir, 'Log Files'));
end

function S = step_readout(S, opts)
% Read-out call 2 (or the only call in 'full' mode): center, window, sites, tables, figure. Then the
% automatic decisions (c) and (f): the center is accepted when the sanity checks pass, else it falls
% back to the Protocol 1 peak (tables shifted there, flagged); the 'center' option overrides (logged).
out_dir = fullfile(S.exp_folder, 'readout'); if ~isfolder(out_dir), mkdir(out_dir); end
tag = S.cell;
method = call_readout(2, S, opts, out_dir);
R = read_readout(out_dir, tag);
S.out.readout = struct('dir', out_dir, 'method', method, 'summary', R.summary, 'tables', R.files);
say(S, 'Read-out (%s): %s', method, summary_line(R.summary));
if ~isempty(R.summary.warnings), say(S, 'FLAG read-out warnings: %s', strjoin(R.summary.warnings, ' | ')); end
if ~isempty(R.summary.placement)
    S = decide(S, 'readout', 'placement_pd', 'displayed from the read-out (rule applied there)', R.summary.placement, ...
        sprintf('placement PD %s deg (rule %s; snapped %s deg)', num2str(getf(R.summary.placement, 'pd_placement', [])), ...
        getf(R.summary.placement, 'pd_rule', ''), num2str(getf(R.summary.placement, 'pd_deg', []))), false, ...
        tern(getf(R.summary.placement, 'pd_flag_gt5', false), '28/56 disagree > 5 deg', ''));
end
show_readout_figure(opts, fullfile(out_dir, [tag '_locations.png']));

% ---- decision (c)/(f): the center
center = R.summary.center_col_row(:).';
p1xy = [];
if isfield(S.out, 'p1_peak') && isfield(S.out.p1_peak, 'xy') && ~isempty(S.out.p1_peak.xy), p1xy = double(S.out.p1_peak.xy(:).'); end
if isempty(p1xy) && ~isempty(R.summary.p1_peak_col_row), p1xy = R.summary.p1_peak_col_row(:).'; end
checks = struct();
inputs = struct('readout_center', center, 'p1_peak', p1xy, 'centroid_offset_px', R.summary.centroid_offset, ...
    'center_peak_mV', R.summary.center_peak_mV, 'profiles_unimodal', R.summary.profiles_unimodal, ...
    'thresholds', struct('center_max_dist_px', opts.center_max_dist_px, 'center_max_offset_px', opts.center_max_offset_px, 'peak_floor_mV', opts.peak_floor_mV));
if ~isempty(center) && ~isempty(p1xy)
    checks.dist_to_p1_px = norm(center - p1xy); checks.near_p1 = checks.dist_to_p1_px <= opts.center_max_dist_px;
else
    checks.dist_to_p1_px = NaN; checks.near_p1 = 'n/a';
end
if ~isempty(R.summary.centroid_offset)
    co = R.summary.centroid_offset;
    if isstruct(co), co = [getf(co, 's', NaN), getf(co, 'v', NaN)]; end
    if all(isfinite(co)), checks.inside_bar_set = all(abs(co) <= opts.center_max_offset_px); else, checks.inside_bar_set = 'n/a'; end
else
    checks.inside_bar_set = 'n/a';
end
if ~isempty(R.summary.center_peak_mV), checks.peak_above_floor = R.summary.center_peak_mV >= opts.peak_floor_mV; else, checks.peak_above_floor = 'n/a'; end
if ~isempty(R.summary.profiles_unimodal), checks.unimodal = all(logical(R.summary.profiles_unimodal)); else, checks.unimodal = 'n/a'; end
failed = {};
for c = {'near_p1', 'inside_bar_set', 'peak_above_floor', 'unimodal'}
    v = checks.(c{1});
    if islogical(v) && ~v, failed{end+1} = c{1}; end %#ok<AGROW>
end
no_pd = isfield(S.out, 'readout1') && isfield(S.out.readout1, 'usable_pd') && ~S.out.readout1.usable_pd;
shift = [0 0]; source = 'read-out M6 center'; fallback = false; flags = '';
override = opts.center;
if ~isempty(override)
    if ischar(override) && strcmp(override, 'p1_peak'), target = p1xy; src = 'override: Protocol 1 peak';
    else, target = double(override(:).'); src = sprintf('override: operator center %s', mat2str(target)); end
    if isempty(target), error('run_session:stop', 'center override requested but no Protocol 1 peak position is known'); end
    shift = round(target - center); source = src; flags = 'override';
elseif no_pd
    if isempty(p1xy), error('run_session:stop', 'no usable PD and no Protocol 1 peak position: cannot place the sites'); end
    shift = round(p1xy - center); source = 'Protocol 1 peak (no usable PD after stage 1)'; fallback = true; flags = 'no_pd_fallback';
elseif ~isempty(failed)
    if isempty(p1xy)
        flags = ['checks failed (' strjoin(failed, ', ') ') but no Protocol 1 peak position: read-out center kept'];
        fallback = true;
    else
        shift = round(p1xy - center); source = 'Protocol 1 peak (sanity check failed)'; fallback = true;
        flags = ['failed: ' strjoin(failed, ', ')];
    end
end
S.center_shift_px = shift;
S.out.readout.checks = checks; S.out.readout.center_source = source; S.out.readout.center_shift_px = shift;
inputs.checks = checks;
S = decide(S, 'readout', 'center', ['accept the read-out center when: within center_max_dist_px of the Protocol 1 peak, |M6 offset| <= center_max_offset_px (inside the bar set), ' ...
    'center peak >= peak_floor_mV, profiles unimodal (when reported); else the Protocol 1 peak; the center option overrides'], ...
    inputs, sprintf('%s; tables shifted by [%d %d] px', source, shift(1), shift(2)), fallback, flags);
end

function S = step_build(S, opts)
if ~isfield(S.out, 'readout') || ~isfield(S.out.readout, 'dir'), error('run_session:skip', 'no read-out'); end
rd = S.out.readout.dir; tag = S.cell; shift = S.center_shift_px;
common = struct('cell', tag, 'side', S.side, 'ao_channel', opts.ao_channel, 'pretrial_s', opts.pretrial_s, ...
    'ao_zero_tail_s', opts.ao_zero_tail_s, 'center_shift_px', shift, 'readout_dir', rd, 'dry_run', opts.dry_run, ...
    'center_source', getf(S.out.readout, 'center_source', ''));
p3b = build_prot3b_package(rd, fullfile(S.exp_folder, 'p3b'), common);
p4 = build_prot4_package(rd, fullfile(S.exp_folder, 'p4'), common);
S.out.build = struct('p3b', p3b, 'p4', p4);
S.steps(strcmp({S.steps.id}, 'p3b')).est_min = p3b.programmed_s / 60;
S.steps(strcmp({S.steps.id}, 'p4')).est_min = p4.programmed_s / 60;
say(S, 'Prot3b: %d frames, %d blocks, %.1f min; Protocol 4: %d frames, %d flashes, %.1f min.%s', ...
    p3b.n_frames, p3b.n_blocks, p3b.programmed_s/60, p4.n_frames, p4.n_flashes, p4.programmed_s/60, ...
    tern(p3b.pat_written && p4.pat_written, '', '  (.pat files NOT written: no G4 pattern encoder on this machine)'));
if ~isempty(p3b.warnings), say(S, 'FLAG Prot3b build: %s', strjoin(p3b.warnings, ' | ')); end
if ~isempty(p4.warnings), say(S, 'FLAG Protocol 4 build: %s', strjoin(p4.warnings, ' | ')); end
end

function S = step_run_package(S, opts, which)
if ~isfield(S.out, 'build') || ~isfield(S.out.build, which), error('run_session:skip', 'package not built'); end
pkg = S.out.build.(which).folder;
rows = rows_from_package(pkg);
say(S, '%s: %d rows, %.1f min; AO%d carries the current command (0.25 V per pA) added to the fixed holding current.', ...
    which, numel(rows), sum([rows.duration_s])/60, opts.ao_channel);
status = run_rows(opts, pkg, rows, opts.ao_channel, S);
S.out.(which) = struct('folder', pkg, 'rows', numel(rows), 'state', status.state, 'elapsed_s', status.elapsed_s, ...
    'log_folder', fullfile(pkg, 'Log Files'));
end

% =========================================================================================
% read-out calls (the files are the contract; rapid_readout.m or the Python script produce them)
% =========================================================================================
function method = call_readout(call, S, opts, out_dir)
tag = S.cell;
if ~isempty(opts.readout_dir)
    copy_readout_tables(opts.readout_dir, out_dir, tag); method = 'tables given'; return
end
if ~isempty(opts.readout_fn)
    opts.readout_fn(call, S, out_dir); method = 'readout_fn'; return
end
if strcmp(S.p2_mode, 'full'), p2 = S.out.p2.folder; else, p2 = S.out.p2a.folder; end
if exist('rapid_readout', 'file') == 2
    % src/analysis/protocol2/rapid_readout.m, interface of 26 Oct 2026 (evening): call 1 on the stage-1 folder
    % writes <tag>_stage1.json and the stage-2 tables; call 2 adds 'Stage2', <session>/p2b and takes the bar
    % flashes from that log; 'P1Peak' puts the Protocol 1 peak into other_centers and the distance check.
    method = 'rapid_readout.m';
    extra = {'Cell', tag, 'Out', out_dir, 'Figure', opts.show_figures};
    if isfield(S.out, 'p1_peak') && isfield(S.out.p1_peak, 'xy') && ~isempty(S.out.p1_peak.xy)
        extra = [extra, {'P1Peak', double(S.out.p1_peak.xy(:).')}];
    end
    if call == 2 && strcmp(S.p2_mode, 'two_stage')
        extra = [extra, {'Stage2', S.out.p2b.folder}];
    end
    rapid_readout(p2, extra{:});
    return
end
script = fullfile(opts.analysis_repo, 'tools', 'rapid_prot2_fast.py');
if ~isfile(script), script = fullfile(opts.analysis_repo, 'tools', 'rapid_prot2_locations.py'); end
method = ['python3 ' script];
cmd = sprintf('python3 "%s" --session "%s" --cell %s --out "%s"', script, p2, tag, out_dir);
say(S, '%s', cmd);
[rc, txt] = system(cmd);
write_text(fullfile(out_dir, sprintf('readout_call%d_stdout.txt', call)), txt);
if rc ~= 0, error('run_session:readout', 'Read-out script failed (%d): %s', rc, txt); end
end

function show_readout_figure(opts, png)
if opts.show_figures && isfile(png) && usejava('desktop')
    try, figure('Name', 'Rapid read-out', 'NumberTitle', 'off'); imshow(imread(png)); drawnow; catch, end
end
end

function copy_readout_tables(src, dst, tag)
n = 0;
for f = dir(fullfile(src, [tag '_*'])).'
    copyfile(fullfile(f.folder, f.name), fullfile(dst, f.name)); n = n + 1;
end
if n == 0, error('run_session:readout', 'No %s_* read-out files in %s', tag, src); end
end

function R = read_readout(out_dir, tag)
need = {'prot3b_frames', 'prot3b_blocks', 'prot4_frames', 'prot4_sequence'};
R = struct(); R.files = struct();
for k = 1:numel(need)
    f = fullfile(out_dir, sprintf('%s_%s.csv', tag, need{k}));
    if ~isfile(f), error('run_session:readout', 'Missing read-out table %s', f); end
    R.files.(need{k}) = f;
end
R.summary = struct('pd_deg', [], 'pd_vector_deg', [], 'pd_placement', [], 'pd_rule', '', 'placement', [], 'center_col_row', [], ...
    'static_vm_mV', [], 'warnings', {{}}, 'n_sites', [], 'prot3b_min', [], 'prot4_min', [], 'centroid_offset', [], ...
    'center_peak_mV', [], 'profiles_unimodal', [], 'p1_peak_col_row', [], 'sequencer_inputs', struct());
js = fullfile(out_dir, [tag '_locations.json']);
if isfile(js)
    L = jsondecode(fileread(js));
    R.files.locations = js;
    R.summary.pd_deg = getf(L, 'pd_deg', []);
    R.summary.pd_vector_deg = getf(L, 'pd_vector_deg', getf(L, 'pd_vector_28', []));
    R.summary.pd_placement = getf(L, 'pd_placement', []); R.summary.pd_rule = getf(L, 'pd_rule', '');
    R.summary.placement = getf(L, 'placement', []);
    R.summary.center_col_row = getf(L, 'rf_center_col_row_1based', []); R.summary.static_vm_mV = getf(L, 'static_vm_mV', []);
    si = getf(L, 'sequencer_inputs', struct());
    % placement-frame centroid offset when the read-out gives it (sets actually shown), else the TUTL-axes one
    R.summary.centroid_offset = getf(si, 'centroid_offset_placement_px', getf(L, 'centroid_offset_px', []));
    pu = getf(si, 'profiles_unimodal', getf(L, 'profiles_unimodal', []));
    if ischar(pu) || isstring(pu), pu = []; end               % "n/a": no rule the read-out can justify
    R.summary.profiles_unimodal = pu;
    R.summary.center_peak_mV = getf(si, 'center_expected_ON_80ms_peak_mV', []);
    R.summary.sequencer_inputs = si;
    oc = getf(L, 'other_centers', struct());
    R.summary.p1_peak_col_row = getf(oc, 'protocol1_ON_peak', []);
    w = getf(L, 'site_warnings', {}); if isnumeric(w), w = {}; end; if ischar(w), w = {w}; end
    if isstruct(w), w = arrayfun(@(x) jsonencode(x), w, 'UniformOutput', false); end
    R.summary.warnings = cellstr(w(:).');
    if isfield(L, 'sites') && ~isempty(L.sites)
        R.summary.n_sites = numel(L.sites);
        s1 = L.sites(1); if iscell(L.sites), s1 = L.sites{1}; end
        if isempty(R.summary.center_peak_mV) && isfield(s1, 'expected_ON_80ms_peak_mV') && strcmp(getf(s1, 'region', ''), 'center')
            R.summary.center_peak_mV = s1.expected_ON_80ms_peak_mV;
        end
    end
    if isfield(L, 'prot3b') && isfield(L.prot3b, 'total_s'), R.summary.prot3b_min = L.prot3b.total_s / 60; end
    if isfield(L, 'protocol4') && isfield(L.protocol4, 'total_s'), R.summary.prot4_min = L.protocol4.total_s / 60; end
end
end

function s = summary_line(sm)
parts = {};
if ~isempty(sm.pd_deg), parts{end+1} = sprintf('PD %g deg%s%s', sm.pd_deg, tern(isempty(sm.pd_vector_deg), '', sprintf(' (vector %.1f)', sm.pd_vector_deg)), tern(isempty(sm.pd_placement), '', sprintf(', placement %g', sm.pd_placement))); end
if ~isempty(sm.center_col_row), parts{end+1} = sprintf('center col/row %s', mat2str(sm.center_col_row(:).', 4)); end
if ~isempty(sm.static_vm_mV), parts{end+1} = sprintf('static Vm %.1f mV', sm.static_vm_mV); end
if ~isempty(sm.center_peak_mV), parts{end+1} = sprintf('center peak %.1f mV', sm.center_peak_mV); end
if ~isempty(sm.n_sites), parts{end+1} = sprintf('%d sites', sm.n_sites); end
if ~isempty(sm.prot3b_min), parts{end+1} = sprintf('Prot3b %.1f min', sm.prot3b_min); end
if ~isempty(sm.prot4_min), parts{end+1} = sprintf('Protocol 4 %.1f min', sm.prot4_min); end
if isempty(parts), s = 'tables present (no locations.json summary)'; else, s = strjoin(parts, ', '); end
end

% =========================================================================================
% running rows on the controller (modeled on run_protocol3's loop and cleanup guard)
% =========================================================================================
function status = run_rows(opts, pkg, rows, ao_channel, S)
rep_dir = fullfile(pkg, 'Reports'); if ~isfolder(rep_dir), mkdir(rep_dir); end
logdir = fullfile(pkg, 'Log Files');
if isfolder(logdir) && numel(dir(logdir)) > 2
    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    movefile(logdir, fullfile(pkg, ['Log Files_aborted_' stamp]));
    say(S, 'Existing logs moved to Log Files_aborted_%s', stamp);
end
status = struct('state', 'starting', 'package', pkg, 'started', iso_now(), 'finished', '', 'rows_total', numel(rows), ...
    'rows_completed', 0, 'rows', {struct('index', num2cell(1:numel(rows)), 'label', {rows.label}, 'pattern_id', {rows.pattern_id}, ...
    'pfn_id', {rows.pfn_id}, 'afn_id', {rows.afn_id}, 'duration_s', {rows.duration_s}, 'started', '', 'elapsed_s', [], 'result', 'pending')}, ...
    'cleanup', struct('step', {}, 'ok', {}, 'detail', {}), 'error', '', 'elapsed_s', []);
status_file = fullfile(rep_dir, 'run_status.json');
write_json(status_file, status);
ctlr = opts.controller_factory();
ctx = containers.Map('KeyType', 'char', 'ValueType', 'any');
ctx('ctlr') = ctlr; ctx('opened') = false; ctx('log') = false; ctx('done') = false; ctx('ao') = ao_channel;
ctx('cleanup_ids') = cleanup_ids_of(pkg, rows);
guard = onCleanup(@() controller_cleanup(ctx)); %#ok<NASGU>
t_run = tic;
try
    ctlr.open(true); ctx('opened') = true;
    checked('setRootDirectory', ctlr.setRootDirectory(pkg));
    if ~isempty(ao_channel), checked('setActiveAOChannels', ctlr.setActiveAOChannels(ao_channel)); end
    checked('setActiveAIChannels', ctlr.setActiveAIChannels(opts.ai_channels));
    ctx('log') = true;
    checked('startLog', ctlr.startLog());
    status.state = 'running'; write_json(status_file, status);
    for k = 1:numel(rows)
        r = rows(k);
        fprintf('    row %3d/%d  %-26s pat %d  pfn %d%s  %6.1f s\n', k, numel(rows), r.label, r.pattern_id, r.pfn_id, ...
            tern(isempty(ao_channel), '', sprintf('  afn %d', r.afn_id)), r.duration_s);
        status.rows(k).started = iso_now(); t_row = tic;
        checked('setControlMode', ctlr.setControlMode(1));
        checked('setPatternID', ctlr.setPatternID(r.pattern_id));
        checked('setPatternFunctionID', ctlr.setPatternFunctionID(r.pfn_id));
        ctlr.setPositionX(1);
        if ~isempty(ao_channel), checked('setAOFunctionID', ctlr.setAOFunctionID(ao_channel, r.afn_id)); end
        checked('startDisplay', ctlr.startDisplay(r.deciseconds, true));
        status.rows(k).elapsed_s = toc(t_row); status.rows(k).result = 'completed'; status.rows_completed = k;
        write_json(status_file, status);
    end
    ctx('done') = true;
    status.cleanup = controller_cleanup(ctx);
    status.state = 'completed';
catch err
    status.error = sprintf('%s: %s', err.identifier, err.message);
    status.cleanup = controller_cleanup(ctx);
    status.state = 'failed'; status.finished = iso_now(); status.elapsed_s = toc(t_run);
    write_json(status_file, status);
    rethrow(err);
end
status.finished = iso_now(); status.elapsed_s = toc(t_run);
write_json(status_file, status);
end

function steps = controller_cleanup(ctx)
steps = struct('step', {}, 'ok', {}, 'detail', {});
if isKey(ctx, 'finished') && ctx('finished'), return; end
ctx('finished') = true;
ctlr = ctx('ctlr');
if isempty(ctlr), return; end
if ctx('opened')
    steps(end+1) = guarded('stopDisplay', @() ctlr.stopDisplay());
    if isKey(ctx, 'cleanup_ids') && ~isempty(ctx('cleanup_ids'))
        c = ctx('cleanup_ids');
        steps(end+1) = guarded('setControlMode(1)', @() ctlr.setControlMode(1));
        steps(end+1) = guarded('setPatternID(gray)', @() ctlr.setPatternID(c.pattern_id));
        steps(end+1) = guarded('setPatternFunctionID(gray)', @() ctlr.setPatternFunctionID(c.gray_pfn_id));
        steps(end+1) = guarded('setPositionX(1)', @() noret(@() ctlr.setPositionX(1)));
        if ~isempty(ctx('ao')) && c.zero_afn_id > 0
            steps(end+1) = guarded('setAOFunctionID(zero)', @() ctlr.setAOFunctionID(ctx('ao'), c.zero_afn_id));
        end
        steps(end+1) = guarded('startDisplay(5) flush', @() ctlr.startDisplay(5, true));
    end
    if ~isempty(ctx('ao'))
        steps(end+1) = guarded('setActiveAOChannels([]) (AO off)', @() ctlr.setActiveAOChannels(zeros(1, 0)));
    end
    if ctx('log'), steps(end+1) = guarded('stopLog', @() ctlr.stopLog('showTimeoutMessage', true)); end
end
steps(end+1) = guarded('close', @() noret(@() ctlr.close()));
if ~ctx('done')
    fprintf(2, 'run_session: controller cleanup after an interruption (%d steps, %d ok)\n', numel(steps), sum([steps.ok]));
end
end

function s = guarded(label, fn)
s = struct('step', label, 'ok', false, 'detail', '');
try
    rv = fn();
    if (islogical(rv) || isnumeric(rv)) && isscalar(rv), s.ok = logical(rv); s.detail = sprintf('returned %d', rv);
    else, s.ok = true; s.detail = 'no return value'; end
catch e
    s.detail = ['error: ' e.message];
end
end

function rv = noret(fn)
fn(); rv = [];
end

function checked(label, rtn)
if ~(islogical(rtn) || isnumeric(rtn)) || ~isscalar(rtn) || ~rtn
    error('run_session:controller', 'Controller command %s did not confirm; run stopped.', label);
end
end

function c = cleanup_ids_of(pkg, rows)
f = fullfile(pkg, 'Reports', 'rows.json');
c = struct('pattern_id', rows(1).pattern_id, 'gray_pfn_id', rows(1).pfn_id, 'zero_afn_id', 0);
if isfile(f)
    j = jsondecode(fileread(f));
    if isfield(j, 'cleanup'), c = j.cleanup; end
end
end

% =========================================================================================
% Protocol 1 and Protocol 2 packages
% =========================================================================================
function copy_package(src, dst)
if ~isfolder(dst), mkdir(dst); end
for d = {'Patterns', 'Functions', 'Analog Output Functions'}
    if isfolder(fullfile(src, d{1})), copyfile(fullfile(src, d{1}), fullfile(dst, d{1})); end
end
for f = [dir(fullfile(src, '*.g4p')); dir(fullfile(src, 'currentExp.mat'))].'
    copyfile(fullfile(f.folder, f.name), fullfile(dst, f.name));
end
end

function rows = rows_from_g4p(pkg)
d = dir(fullfile(pkg, '*.g4p'));
if ~isscalar(d), error('run_session:p1', 'Expected one .g4p in %s', pkg); end
q = load(fullfile(pkg, d.name), '-mat'); p = q.exp_parameters;
c = load(fullfile(pkg, 'currentExp.mat')); cx = c.currentExp;
pat_names = text_names(cx.pattern.pattNames);
fn_names = text_names(cx.function.functionName);
if p.is_randomized, error('run_session:p1', 'Randomized .g4p not supported'); end
raw = {}; labels = {};
if ~isempty(p.pretrial{12}) && p.pretrial{12} > 0, raw(end+1, :) = p.pretrial; labels{end+1} = 'pretrial'; end
for r = 1:p.repetitions
    for b = 1:size(p.block_trials, 1)
        raw(end+1, :) = p.block_trials(b, :); labels{end+1} = sprintf('rep %d block %d', r, b); %#ok<AGROW>
    end
end
if ~isempty(p.posttrial{12}) && p.posttrial{12} > 0, raw(end+1, :) = p.posttrial; labels{end+1} = 'posttrial'; end
rows = struct('label', {}, 'pattern_id', {}, 'pfn_id', {}, 'afn_id', {}, 'duration_s', {}, 'deciseconds', {});
for k = 1:size(raw, 1)
    pid = find(strcmp(pat_names, raw{k, 2})); fid = find(strcmp(fn_names, raw{k, 3}));
    if ~isscalar(pid) || ~isscalar(fid)
        error('run_session:p1', 'Row %d (%s): pattern/function not found once in currentExp', k, labels{k});
    end
    rows(k) = struct('label', labels{k}, 'pattern_id', pid, 'pfn_id', fid, 'afn_id', 0, ...
        'duration_s', raw{k, 12}, 'deciseconds', deciseconds(raw{k, 12}, labels{k}));
end
end

function names = text_names(c)
c = c(:);
c = c(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), c));
names = regexprep(cellstr(c), '\.mat$', '');
end

function [xy, polarity] = frame_to_xy(pkg, frame)
% Center (col, row) and polarity of a Protocol 1 frame, from the package's 6-px pattern (as
% patt_frame_to_coord does from the rig's fixed path).
d = dir(fullfile(pkg, 'Patterns', '0002_*.mat'));
if isempty(d), d = dir(fullfile(pkg, 'Patterns', '*.mat')); end
if isempty(d), xy = []; polarity = ''; return; end
q = load(fullfile(d(1).folder, d(1).name)); Pats = q.pattern.Pats;
if frame > size(Pats, 3), error('run_session:p1', 'Peak frame %d exceeds the %d frames of %s', frame, size(Pats, 3), d(1).name); end
f = Pats(:, :, frame); bkg = mode(f(:));
[r, c] = find(f ~= bkg);
if isempty(r), xy = []; polarity = 'background'; return; end
xy = [median(c), median(r)];
polarity = tern(max(f(:)) > bkg, 'ON', 'OFF');
end

function md = p2_metadata(S)
md = struct('Frame', S.out.p1_peak.frame, 'Side', S.side, 'Age', S.age, 'Strain', S.strain, ...
    'Date', char(datetime('now', 'Format', 'yyyy_MM_dd')), 'Time', char(datetime('now', 'Format', 'HH-mm-ss')), ...
    'NotesStart', sprintf('run_session %s', S.cell));
end

function o = dry_p2_placeholder(S, id, what)
folder = fullfile(S.exp_folder, id);
if ~isfolder(folder), mkdir(folder); end
write_text(fullfile(folder, 'DRY_RUN.txt'), sprintf('Dry run: %s was not generated or run (needs the rig paths and G4_Display_Tools). Live: %s.\n', id, what));
o = struct('folder', folder, 'mode', 'dry_run', 'state', 'dry_run');
say(S, 'Dry run: %s skipped (placeholder %s).', id, folder);
end

function [exp_folder, pattern_order, func_order, trial_dur] = generate_protocol2_package(peak_frame, side, metadata)
% Mirrors generate_protocol2 (nested_RF_stimulus, 26 Oct 2026) without its dialog, the run and the
% Google-sheet export. Keep in step with that file. Needs the rig (G4_Display_Tools, C:\matlabroot paths).
px_intensity = [4, 0, 15]; px_crop_flash = 30; px_crop_bar = 30; n_flank = 5; n_reps = 3;
[x, y, on_off] = patt_frame_to_coord(peak_frame, side);
if x < 17 + px_crop_flash/2 || x > 192 - px_crop_flash/2, warning('run_session:p2', 'x %d is near the screen edge', x); end
if y < 1 + px_crop_flash/2 || y > 48 - px_crop_flash/2, warning('run_session:p2', 'y %d is near the screen edge', y); end
exp_path = 'C:\matlabroot\G4_Protocols\nested_RF_protocol2';
exp_name = string(datetime('now', 'TimeZone', 'local', 'Format', 'yyyy_MM_dd_HH_mm'));
exp_folder = create_exp_dir_G4(exp_name, exp_path);
generate_flash_stimulus_xy(x, y, px_intensity, px_crop_flash, 4, on_off, 0.16, 0.44, exp_folder);
generate_bar_stimulus_xy(x, y, px_intensity, px_crop_bar, on_off, exp_folder);
fn_dir = fullfile(exp_folder, 'Functions');
generate_bar_pos_fns(fn_dir, px_crop_bar);
generate_bar_flash_stimulus_xy(x, y, px_intensity, px_crop_bar, on_off, exp_folder, n_flank);
generate_bar_flash_pos_fns(fn_dir, n_flank, n_reps);
[pattern_order, func_order, trial_dur] = create_protocol2(exp_folder, metadata);
end

function P = p2_plan(folder, opts)
% The condition list of a Protocol 2 package (currentExp.mat: pattern_order, func_order, trial_dur) split
% into stage 1 (everything but the bar-flash patterns) and stage 2 (the bar-flash conditions as
% run_protocol2 plays them: patterns max-1 / max, substituted by +1 / +2 on repeats 2 and 3).
c = load(fullfile(folder, 'currentExp.mat'));
if ~isfield(c, 'pattern_order'), error('run_session:p2', '%s has no pattern_order (not a create_protocol2 package)', folder); end
po = double(c.pattern_order(:).'); fo = double(c.func_order(:).'); td = double(c.trial_dur(:).');
n_reps = opts.n_reps_p2;
slow = max(po) - 1; fast = max(po);
if ~isempty(opts.p2_bar_flash_pattern_id), slow = opts.p2_bar_flash_pattern_id; fast = slow + 1; end
is_bar = po == slow | po == fast;
P = struct('pattern_order', po, 'func_order', fo, 'trial_dur', td, 'n_reps', n_reps, 'bar_flash_patterns', [slow fast], ...
    'stage1_idx', find(~is_bar), 'stage2_idx', find(is_bar));
% Protocol 2 packages carry trial durations such as 3.328 s; run_protocol2 plays them as
% ceil(dur*10)-1 deciseconds, and so does the sequencer for these rows (never for the packages it builds).
rows = struct('label', {}, 'pattern_id', {}, 'pfn_id', {}, 'afn_id', {}, 'duration_s', {}, 'deciseconds', {});
for r = 1:n_reps
    for i = P.stage1_idx
        rows(end+1) = row_p2(sprintf('rep %d cond %d', r, i), po(i), fo(i), td(i)); %#ok<AGROW>
    end
end
P.stage1_rows = rows;
rows = struct('label', {}, 'pattern_id', {}, 'pfn_id', {}, 'afn_id', {}, 'duration_s', {}, 'deciseconds', {});
for r = 1:n_reps
    for i = P.stage2_idx
        pid = po(i) + (r - 1) * (r <= 3);         % run_protocol2: reps 2 and 3 use pattern +1 / +2
        rows(end+1) = row_p2(sprintf('rep %d bars cond %d', r, i), pid, fo(i), td(i)); %#ok<AGROW>
    end
end
P.stage2_rows_all88 = rows;
end

function r = row_p2(label, pattern_id, pfn_id, duration_s)
ds = max(1, ceil(duration_s * 10) - 1);            % run_protocol2's startDisplay(ceil(trial_t*10)-1)
r = struct('label', label, 'pattern_id', pattern_id, 'pfn_id', pfn_id, 'afn_id', 0, 'duration_s', duration_s, 'deciseconds', ds);
end

function B = build_stage2_package(p2a, seq_csv, out_dir, P, opts, tag)
% Stage-2 package: the bar-flash pattern copied from the stage-1 package as pattern 1 (its name keeps
% "FLASHES"), one position function per repeat and class from the stage-2 sequence table (repeat, cls,
% order, bar_pattern_frame_1based, flash_ms, gray_ms), 3 s gray first; 2 s gray pretrial and posttrial.
if isempty(seq_csv) || ~isfile(seq_csv), error('run_session:p2', 'stage-2 sequence table missing (%s)', seq_csv); end
Q = readtable(seq_csv, 'TextType', 'string', 'Delimiter', ',');
pat_id = P.bar_flash_patterns(1);
src = dir(fullfile(p2a, 'Patterns', sprintf('%04d_*.mat', pat_id)));
if ~isscalar(src), error('run_session:p2', 'bar-flash pattern %04d_*.mat not found once in %s', pat_id, fullfile(p2a, 'Patterns')); end
dirs = new_package_dir(out_dir);
copyfile(fullfile(src.folder, src.name), fullfile(dirs.pat, ['0001_' src.name(6:end)]));
pat_bin = fullfile(p2a, 'Patterns', sprintf('pat%04d.pat', pat_id));
if isfile(pat_bin), copyfile(pat_bin, fullfile(dirs.pat, 'pat0001.pat')); end
q = load(fullfile(src.folder, src.name)); n_frames = size(q.pattern.Pats, 3);
frames_all = Q.bar_pattern_frame_1based.';
if any(frames_all < 2 | frames_all > n_frames), error('run_session:p2', 'stage-2 frames outside 2..%d of the bar pattern', n_frames); end
reps_in_table = unique(Q.repeat).';
n_reps = opts.n_reps_stage2;
settle_s = 3;
pfn(1) = write_function(pfn_sections(opts.pretrial_s, 1, [], [], [], n_frames), 'pfn', 1, dirs.fn, sprintf('gray_%gs', opts.pretrial_s), n_frames);
rows = struct('label', {}, 'pattern_id', {}, 'pfn_id', {}, 'afn_id', {}, 'duration_s', {}, 'deciseconds', {});
rows(1) = row('pretrial', 1, 1, 0, opts.pretrial_s);
rng(uint32(mod(sum(double(tag)) * 7919, 2^31)), 'twister');   % reproducible extra repeats, seeded from the cell tag
% The read-out classifies stage-2 trials by the position-function file names ("flashes_80ms" / "flashes_14ms",
% the Protocol 2 generator's names) and finds the bar geometry in the one pattern named "FLASHES": one
% function per repeat and class, 3 s gray first, then the table's order of that class.
eh = {'repeat', 'cls', 'order', 'pattern_frame', 'flash_ms', 'gray_ms', 'position_function', 'flash_onset_s', 'flash_onset_programmed_s'};
er = {}; t = opts.pretrial_s; n_flash = 0; fid_next = 2;
for r = 1:n_reps
    src_rep = reps_in_table(mod(r - 1, numel(reps_in_table)) + 1);
    T = Q(Q.repeat == src_rep, :);
    [~, first] = unique(T.cls, 'stable');
    cls_order = T.cls(sort(first));
    for ci = 1:numel(cls_order)
        Tc = sortrows(T(T.cls == cls_order(ci), :), 'order');
        if r > numel(reps_in_table), Tc = Tc(randperm(height(Tc)), :); end   % more repeats than the table: fresh order
        frames = Tc.bar_pattern_frame_1based.'; fl = Tc.flash_ms.'/1000; gr = Tc.gray_ms.'/1000;
        func = pfn_sections(settle_s, 1, frames, fl, gr, n_frames);
        dur = r6(numel(func)/500);
        name = sprintf('bar_flashes_%dms_%dms_%dbars_rep%d', Tc.flash_ms(1), Tc.gray_ms(1), numel(frames), r);
        pfn(fid_next) = write_function(func, 'pfn', fid_next, dirs.fn, name, n_frames); %#ok<AGROW>
        rows(end+1) = row(sprintf('stage 2 rep %d %s (%d bars)', r, char(cls_order(ci)), numel(frames)), 1, fid_next, 0, dur); %#ok<AGROW>
        onset = settle_s;
        for e = 1:numel(frames)
            er(end+1, :) = {r, char(Tc.cls(e)), Tc.order(e), frames(e), Tc.flash_ms(e), Tc.gray_ms(e), pfn(fid_next).stem, r6(onset), r6(t + onset)}; %#ok<AGROW>
            onset = onset + fl(e) + gr(e);
        end
        t = t + dur; n_flash = n_flash + numel(frames); fid_next = fid_next + 1;
    end
end
rows(end+1) = row('posttrial', 1, 1, 0, opts.pretrial_s);
programmed = sum([rows.duration_s]);
pat = struct('id', 1, 'stem', ['0001_' src.name(6:end-4)], 'x_num', n_frames, 'y_num', 1, 'gs_val', 4, 'arena_pitch', 0);
write_current_exp(out_dir, pat, pfn, []);
write_json(fullfile(dirs.rep, 'rows.json'), struct('rows', rows, 'cleanup', struct('pattern_id', 1, 'gray_pfn_id', 1, 'zero_afn_id', 0), 'programmed_s', programmed));
write_csv(fullfile(dirs.rep, 'event_ledger.csv'), eh, er);
cfg = struct('schema', 'run_session/2', 'protocol', 'Protocol 2 stage 2', 'cell', tag, 'generated', iso_now(), 'source_sequence', seq_csv, ...
    'source_sequence_sha256', sha256_file(seq_csv), 'bar_pattern_from', fullfile(src.folder, src.name), 'bar_pattern_id_in_stage1', pat_id, ...
    'n_bars', numel(unique(frames_all)), 'repeats', n_reps, 'repeats_in_table', numel(reps_in_table), 'settle_s', settle_s, ...
    'programmed_duration_s', programmed, 'flashes', n_flash);
write_json(fullfile(dirs.rep, 'config.json'), cfg);
B = struct('folder', out_dir, 'n_bars', numel(unique(frames_all)), 'n_reps', n_reps, 'n_flashes', n_flash, 'programmed_s', programmed, 'pattern_id', pat_id);
end

% =========================================================================================
% Prot3b package from the read-out tables
% =========================================================================================
function out = build_prot3b_package(readout_dir, out_dir, c)
tag = c.cell;
F = readtable(fullfile(readout_dir, [tag '_prot3b_frames.csv']), 'TextType', 'string', 'Delimiter', ',');
B = readtable(fullfile(readout_dir, [tag '_prot3b_blocks.csv']), 'TextType', 'string', 'Delimiter', ',');
[Pats, G, warn] = pats_from_frame_table(F, c.center_shift_px);
n_frames = size(Pats, 3);
dirs = new_package_dir(out_dir);
pat = write_pattern(Pats, dirs.pat, 1, sprintf('prot3b_%dframes_4px_sites_mi9', n_frames), struct('source', 'prot3b_frames.csv', 'cell', tag));
gray = pfn_sections(c.pretrial_s, 1, [], [], [], n_frames);
pfn(1) = write_function(gray, 'pfn', 1, dirs.fn, sprintf('gray_%gs', c.pretrial_s), n_frames);
nb = height(B);
for b = 1:nb
    frames = sscanf(char(B.frames(b)), '%d').';
    if numel(frames) ~= B.n_flashes(b), error('run_session:build', 'Block %d: %d frames listed, n_flashes %d', b, numel(frames), B.n_flashes(b)); end
    if any(frames < 2 | frames > n_frames), error('run_session:build', 'Block %d references a frame outside 2..%d', b, n_frames); end
    func = pfn_sections(B.settle_s(b), 1, frames, B.flash_ms(b)/1000, B.gray_ms(b)/1000, n_frames);
    if abs(numel(func)/500 - B.duration_s(b)) > 1e-6
        error('run_session:build', 'Block %d: function lasts %g s, table says %g s', b, numel(func)/500, B.duration_s(b));
    end
    pfn(1+b) = write_function(func, 'pfn', 1+b, dirs.fn, sprintf('prot3b_block%02d_%s_%spA', b, char(B.phase(b)), tok(B.current_pA(b))), n_frames); %#ok<AGROW>
end
afn(1) = write_function(zeros(1, round(c.pretrial_s*1000)), 'afn', 1, dirs.ao, sprintf('zero_%gs', c.pretrial_s), n_frames);
keyv = [B.duration_s, B.ao_volts];
[uk, ~, ik] = unique(keyv, 'rows', 'stable');
afn_of_block = zeros(nb, 1);
for u = 1:size(uk, 1)
    dur = uk(u, 1); v = uk(u, 2);
    func = [v*ones(1, round((dur - c.ao_zero_tail_s)*1000)), zeros(1, round(c.ao_zero_tail_s*1000))];
    pa = B.current_pA(find(ik == u, 1));
    afn(1+u) = write_function(func, 'afn', 1+u, dirs.ao, sprintf('offset_%spA_%sV_%gs', tok(pa), tok(v), dur), n_frames); %#ok<AGROW>
    afn_of_block(ik == u) = 1 + u;
end
rows = struct('label', {}, 'pattern_id', {}, 'pfn_id', {}, 'afn_id', {}, 'duration_s', {}, 'deciseconds', {});
rows(1) = row('pretrial', 1, 1, 1, c.pretrial_s);
for b = 1:nb, rows(1+b) = row(sprintf('block %02d %s %+g pA', b, char(B.phase(b)), B.current_pA(b)), 1, 1+b, afn_of_block(b), B.duration_s(b)); end
rows(nb+2) = row('posttrial', 1, 1, 1, c.pretrial_s);
programmed = sum([rows.duration_s]);
write_current_exp(out_dir, pat, pfn, afn);
write_json(fullfile(dirs.rep, 'rows.json'), struct('rows', rows, 'cleanup', struct('pattern_id', 1, 'gray_pfn_id', 1, 'zero_afn_id', 1), 'programmed_s', programmed));
bh = {'block_index', 'phase', 'requested_current_pA', 'ao_volts', 'settle_s', 'n_flashes', 'flash_ms', 'gray_ms', ...
    'duration_s', 'programmed_start_s', 'position_function', 'ao_function', 'pattern', 'frames'};
br = cell(nb, numel(bh)); eh = {'block_index', 'phase', 'requested_current_pA', 'event_index_in_block', 'pattern_frame', 'item', ...
    'polarity', 'size_px', 'region', 's_px', 'v_px', 'square_row_start', 'square_row_end', 'square_col_start', 'square_col_end', ...
    'square_center_x', 'square_center_y', 'flash_onset_s', 'flash_onset_programmed_s', 'flash_duration_s', 'gray_after_s'};
er = {}; t = c.pretrial_s;
for b = 1:nb
    frames = sscanf(char(B.frames(b)), '%d').';
    br(b, :) = {b, char(B.phase(b)), B.current_pA(b), B.ao_volts(b), B.settle_s(b), B.n_flashes(b), B.flash_ms(b), B.gray_ms(b), ...
        B.duration_s(b), r6(t), pfn(1+b).stem, afn(afn_of_block(b)).stem, pat.stem, char(B.frames(b))};
    for e = 1:numel(frames)
        g = G(frames(e)); onset = r6(B.settle_s(b) + (e-1)*(B.flash_ms(b) + B.gray_ms(b))/1000);
        er(end+1, :) = {b, char(B.phase(b)), B.current_pA(b), e, frames(e), g.item, g.polarity, g.size_px, g.region, g.s_px, g.v_px, ...
            g.row_start, g.row_end, g.col_start, g.col_end, g.col_center, g.row_center, onset, r6(t + onset), B.flash_ms(b)/1000, B.gray_ms(b)/1000}; %#ok<AGROW>
    end
    t = t + B.duration_s(b);
end
write_csv(fullfile(dirs.rep, 'block_ledger.csv'), bh, br);
write_csv(fullfile(dirs.rep, 'event_ledger.csv'), eh, er);
cfg = package_config(c, 'Prot3b', readout_dir, {[tag '_prot3b_frames.csv'], [tag '_prot3b_blocks.csv']}, n_frames, programmed, warn);
cfg.blocks = nb; cfg.ao_levels_pA = unique(B.current_pA).'; cfg.ao_volts_per_pA = 0.25;
cfg.ao_shape = sprintf('level from t = 0 to %g s before the block end, then 0 (settle %g s at the level before the first flash)', c.ao_zero_tail_s, B.settle_s(1));
cfg.frames = 'frame 1 background; 2..12 sites OFF; 13..23 sites ON; 24..26 Mi9 series 3-px OFF (as the read-out table)';
write_json(fullfile(dirs.rep, 'config.json'), cfg);
out = struct('folder', out_dir, 'n_frames', n_frames, 'n_blocks', nb, 'n_events', size(er, 1), 'programmed_s', programmed, ...
    'pat_written', pat.pat_written, 'warnings', {warn}, 'rows', numel(rows));
end

% =========================================================================================
% Protocol 4 package from the read-out tables
% =========================================================================================
function out = build_prot4_package(readout_dir, out_dir, c)
tag = c.cell;
F = readtable(fullfile(readout_dir, [tag '_prot4_frames.csv']), 'TextType', 'string', 'Delimiter', ',');
Q = readtable(fullfile(readout_dir, [tag '_prot4_sequence.csv']), 'TextType', 'string', 'Delimiter', ',');
[Pats, G, warn] = pats_from_frame_table(F, c.center_shift_px);
n_frames = size(Pats, 3);
dirs = new_package_dir(out_dir);
pat = write_pattern(Pats, dirs.pat, 1, sprintf('prot4_%dframes_lattice_core_line_off', n_frames), struct('source', 'prot4_frames.csv', 'cell', tag));
settle_s = 3;
frames = Q.frame.'; flash_s = Q.flash_ms.'/1000; gray_s = Q.gray_ms.'/1000;
if any(frames < 2 | frames > n_frames), error('run_session:build', 'Protocol 4 sequence references a frame outside 2..%d', n_frames); end
func = pfn_sections(settle_s, 1, frames, flash_s, gray_s, n_frames);
dur = r6(numel(func)/500);
t_tab = Q.t_start_s.'; t_cum = [0, cumsum(flash_s(1:end-1) + gray_s(1:end-1))];
if max(abs(t_tab - t_cum)) > 0.051, warn{end+1} = sprintf('t_start_s in the sequence table deviates from the cumulative durations by up to %.3f s', max(abs(t_tab - t_cum))); end
gray = pfn_sections(c.pretrial_s, 1, [], [], [], n_frames);
pfn(1) = write_function(gray, 'pfn', 1, dirs.fn, sprintf('gray_%gs', c.pretrial_s), n_frames);
pfn(2) = write_function(func, 'pfn', 2, dirs.fn, sprintf('prot4_sequence_%dflashes_%gs', numel(frames), dur), n_frames);
afn(1) = write_function(zeros(1, round(c.pretrial_s*1000)), 'afn', 1, dirs.ao, sprintf('zero_%gs', c.pretrial_s), n_frames);
afn(2) = write_function(zeros(1, round(dur*1000)), 'afn', 2, dirs.ao, sprintf('zero_%gs', dur), n_frames);
rows = struct('label', {}, 'pattern_id', {}, 'pfn_id', {}, 'afn_id', {}, 'duration_s', {}, 'deciseconds', {});
rows(1) = row('pretrial', 1, 1, 1, c.pretrial_s);
rows(2) = row('protocol 4 sequence', 1, 2, 2, dur);
rows(3) = row('posttrial', 1, 1, 1, c.pretrial_s);
programmed = sum([rows.duration_s]);
write_current_exp(out_dir, pat, pfn, afn);
write_json(fullfile(dirs.rep, 'rows.json'), struct('rows', rows, 'cleanup', struct('pattern_id', 1, 'gray_pfn_id', 1, 'zero_afn_id', 1), 'programmed_s', programmed));
eh = {'event', 'pass', 'order', 'pattern_frame', 'item', 'polarity', 'size_px', 's_px', 'v_px', 'square_row_start', 'square_row_end', ...
    'square_col_start', 'square_col_end', 'square_center_x', 'square_center_y', 'flash_onset_s', 'flash_onset_programmed_s', 'flash_duration_s', 'gray_after_s'};
n = numel(frames); er = cell(n, numel(eh));
for e = 1:n
    g = G(frames(e)); onset = r6(settle_s + t_cum(e));
    er(e, :) = {e, Q.pass(e), Q.order(e), frames(e), g.item, g.polarity, g.size_px, g.s_px, g.v_px, g.row_start, g.row_end, g.col_start, g.col_end, ...
        g.col_center, g.row_center, onset, r6(c.pretrial_s + onset), flash_s(e), gray_s(e)};
end
write_csv(fullfile(dirs.rep, 'event_ledger.csv'), eh, er);
write_csv(fullfile(dirs.rep, 'block_ledger.csv'), {'block_index', 'phase', 'requested_current_pA', 'n_flashes', 'duration_s', 'programmed_start_s', 'position_function', 'ao_function', 'pattern'}, ...
    {1, 'protocol4', 0, n, dur, c.pretrial_s, pfn(2).stem, afn(2).stem, pat.stem});
items = unique(Q.item); counts = struct();
for k = 1:numel(items), counts.(matlab.lang.makeValidName(char(items(k)))) = sum(Q.item == items(k)); end
cfg = package_config(c, 'Protocol 4', readout_dir, {[tag '_prot4_frames.csv'], [tag '_prot4_sequence.csv']}, n_frames, programmed, warn);
cfg.passes = max(Q.pass); cfg.flashes = n; cfg.flashes_per_item = counts; cfg.settle_s = settle_s;
cfg.periods = 'lattice 160 ms flash + 440 ms gray; core, line and OFF squares 160 ms + 240 ms (per row of the sequence table)';
cfg.ao = sprintf('AO%d held at 0 V for the whole run', c.ao_channel);
write_json(fullfile(dirs.rep, 'config.json'), cfg);
out = struct('folder', out_dir, 'n_frames', n_frames, 'n_flashes', n, 'programmed_s', programmed, 'pat_written', pat.pat_written, ...
    'warnings', {warn}, 'rows', numel(rows));
end

% =========================================================================================
% shared builders
% =========================================================================================
function [Pats, G, warn] = pats_from_frame_table(F, shift)
warn = {};
n = max(F.frame);
if ~isequal(sort(F.frame(:)).', 1:n), error('run_session:build', 'Frame numbers must be 1..N without gaps'); end
bkg = 4; on = 15; off = 0;
Pats = bkg * ones(48, 192, n);
G = repmat(struct('frame', 0, 'item', 'background', 'polarity', 'gray', 'size_px', 0, 'region', '', 's_px', NaN, 'v_px', NaN, ...
    'row_start', NaN, 'row_end', NaN, 'col_start', NaN, 'col_end', NaN, 'col_center', NaN, 'row_center', NaN), 1, n);
has_region = ismember('region', F.Properties.VariableNames);
for k = 1:height(F)
    fr = F.frame(k); g = G(fr); g.frame = fr; g.item = txt(F.item(k)); g.polarity = txt(F.polarity(k)); g.size_px = F.size_px(k);
    if has_region, g.region = txt(F.region(k)); end
    if ismember('s_px', F.Properties.VariableNames), g.s_px = F.s_px(k); g.v_px = F.v_px(k); end
    if strcmp(g.item, 'background') || isnan(F.col_start(k)), G(fr) = g; continue; end
    c0 = F.col_start(k) + shift(1); r0 = F.row_start(k) + shift(2); sz = F.size_px(k);
    c1 = c0 + sz - 1; r1 = r0 + sz - 1;
    g.col_start = c0; g.col_end = c1; g.row_start = r0; g.row_end = r1; g.col_center = (c0 + c1)/2; g.row_center = (r0 + r1)/2;
    cc = max(c0, 17):min(c1, 192); rr = max(r0, 1):min(r1, 48);
    if numel(cc) < sz || numel(rr) < sz
        warn{end+1} = sprintf('frame %d (%s %s) %s screen after shift [%d %d]', fr, g.item, g.polarity, tern(isempty(cc) || isempty(rr), 'off', 'clipped at'), shift(1), shift(2)); %#ok<AGROW>
    end
    if ~isempty(cc) && ~isempty(rr)
        Pats(rr, cc, fr) = tern(strcmpi(g.polarity, 'ON'), on, off);
    end
    G(fr) = g;
end
end

function func = pfn_sections(settle_s, bkg_frame, frames, flash_s, gray_s, n_frames) %#ok<INUSD>
q = 0.002;
chk = @(x, nm) assert(all(abs(x/q - round(x/q)) < 1e-9), 'run_session:build', '%s must be whole multiples of 2 ms', nm);
chk(settle_s, 'settle');
func = bkg_frame * ones(1, round(settle_s/q));
if isempty(frames), return; end
if isscalar(flash_s), flash_s = repmat(flash_s, size(frames)); end
if isscalar(gray_s), gray_s = repmat(gray_s, size(frames)); end
chk(flash_s, 'flash'); chk(gray_s, 'gray');
nf = round(flash_s/q); ng = round(gray_s/q);
parts = cell(1, 2*numel(frames));
for k = 1:numel(frames)
    parts{2*k-1} = frames(k) * ones(1, nf(k)); parts{2*k} = bkg_frame * ones(1, ng(k));
end
func = [func, parts{:}];
end

function d = new_package_dir(out_dir)
if isfolder(out_dir) && ~isempty(dir(fullfile(out_dir, 'Patterns', '*.mat')))
    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    movefile(out_dir, [out_dir '_replaced_' stamp]);
end
d = struct('root', out_dir, 'pat', fullfile(out_dir, 'Patterns'), 'fn', fullfile(out_dir, 'Functions'), ...
    'ao', fullfile(out_dir, 'Analog Output Functions'), 'rep', fullfile(out_dir, 'Reports'));
for f = {'root', 'pat', 'fn', 'ao', 'rep'}, if ~isfolder(d.(f{1})), mkdir(d.(f{1})); end, end
end

function info = write_pattern(Pats, save_dir, id, name, extra)
param = struct('ID', id, 'stretch', zeros(size(Pats, 3), 1), 'gs_val', 4, 'arena_pitch', 0, 'checker_layout', 0, ...
    'px_rng', [1 48 17 192], 'generated_by', 'run_session', 'generated', iso_now());
for f = fieldnames(extra).', param.(f{1}) = extra.(f{1}); end
stem = sprintf('%04d_%s_G4', id, name);
info = struct('id', id, 'stem', stem, 'mat', fullfile(save_dir, [stem '.mat']), 'pat', fullfile(save_dir, sprintf('pat%04d.pat', id)), ...
    'pat_written', false, 'encoder', 'none', 'x_num', size(Pats, 3), 'y_num', 1, 'gs_val', 4, 'arena_pitch', 0);
if exist('save_pattern_G4', 'file') == 2
    save_pattern_G4(Pats, param, save_dir, name);
    info.pat_written = isfile(info.pat); info.encoder = 'save_pattern_G4';
    return
end
pattern = struct('Pats', Pats, 'x_num', size(Pats, 3), 'y_num', 1, 'gs_val', 4, 'stretch', param.stretch, 'param', param);
if exist('maDisplayTools', 'class') == 8 || exist('maDisplayTools', 'file') == 2
    try
        enc = pattern; enc.gs_val = 16;
        pattern.data = maDisplayTools.make_pattern_vector_g4(enc);
        fid = fopen(info.pat, 'w'); fwrite(fid, pattern.data); fclose(fid);
        info.pat_written = true; info.encoder = 'maDisplayTools.make_pattern_vector_g4';
    catch err
        warning('run_session:pattern', 'maDisplayTools encoder failed (%s); .pat not written', err.message);
    end
end
save(info.mat, 'pattern', '-v7');
end

function info = write_function(func, type, id, save_dir, name, n_frames)
stem = sprintf('%04d_%s_G4', id, name);
param = struct('func', func, 'ID', id, 'frames', n_frames, 'gs_val', 4, 'type', type, 'dur', numel(func) / tern(strcmp(type, 'pfn'), 500, 1000), ...
    'generated_by', 'run_session');
if strcmp(type, 'pfn')
    data = int16(func - 1); prefix = 'fun';
else
    assert(max(func) <= 10 && min(func) >= -10, 'run_session:build', 'AO function exceeds +/-10 V');
    data = int16(adconvert(func)); prefix = 'ao';
end
funcname = sprintf('%s%04d', prefix, id);
header = zeros(1, 512);
header(1:4) = double(typecast(uint32(numel(data)*2), 'uint8'));
header(5) = numel(funcname); header(6:5+numel(funcname)) = double(funcname);
bytes = [uint8(header), typecast(data(:).', 'uint8')];
param.size = numel(data) * 2;
bin = fullfile(save_dir, [funcname '.' type]);
fid = fopen(bin, 'w'); fwrite(fid, bytes, 'uint8'); fclose(fid);
if strcmp(type, 'pfn'), pfnparam = param; save(fullfile(save_dir, [stem '.mat']), 'pfnparam', '-v7');
else, afnparam = param; save(fullfile(save_dir, [stem '.mat']), 'afnparam', '-v7'); end
info = struct('id', id, 'stem', stem, 'bin', bin, 'size', param.size, 'dur', param.dur, 'type', type);
end

function code = adconvert(volts)
if exist('ADConvert', 'file') == 2, code = double(ADConvert(volts)); return; end
a = (double(volts) + 10) / 20 * 65535 - 32768;
code = min(max(sign(a) .* floor(abs(a) + 0.5), -32768), 32767);
end

function write_current_exp(out_dir, pat, pfn, afn)
x = struct();
x.pattern.pattNames = {[pat.stem '.mat']}; x.pattern.patternList = {sprintf('pat%04d.pat', pat.id)};
x.pattern.x_num = pat.x_num; x.pattern.y_num = pat.y_num; x.pattern.gs_val = pat.gs_val; x.pattern.arena_pitch = pat.arena_pitch;
x.pattern.num_patterns = 1;
x.function.functionName = arrayfun(@(p) [p.stem '.mat'], pfn, 'UniformOutput', false);
x.function.functionList = arrayfun(@(p) sprintf('fun%04d.pfn', p.id), pfn, 'UniformOutput', false);
x.function.functionSize = num2cell([pfn.size]); x.function.numFunc = numel(pfn);
if ~isempty(afn)
    x.aoFunction.aoFunctionName = arrayfun(@(p) [p.stem '.mat'], afn, 'UniformOutput', false);
    x.aoFunction.aoFunctionList = arrayfun(@(p) sprintf('ao%04d.afn', p.id), afn, 'UniformOutput', false);
    x.aoFunction.aoFunctionSize = num2cell([afn.size]); x.aoFunction.numaoFunc = numel(afn);
end
x.totalDuration = sum([pfn.dur]);
currentExp = x; %#ok<NASGU>
save(fullfile(out_dir, 'currentExp.mat'), 'currentExp', '-v7');
end

function rows = rows_from_package(pkg)
j = jsondecode(fileread(fullfile(pkg, 'Reports', 'rows.json')));
rows = j.rows(:).';
for k = 1:numel(rows)
    rows(k).deciseconds = deciseconds(rows(k).duration_s, rows(k).label);
end
end

function r = row(label, pattern_id, pfn_id, afn_id, duration_s)
r = struct('label', label, 'pattern_id', pattern_id, 'pfn_id', pfn_id, 'afn_id', afn_id, 'duration_s', duration_s, ...
    'deciseconds', deciseconds(duration_s, label));
end

function ds = deciseconds(dur, label)
raw = double(dur) * 10; ds = round(raw);
if abs(raw - ds) > 1e-6 || ds < 1 || ds > 65535
    error('run_session:duration', '%s: %g s is not a whole number of deciseconds in 1..65535', label, dur);
end
end

function cfg = package_config(c, name, readout_dir, tables, n_frames, programmed, warn)
cfg = struct('schema', 'run_session/2', 'protocol', name, 'cell', c.cell, 'side', c.side, 'generated', iso_now(), ...
    'readout_dir', readout_dir, 'source_tables', {tables}, 'center_shift_px', c.center_shift_px, 'center_source', getf(c, 'center_source', ''), ...
    'pattern_frames', n_frames, 'programmed_duration_s', programmed, 'programmed_duration_min', r6(programmed/60), 'pretrial_s', c.pretrial_s, ...
    'ao_channel', c.ao_channel, 'build_warnings', {warn}, 'dry_run', c.dry_run, ...
    'conventions', 'arena col/row 1-based, row 1 at the top; background 4, ON 15, OFF 0; s + = toward the PD (trailing/proximal); v + = 90 deg clockwise from PD');
for k = 1:numel(tables)
    f = fullfile(readout_dir, tables{k});
    if isfile(f), cfg.source_sha256.(matlab.lang.makeValidName(tables{k})) = sha256_file(f); end
end
end

% =========================================================================================
% Host configuration (live only)
% =========================================================================================
function S = apply_host_config_once(S, opts)
if isfield(S.out.setup, 'host_config_applied') && S.out.setup.host_config_applied, return; end
f = opts.host_config_file;
if isempty(f) && exist('G4_Protocol_Designer_Settings', 'file') == 2
    s = G4_Protocol_Designer_Settings(); f = char(s.Configuration_Filepath);
end
if isempty(f) || ~isfile(f)
    say(S, 'Host configuration not edited (no file given / found); the Host keeps its current ADC rates.');
    S.out.setup.host_config_applied = false; return
end
txt = fileread(f);
if contains(txt, sprintf('\r\n')), eol = sprintf('\r\n'); else, eol = newline; end
lines = strsplit(txt, eol, 'CollapseDelimiters', false);
targets = {'^\s*Number of Rows\s*=', 'Number of Rows = 3'};
for a = 0:3, targets(end+1, :) = {sprintf('^\\s*ADC%d Rate \\(Hz\\)\\s*=', a), sprintf('ADC%d Rate (Hz) = %d', a, opts.adc_rates_hz(a+1))}; end %#ok<AGROW>
for k = 1:size(targets, 1)
    hit = find(~cellfun(@isempty, regexp(lines, targets{k, 1}, 'once')));
    if ~isscalar(hit), error('run_session:hostcfg', 'Host configuration %s: expected one line matching %s', f, targets{k, 1}); end
    lines{hit} = targets{k, 2};
end
backup = fullfile(S.exp_folder, ['host_config_backup_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.ini']);
copyfile(f, backup);
write_text(f, strjoin(lines, eol));
S.out.setup = struct('host_config_applied', true, 'host_config_file', f, 'backup', backup, 'adc_rates_hz', opts.adc_rates_hz);
say(S, 'Host configuration %s: ADC rates set to %s (backup %s).', f, mat2str(opts.adc_rates_hz), backup);
end

% =========================================================================================
% console, files, small helpers
% =========================================================================================
function print_checklist(S)
fprintf('\n  Session %s  (%s)%s  p2_mode %s\n', S.cell, S.exp_folder, tern(S.dry_run, '  [DRY RUN]', ''), S.p2_mode);
fprintf('  %-3s %-9s %-11s %-8s %8s %7s  %s\n', '#', 'step', 'status', 'start', 'elapsed', 'est', 'note');
tot = 0;
for k = 1:numel(S.steps)
    st = S.steps(k);
    t0 = ''; if ~isempty(st.t_start), t0 = st.t_start(12:19); end
    el = ''; if st.elapsed_s > 0, el = sprintf('%.1f min', st.elapsed_s/60); tot = tot + st.elapsed_s; end
    mark = '.';
    switch st.status
        case 'done', mark = 'x';
        case 'running', mark = '>';
        case {'failed', 'interrupted'}, mark = '!';
        case 'skipped', mark = '-';
    end
    fprintf('  %s%-2d %-9s %-11s %-8s %8s %5.1f    %s\n', mark, k, st.id, st.status, t0, el, st.est_min, st.note);
end
fprintf('  elapsed so far %.1f min\n\n', tot/60);
end

function say(S, fmt, varargin)
msg = sprintf(fmt, varargin{:});
line = sprintf('[%s] %s', char(datetime('now', 'Format', 'HH:mm:ss')), msg);
fprintf('%s\n', line);
try
    fid = fopen(fullfile(S.exp_folder, 'session_log.txt'), 'a'); fprintf(fid, '%s %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd')), line); fclose(fid);
catch
end
end

function r = ask(opts, prompt, default)
% Used only for the operator's physical hand-offs (holding current; Conductor results folder).
if opts.auto_accept
    fprintf('%s[auto: %s]\n', prompt, tern(isempty(default), 'Enter', default));
    r = default; return
end
r = input(prompt, 's');
if isempty(r), r = default; end
end

function write_csv(pathname, header, rows)
fid = fopen(pathname, 'w'); if fid < 0, error('run_session:io', 'Could not write %s', pathname); end
c = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strjoin(header, ','));
for i = 1:size(rows, 1)
    parts = cell(1, size(rows, 2));
    for j = 1:size(rows, 2)
        v = rows{i, j};
        if ischar(v) || isstring(v), parts{j} = char(v);
        elseif islogical(v), parts{j} = lower(mat2str(v));
        elseif isempty(v) || (isscalar(v) && isnan(v)), parts{j} = '';
        elseif v == fix(v) && abs(v) < 1e15, parts{j} = sprintf('%d', v);
        else, parts{j} = sprintf('%.15g', v); end
        if strcmp(parts{j}, '-0'), parts{j} = '0'; end
        if any(parts{j} == ',' | parts{j} == '"' | parts{j} == ' ')
            parts{j} = ['"' strrep(parts{j}, '"', '""') '"'];
        end
    end
    fprintf(fid, '%s\n', strjoin(parts, ','));
end
end

function write_json(pathname, s)
write_text(pathname, [jsonencode(s, 'PrettyPrint', true) newline]);
end

function write_text(pathname, text)
fid = fopen(pathname, 'w'); if fid < 0, error('run_session:io', 'Could not write %s', pathname); end
fwrite(fid, unicode2native(char(text), 'UTF-8'), 'uint8'); fclose(fid);
end

function h = sha256_file(pathname)
fid = fopen(pathname, 'r'); bytes = fread(fid, inf, '*uint8'); fclose(fid);
md = java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end

function v = getf(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

function t = tok(x)
t = strrep(strrep(sprintf('%g', x), '-', 'minus'), '.', 'p');
if x > 0, t = ['plus' t]; end
end

function x = r6(x)
x = round(x * 1e6) / 1e6;
end

function s = txt(x)
if isstring(x) && (isempty(x) || ismissing(x)), s = ''; else, s = char(x); end
end

function t = iso_now()
t = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end

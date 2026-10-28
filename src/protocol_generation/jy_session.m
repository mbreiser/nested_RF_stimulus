function S = jy_session(tag, side, varargin)
%JY_SESSION  One cell's session on the rig, from one command: checks, the arena viewer, run_session (28 Oct 2026).
%
%   In MATLAB on the rig PC, after break-in (current clamp, 0 pA):
%       addpath('C:\matlabroot\GitHub\nested_RF_stimulus_session_ui\src\protocol_generation')
%       jy_session('JY1028a', 'L')
%   The same command again RESUMES the session after a FLAG, an error or Ctrl-C (run_session continues
%   from the first step that is not done). When the cell is lost:
%       jy_session('JY1028a', 'L', 'abandon', true)
%   Before the first cell, to check the PC without running anything:
%       jy_session('JY1028a', 'L', 'check_only', true)
%
%   What it does:
%     1. puts this repository's src first on the path (restored when it returns or is stopped) and checks
%        that run_session, the Protocol 2 generators and the read-outs come from it, that the real
%        PanelsController and G4_Display_Tools are there, the Host configuration (ADC1 and ADC2 at 10 kHz,
%        3 rows; it is checked, not edited) and the free disk space;
%     2. the session folder <root>\<cell tag>; a finished session is never re-run (choose a new cell tag);
%     3. the arena viewer (arena_viewer) in a second MATLAB, unless one is already running for this folder;
%        it closes itself 30 s after the session finishes, or when its window is closed;
%     4. run_session with the default two-stage sequence: dark30 (then the one prompt: type the holding
%        current), p1, p1_peak (rapid_p1_peak), p2a (the package generated from the peak as
%        generate_protocol2 makes it, in C:\matlabroot\G4_Protocols\nested_RF_protocol2\data), readout1,
%        p2b, readout, build, p3b, p4, the arena dark at the end. About 40 min.
%
%   Options (name, value): 'root' 'D:\nested_RF_sessions', 'p1_package' (default: the 6-px LHS package
%   for side 'L'; required for 'R'), 'viewer' true, 'check_only' false, 'branch'
%   'feat/session-progress-plots' ('' skips the check). Any other option goes to run_session
%   ('abandon', 'center', 'peak_frame', 'age', 'strain', ...; see help run_session).

o = struct('root', 'D:\nested_RF_sessions', 'p1_package', '', 'viewer', true, 'check_only', false, ...
    'branch', 'feat/session-progress-plots');
pass = {};
for k = 1:2:numel(varargin)
    if isfield(o, varargin{k}), o.(varargin{k}) = varargin{k+1}; else, pass = [pass, varargin(k:k+1)]; end %#ok<AGROW>
end
tag = char(tag); side = upper(char(side));
assert(~isempty(regexp(tag, '^[A-Za-z0-9_-]+$', 'once')), 'jy_session: cell tag %s: letters, digits, _ and - only', tag);
assert(ismember(side, {'L', 'R'}), 'jy_session: side must be L or R');
if isempty(o.p1_package)
    assert(strcmp(side, 'L'), 'jy_session: no default Protocol 1 package for side R; give ''p1_package''');
    o.p1_package = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\bkg4\LHS\protocol1_10kHz_4reps_6pxONLY_LHS_2sbkg_200msfl_150msint_09-18-26_12-59-40';
end

%% 1. path and checks
src = fileparts(fileparts(mfilename('fullpath')));
repo = fileparts(src);
old_path = path;
restore = onCleanup(@() path(old_path)); %#ok<NASGU>
addpath(genpath(src));
if ~isempty(o.branch)
    g = fullfile(repo, '.git');
    if isfile(g), g = strtrim(erase(fileread(g), 'gitdir:')); end
    head = strtrim(erase(fileread(fullfile(g, 'HEAD')), 'ref: refs/heads/'));
    assert(strcmp(head, o.branch), 'jy_session: %s is on %s, not %s', repo, head, o.branch);
end
for f = {'run_session', 'run_protocol2', 'create_protocol2', 'generate_bar_pos_fns', 'generate_bar_flash_pos_fns', ...
        'rapid_p1_peak', 'rapid_readout', 'arena_viewer'}
    assert(startsWith(which(f{1}), src), 'jy_session: %s comes from %s, not this repository', f{1}, which(f{1}));
end
pc = which('PanelsController');
assert(~isempty(pc) && ~contains(lower(pc), 'mock'), 'jy_session: the real PanelsController is not on the path (%s)', pc);
assert(~isempty(which('create_exp_dir_G4')) && ~isempty(which('G4_Function_Generator')), 'jy_session: G4_Display_Tools is not on the path');
s = G4_Protocol_Designer_Settings(); txt = fileread(char(s.Configuration_Filepath));
val = @(key) str2double(regexp(txt, [regexprep(key, '([()])', '\\$1') '\s*=\s*(\d+)'], 'tokens', 'once'));
rates = arrayfun(@(a) val(sprintf('ADC%d Rate (Hz)', a)), 0:3);
assert(isequal(rates, [0 10000 10000 0]) && val('Number of Rows') == 3, ...
    'jy_session: Host configuration %s has ADC rates %s and %g rows; expected [0 10000 10000 0] and 3', char(s.Configuration_Filepath), mat2str(rates), val('Number of Rows'));
assert(isfolder(o.p1_package) && ~isempty(dir(fullfile(o.p1_package, '*.g4p'))), 'jy_session: no Protocol 1 package (.g4p) in %s', o.p1_package);
if ~isfolder(o.root), mkdir(o.root); end
free_gb = @(d) java.io.File(d).getUsableSpace() / 1e9;
assert(free_gb(o.root) > 5, 'jy_session: less than 5 GB free for %s', o.root);
assert(free_gb('C:\matlabroot\G4_Protocols') > 2, 'jy_session: less than 2 GB free on C: (the Protocol 2 package goes there)');

%% 2. session folder
exp = fullfile(o.root, tag);
state = fullfile(exp, 'session_state.json');
resume = isfile(state);
if resume
    St = jsondecode(fileread(state));
    assert(~isfield(St, 'finished') || isempty(St.finished), 'jy_session: session %s finished at %s; use a new cell tag', exp, char(string(St.finished)));
end
fprintf('\njy_session %s, side %s: %s session %s\n', tag, side, tern(resume, 'RESUMING', 'new'), exp);
fprintf('  code %s (%s); Host config ADC %s, 3 rows; %.0f GB free\n', repo, o.branch, mat2str(rates), free_gb(o.root));
if o.check_only
    fprintf('  check_only: all checks passed; nothing run.\n\n'); S = []; return
end
if ~isfolder(exp), mkdir(exp); end
diary(fullfile(exp, ['jy_session_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.log']));
stop_diary = onCleanup(@() diary('off')); %#ok<NASGU>

%% 3. the viewer, in a second MATLAB (it only reads files)
hb = fullfile(exp, 'viewer_alive.txt');
d = dir(hb);
running = ~isempty(d) && (now - d.datenum) * 86400 < 20; %#ok<TNOW1>
if o.viewer && ~running
    r = sprintf(['addpath(genpath(''%s'')); V = arena_viewer(''%s''); ' ...
        'while ishghandle(V.figure), fid = fopen(''%s'', ''w''); fprintf(fid, ''%%s'', datestr(now)); fclose(fid); pause(5); ' ...
        'try, s = jsondecode(fileread(''%s'')); if ~isempty(s.finished), pause(30); break; end, catch, end, end; ' ...
        'try, V.stop(); catch, end; exit'], src, exp, hb, state);
    system(sprintf('start "arena viewer" "%s" -nodesktop -nosplash -logfile "%s" -r "%s"', ...
        fullfile(matlabroot, 'bin', 'win64', 'MATLAB.exe'), fullfile(exp, 'arena_viewer.log'), r));
    fprintf('  arena viewer starting in a second MATLAB (about 20 s)\n');
elseif running
    fprintf('  arena viewer already running for this session\n');
end

%% 4. the session
opts = [{'cell', tag, 'side', side, 'p1_package', o.p1_package, ...
    'host_config_file', fullfile(exp, 'not edited (checked by jy_session)'), 'show_figures', true, 'show_status', true}, pass];
S = run_session(exp, opts{:});
fprintf('\nSession folder %s\n', exp);
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end

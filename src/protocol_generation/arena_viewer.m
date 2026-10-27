function V = arena_viewer(exp_folder, varargin)
%ARENA_VIEWER Live preview of what the G4 arena should be showing while run_session runs.
%   Run it in a second MATLAB, beside the one that runs run_session (that MATLAB cannot draw during a
%   row: PanelsController.startDisplay waits for the arena without a pause):
%       addpath(genpath('<repo>\src'))
%       arena_viewer('D:\...\JY1008')          % the exp_folder given to run_session
%   or from a command prompt (start it before or after run_session; it waits for the file):
%       matlab -nodesktop -r "addpath(genpath('<repo>\src')); arena_viewer('D:\...\JY1008')"
%   Cost on the rig PC: about 12 % of one core and 1.2 GB (mostly MATLAB itself), measured 27 Oct 2026.
%   It reads EXP_FOLDER\now_playing.json every 0.1 s (run_session rewrites it before every row) and, when
%   the row changes, that package's pattern (Patterns\<id>_*.mat, pattern.Pats) and position function
%   (Functions\<id>_*.mat, pfnparam.func at 500 Hz, 1-based frame numbers). It draws the frame the function
%   gives at the time elapsed in the row (since started_posix_ms; 'started', whole seconds, when a file has
%   no started_posix_ms, shown with a ~): 48 x 192 pixels, gray levels 0..15 (gs_val 4) as the arena shows
%   them; columns 1..16 are shaded (not on screen, the arena shows 17..192). The header gives the row k/n,
%   its label, the time in the row and the frame number; the strip below plots the row's function with a
%   cursor (flashes shorter than the 0.1 s refresh show there even when the image misses them; rows over
%   30 s show 20 s pages).
%   Between rows and during steps without rows it says what the session is doing (session_state.json).
%   Protocol 2 in 'full' mode runs through run_protocol2, which writes no now_playing.json: no preview.
%   It only reads files: it never opens the controller or the G4 Host. Close the window to stop it.
%
%   Options (name, value): 'period_s' 0.1 (refresh), 'lead_s' 0 (added to the elapsed time; the row's
%   start time is written a few controller calls before startDisplay), 'colors' 'gray' | 'green',
%   'log' '' (a CSV of every frame drawn, for tests), 'position' figure position in pixels.
%   V = arena_viewer(...) returns the figure, the timer and V.stop().
%   27 October 2026.

o = struct('period_s', 0.1, 'lead_s', 0, 'colors', 'gray', 'log', '', 'position', [80 80 1200 470]);
for k = 1:2:numel(varargin)
    assert(isfield(o, varargin{k}), 'arena_viewer:option', 'Unknown option %s', varargin{k});
    o.(varargin{k}) = varargin{k+1};
end
exp_folder = char(exp_folder);
np_file = fullfile(exp_folder, 'now_playing.json');
state_file = fullfile(exp_folder, 'session_state.json');
SR = 500;                                               % position function samples per second
cache = struct('file', {}, 'pats', {});                 % the last few patterns, uint8
row = [];                                               % the row being shown (now_playing + loaded files)
last_text = ''; t_state = -inf; step_msg = ''; last_tick_ms = 0;
w = reshape(1:48*192, 48, 192);                         % weights of the log's frame checksum
logfid = -1;
if ~isempty(o.log)
    logfid = fopen(o.log, 'a');
    fprintf(logfid, 'now_ms,package,row,t_row_s,sample,frame,pix_sum,pix_wsum\n');
end

levels = linspace(0, 1, 16).';
if strcmpi(o.colors, 'green'), cmap = [0*levels levels 0*levels]; else, cmap = repmat(levels, 1, 3); end
fig = figure('Name', ['Arena viewer - ' exp_folder], 'NumberTitle', 'off', 'Color', [0.12 0.12 0.12], 'MenuBar', 'none', ...
    'ToolBar', 'none', 'Position', o.position, 'Colormap', cmap, 'CloseRequestFcn', @(~, ~) close_viewer());
fg = [0.92 0.92 0.92];
hd = axes(fig, 'Position', [0.03 0.83 0.94 0.15], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
h.title = text(hd, 0, 0.78, 'waiting for now_playing.json', 'Color', fg, 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'none');
h.sub = text(hd, 0, 0.25, exp_folder, 'Color', [0.7 0.7 0.7], 'FontSize', 11, 'Interpreter', 'none');
h.time = text(hd, 1, 0.78, '', 'Color', fg, 'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'FontName', 'Consolas');
h.frame = text(hd, 1, 0.25, '', 'Color', [0.7 0.7 0.7], 'FontSize', 11, 'HorizontalAlignment', 'right', 'FontName', 'Consolas');
ax = axes(fig, 'Position', [0.03 0.33 0.94 0.47]);
h.img = image(ax, zeros(48, 192, 'uint8'), 'CDataMapping', 'direct');
axis(ax, 'image'); hold(ax, 'on');                     % axes styles after image(), which resets them
set(ax, 'Color', 'k', 'XColor', fg, 'YColor', fg, 'TickDir', 'out', 'Layer', 'top', 'XTick', [1 17 48 96 144 192], 'YTick', [1 24 48], 'FontSize', 9);
patch(ax, [0.5 16.5 16.5 0.5], [0.5 0.5 48.5 48.5], [0.8 0.15 0.15], 'FaceAlpha', 0.45, 'EdgeColor', 'none');
text(ax, 8.5, 24.5, 'not on screen', 'Rotation', 90, 'Color', fg, 'FontSize', 9, 'HorizontalAlignment', 'center');
h.msg = text(ax, 104.5, 24.5, '', 'Color', [1 0.8 0.3], 'FontSize', 13, 'HorizontalAlignment', 'center', 'Interpreter', 'none');
tl = axes(fig, 'Position', [0.045 0.08 0.925 0.17]);
h.fn = stairs(tl, NaN, NaN, 'Color', [0.45 0.75 1]); hold(tl, 'on');
h.cursor = plot(tl, [0 0], [0 1], 'Color', [1 0.8 0.3], 'LineWidth', 1.5);
set(tl, 'Color', [0.05 0.05 0.05], 'XColor', fg, 'YColor', fg, 'FontSize', 9, 'Box', 'on');
xlabel(tl, 'time in row (s)'); ylabel(tl, 'frame');

tm = timer('Name', 'arena_viewer', 'ExecutionMode', 'fixedRate', 'Period', o.period_s, 'BusyMode', 'drop', 'TimerFcn', @(~, ~) tick());
start(tm);
V = struct('figure', fig, 'timer', tm, 'stop', @close_viewer);

    function tick()
        try
            txt = read_text(np_file);
            if ~isempty(txt) && ~strcmp(txt, last_text)      % a new row (each row's file differs in its start time)
                try
                    np = jsondecode(txt);                % fails on a half-written file: try again next tick
                    last_text = txt;
                    load_row(np);
                catch
                end
            end
            now_ms = posixtime(datetime('now', 'TimeZone', 'UTC')) * 1000;   % after a pattern load, which can take a while
            last_tick_ms = now_ms;
            draw(now_ms);
        catch err
            h.msg.String = ['viewer error: ' err.message];
        end
        drawnow limitrate
    end

    function load_row(np)
        r = struct('np', np, 'pkg', resolve_pkg(np.package), 'pats', [], 'func', [], 'err', '', 'start_ms', NaN, 'approx', false);
        if isfield(np, 'started_posix_ms') && ~isempty(np.started_posix_ms)
            r.start_ms = double(np.started_posix_ms);
        else                                    % an older run_session: whole seconds, local time
            r.start_ms = posixtime(datetime(np.started, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss', 'TimeZone', 'local')) * 1000;
            r.approx = true;
        end
        try
            r.pats = load_pattern(find_file(fullfile(r.pkg, 'Patterns'), np.pattern_id));
            F = load(find_file(fullfile(r.pkg, 'Functions'), np.pfn_id), 'pfnparam');
            r.func = double(F.pfnparam.func(:).');
        catch err
            r.err = err.message;
        end
        row = r;
        n = numel(r.func);
        set(h.fn, 'XData', (0:n-1) / SR, 'YData', r.func);
        set(tl, 'XLim', [0 max([np.duration_s, n / SR, 0.1])], 'YLim', [0 max([r.func, 1]) + 1]);
        [~, pname] = fileparts(r.pkg);
        h.title.String = sprintf('row %d/%d   %s', np.row, np.n_rows, np.label);
        h.sub.String = sprintf('package %s   pattern %d (%d frames)   pfn %d   afn %d', pname, np.pattern_id, size(r.pats, 3), np.pfn_id, np.afn_id);
    end

    function draw(now_ms)
        if isempty(row), idle('Waiting for the first row (no now_playing.json yet).'); return, end
        if ~isempty(row.err), idle(['Cannot show this row: ' row.err]); return, end
        t = (now_ms - row.start_ms) / 1000 + o.lead_s;
        dur = row.np.duration_s;
        if t > dur + 0.5
            idle(sprintf('Row %d/%d ended %s ago.', row.np.row, row.np.n_rows, dur_text(t - dur)));
            return
        end
        t = min(max(t, 0), dur);
        s = min(floor(t * SR), round(dur * SR) - 1);    % the last 0.5 s: the row's last frame, as the arena holds it
        s = mod(s, numel(row.func)) + 1;                % the controller repeats a function shorter than the row
        f = row.func(s);
        im = row.pats(:, :, f);
        h.img.CData = im;
        h.msg.String = ''; h.title.Color = fg;
        h.time.String = sprintf('%s%6.2f / %.2f s', tern(row.approx, '~', ''), t, dur);
        h.frame.String = sprintf('frame %d of %d   sample %d', f, size(row.pats, 3), s);
        set(h.cursor, 'XData', [t t], 'YData', tl.YLim);
        if dur > 30                                     % long rows: 20 s pages (redrawing the whole function costs CPU)
            x0 = min(floor(t / 20) * 20, dur - 20);
            if tl.XLim(1) ~= x0, tl.XLim = x0 + [0 20]; end
        end
        if logfid > 0
            [~, pname] = fileparts(row.pkg);
            fprintf(logfid, '%.0f,%s,%d,%.4f,%d,%d,%d,%d\n', now_ms, pname, row.np.row, t, s, f, sum(im(:)), sum(double(im(:)) .* w(:)));
        end
    end

    function idle(msg)
        if now_ms_since(t_state) > 1000                  % session_state.json can be large: read it once a second
            t_state = last_tick_ms; step_msg = step_text();
        end
        h.img.CData = zeros(48, 192, 'uint8'); h.title.Color = [0.55 0.55 0.55];
        h.msg.String = sprintf('%s\n%s', msg, step_msg);
        h.time.String = ''; h.frame.String = '';
        set(h.cursor, 'XData', [NaN NaN]);
    end

    function m = step_text()
        m = '';
        try
            S = jsondecode(read_text(state_file));
            st = S.steps(strcmp({S.steps.status}, 'running'));
            if isempty(st)
                if ~isempty(S.steps) && all(ismember({S.steps.status}, {'done', 'skipped'})), m = 'Session finished.'; end
                return
            end
            switch st(1).id
                case 'p2'
                    m = sprintf(['Step p2 (Protocol 2, full mode) is running: run_protocol2 writes no now_playing.json,\n' ...
                        'so there is no preview for this step. Started %s.'], st(1).t_start);
                case {'p1', 'p2a', 'p2b', 'p3b', 'p4'}
                    m = sprintf('Step %s (%s) is running: waiting for its next row.', st(1).id, st(1).label);
                otherwise
                    m = sprintf('Step %s (%s) is running: no stimulus rows.', st(1).id, st(1).label);
            end
        catch
        end
    end

    function d = now_ms_since(t0)
        d = last_tick_ms - t0;
    end

    function pkg = resolve_pkg(p)
        pkg = char(p);
        if ~isfolder(pkg)                               % the session folder was moved or is seen from another PC
            [~, nm] = fileparts(pkg);
            if isfolder(fullfile(exp_folder, nm)), pkg = fullfile(exp_folder, nm); end
        end
    end

    function P = load_pattern(file)
        j = find(strcmp({cache.file}, file), 1);
        if ~isempty(j), P = cache(j).pats; cache = [cache(j) cache(setdiff(1:numel(cache), j))]; return, end
        S = load(file, 'pattern');
        P = uint8(S.pattern.Pats(:, :, :, 1));          % y position 1: run_session sets position X only
        cache = [struct('file', file, 'pats', P) cache(1:min(end, 5))];
    end

    function close_viewer()
        if isvalid(tm), stop(tm); delete(tm); end
        if logfid > 0, fclose(logfid); logfid = -1; end
        if isvalid(fig), delete(fig); end
    end
end

function file = find_file(folder, id)
d = dir(fullfile(folder, sprintf('%04d_*.mat', id)));
if isempty(d), error('arena_viewer:file', 'No %04d_*.mat in %s', id, folder); end
file = fullfile(d(1).folder, d(1).name);
end

function txt = read_text(file)
txt = '';
fid = fopen(file, 'r');
if fid < 0, return, end
txt = native2unicode(fread(fid, inf, '*uint8').', 'UTF-8');
fclose(fid);
end

function s = dur_text(sec)
if sec < 60, s = sprintf('%.0f s', sec); else, s = sprintf('%.1f min', sec / 60); end
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end

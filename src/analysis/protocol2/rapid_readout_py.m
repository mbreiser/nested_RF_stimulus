function R = rapid_readout_py(session_folder, varargin)
% RAPID_READOUT_PY  Stopgap: run the Python rapid read-out from MATLAB and draw its placement figure in MATLAB.
%
%   R = RAPID_READOUT_PY(SESSION_FOLDER) calls tools/rapid_prot2_fast.py of the
%   ephys-standard-analysis repository on the session folder (it reads the raw
%   TDMS log directly, about 1 s), loads the JSON and CSV tables it writes
%   (<cell>_locations.json, _sites.csv, _lattice.csv, _prot4_frames.csv, ...),
%   prints its report and draws the arena-map panel of the placement figure
%   with PLOT_RAPID_READOUT. The Python figure (<cell>_locations.png, with the
%   two axis profiles) is shown in a second window when it exists.
%
%   Use RAPID_READOUT (pure MATLAB, same tables) when it is available; this
%   wrapper exists in case that one is not ready at the rig.
%
%   Options (name/value):
%     'Repo'    path of the ephys-standard-analysis checkout (default: the
%               environment variable EPHYS_STANDARD_ANALYSIS, else a sibling
%               folder of nested_RF_stimulus named ephys-standard-analysis)
%     'Python'  interpreter (default 'python3'; on the Windows rig PC
%               usually 'python' or 'py -3')
%     'Cell'    cell label (default: the session folder name)
%     'Out'     output folder (default: <session_folder>/rapid_readout;
%               refused under /Volumes/reiserlab)
%     'Metrics' metrics.csv to validate against (passed through as --metrics)
%     'Figure'  draw the MATLAB figure (default true)
%
%   Rig PC requirements for the Python path: Python 3.9 or later with numpy,
%   scipy (scipy.io, for the pattern .mat) and matplotlib (the figure); a
%   checkout of mbreiser/ephys-standard-analysis. No nptdms, no MATLAB toolbox.
%
%   See also RAPID_READOUT, PLOT_RAPID_READOUT, FAST_TDMS_FOLDER2STRUCT.

    here = fileparts(mfilename('fullpath'));
    addpath(here, fullfile(here, 'pipeline'));
    session_folder = char(session_folder);
    if any(session_folder(end) == '/\'), session_folder = session_folder(1:end-1); end
    [~, sid] = fileparts(session_folder);
    repo_default = getenv('EPHYS_STANDARD_ANALYSIS');
    if isempty(repo_default)
        repo_default = fullfile(here, '..', '..', '..', '..', 'ephys-standard-analysis');
    end
    opts = struct('Repo', repo_default, 'Python', 'python3', 'Cell', sid, 'Out', fullfile(session_folder, 'rapid_readout'), 'Metrics', '', 'Figure', true);
    for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
    if startsWith(opts.Out, '/Volumes/reiserlab')
        error('rapid_readout_py:ReadOnlyShare', 'the Janelia share is read-only: give an ''Out'' folder on the local disk');
    end
    script = fullfile(opts.Repo, 'tools', 'rapid_prot2_fast.py');
    assert(isfile(script), 'rapid_readout_py:NoScript', 'rapid_prot2_fast.py not found at %s (set ''Repo'' or EPHYS_STANDARD_ANALYSIS)', script);
    if ~isfolder(opts.Out), mkdir(opts.Out); end

    %% run the Python read-out
    cmd = sprintf('%s "%s" --session "%s" --cell "%s" --out "%s"', opts.Python, script, session_folder, opts.Cell, opts.Out);
    if ~isempty(opts.Metrics), cmd = sprintf('%s --metrics "%s"', cmd, opts.Metrics); end
    t0 = tic;
    [status, txt] = system(cmd);
    t_py = toc(t0);
    fprintf('%s', txt);
    if status ~= 0
        error('rapid_readout_py:PythonFailed', 'the Python read-out failed (exit %d):\n%s', status, txt);
    end
    fprintf('python read-out: %.1f s\n', t_py);

    %% load its outputs
    tag = opts.Cell;
    J = jsondecode(fileread(fullfile(opts.Out, [tag '_locations.json'])));
    sites = table2struct(readtable(fullfile(opts.Out, [tag '_sites.csv']), 'TextType', 'char'))';
    kept = table2struct(readtable(fullfile(opts.Out, [tag '_lattice.csv']), 'TextType', 'char'))';
    p4 = table2struct(readtable(fullfile(opts.Out, [tag '_prot4_frames.csv']), 'TextType', 'char'))';
    R = struct('cell', tag, 'sid', sid, 'json', J, 'sites', sites, 'lattice', kept, 'prot4_frames', p4, 'python_stdout', txt, 'python_seconds', t_py, 'out_dir', opts.Out, 'figure', []);

    %% figure (arena map; the profiles live only in the Python PNG)
    if ~opts.Figure, return; end
    if isfield(J, 'pd_placement'), pd = J.pd_placement; else, pd = J.pd_deg; end        % placement PD (revision of 26 Oct 15:05) or the snapped PD
    e_pd = [cosd(pd), sind(pd)]; e_or = [cosd(pd - 90), sind(pd - 90)];
    center = J.rf_center_col_row_1based(:)';
    [cx, cy, go] = bar_geometry(session_folder);
    if isfield(J, 'bar_sets'), g_pd = J.bar_sets.g_perp + 1; g_or = J.bar_sets.g_par + 1;
    else, [~, g_pd] = min(abs(mod(go - (pd + 90) + 90, 180) - 90)); [~, g_or] = min(abs(mod(go - pd + 90, 180) - 90)); end
    P = struct('cell', tag, 'pd', pd, 'pd_vector', NaN, 'center', center, 'center_xy', [center(1) - 1, -(center(2) - 1)], 'e_pd', e_pd, 'e_or', e_or, ...
               'bar_cx', cx + 1, 'bar_cy', cy + 1, 'g_pd', g_pd - 1, 'g_or', g_or - 1, 'lattice', kept, 'sites', sites, 'p4_frames', p4, ...
               'window', struct('s_lead', J.window_px.s_lead, 's_trail', J.window_px.s_trail, 'v_half', J.window_px.v_half));
    if isfield(J, 'pd_vector_28'), P.pd_vector = J.pd_vector_28; P.pd_vector_56 = J.pd_vector_56; P.pd_snapped = J.pd_deg;
    elseif isfield(J, 'pd_vector_deg'), P.pd_vector = J.pd_vector_deg; end
    R.figure = plot_rapid_readout(P);
    png = fullfile(opts.Out, [tag '_locations.png']);
    if isfile(png)
        f2 = figure('Name', ['Python figure ' tag], 'Color', 'w', 'Position', [60 60 1500 500]);
        imshow(imread(png), 'Parent', axes(f2), 'Border', 'tight');
    end
end

function [cx, cy, go] = bar_geometry(sdir)
% Centers (0-based) and orientations of the 88 bars of the FLASHES pattern, as rapid_readout computes them.
    d = dir(fullfile(sdir, 'Patterns', '*FLASHES*.mat'));
    assert(numel(d) == 1, 'rapid_readout_py:Flashes', 'expected one FLASHES pattern, found %d', numel(d));
    pz = load(fullfile(d(1).folder, d(1).name), 'pattern'); Pm = single(pz.pattern.Pats);
    if max(Pm(:)) > 4, W = min(max((Pm - 4) / 11, 0), 1); else, W = min(max((4 - Pm) / 4, 0), 1); end
    [H, Wd, ~] = size(W);
    fp = double(reshape(W(:, :, 2:89), H * Wd, 88))';
    [rr, cc] = ndgrid(0:H-1, 0:Wd-1); X = cc(:)'; Y = rr(:)';
    a = sum(fp, 2); cx = (fp * X')' ./ a'; cy = (fp * Y')' ./ a';
    go = zeros(1, 8);
    for g = 0:7
        pts = [];
        for k = 1:11, m = fp(g * 11 + k, :) > 0.5; pts = [pts; X(m)', -Y(m)']; end %#ok<AGROW>
        pts = pts - mean(pts, 1); [~, ~, V] = svd(pts, 'econ'); u = V(:, 1);
        go(g + 1) = mod(atan2d(u(2), u(1)), 180);
    end
end

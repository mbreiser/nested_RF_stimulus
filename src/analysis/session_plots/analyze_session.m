function R = analyze_session(run_dirs, varargin)
%ANALYZE_SESSION  The per-protocol figures and metrics of one or more RUN_SESSION folders, and their comparison.
%
%   R = ANALYZE_SESSION(RUN_DIR) for a session folder written by RUN_SESSION (Auto_T4_v1): draws the
%   Prot3b (ANALYZE_PROT3B) and Protocol 4 (ANALYZE_P4) figures from the recorded p3b/ and p4/ logs, and
%   collects the session's numbers in one table row: Vm per step (the step overviews), the Protocol 1 ON
%   peak, the placement PD with L and DSI (read-out 1), the center 80 ms peak (read-out 2), R_gray and the
%   a-priori dR read-outs, the Protocol 4 map numbers and the QC counts. Figures already made by the
%   session (Protocol 1 maps, stage 1, read-out 2 locations, step overviews) are listed, not redrawn.
%   R = ANALYZE_SESSION({RUN_DIR1, RUN_DIR2, ...}) does the same for each and adds a comparison figure
%   (the numbers side by side, Vm per step against the time of day: e.g. two runs on the same cell).
%
%   Options: 'Out' (folder of the PNGs and <cell>_session_metrics.json; default <run_dir>/plots, the
%   comparison in the first run's), 'Figure' (false: draw off screen), 'Labels' (cellstr; default the
%   folder names). Folders renamed after the run are handled: paths in session_state.json under its
%   exp_folder are moved to RUN_DIR.
%   R: .metrics (table, one row per session), .p3b / .p4 (the analysis structs), .pngs.
%   28 October 2026.

opts = struct('Out', '', 'Figure', false, 'Labels', {{}});
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
run_dirs = cellstr(run_dirs);
n = numel(run_dirs);
if isempty(opts.Labels), opts.Labels = cellfun(@(d) last_part(d), run_dirs, 'UniformOutput', false); end
M = table(); R = struct('p3b', {cell(1, n)}, 'p4', {cell(1, n)}, 'pngs', {{}}, 'steps', {cell(1, n)});
for i = 1:n
    rd = char(run_dirs{i});
    S = jsondecode(fileread(fullfile(rd, 'session_state.json'))); if iscell(S.steps), S.steps = [S.steps{:}]; end
    mv = @(p) relocate(p, S.exp_folder, rd);
    out = opts.Out; if isempty(out), out = fullfile(rd, 'plots'); end
    if ~isfolder(out), mkdir(out); end
    lab = opts.Labels{i};
    m = struct('label', string(lab), 'cell', string(S.cell), 'started', string(S.created));
    % Vm per step (overviews)
    steps = {'dark30', 'p1', 'p2a', 'p2b', 'p3b', 'p4'}; st = table();
    for s = steps
        if isfield(S.out, s{1}) && isfield(S.out.(s{1}), 'overview')
            o = S.out.(s{1}).overview; k = find(strcmp({S.steps.id}, s{1}), 1);
            st = [st; table(string(s{1}), datetime(S.steps(k).t_start, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss'), o.vm_median_mV, o.vm_drift_mV, ...
                'VariableNames', {'step', 't_start', 'vm_median_mV', 'vm_drift_mV'})]; %#ok<AGROW>
            m.(['vm_' s{1} '_mV']) = o.vm_median_mV;
        else
            m.(['vm_' s{1} '_mV']) = NaN;
        end
    end
    R.steps{i} = st;
    % Protocol 1, read-outs
    m.p1_peak_mV = getf(getf(S.out, 'p1_peak', struct()), 'peak_mV', NaN);
    j1 = fullfile(rd, 'readout', [S.cell '_stage1.json']); J1 = struct();
    if isfile(j1), J1 = jsondecode(fileread(j1)); end
    m.pd_placement_deg = getf(J1, 'pd_placement', NaN); m.L28 = getf(J1, 'vector_length_28', NaN); m.L56 = getf(J1, 'vector_length_56', NaN);
    m.dsi28 = getf(J1, 'dsi_28', NaN); m.dsi56 = getf(J1, 'dsi_56', NaN);
    jl = fullfile(rd, 'readout', [S.cell '_locations.json']); m.center_80ms_peak_mV = NaN; m.static_vm_mV = NaN;
    if isfile(jl)
        JL = jsondecode(fileread(jl)); m.static_vm_mV = getf(JL, 'static_vm_mV', NaN);
        if isfield(JL, 'sites') && ~isempty(JL.sites), s1 = JL.sites; if iscell(s1), s1 = s1{1}; end, m.center_80ms_peak_mV = getf(s1(1), 'expected_ON_80ms_peak_mV', NaN); end
    end
    % Prot3b
    f3 = fullfile(rd, 'p3b');
    if has_log(f3)
        png = fullfile(out, sprintf('%s_p3b_analysis.png', lab));
        A3 = analyze_prot3b(f3, 'Cell', sprintf('%s (%s)', lab, S.cell), 'Png', png, 'Figure', opts.Figure, 'Name', ['Prot3b ' lab]);
        R.p3b{i} = rmfield(A3, 'figure'); R.pngs{end+1} = png;
        if ~opts.Figure, close(A3.figure); end
        m.R_gray_GOhm = A3.summary.R_gray_GOhm; m.pairs_kept = A3.summary.pairs_kept;
        m.dR_ON_trailing = A3.summary.ON_trailing_dR_300_500; m.dR_OFF_center = A3.summary.OFF_center_dR_100_200; m.dR_ON_center = A3.summary.ON_center_dR_150_250;
        x = A3.summary.mi9_dR_off_s0; if isempty(x), x = NaN; end, m.mi9_dR_off_s0 = x;
    else
        [m.R_gray_GOhm, m.pairs_kept, m.dR_ON_trailing, m.dR_OFF_center, m.dR_ON_center, m.mi9_dR_off_s0] = deal(NaN);
    end
    % Protocol 4
    f4 = fullfile(rd, 'p4');
    if has_log(f4)
        png = fullfile(out, sprintf('%s_p4_analysis.png', lab));
        png2 = fullfile(out, sprintf('%s_p4_detail.png', lab));
        A4 = analyze_p4(f4, 'Cell', sprintf('%s (%s)', lab, S.cell), 'Png', png, 'DetailPng', png2, 'Figure', opts.Figure, 'Name', ['Protocol 4 ' lab]);
        R.p4{i} = rmfield(A4, {'figure', 'detail_figure'}); R.pngs = [R.pngs, {png, png2}];
        if ~opts.Figure, close([A4.figure, A4.detail_figure]); end
        m.p4_lattice_ON_peak_mV = A4.summary.lattice_ON_peak_mV; m.p4_gauss_sigma_s_px = A4.summary.gauss_sigma_s_px; m.p4_gauss_sigma_v_px = A4.summary.gauss_sigma_v_px;
        m.p4_trailing_extent_px = A4.summary.trailing_extent_px; m.p4_late_min_mV = A4.summary.lattice_ON_late_min_mV; m.p4_passes_kept = A4.summary.passes_kept;
    else
        [m.p4_lattice_ON_peak_mV, m.p4_gauss_sigma_s_px, m.p4_gauss_sigma_v_px, m.p4_trailing_extent_px, m.p4_late_min_mV, m.p4_passes_kept] = deal(NaN);
    end
    for f = {'_p1_maps.png', '_stage1.png', '_locations.png'}
        p = fullfile(rd, 'readout', [S.cell f{1}]); if isfile(p), R.pngs{end+1} = p; end
    end
    write_json(fullfile(out, sprintf('%s_session_metrics.json', lab)), m);
    M = [M; struct2table(m, 'AsArray', true)]; %#ok<AGROW>
end
R.metrics = M;
if n > 1
    out = opts.Out; if isempty(out), out = fullfile(char(run_dirs{1}), 'plots'); end
    png = fullfile(out, sprintf('compare_%s.png', strjoin(opts.Labels, '_vs_')));
    fig = compare_figure(M, R.steps, opts); exportgraphics(fig, png, 'Resolution', 120); R.pngs{end+1} = png;
    if ~opts.Figure, close(fig); end
end
end

% ------------------------------------------------------------------------------------------------
function fig = compare_figure(M, steps, opts)
n = height(M); col = lines(max(n, 2)); col(1:2, :) = [0.1 0.3 0.75; 0.85 0.35 0.1];
fig = session_figure('Session comparison', opts.Figure, [80 80 1500 820]);
tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Session comparison: %s', strjoin(cellstr(M.label), '  vs  ')), 'FontWeight', 'bold', 'Interpreter', 'none');
% Vm per step against the time of day
ax = nexttile(tl, [1 2]); hold(ax, 'on');
for i = 1:n
    s = steps{i}; if isempty(s), continue; end
    plot(ax, s.t_start, s.vm_median_mV, '-o', 'Color', col(i, :), 'MarkerFaceColor', col(i, :), 'LineWidth', 1.5, 'DisplayName', M.label(i));
    text(ax, s.t_start, s.vm_median_mV + 1.2, cellstr(s.step), 'FontSize', 7, 'Color', col(i, :), 'HorizontalAlignment', 'center');
end
ylabel(ax, 'median Vm per step (mV)'); xlabel(ax, 'time of day'); legend(ax, 'Location', 'best', 'Interpreter', 'none'); box(ax, 'on'); grid(ax, 'on');
title(ax, 'Vm through each session (step overviews)');
panels = { ...
    'responses (mV)', {'p1_peak_mV', 'Protocol 1 ON peak'; 'center_80ms_peak_mV', 'center bar 80 ms'; 'p4_lattice_ON_peak_mV', '4-px lattice ON'; 'p4_late_min_mV', '4-px late min'}; ...
    'direction', {'L28', 'L 28'; 'L56', 'L 56'; 'dsi28', 'DSI 28'; 'dsi56', 'DSI 56'}; ...
    'resistance (G\Omega)', {'R_gray_GOhm', 'R_{gray}'; 'dR_ON_trailing', '\DeltaR ON trail'; 'dR_OFF_center', '\DeltaR OFF ctr'; 'dR_ON_center', '\DeltaR ON ctr'; 'mi9_dR_off_s0', '\DeltaR_{OFF} Mi9'}; ...
    'RF size (px)', {'p4_gauss_sigma_s_px', '\sigma_s'; 'p4_gauss_sigma_v_px', '\sigma_v'; 'p4_trailing_extent_px', 'trailing hyperpol.'}; ...
    'QC', {'pairs_kept', 'Prot3b pairs kept (/8)'; 'p4_passes_kept', 'P4 passes kept (/6)'}};
for j = 1:size(panels, 1)
    ax = nexttile(tl); hold(ax, 'on');
    f = panels{j, 2}; Y = nan(size(f, 1), n);
    for a = 1:size(f, 1), Y(a, :) = M.(f{a, 1}).'; end
    b = bar(ax, Y, 'grouped', 'EdgeColor', 'none');
    for i = 1:n, b(i).FaceColor = col(i, :); b(i).DisplayName = M.label(i); end
    for i = 1:n
        xt = b(i).XEndPoints; yt = b(i).YEndPoints;
        text(ax, xt, yt, compose('%.2g', Y(:, i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', tern(all(yt >= 0), 'bottom', 'middle'), 'FontSize', 7);
    end
    xticks(ax, 1:size(f, 1)); xticklabels(ax, f(:, 2)); xtickangle(ax, 20); yline(ax, 0, '-', 'HandleVisibility', 'off');
    title(ax, panels{j, 1}); box(ax, 'on');
    if j == 1, legend(ax, 'Location', 'best', 'Interpreter', 'none', 'FontSize', 7); end
end
ax = nexttile(tl); axis(ax, 'off');
txt = {'\bfpath into each number\rm', 'Vm: step overview medians', 'P1: rapid\_p1\_peak (6-px ON)', 'L, DSI: read-out 1 (stage 1)', ...
    'center 80 ms: read-out 2', 'R, \DeltaR: analyze\_prot3b (kept pairs)', 'P4: analyze\_p4 (kept passes)'};
text(ax, 0, 1, txt, 'VerticalAlignment', 'top', 'FontSize', 8);
drawnow;
end

function p = relocate(p, old_root, new_root)
p = char(p); if isfolder(p) || isfile(p), return; end
if startsWith(lower(p), lower(char(old_root))), p = [char(new_root) p(numel(char(old_root)) + 1:end)]; end
end

function tf = has_log(f)
tf = isfolder(f) && ~isempty(dir(fullfile(f, 'Log Files', '**', '*.tdms'))) && isfile(fullfile(f, 'Reports', 'event_ledger.csv'));
end

function s = last_part(d)
[~, s] = fileparts(char(d)); if isempty(s), [~, s] = fileparts(fileparts(char(d))); end
end

function v = getf(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

function write_json(f, x)
fid = fopen(f, 'w'); fprintf(fid, '%s', jsonencode(x, 'PrettyPrint', true)); fclose(fid);
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end

function html_file = session_results_page(run_dirs, out_dir, varargin)
%SESSION_RESULTS_PAGE  One HTML page of the standard result figures of one or more RUN_SESSION folders.
%
%   F = SESSION_RESULTS_PAGE(RUN_DIRS, OUT_DIR) runs ANALYZE_SESSION on the session folders (the Prot3b and
%   Protocol 4 figures, the metrics and, for several folders, the comparison figure), copies every figure
%   into OUT_DIR/img and writes OUT_DIR/index.html:
%     numbers      the main read-outs of each session side by side (and run 2 - run 1 for two sessions)
%     comparison   ANALYZE_SESSION's comparison figure
%     protocols    one section per protocol, the sessions in columns: recording quality (dark, step Vm),
%                  Protocol 1 (map), Protocol 2 stage 1 (direction), Protocol 2 read-out 2 (center, bar
%                  profiles), Prot3b (input resistance), Protocol 4 (map); step overviews folded under each
%   Images open full size on click. Options: 'Labels' (cellstr, default the folder names), 'Title',
%   'Notes' (an HTML fragment file, shown first as its own section, e.g. the assessment against the
%   protocol's objectives), 'Embed' (true: the figures are written into index.html as data URIs, so the
%   page is one self-contained file that previews and mails intact; false: linked from img/).
%   The session folders are only read. 28 October 2026.

opts = struct('Labels', {{}}, 'Title', 'Auto_T4_v1 results', 'Notes', '', 'Embed', true);
for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
run_dirs = cellstr(run_dirs); n = numel(run_dirs);
if isempty(opts.Labels), opts.Labels = cellfun(@last_part, run_dirs, 'UniformOutput', false); end
img = fullfile(out_dir, 'img'); if ~isfolder(img), mkdir(img); end
R = analyze_session(run_dirs, 'Out', img, 'Labels', opts.Labels);
M = R.metrics;

% figures per session, by key
keys = {'p1_maps', 'stage1', 'locations', 'p3b', 'p4', 'p4_detail', 'ov_dark30', 'ov_p1', 'ov_p2a', 'ov_p2b', 'ov_p3b', 'ov_p4'};
F = repmat({''}, n, numel(keys));
for i = 1:n
    rd = char(run_dirs{i}); lab = opts.Labels{i}; S = jsondecode(fileread(fullfile(rd, 'session_state.json')));
    src = {fullfile(rd, 'readout', [S.cell '_p1_maps.png']), fullfile(rd, 'readout', [S.cell '_stage1.png']), fullfile(rd, 'readout', [S.cell '_locations.png']), ...
           fullfile(img, [lab '_p3b_analysis.png']), fullfile(img, [lab '_p4_analysis.png']), fullfile(img, [lab '_p4_detail.png'])};
    for s = {'dark30', 'p1', 'p2a', 'p2b', 'p3b', 'p4'}
        d = dir(fullfile(rd, 'plots', sprintf('*_%s_overview.png', s{1}))); if isempty(d), src{end+1} = ''; else, src{end+1} = fullfile(d(1).folder, d(1).name); end %#ok<AGROW>
    end
    for k = 1:numel(keys)
        if isempty(src{k}) || ~isfile(src{k}), continue; end
        dst = fullfile(img, sprintf('%s_%s.png', lab, keys{k}));
        if ~strcmp(src{k}, dst), copyfile(src{k}, dst); end
        F{i, k} = ['img/' sprintf('%s_%s.png', lab, keys{k})];
    end
end
cmp = ''; if n > 1, [~, nm, ex] = fileparts(R.pngs{end}); cmp = ['img/' nm ex]; end

% numbers
rows = { ...
    'Recording', '', ''; ...
    'Vm, dark (break-in, 0 pA)', 'vm_dark30_mV', '%.1f mV'; 'Vm, Protocol 1', 'vm_p1_mV', '%.1f mV'; 'Vm, Protocol 2 stage 1', 'vm_p2a_mV', '%.1f mV'; ...
    'Vm, Prot3b', 'vm_p3b_mV', '%.1f mV'; 'Vm, Protocol 4', 'vm_p4_mV', '%.1f mV'; 'R_gray (Prot3b, pre-flash)', 'R_gray_GOhm', '%.2f G&Omega;'; ...
    'Protocol 1', '', ''; 'ON peak (6-px)', 'p1_peak_mV', '%.1f mV'; ...
    'Protocol 2', '', ''; 'placement PD', 'pd_placement_deg', '%.1f&deg;'; 'L 28 / 56 deg/s', {'L28', 'L56'}, '%.2f / %.2f'; 'DSI 28 / 56 deg/s', {'dsi28', 'dsi56'}, '%.2f / %.2f'; ...
    'center bar, 80 ms peak', 'center_80ms_peak_mV', '%.1f mV'; ...
    'Prot3b', '', ''; 'pairs kept (of 8)', 'pairs_kept', '%d'; '&Delta;R ON trailing 300-500 ms', 'dR_ON_trailing', '%.2f G&Omega;'; ...
    '&Delta;R OFF center 100-200 ms', 'dR_OFF_center', '%.2f G&Omega;'; '&Delta;R ON center 150-250 ms', 'dR_ON_center', '%.2f G&Omega;'; '&Delta;R<sub>OFF</sub> Mi9 column', 'mi9_dR_off_s0', '%.2f G&Omega;'; ...
    'Protocol 4', '', ''; 'passes kept (of 6)', 'p4_passes_kept', '%d'; '4-px lattice ON peak', 'p4_lattice_ON_peak_mV', '%.1f mV'; ...
    'Gaussian &sigma; s / v', {'p4_gauss_sigma_s_px', 'p4_gauss_sigma_v_px'}, '%.1f / %.1f px'; 'late (150-500 ms) minimum', 'p4_late_min_mV', '%.1f mV'; ...
    'trailing hyperpolarization to', 'p4_trailing_extent_px', 's = %+.1f px'};
tb = sprintf('<table class="num"><thead><tr><th></th>%s%s</tr></thead><tbody>', strjoin(cellfun(@(l) sprintf('<th>%s</th>', esc(l)), opts.Labels, 'UniformOutput', false), ''), ...
    tern(n == 2, sprintf('<th>%s &minus; %s</th>', esc(opts.Labels{2}), esc(opts.Labels{1})), ''));
for r = 1:size(rows, 1)
    if isempty(rows{r, 2}), tb = [tb sprintf('<tr class="grp"><td colspan="%d">%s</td></tr>', n + 1 + (n == 2), rows{r, 1})]; continue; end %#ok<AGROW>
    f = cellstr(rows{r, 2}); cells = '';
    for i = 1:n, v = cellfun(@(x) M.(x)(i), f); cells = [cells '<td>' sprintf(rows{r, 3}, v) '</td>']; end %#ok<AGROW>
    if n == 2
        v = cellfun(@(x) M.(x)(2) - M.(x)(1), f); fm = strrep(strrep(rows{r, 3}, '%.', '%+.'), '%d', '%+d'); fm = strrep(fm, '%++', '%+');
        cells = [cells '<td class="d">' sprintf(fm, v) '</td>']; %#ok<AGROW>
    end
    tb = [tb sprintf('<tr><td>%s</td>%s</tr>', rows{r, 1}, cells)]; %#ok<AGROW>
end
tb = [tb '</tbody></table>'];

sections = { ...
    'Recording quality', 'Dark period at break-in and Vm through each step (min/max envelope, per-row medians, frames, current).', {'ov_dark30', 'ov_p4'}, {'ov_p1', 'ov_p2a', 'ov_p2b', 'ov_p3b'}; ...
    'Protocol 1: coarse map', '6-px ON and OFF maps (dV 50-190 ms); + the ON peak that aims Protocol 2.', {'p1_maps'}, {'ov_p1'}; ...
    'Protocol 2 stage 1: direction', '28 and 56 deg/s tuning, vector PDs, placement PD; the 28 deg/s trace per direction.', {'stage1'}, {'ov_p2a'}; ...
    'Protocol 2 stage 2 / read-out 2: center', 'The two placement-frame bar sets (80 ms peak and post-flash minimum), the RF center, the Prot3b sites and the Protocol 4 window.', {'locations'}, {'ov_p2b'}; ...
    'Prot3b: input resistance', '&Delta;R(t) = R(t) &minus; R<sub>gray</sub> from the &plusmn;0.5 pA pairs per region; 0 pA responses; the Mi9 holding series; pair QC on the Vm record.', {'p3b'}, {'ov_p3b'}; ...
    'Protocol 4: combined map', ['4-px lattice ON/OFF early and late, Gaussian fit, trailing hyperpolarization, 3-px core, 2-px line scans, OFF series, pass QC. ' ...
        'Below it, the detail figure: position &times; time images, maps on an asymmetric scale (blue saturates at each panel''s minimum), and the three square sizes on common profiles.'], {'p4', 'p4_detail'}, {'ov_p4'}};
body = '';
for s = 1:size(sections, 1)
    body = [body sprintf('<section id="s%d"><h2>%s</h2><p class="cap">%s</p>', s, sections{s, 1}, sections{s, 2})]; %#ok<AGROW>
    body = [body grid_html(sections{s, 3}, F, keys, opts.Labels)]; %#ok<AGROW>
    ex = sections{s, 4};
    if ~isempty(ex) && s > 1
        body = [body '<details><summary>step overview' tern(numel(ex) > 1, 's', '') '</summary>' grid_html(ex, F, keys, opts.Labels) '</details>']; %#ok<AGROW>
    elseif s == 1
        body = [body '<details><summary>the other step overviews</summary>' grid_html(ex, F, keys, opts.Labels) '</details>']; %#ok<AGROW>
    end
    body = [body '</section>']; %#ok<AGROW>
end
nav = strjoin(arrayfun(@(s) sprintf('<a href="#s%d">%s</a>', s, strtok(sections{s, 1}, ':')), 1:size(sections, 1), 'UniformOutput', false), ' ');
notes = '';
if ~isempty(opts.Notes)
    notes = sprintf('<section id="notes" class="notes">%s</section>', fileread(opts.Notes));
    nav = ['<a href="#notes">Objectives</a> ' nav];
end
cmp_html = ''; if ~isempty(cmp), cmp_html = sprintf('<section id="cmp"><h2>Comparison</h2><p class="cap">The read-outs side by side; Vm per step against the time of day.</p><img class="wide" src="%s"></section>', cmp); end
css = ['body{font:14px/1.45 -apple-system,Segoe UI,Arial,sans-serif;margin:0;background:#f6f7f9;color:#1d2330}' ...
    'header{background:#1d2330;color:#fff;padding:14px 24px;position:sticky;top:0;z-index:2}header h1{margin:0;font-size:19px}' ...
    'header .sub{opacity:.75;font-size:12px}nav a{color:#9cc3ff;margin-right:14px;font-size:13px;text-decoration:none}' ...
    'main{padding:12px 24px 40px;max-width:1900px}section{background:#fff;border:1px solid #e3e6eb;border-radius:8px;padding:12px 16px;margin:14px 0}' ...
    'h2{font-size:16px;margin:2px 0 4px}.cap{margin:0 0 10px;color:#5a6272;font-size:13px}' ...
    '.grid{display:grid;gap:10px}.cell h3{font-size:13px;margin:0 0 4px;color:#33415c}.cell img,img.wide{width:100%;border:1px solid #eceef2;border-radius:4px;background:#fff}' ...
    '.missing{color:#a0a6b1;font-style:italic;padding:30px;text-align:center;border:1px dashed #d5d9e0;border-radius:4px}' ...
    'details{margin-top:8px}summary{cursor:pointer;color:#3d5a99;font-size:13px}' ...
    'table.num{border-collapse:collapse;font-size:13px;min-width:520px}table.num th,table.num td{padding:3px 12px;text-align:right;border-bottom:1px solid #eef0f3}' ...
    'table.num td:first-child,table.num th:first-child{text-align:left}table.num tr.grp td{font-weight:600;background:#f1f4f9;text-align:left;color:#33415c}' ...
    'td.d{color:#6b4ca8}.top{display:grid;grid-template-columns:minmax(520px,640px) 1fr;gap:14px;align-items:start}@media(max-width:1100px){.top{grid-template-columns:1fr}}' ...
    '.cell img,img.wide{cursor:zoom-in}#lb{display:none;position:fixed;inset:0;background:rgba(20,24,32,.88);z-index:9;align-items:center;justify-content:center;cursor:zoom-out}' ...
    '#lb img{max-width:97vw;max-height:95vh;background:#fff;border-radius:4px}' ...
    '.notes table{border-collapse:collapse;font-size:13px;width:100%;table-layout:fixed}.notes th,.notes td{padding:5px 9px;border-bottom:1px solid #eef0f3;text-align:left;vertical-align:top}' ...
    '.notes thead th:nth-child(1){width:19%}.notes thead th:nth-child(2){width:17%}.notes thead th:nth-child(3){width:44%}.notes thead th:nth-child(4){width:20%}' ...
    '.notes th{background:#f1f4f9;color:#33415c}.notes td.v{font-weight:600}.ok{color:#1d7a3a}.part{color:#a86400}.no{color:#b3261e}.open{color:#5a6272}' ...
    '.notes h3{font-size:14px;margin:14px 0 4px}.notes p,.notes li{font-size:13px}'];
lb = ['<div id="lb" onclick="this.style.display=''none''"><img alt=""></div><script>' ...
    'document.querySelectorAll(''.cell img,img.wide'').forEach(function(i){i.onclick=function(){var b=document.getElementById(''lb'');' ...
    'b.firstChild.src=i.src;b.style.display=''flex'';};});</script>'];
html = sprintf(['<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>%s</title><style>%s</style></head><body>' ...
    '<header><h1>%s</h1><div class="sub">%s &middot; generated %s by session_results_page.m</div><nav>%s <a href="#cmp">Comparison</a></nav></header><main>' ...
    '%s<div class="top"><section><h2>Numbers</h2>%s</section>%s</div>%s</main>%s</body></html>'], ...
    esc(opts.Title), css, esc(opts.Title), esc(strjoin(cellfun(@(d) char(d), run_dirs, 'UniformOutput', false), '  |  ')), char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), nav, notes, tb, cmp_html, body, lb);
if opts.Embed                                                             % figures as data URIs: one self-contained file
    for p = unique(regexp(html, 'img/[^"]+\.png', 'match'))
        fid = fopen(fullfile(out_dir, p{1}), 'r'); b = fread(fid, Inf, '*uint8'); fclose(fid);
        html = strrep(html, ['"' p{1} '"'], ['"data:image/png;base64,' matlab.net.base64encode(b) '"']);
    end
end
html_file = fullfile(out_dir, 'index.html');
fid = fopen(html_file, 'w', 'n', 'UTF-8'); fprintf(fid, '%s', html); fclose(fid);
end

function h = grid_html(ks, F, keys, labels)
n = size(F, 1); h = '';
for k = ks
    j = find(strcmp(keys, k{1}));
    h = [h sprintf('<div class="grid" style="grid-template-columns:repeat(%d,1fr)">', n)]; %#ok<AGROW>
    for i = 1:n
        if isempty(F{i, j}), im = '<div class="missing">not recorded</div>'; else, im = sprintf('<img src="%s">', F{i, j}); end
        h = [h sprintf('<div class="cell"><h3>%s</h3>%s</div>', esc(labels{i}), im)]; %#ok<AGROW>
    end
    h = [h '</div>']; %#ok<AGROW>
end
end

function s = esc(s)
s = strrep(strrep(strrep(char(s), '&', '&amp;'), '<', '&lt;'), '>', '&gt;');
end

function s = last_part(d)
[~, s] = fileparts(char(d)); if isempty(s), [~, s] = fileparts(fileparts(char(d))); end
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end

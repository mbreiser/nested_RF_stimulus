function R = rapid_p1_peak(p1_folder, varargin)
% RAPID_P1_PEAK  Peak of the Protocol 1 6-px square map (ON and OFF) from the raw TDMS log, in about a second.
%
%   R = RAPID_P1_PEAK(P1_FOLDER) reads the Protocol 1 log with
%   FAST_TDMS_FOLDER2STRUCT and computes the ON and OFF maps exactly as
%   datasets/20261022_jy_rnai_rf/protocol1_maps.py (ephys-standard-analysis;
%   itself equal to qc.py's Protocol 1 maps): for every flash, taken from the
%   frame log (each change of the frame index to a nonzero frame), dV = mean
%   Vm over 50..190 ms minus mean Vm over -80..-20 ms from the frame onset;
%   the map is the mean over the repeats per frame; the single-presentation
%   noise is the RMS over frames of the SD across repeats; SNR vs the mean
%   map = peak / (noise / sqrt(n_repeats)). The peak is the frame with the
%   largest dV per polarity (ON frames: squares brighter than the background
%   4; OFF: darker).
%
%   R fields: frame (ON peak frame, 1-based, as the Protocol 1 'Frame' noted in
%   the experiment log), xy (its center, 1-based [col row]), peak_mV, polarity
%   ('ON'), frame12 = [ON OFF] peak frames, peak12_mV = [ON OFF] peak dV, ON and
%   OFF (frame_1based, xy, peak_mV, neg_mV, snr_single, snr_mean,
%   noise_single_mV), rows (the map table), n_flashes, n_repeats, layout,
%   pattern, timing, out_dir. Files in 'Out': <cell>_p1_summary.json,
%   <cell>_p1_rf_map.csv (columns as protocol1_maps.py: frame_1based,
%   polarity, x_center_1based, y_center_1based, n_repeats, dV_mean_mV,
%   dV_sd_across_repeats_mV) and <cell>_p1_maps.png.
%
%   Options: 'Out' (default <p1_folder>/rapid_readout), 'Cell' (default: the
%   folder name), 'Pattern' (path of the 6-px square pattern .mat; default: the
%   first *square*.mat under Patterns/ of P1_FOLDER, its parent or
%   grandparent, where the Protocol 1 root keeps it), 'Figure' (true),
%   'SavePng' (true).
%
%   Both channel layouts are handled (historic ADC0 frame index / ADC1 Vm;
%   October 2026 ADC1 Vm / ADC2 current); the frame onsets come from the frame
%   log in both. When the log holds several patterns (e.g. a Protocol 2 session
%   used for testing), only the flashes shown while the square pattern was the
%   active Set Pattern ID are used; on a Protocol 1 session this changes
%   nothing. Validated against the reference arithmetic on 25 TUTL sessions
%   and against p1_summary.csv on five October 2026 Protocol 1 folders (peak
%   frame, x, y identical; mV within the CSV's rounding); record in
%   docs/notes/MATLAB_FAST_READOUT_20261026.md.
%
%   See also RAPID_READOUT, FAST_TDMS_FOLDER2STRUCT.

    t_all = tic;
    here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here, 'pipeline'));
    p1_folder = char(p1_folder); if any(p1_folder(end) == '/\'), p1_folder = p1_folder(1:end-1); end
    [~, tag] = fileparts(p1_folder);
    opts = struct('Out', fullfile(p1_folder, 'rapid_readout'), 'Cell', tag, 'Pattern', '', 'Figure', true, 'SavePng', true);
    for k = 1:2:numel(varargin), opts.(varargin{k}) = varargin{k+1}; end
    if startsWith(opts.Out, '/Volumes/reiserlab'), error('rapid_p1_peak:ReadOnlyShare', 'the Janelia share is read-only: give an ''Out'' folder on the local disk'); end
    BG = 4; FS = 10000; BASE_MS = [-80 -20]; RESP_MS = [50 190];

    %% pattern
    pat_path = opts.Pattern;
    if isempty(pat_path)
        cands = {p1_folder, fileparts(p1_folder), fileparts(fileparts(p1_folder))};
        for c = 1:numel(cands)
            d = dir(fullfile(cands{c}, 'Patterns', '*square*.mat'));
            if ~isempty(d), [~, i] = sort({d.name}); pat_path = fullfile(d(i(1)).folder, d(i(1)).name); break; end
        end
    end
    assert(~isempty(pat_path) && isfile(pat_path), 'rapid_p1_peak:NoPattern', 'no square pattern .mat found; give ''Pattern''');
    pz = load(pat_path, 'pattern'); Pats = double(pz.pattern.Pats); n_frames = size(Pats, 3);
    [~, pat_name] = fileparts(pat_path); pat_id = str2double(pat_name(1:4));

    %% log
    t0 = tic; [Log, info] = fast_tdms_folder2struct(p1_folder); t_read = toc(t0);
    t0 = tic;
    av = Log.ADC.Volts; at = double(Log.ADC.Time);
    if mean(abs(av(1, :) - round(av(1, :))) < 1e-3) > 0.9, ch = 2; layout = 'ADC0 frames / ADC1 Vm'; else, ch = 1; layout = 'ADC1 Vm / ADC2 current + frame log'; end
    vm = av(ch, :) * 10; t = at(ch, :); n = numel(vm);
    fx = double(Log.Frames.Position(1, :)); ft = double(Log.Frames.Time(1, :));
    edges = [1, find(diff(fx) ~= 0) + 1];
    % active pattern ID at each edge (restricts to the square pattern when the log holds several patterns)
    names = Log.Commands.Name; data = Log.Commands.Data; ct = double(Log.Commands.Time);
    isp = strcmp(names, 'Set Pattern ID'); pid_t = ct(isp); pid_v = cellfun(@(d) hex2dec(strtok(d, '/')), data(isp));
    n_patterns = numel(unique(pid_v));
    if ~isempty(pid_t)
        ip = discretize(ft(edges), [pid_t, Inf]); active = nan(size(edges)); ok = ~isnan(ip); active(ok) = pid_v(ip(ok));
    else
        active = nan(size(edges));
    end
    fi = fx(edges);
    k0 = discretize(ft(edges), [t, Inf]); k0(isnan(k0)) = 0;                       % count of ticks <= onset time
    hit = k0 >= 1; hit(hit) = t(k0(hit)) == ft(edges(hit)); k0(hit) = k0(hit) - 1;   % searchsorted 'left': count of ticks < onset
    use = fi ~= 0 & fi < n_frames & (n_patterns <= 1 | isnan(active) | active == pat_id);
    n_skipped_other = sum(fi ~= 0 & ~use);
    fi = fi(use); k0 = k0(use);
    dv = nan(size(fi));
    for j = 1:numel(fi)
        dv(j) = wmean(vm, k0(j), RESP_MS, FS, n) - wmean(vm, k0(j), BASE_MS, FS, n);
    end
    frames = unique(fi); rows = struct([]);
    for j = 1:numel(frames)
        f = frames(j); img = Pats(:, :, f + 1); m = img ~= BG; [rr, cc] = find(m); d = dv(fi == f);
        if any(img(m) > BG), pol = 'ON'; else, pol = 'OFF'; end
        if numel(d) > 1, sd = std(d); else, sd = NaN; end
        row = struct('frame_1based', f + 1, 'polarity', pol, 'x_center_1based', mean(cc - 1) + 1, 'y_center_1based', mean(rr - 1) + 1, ...
                     'n_repeats', numel(d), 'dV_mean_mV', mean(d), 'dV_sd_across_repeats_mV', sd);
        if isempty(rows), rows = row; else, rows(end+1) = row; end %#ok<AGROW>
    end
    nrep = max([rows.n_repeats]);
    R = struct('cell', opts.Cell, 'p1_folder', p1_folder, 'pattern', pat_path, 'layout', layout, 'n_flashes', numel(fi), 'n_repeats', nrep, 'n_skipped_other_pattern', n_skipped_other);
    for pol = {'ON', 'OFF'}
        rs = rows(strcmp({rows.polarity}, pol{1}));
        if isempty(rs), R.(pol{1}) = struct('frame_1based', NaN, 'xy', [NaN NaN], 'peak_mV', NaN, 'neg_mV', NaN, 'snr_single', NaN, 'snr_mean', NaN, 'noise_single_mV', NaN); continue; end
        a = [rs.dV_mean_mV]; sd = [rs.dV_sd_across_repeats_mV]; sig = sqrt(mean(sd(~isnan(sd)) .^ 2));
        [~, ip] = max(a); [~, ineg] = min(a);
        R.(pol{1}) = struct('frame_1based', rs(ip).frame_1based, 'xy', [rs(ip).x_center_1based, rs(ip).y_center_1based], 'peak_mV', a(ip), 'neg_mV', a(ineg), ...
                            'snr_single', a(ip) / sig, 'snr_mean', a(ip) / (sig / sqrt(nrep)), 'noise_single_mV', sig);
    end
    R.frame = R.ON.frame_1based; R.xy = R.ON.xy; R.peak_mV = R.ON.peak_mV; R.polarity = 'ON';
    R.frame12 = [R.ON.frame_1based, R.OFF.frame_1based]; R.peak12_mV = [R.ON.peak_mV, R.OFF.peak_mV];
    R.rows = rows; t_an = toc(t0);

    %% print, files, figure
    if ~isfolder(opts.Out), mkdir(opts.Out); end
    fprintf('Protocol 1 %s: %s; %d flashes, %d repeats, %d frames of %s; layout %s%s\n', opts.Cell, p1_folder, R.n_flashes, nrep, numel(frames), pat_name, layout, ...
        iff(n_skipped_other > 0, sprintf(' (%d flashes of other patterns skipped)', n_skipped_other), ''));
    for pol = {'ON', 'OFF'}
        q = R.(pol{1});
        fprintf('   %-3s peak frame %3d at col %.1f row %.1f: %.2f mV (most negative %.2f mV; noise %.2f mV single, SNR %.1f single / %.1f mean of %d)\n', ...
            pol{1}, q.frame_1based, q.xy(1), q.xy(2), q.peak_mV, q.neg_mV, q.noise_single_mV, q.snr_single, q.snr_mean, nrep);
    end
    out = R; out = rmfield(out, 'rows');
    fid = fopen(fullfile(opts.Out, [opts.Cell '_p1_summary.json']), 'w'); fprintf(fid, '%s', jsonencode(out, 'PrettyPrint', true)); fclose(fid);
    fid = fopen(fullfile(opts.Out, [opts.Cell '_p1_rf_map.csv']), 'w');
    fprintf(fid, 'frame_1based,polarity,x_center_1based,y_center_1based,n_repeats,dV_mean_mV,dV_sd_across_repeats_mV\n');
    for j = 1:numel(rows)
        r = rows(j); fprintf(fid, '%d,%s,%s,%s,%d,%s,%s\n', r.frame_1based, r.polarity, num(r.x_center_1based), num(r.y_center_1based), r.n_repeats, num(r.dV_mean_mV), num(r.dV_sd_across_repeats_mV));
    end
    fclose(fid);
    R.figure = [];
    if opts.Figure
        fig = figure('Color', 'w', 'Position', [60 60 1000 380], 'Name', ['Protocol 1 ' opts.Cell]);
        tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
        lim = max(abs([rows.dV_mean_mV]));
        for pol = {'ON', 'OFF'}
            rs = rows(strcmp({rows.polarity}, pol{1})); ax = nexttile(tl);
            if isempty(rs), title(ax, [pol{1} ': no frames']); continue; end
            xs = unique([rs.x_center_1based]); ys = unique([rs.y_center_1based]); M = nan(numel(ys), numel(xs));
            for j = 1:numel(rs), M(ys == rs(j).y_center_1based, xs == rs(j).x_center_1based) = rs(j).dV_mean_mV; end
            imagesc(ax, xs, ys, M, [-lim lim]); colormap(ax, redblue_map()); axis(ax, 'image'); set(ax, 'YDir', 'reverse'); hold(ax, 'on');
            q = R.(pol{1}); plot(ax, q.xy(1), q.xy(2), 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
            title(ax, sprintf('%s  %s peak %.1f mV (frame %d, SNR %.1f)', strrep(opts.Cell, '_', '\_'), pol{1}, q.peak_mV, q.frame_1based, q.snr_mean), 'FontSize', 9);
            xlabel(ax, 'arena column'); ylabel(ax, 'row');
        end
        cb = colorbar(ax); cb.Label.String = 'dV 50-190 ms (mV)';
        R.figure = fig;
        if opts.SavePng, exportgraphics(fig, fullfile(opts.Out, [opts.Cell '_p1_maps.png']), 'Resolution', 100); end
    end
    R.timing = struct('read_s', t_read, 'analysis_s', t_an, 'total_s', toc(t_all), 'read_info', info); R.out_dir = opts.Out;
    fprintf('written: %s   (total %.1f s: read %.2f, analysis %.2f)\n', fullfile(opts.Out, [opts.Cell '_{p1_summary.json,p1_rf_map.csv,p1_maps.png}']), R.timing.total_s, t_read, t_an);
end

function y = wmean(v, k, win, FS, n)
% mean of v over [k + win(1) ms, k + win(2) ms) with k the 0-based onset sample (protocol1_maps.wmean)
    a = k + round(win(1) * FS / 1000); b = k + round(win(2) * FS / 1000);
    if 0 <= a && a < b && b <= n, y = mean(v(a + 1:b)); else, y = NaN; end
end

function s = num(x)
    if isnan(x), s = ''; elseif x == round(x), s = sprintf('%.1f', x); else, s = sprintf('%.15g', x); end
end

function s = iff(c, a, b)
    if c, s = a; else, s = b; end
end

function cm = redblue_map()
    n = 128; r = [linspace(0.02, 1, n)'; ones(n, 1)]; g = [linspace(0.19, 1, n)'; linspace(1, 0.02, n)']; b = [ones(n, 1); linspace(1, 0.19, n)'];
    cm = [r g b];
end

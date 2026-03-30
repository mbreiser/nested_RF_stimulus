% GENERATE_FG_M2_VS_M5  Classic + Gruntman panels for both M2 and M5 alignment.
%
%   Generates 8 panels total:
%     Classic F/G (per-position raw mV amplitudes, Wilcoxon p<0.05):
%       fig_ds_panel_F_amplitudes_ON_M2  / _M5
%       fig_ds_panel_G_amplitudes_OFF_M2 / _M5
%     Gruntman F2/G2 (normalized profiles, sliding-mean Wilcoxon p<0.05):
%       fig_ds_panel_F2_gruntman_ON_M2   / _M5
%       fig_ds_panel_G2_gruntman_OFF_M2  / _M5
%
%   Asterisks use raw (uncorrected) p < 0.05 threshold throughout.
%
%   Uses 25-cell late dataset (batch_results.mat).
%
%   Outputs (in manuscript_figures/):
%     8 panel PDFs + PNGs
%     fig_ds_fg_m2_vs_m5_stats.txt
%
%   Usage:
%     run('scripts/generate_fg_m2_vs_m5.m')

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
out_dir   = fullfile(data_root, 'manuscript_figures');
if ~isfolder(out_dir), mkdir(out_dir); end

%% Load batch results
fprintf('Loading batch results...\n');
S = load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
results = S.results;
fprintf('  %d cells loaded\n', numel(results));

% Group masks
on_mask   = [results.is_on];
off_mask  = ~on_mask;
ctrl_mask = ~[results.is_ttl];
ttl_mask  = [results.is_ttl];

fprintf('  ON ctrl=%d, ON tutl-=%d, OFF ctrl=%d, OFF tutl-=%d\n', ...
    sum(on_mask & ctrl_mask), sum(on_mask & ttl_mask), ...
    sum(off_mask & ctrl_mask), sum(off_mask & ttl_mask));

%% Manuscript defaults
FONT_NAME    = 'Helvetica';
FONT_AX      = 7;
FONT_LABEL   = 8;
FONT_TITLE   = 9;
FONT_STAT    = 6;
FONT_LEGEND  = 6;

COL_CTRL      = [0 0 0];
COL_CTRL_FILL = [0.65 0.65 0.65];
COL_TTL       = [0.85 0 0];
COL_TTL_FILL  = [1.0 0.7 0.7];
ALPHA_CTRL    = 0.30;
ALPHA_TTL     = 0.20;

STAT_COLOR   = [0.15 0.15 0.6];
TTL_TEX      = '{\ittutl}^{-}';

% Response windows
STIM_ON   = 5001;
STIM_OFF  = 5801;
RESP_END  = 6551;

% Classic F/G
DEP_PCTILE    = 99.5;
HYP_PCTILE    = 0.5;
REJECT_THRESH = 0.5;  % mV

% Gruntman thresholds
DEP_SD_MULT  = 3;
HYP_SD_MULT  = 2;

% Alignment configurations
align_configs = struct( ...
    'name',  {'M2', 'M5'}, ...
    'field', {'pd_flash_peak_aligned', 'pd_flash_m5_aligned'}, ...
    'center_label', {'0 (M2 peak)', '0 (M5 centroid)'});

positions = -5:5;

%% Open stats file
stats_file = fullfile(out_dir, 'fig_ds_fg_m2_vs_m5_stats.txt');
fid = fopen(stats_file, 'w');
fprintf(fid, 'Classic F/G + Gruntman F2/G2 — M2 vs M5 Alignment Comparison\n');
fprintf(fid, 'All asterisks use raw (uncorrected) p < 0.05\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

%% Generate all panels
fprintf('\n=== Generating F/G panels for M2 and M5 alignment ===\n');

for ai = 1:numel(align_configs)
    aname = align_configs(ai).name;
    afield = align_configs(ai).field;
    acenter = align_configs(ai).center_label;

    fprintf(fid, '================================================================\n');
    fprintf(fid, '  ALIGNMENT: %s (%s)\n', aname, afield);
    fprintf(fid, '================================================================\n\n');

    for panel_idx = 1:2
        if panel_idx == 1
            type_str = 'ON'; type_long = 'T4 (ON)';
            label_classic = 'F'; label_gruntman = 'F2';
            mask = on_mask;
        else
            type_str = 'OFF'; type_long = 'T5 (OFF)';
            label_classic = 'G'; label_gruntman = 'G2';
            mask = off_mask;
        end

        ctrl_sel = mask & ctrl_mask;
        ttl_sel  = mask & ttl_mask;
        nc = sum(ctrl_sel);
        nt = sum(ttl_sel);

        fprintf('%s %s: %s — ctrl n=%d, tutl- n=%d\n', aname, label_classic, type_long, nc, nt);

        % =============================================================
        %  Classic F/G: per-position raw mV amplitudes
        % =============================================================
        [ctrl_dep, ctrl_hyp] = extract_robust_amps(results(ctrl_sel), ...
            afield, [STIM_ON RESP_END], [STIM_ON Inf], DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
        [ttl_dep, ttl_hyp]   = extract_robust_amps(results(ttl_sel), ...
            afield, [STIM_ON RESP_END], [STIM_ON Inf], DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

        % Figure
        fig = figure('Units', 'centimeters', 'Position', [2 2 9 6], ...
            'Color', 'w', 'PaperUnits', 'centimeters', ...
            'PaperSize', [9 6], 'PaperPosition', [0 0 9 6]);
        set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultTextFontName', FONT_NAME);
        tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'tight');

        % Top: depolarization
        ax1 = nexttile(tl); hold(ax1, 'on');
        dep_stats = draw_amplitude_panel(ax1, positions, ctrl_dep, ttl_dep, ...
            COL_CTRL, COL_TTL, STAT_COLOR, FONT_AX, FONT_STAT);
        xline(ax1, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
        ylabel(ax1, sprintf('Depol. (%.0f%% pctile, mV)', DEP_PCTILE), 'FontSize', FONT_LABEL);
        set(ax1, 'FontSize', FONT_AX, 'TickDir', 'out');
        xlim(ax1, [-5.5 5.5]);
        box(ax1, 'off');
        set(ax1, 'XTickLabel', []);

        % Legend
        h1 = plot(ax1, NaN, NaN, 'o', 'Color', COL_CTRL, 'MarkerFaceColor', COL_CTRL, 'MarkerSize', 3);
        h2 = plot(ax1, NaN, NaN, 'o', 'Color', COL_TTL, 'MarkerFaceColor', COL_TTL, 'MarkerSize', 3);
        legend(ax1, [h1 h2], {sprintf('ctrl (n=%d)', nc), sprintf('%s (n=%d)', TTL_TEX, nt)}, ...
            'FontSize', FONT_STAT, 'Location', 'northeast', 'Box', 'off', 'Interpreter', 'tex');

        % Bottom: hyperpolarization
        ax2 = nexttile(tl); hold(ax2, 'on');
        hyp_stats = draw_amplitude_panel(ax2, positions, ctrl_hyp, ttl_hyp, ...
            COL_CTRL, COL_TTL, STAT_COLOR, FONT_AX, FONT_STAT);
        xline(ax2, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
        ylabel(ax2, sprintf('Hyperpol. (%.0f%% pctile, mV)', HYP_PCTILE), 'FontSize', FONT_LABEL);
        xlabel(ax2, sprintf('Position relative to %s center', aname), 'FontSize', FONT_LABEL);
        set(ax2, 'FontSize', FONT_AX, 'TickDir', 'out');
        xlim(ax2, [-5.5 5.5]);
        box(ax2, 'off');

        title(tl, sprintf('%s — Per-Position Amplitudes (%s-aligned)', type_long, aname), ...
            'FontSize', FONT_TITLE, 'FontWeight', 'bold', 'FontName', FONT_NAME);

        export_panel(fig, out_dir, sprintf('fig_ds_panel_%s_amplitudes_%s_%s', ...
            label_classic, type_str, aname));
        close(fig);

        % Write classic stats
        fprintf(fid, '--- Panel %s (%s): %s Per-Position Amplitudes ---\n', ...
            label_classic, aname, type_long);
        fprintf(fid, 'ctrl n=%d, tutl- n=%d\n\n', nc, nt);
        write_amp_stats(fid, 'Depolarization', positions, ctrl_dep, ttl_dep, dep_stats);
        write_amp_stats(fid, 'Hyperpolarization', positions, ctrl_hyp, ttl_hyp, hyp_stats);

        % =============================================================
        %  Gruntman F2/G2: normalized profiles with sliding-mean stats
        % =============================================================
        fprintf('%s %s: %s Gruntman — ctrl n=%d, tutl- n=%d\n', ...
            aname, label_gruntman, type_long, nc, nt);

        % Extract with SD thresholding
        [ctrl_dep_raw, ctrl_hyp_raw, ctrl_dep_z, ctrl_hyp_z] = ...
            extract_position_profiles(results(ctrl_sel), afield, ...
            STIM_ON, RESP_END, DEP_PCTILE, HYP_PCTILE, DEP_SD_MULT, HYP_SD_MULT);
        [ttl_dep_raw, ttl_hyp_raw, ttl_dep_z, ttl_hyp_z] = ...
            extract_position_profiles(results(ttl_sel), afield, ...
            STIM_ON, RESP_END, DEP_PCTILE, HYP_PCTILE, DEP_SD_MULT, HYP_SD_MULT);

        fprintf('  SD zeroed: ctrl dep=%d hyp=%d, TTL dep=%d hyp=%d\n', ...
            ctrl_dep_z, ctrl_hyp_z, ttl_dep_z, ttl_hyp_z);

        % Normalize
        ctrl_dep_norm = normalize_profiles(ctrl_dep_raw, 'depol');
        ctrl_hyp_norm = normalize_profiles(ctrl_hyp_raw, 'hyperpol');
        ttl_dep_norm  = normalize_profiles(ttl_dep_raw, 'depol');
        ttl_hyp_norm  = normalize_profiles(ttl_hyp_raw, 'hyperpol');

        % Mean +/- SEM
        [cd_m, cd_s, cd_n] = mean_sem(ctrl_dep_norm);
        [ch_m, ch_s, ch_n] = mean_sem(ctrl_hyp_norm);
        [td_m, td_s, td_n] = mean_sem(ttl_dep_norm);
        [th_m, th_s, th_n] = mean_sem(ttl_hyp_norm);

        % Sliding-mean Wilcoxon (raw p, no FDR)
        dep_sm = sliding_mean_wilcoxon_raw(ctrl_dep_norm, ttl_dep_norm);
        hyp_sm = sliding_mean_wilcoxon_raw(ctrl_hyp_norm, ttl_hyp_norm);

        % Figure
        fig = figure('Units', 'centimeters', 'Position', [2 2 7 5.5], ...
            'Color', 'w', 'PaperUnits', 'centimeters', ...
            'PaperSize', [7 5.5], 'PaperPosition', [0 0 7 5.5]);
        set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultTextFontName', FONT_NAME);

        ax = axes(fig, 'Units', 'centimeters', 'Position', [1.0 0.9 5.5 3.8]);
        hold(ax, 'on');

        % Shaded SEM + mean lines
        draw_shaded_line(ax, positions, cd_m, cd_s, COL_CTRL, COL_CTRL_FILL, ALPHA_CTRL, '-', 1.5);
        draw_shaded_line(ax, positions, ch_m, ch_s, COL_CTRL, COL_CTRL_FILL, ALPHA_CTRL, '-', 1.5);
        draw_shaded_line(ax, positions, td_m, td_s, COL_TTL, COL_TTL_FILL, ALPHA_TTL, '-', 1.2);
        draw_shaded_line(ax, positions, th_m, th_s, COL_TTL, COL_TTL_FILL, ALPHA_TTL, '-', 1.2);

        % Markers
        plot(ax, positions, cd_m, 'o', 'Color', COL_CTRL, 'MarkerFaceColor', COL_CTRL, 'MarkerSize', 3.5, 'LineWidth', 0.5);
        plot(ax, positions, ch_m, 'o', 'Color', COL_CTRL, 'MarkerFaceColor', COL_CTRL, 'MarkerSize', 3.5, 'LineWidth', 0.5);
        plot(ax, positions, td_m, 'o', 'Color', COL_TTL, 'MarkerFaceColor', COL_TTL, 'MarkerSize', 3.5, 'LineWidth', 0.5);
        plot(ax, positions, th_m, 'o', 'Color', COL_TTL, 'MarkerFaceColor', COL_TTL, 'MarkerSize', 3.5, 'LineWidth', 0.5);

        % Reference lines
        yline(ax, 0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
        xline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);

        % Axes formatting
        set(ax, 'FontSize', FONT_AX, 'TickDir', 'out', 'Box', 'off');
        xlim(ax, [-5.5 5.5]);
        ylim(ax, [-1.15 1.15]);
        set(ax, 'XTick', -5:5, 'YTick', -1:0.5:1);
        xlabel(ax, sprintf('Stimulus position (%s-aligned)', aname), 'FontSize', FONT_LABEL);
        ylabel(ax, 'Normalized response', 'FontSize', FONT_LABEL);

        % Custom x-tick labels
        xtl = cell(1, 11);
        for i = 1:11
            if i == 1,       xtl{i} = 'ND';
            elseif i == 6,   xtl{i} = '0';
            elseif i == 11,  xtl{i} = 'PD';
            else,             xtl{i} = num2str(i - 6);
            end
        end
        set(ax, 'XTickLabel', xtl);

        % Legend
        h_c = plot(ax, NaN, NaN, '-o', 'Color', COL_CTRL, 'MarkerFaceColor', COL_CTRL, 'MarkerSize', 3, 'LineWidth', 1.2);
        h_t = plot(ax, NaN, NaN, '-o', 'Color', COL_TTL, 'MarkerFaceColor', COL_TTL, 'MarkerSize', 3, 'LineWidth', 1.0);
        leg = legend(ax, [h_c h_t], ...
            {sprintf('ctrl (n=%d)', nc), sprintf('%s (n=%d)', TTL_TEX, nt)}, ...
            'FontSize', FONT_LEGEND, 'Location', 'northeast', 'Box', 'off', 'Interpreter', 'tex');
        leg.ItemTokenSize = [10 10];

        title(ax, sprintf('%s — %s-aligned', type_long, aname), ...
            'FontSize', FONT_TITLE, 'FontWeight', 'bold');

        % Asterisks (raw p<0.05, depol above, hyperpol below)
        for i = 1:numel(dep_sm)
            if dep_sm(i).sig
                center = dep_sm(i).center;
                ci = center + 6;
                y_star = max([cd_m(ci) + cd_s(ci), td_m(ci) + td_s(ci)]) + 0.05;
                text(ax, center, y_star, sig_stars(dep_sm(i).p), ...
                    'FontSize', FONT_STAT, 'Color', STAT_COLOR, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            end
        end
        for i = 1:numel(hyp_sm)
            if hyp_sm(i).sig
                center = hyp_sm(i).center;
                ci = center + 6;
                y_star = min([ch_m(ci) - ch_s(ci), th_m(ci) - th_s(ci)]) - 0.05;
                text(ax, center, y_star, sig_stars(hyp_sm(i).p), ...
                    'FontSize', FONT_STAT, 'Color', STAT_COLOR, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            end
        end

        export_panel(fig, out_dir, sprintf('fig_ds_panel_%s_gruntman_%s_%s', ...
            label_gruntman, type_str, aname));
        close(fig);

        % Write Gruntman stats
        fprintf(fid, '--- Panel %s (%s): %s Gruntman-style ---\n', ...
            label_gruntman, aname, type_long);
        fprintf(fid, 'ctrl n=%d, tutl- n=%d\n', nc, nt);
        fprintf(fid, 'SD zeroed: ctrl dep=%d hyp=%d, TTL dep=%d hyp=%d\n\n', ...
            ctrl_dep_z, ctrl_hyp_z, ttl_dep_z, ttl_hyp_z);

        write_descriptive_stats(fid, 'Depol (normalized)', positions, cd_m, cd_s, cd_n, td_m, td_s, td_n);
        write_descriptive_stats(fid, 'Hyperpol (normalized)', positions, ch_m, ch_s, ch_n, th_m, th_s, th_n);
        write_sliding_mean_stats(fid, 'Depol sliding-mean (3-pos, raw p)', dep_sm);
        write_sliding_mean_stats(fid, 'Hyperpol sliding-mean (3-pos, raw p)', hyp_sm);
    end
end

fclose(fid);
fprintf('\nSaved stats: %s\n', stats_file);
fprintf('=== Done ===\n');


%% ========================= Local Functions ===============================

function [dep_mat, hyp_mat] = extract_robust_amps(results_sub, trace_field, ...
    dep_win, hyp_win, dep_pct, hyp_pct, reject_thresh)
% Extract per-cell robust amplitude metrics at each of 11 positions.
    n = numel(results_sub);
    dep_mat = NaN(n, 11);
    hyp_mat = NaN(n, 11);

    for k = 1:n
        if ~isfield(results_sub(k), trace_field), continue; end
        traces = results_sub(k).(trace_field);
        if isempty(traces), continue; end
        n_samp = size(traces, 2);

        for pos = 1:11
            row = traces(pos, :);
            if all(isnan(row)), continue; end

            stim_start = dep_win(1);
            if stim_start > n_samp, continue; end
            mean_resp = mean(row(stim_start:end), 'omitnan');
            if abs(mean_resp) < reject_thresh, continue; end

            % Depolarization
            de = min(dep_win(2), n_samp);
            if ~isinf(de)
                dep_mat(k, pos) = prctile(row(dep_win(1):de), dep_pct);
            else
                dep_mat(k, pos) = prctile(row(dep_win(1):end), dep_pct);
            end

            % Hyperpolarization
            he = hyp_win(2);
            if isinf(he), he = n_samp; end
            he = min(he, n_samp);
            hyp_mat(k, pos) = prctile(row(hyp_win(1):he), hyp_pct);
        end
    end
end


function stats = draw_amplitude_panel(ax, positions, ctrl_data, ttl_data, ...
    col_c, col_t, col_stat, font_ax, font_stat)
% Draw jittered dots + mean/SEM + Wilcoxon asterisks (raw p<0.05).
    n_pos = numel(positions);
    stats = struct('pos', num2cell(positions), 'n_ctrl', 0, 'n_ttl', 0, ...
        'mean_ctrl', NaN, 'mean_ttl', NaN, 'p', NaN, 'sig', false);

    for i = 1:n_pos
        x = positions(i);

        cv = ctrl_data(:, i); cv = cv(~isnan(cv));
        tv = ttl_data(:, i);  tv = tv(~isnan(tv));
        stats(i).n_ctrl = numel(cv);
        stats(i).n_ttl  = numel(tv);
        stats(i).mean_ctrl = mean(cv);
        stats(i).mean_ttl  = mean(tv);

        % Jittered dots
        if ~isempty(cv)
            jx = x - 0.12 + 0.08 * (rand(size(cv)) - 0.5);
            scatter(ax, jx, cv, 12, col_c, 'filled', 'MarkerFaceAlpha', 0.45);
        end
        if ~isempty(tv)
            jx = x + 0.12 + 0.08 * (rand(size(tv)) - 0.5);
            scatter(ax, jx, tv, 12, col_t, 'filled', 'MarkerFaceAlpha', 0.45);
        end

        % Mean +/- SEM
        if numel(cv) >= 2
            m = mean(cv); s = std(cv) / sqrt(numel(cv));
            errorbar(ax, x - 0.12, m, s, 'o', 'Color', col_c, ...
                'MarkerFaceColor', col_c, 'MarkerSize', 3.5, ...
                'LineWidth', 1.0, 'CapSize', 3);
        end
        if numel(tv) >= 2
            m = mean(tv); s = std(tv) / sqrt(numel(tv));
            errorbar(ax, x + 0.12, m, s, 'o', 'Color', col_t, ...
                'MarkerFaceColor', col_t, 'MarkerSize', 3.5, ...
                'LineWidth', 1.0, 'CapSize', 3);
        end

        % Wilcoxon rank-sum (raw p<0.05)
        if numel(cv) >= 2 && numel(tv) >= 2
            p = ranksum(cv, tv);
            stats(i).p = p;
            stats(i).sig = p < 0.05;
            if p < 0.05
                yl = ylim(ax);
                y_star = max([cv; tv]) + 0.06 * diff(yl);
                text(ax, x, y_star, sig_stars(p), 'FontSize', font_stat, ...
                    'HorizontalAlignment', 'center', 'Color', col_stat, ...
                    'FontWeight', 'bold');
            end
        end
    end
end


function [dep_mat, hyp_mat, n_dep_zeroed, n_hyp_zeroed] = extract_position_profiles( ...
    results_sub, trace_field, stim_on, resp_end, dep_pct, hyp_pct, dep_sd_mult, hyp_sd_mult)
% Gruntman-style extraction with pooled baseline SD thresholding.
    n = numel(results_sub);
    dep_mat = NaN(n, 11);
    hyp_mat = NaN(n, 11);
    n_dep_zeroed = 0;
    n_hyp_zeroed = 0;

    for k = 1:n
        if ~isfield(results_sub(k), trace_field), continue; end
        traces = results_sub(k).(trace_field);
        if isempty(traces), continue; end
        n_samp = size(traces, 2);
        if stim_on > n_samp, continue; end

        % Pooled baseline SD
        baseline_samples = [];
        for pos = 1:11
            row = traces(pos, :);
            if all(isnan(row)), continue; end
            bl = row(1:stim_on-1);
            bl = bl(~isnan(bl));
            baseline_samples = [baseline_samples, bl]; %#ok<AGROW>
        end
        if numel(baseline_samples) < 100, continue; end
        baseline_sd = std(baseline_samples);

        for pos = 1:11
            row = traces(pos, :);
            if all(isnan(row)), continue; end

            de = min(resp_end, n_samp);
            win_dep = row(stim_on:de);
            win_dep = win_dep(~isnan(win_dep));
            if ~isempty(win_dep)
                val = prctile(win_dep, dep_pct);
                if val < dep_sd_mult * baseline_sd
                    dep_mat(k, pos) = 0;
                    n_dep_zeroed = n_dep_zeroed + 1;
                else
                    dep_mat(k, pos) = val;
                end
            end

            win_hyp = row(stim_on:end);
            win_hyp = win_hyp(~isnan(win_hyp));
            if ~isempty(win_hyp)
                val = prctile(win_hyp, hyp_pct);
                if abs(val) < hyp_sd_mult * baseline_sd
                    hyp_mat(k, pos) = 0;
                    n_hyp_zeroed = n_hyp_zeroed + 1;
                else
                    hyp_mat(k, pos) = val;
                end
            end
        end
    end
end


function norm_mat = normalize_profiles(raw_mat, polarity)
% Per-cell normalization: depol -> [0,1], hyperpol -> [-1,0].
    [n, p] = size(raw_mat);
    norm_mat = NaN(n, p);
    for k = 1:n
        profile = raw_mat(k, :);
        valid = ~isnan(profile);
        if sum(valid) < 3, continue; end

        if strcmp(polarity, 'depol')
            profile(profile < 0) = 0;
            pk = max(profile);
            if pk > 0
                norm_mat(k, :) = profile / pk;
            end
        else
            profile(profile > 0) = 0;
            mn = min(profile);
            if mn < 0
                norm_mat(k, :) = profile / abs(mn);
            end
        end
    end
end


function [m, s, n] = mean_sem(data)
% Column-wise mean, SEM, count (ignoring NaN).
    n_pos = size(data, 2);
    m = NaN(1, n_pos);
    s = NaN(1, n_pos);
    n = zeros(1, n_pos);
    for i = 1:n_pos
        col = data(:, i);
        col = col(~isnan(col));
        n(i) = numel(col);
        if n(i) >= 2
            m(i) = mean(col);
            s(i) = std(col) / sqrt(n(i));
        elseif n(i) == 1
            m(i) = col;
            s(i) = 0;
        end
    end
end


function stats = sliding_mean_wilcoxon_raw(ctrl_norm, ttl_norm)
% 3-position sliding-mean Wilcoxon rank-sum with RAW p<0.05 (no FDR).
    centers = -4:4;
    n_c = numel(centers);
    stats = struct('center', num2cell(centers), ...
        'n_ctrl', 0, 'n_ttl', 0, ...
        'mean_ctrl', NaN, 'mean_ttl', NaN, ...
        'p', NaN, 'sig', false);

    for i = 1:n_c
        ci = centers(i) + 6;
        cols = (ci-1):(ci+1);

        ctrl_vals = nanmean(ctrl_norm(:, cols), 2); %#ok<NANMEAN>
        ctrl_vals = ctrl_vals(~isnan(ctrl_vals));
        ttl_vals = nanmean(ttl_norm(:, cols), 2); %#ok<NANMEAN>
        ttl_vals = ttl_vals(~isnan(ttl_vals));

        stats(i).n_ctrl = numel(ctrl_vals);
        stats(i).n_ttl  = numel(ttl_vals);
        if ~isempty(ctrl_vals), stats(i).mean_ctrl = mean(ctrl_vals); end
        if ~isempty(ttl_vals),  stats(i).mean_ttl  = mean(ttl_vals); end

        if numel(ctrl_vals) >= 2 && numel(ttl_vals) >= 2
            stats(i).p = ranksum(ctrl_vals, ttl_vals);
            stats(i).sig = stats(i).p < 0.05;  % raw threshold, no FDR
        end
    end
end


function draw_shaded_line(ax, x, y_mean, y_sem, line_col, fill_col, alpha, style, lw)
% Line with shaded SEM band.
    valid = ~isnan(y_mean) & ~isnan(y_sem);
    if sum(valid) < 2, return; end
    xv = x(valid); yv = y_mean(valid); sv = y_sem(valid);
    fill_x = [xv, fliplr(xv)];
    fill_y = [yv + sv, fliplr(yv - sv)];
    fill(ax, fill_x, fill_y, fill_col, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    plot(ax, xv, yv, style, 'Color', line_col, 'LineWidth', lw);
end


function str = sig_stars(p)
    if p < 0.001,     str = '***';
    elseif p < 0.01,  str = '**';
    else,              str = '*';
    end
end


function export_panel(fig, out_dir, name)
    ts = datestr(now, 'yyyymmdd_HHMM');
    png_ts = fullfile(out_dir, sprintf('%s_%s.png', name, ts));
    exportgraphics(fig, png_ts, 'Resolution', 300);
    fprintf('  Saved: %s\n', png_ts);
end


function write_amp_stats(fid, metric_label, positions, ctrl_data, ttl_data, stats)
    fprintf(fid, '%s:\n', metric_label);
    fprintf(fid, '%-5s  %5s  %5s  %8s  %8s  %8s  %s\n', ...
        'Pos', 'nCtrl', 'nTTL', 'MnCtrl', 'MnTTL', 'p', 'Sig');
    fprintf(fid, '%s\n', repmat('-', 1, 55));
    for i = 1:numel(positions)
        cv = ctrl_data(:, i); cv = cv(~isnan(cv));
        tv = ttl_data(:, i);  tv = tv(~isnan(tv));
        sig_str = '';
        if stats(i).sig, sig_str = sig_stars(stats(i).p); end
        fprintf(fid, '%+3d    %5d  %5d  %8.3f  %8.3f  %8.4f  %s\n', ...
            positions(i), stats(i).n_ctrl, stats(i).n_ttl, ...
            mean(cv), mean(tv), stats(i).p, sig_str);
    end
    fprintf(fid, '\n');
end


function write_descriptive_stats(fid, metric_label, positions, ...
    ctrl_mean, ctrl_sem, ctrl_n, ttl_mean, ttl_sem, ttl_n)
    fprintf(fid, '%s:\n', metric_label);
    fprintf(fid, '%-5s  %5s  %5s  %8s  %8s  %8s  %8s\n', ...
        'Pos', 'nCtrl', 'nTTL', 'MnCtrl', 'SEMctrl', 'MnTTL', 'SEMttl');
    fprintf(fid, '%s\n', repmat('-', 1, 55));
    for i = 1:numel(positions)
        fprintf(fid, '%+3d    %5d  %5d  %8.3f  %8.3f  %8.3f  %8.3f\n', ...
            positions(i), ctrl_n(i), ttl_n(i), ...
            ctrl_mean(i), ctrl_sem(i), ttl_mean(i), ttl_sem(i));
    end
    fprintf(fid, '\n');
end


function write_sliding_mean_stats(fid, metric_label, sm_stats)
    fprintf(fid, '%s:\n', metric_label);
    fprintf(fid, '%-7s  %5s  %5s  %8s  %8s  %8s  %s\n', ...
        'Center', 'nCtrl', 'nTTL', 'MnCtrl', 'MnTTL', 'p', 'Sig');
    fprintf(fid, '%s\n', repmat('-', 1, 55));
    for i = 1:numel(sm_stats)
        sig_str = '';
        if sm_stats(i).sig, sig_str = sig_stars(sm_stats(i).p); end
        fprintf(fid, '%+3d      %5d  %5d  %8.3f  %8.3f  %8.4f  %s\n', ...
            sm_stats(i).center, sm_stats(i).n_ctrl, sm_stats(i).n_ttl, ...
            sm_stats(i).mean_ctrl, sm_stats(i).mean_ttl, ...
            sm_stats(i).p, sig_str);
    end
    fprintf(fid, '\n');
end

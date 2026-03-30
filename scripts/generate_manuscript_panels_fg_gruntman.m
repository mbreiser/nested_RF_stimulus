% GENERATE_MANUSCRIPT_PANELS_FG_GRUNTMAN  Gruntman-style 1D RF profile panels.
%
%   Generates Gruntman et al. (2019) Figure 2C-style panels showing
%   normalized depolarization and hyperpolarization spatial profiles
%   by M2-aligned stimulus position.
%
%   For each cell type (T4/T5), one panel shows:
%     - Depolarization (positive): 99.5th pctile in [5001, 6551]
%     - Hyperpolarization (negative): 0.5th pctile in [5001, end]
%     - SD threshold matching Gruntman methods:
%         depol zeroed if < 3 × baseline SD; hyperpol zeroed if < 2 × SD
%         (pooled baseline SD from samples 1:5000 across all 11 positions)
%     - Per-cell normalization: depol to [0,1], hyperpol to [-1,0]
%     - Ctrl (black) vs tutl- (red), all solid lines
%     - Mean +/- SEM across cells
%
%   Statistics: 3-position sliding-mean Wilcoxon rank-sum (9 centers,
%   positions -4 to +4), BH FDR q=0.05.
%
%   Uses 25-cell late dataset (batch_results.mat).
%
%   Outputs (in manuscript_figures/):
%     fig_ds_panel_F2_gruntman_ON.pdf / .png
%     fig_ds_panel_G2_gruntman_OFF.pdf / .png
%     fig_ds_gruntman_stats.txt
%
%   Usage:
%     run('scripts/generate_manuscript_panels_fg_gruntman.m')

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
out_dir   = fullfile(data_root, 'manuscript_figures');
if ~isfolder(out_dir), mkdir(out_dir); end

%% Load batch results (25-cell late dataset)
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
FONT_LEGEND  = 6;
FONT_STAT    = 6;

% Simplified color scheme: ctrl = black, TTL = red
COL_CTRL      = [0 0 0];             % black
COL_CTRL_FILL = [0.65 0.65 0.65];    % gray fill
COL_TTL       = [0.85 0 0];          % dark red
COL_TTL_FILL  = [1.0 0.7 0.7];       % light red fill

ALPHA_CTRL   = 0.30;
ALPHA_TTL    = 0.20;

% Response windows
STIM_ON   = 5001;
RESP_END  = 6551;  % depol window end
DEP_PCTILE  = 99.5;
HYP_PCTILE  = 0.5;

% Gruntman SD thresholds
DEP_SD_MULT  = 3;   % depol threshold: 3 × baseline SD
HYP_SD_MULT  = 2;   % hyperpol threshold: 2 × baseline SD

STAT_COLOR  = [0.15 0.15 0.6];
TTL_TEX     = '{\ittutl}^{-}';

% Display scaling
HYPOL_SCALE = 2;  % magnify hyperpol display by this factor

%% Generate panels — dual M2/M5 alignment
fprintf('\n=== Generating Gruntman-style 1D RF panels (M2 + M5) ===\n');

% Stats file (timestamped)
ts = datestr(now, 'yyyymmdd_HHMM');
stats_file = fullfile(out_dir, sprintf('fig_ds_gruntman_norm_stats_%s.txt', ts));
fid = fopen(stats_file, 'w');
fprintf(fid, 'Gruntman-style normalized 1D RF profile statistics\n');
fprintf(fid, 'Depol: 99.5th pctile in [%d, %d], threshold = %d × baseline SD\n', ...
    STIM_ON, RESP_END, DEP_SD_MULT);
fprintf(fid, 'Hyperpol: 0.5th pctile in [%d, end], threshold = %d × baseline SD\n', ...
    STIM_ON, HYP_SD_MULT);
fprintf(fid, 'Baseline SD: pooled from samples [1:%d] across all 11 positions\n', STIM_ON-1);
fprintf(fid, 'Depol norm: per-cell / max(depol across positions) -> [0, 1]\n');
fprintf(fid, 'Hyperpol norm: per-cell value / abs(min(hyperpol)) -> [-1, 0]\n');
fprintf(fid, 'Hyperpol display: ×%d magnification\n', HYPOL_SCALE);
fprintf(fid, 'Stats: 3-position sliding-mean Wilcoxon rank-sum, BH FDR q=0.05\n');
fprintf(fid, 'FWHM: per-cell depol profile, Wilcoxon rank-sum\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

positions = -5:5;  % 11 positions

% Alignment configs: M2 and M5
align_configs = struct( ...
    'name',        {'M2', 'M5'}, ...
    'trace_field', {'pd_flash_peak_aligned', 'pd_flash_m5_aligned'}, ...
    'center_label', {'0 (M2 peak)', '0 (M5 centroid)'}, ...
    'xlabel_str',  {'Stimulus position (M2-aligned)', 'Stimulus position (M5-aligned)'});

for ai = 1:numel(align_configs)
    ac = align_configs(ai);

    fprintf(fid, '################ %s alignment ################\n\n', ac.name);

    for panel_idx = 1:2
        if panel_idx == 1
            label = 'F2'; type_str = 'ON'; type_long = 'T4 (ON)';
            mask = on_mask;
        else
            label = 'G2'; type_str = 'OFF'; type_long = 'T5 (OFF)';
            mask = off_mask;
        end

        ctrl_sel = mask & ctrl_mask;
        ttl_sel  = mask & ttl_mask;
        nc = sum(ctrl_sel);
        nt = sum(ttl_sel);

        fprintf('Panel %s (%s): %s — ctrl n=%d, tutl- n=%d\n', label, ac.name, type_long, nc, nt);

        % Extract raw per-cell amplitude profiles with SD thresholding
        [ctrl_dep_raw, ctrl_hyp_raw, ctrl_dep_zeroed, ctrl_hyp_zeroed] = ...
            extract_position_profiles(results(ctrl_sel), ac.trace_field, ...
            STIM_ON, RESP_END, DEP_PCTILE, HYP_PCTILE, DEP_SD_MULT, HYP_SD_MULT);
        [ttl_dep_raw, ttl_hyp_raw, ttl_dep_zeroed, ttl_hyp_zeroed] = ...
            extract_position_profiles(results(ttl_sel), ac.trace_field, ...
            STIM_ON, RESP_END, DEP_PCTILE, HYP_PCTILE, DEP_SD_MULT, HYP_SD_MULT);

        fprintf('  SD threshold zeroed: ctrl depol=%d, ctrl hyperpol=%d, TTL depol=%d, TTL hyperpol=%d\n', ...
            ctrl_dep_zeroed, ctrl_hyp_zeroed, ttl_dep_zeroed, ttl_hyp_zeroed);

        % Normalize per cell: depol -> [0, 1], hyperpol -> [-1, 0]
        ctrl_dep_norm = normalize_profiles(ctrl_dep_raw, 'depol');
        ctrl_hyp_norm = normalize_profiles(ctrl_hyp_raw, 'hyperpol');
        ttl_dep_norm  = normalize_profiles(ttl_dep_raw, 'depol');
        ttl_hyp_norm  = normalize_profiles(ttl_hyp_raw, 'hyperpol');

        % Compute mean +/- SEM (ignoring NaN positions from edge cells)
        [ctrl_dep_mean, ctrl_dep_sem, ctrl_dep_n] = mean_sem(ctrl_dep_norm);
        [ctrl_hyp_mean, ctrl_hyp_sem, ctrl_hyp_n] = mean_sem(ctrl_hyp_norm);
        [ttl_dep_mean, ttl_dep_sem, ttl_dep_n]    = mean_sem(ttl_dep_norm);
        [ttl_hyp_mean, ttl_hyp_sem, ttl_hyp_n]    = mean_sem(ttl_hyp_norm);

        % --- FWHM of depolarization per cell (for statistics) ---
        fwhm_ctrl_cells = nan(nc, 1);
        for ci_f = 1:nc
            fwhm_ctrl_cells(ci_f) = compute_profile_fwhm(positions, ctrl_dep_norm(ci_f, :));
        end
        fwhm_ttl_cells = nan(nt, 1);
        for ci_f = 1:nt
            fwhm_ttl_cells(ci_f) = compute_profile_fwhm(positions, ttl_dep_norm(ci_f, :));
        end
        fwhm_ctrl_v = fwhm_ctrl_cells(~isnan(fwhm_ctrl_cells));
        fwhm_ttl_v  = fwhm_ttl_cells(~isnan(fwhm_ttl_cells));

        % FWHM of mean curves (for drawing on plot)
        [fwhm_ctrl_mean, fwhm_ctrl_L, fwhm_ctrl_R, hm_ctrl] = ...
            compute_profile_fwhm(positions, ctrl_dep_mean);
        [fwhm_ttl_mean, fwhm_ttl_L, fwhm_ttl_R, hm_ttl] = ...
            compute_profile_fwhm(positions, ttl_dep_mean);

        % Wilcoxon rank-sum on per-cell FWHM
        if numel(fwhm_ctrl_v) >= 2 && numel(fwhm_ttl_v) >= 2
            p_fwhm = ranksum(fwhm_ctrl_v, fwhm_ttl_v);
        else
            p_fwhm = NaN;
        end

        % 3-position sliding-mean Wilcoxon rank-sum
        dep_sm_stats = sliding_mean_wilcoxon(ctrl_dep_norm, ttl_dep_norm);
        hyp_sm_stats = sliding_mean_wilcoxon(ctrl_hyp_norm, ttl_hyp_norm);

        % ---- Write stats ----
        fprintf(fid, '=== Panel %s (%s): %s Gruntman-style normalized 1D RF ===\n', ...
            label, ac.name, type_long);
        fprintf(fid, 'ctrl n=%d, tutl- n=%d\n', nc, nt);
        fprintf(fid, 'SD threshold zeroed: ctrl dep=%d, ctrl hyp=%d, TTL dep=%d, TTL hyp=%d\n\n', ...
            ctrl_dep_zeroed, ctrl_hyp_zeroed, ttl_dep_zeroed, ttl_hyp_zeroed);

        write_descriptive_stats(fid, 'Depolarization (normalized)', positions, ...
            ctrl_dep_mean, ctrl_dep_sem, ctrl_dep_n, ...
            ttl_dep_mean, ttl_dep_sem, ttl_dep_n);
        write_descriptive_stats(fid, 'Hyperpolarization (normalized)', positions, ...
            ctrl_hyp_mean, ctrl_hyp_sem, ctrl_hyp_n, ...
            ttl_hyp_mean, ttl_hyp_sem, ttl_hyp_n);

        fprintf(fid, 'Depol FWHM (bar widths, normalized profiles):\n');
        fprintf(fid, '  ctrl: %.2f +/- %.2f (n=%d)\n', ...
            mean(fwhm_ctrl_v), std(fwhm_ctrl_v)/sqrt(numel(fwhm_ctrl_v)), numel(fwhm_ctrl_v));
        fprintf(fid, '  TTL:  %.2f +/- %.2f (n=%d)\n', ...
            mean(fwhm_ttl_v), std(fwhm_ttl_v)/sqrt(numel(fwhm_ttl_v)), numel(fwhm_ttl_v));
        fprintf(fid, '  Mean curve FWHM: ctrl=%.2f, TTL=%.2f bw\n', ...
            fwhm_ctrl_mean, fwhm_ttl_mean);
        fprintf(fid, '  ctrl vs TTL: p = %.4f (Wilcoxon rank-sum)\n\n', p_fwhm);

        write_sliding_mean_stats(fid, 'Depol sliding-mean Wilcoxon (3-pos window)', dep_sm_stats);
        write_sliding_mean_stats(fid, 'Hyperpol sliding-mean Wilcoxon (3-pos window)', hyp_sm_stats);

        % ---- Figure ----
        fig = figure('Units', 'centimeters', 'Position', [2 2 8 7], ...
            'Color', 'w', 'PaperUnits', 'centimeters', ...
            'PaperSize', [8 7], 'PaperPosition', [0 0 8 7]);
        set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultTextFontName', FONT_NAME);

        ax = axes(fig, 'Position', [0.15 0.18 0.78 0.72]);
        hold(ax, 'on');

        % Scale hyperpol for display
        ctrl_hyp_mean_d = ctrl_hyp_mean * HYPOL_SCALE;
        ctrl_hyp_sem_d  = ctrl_hyp_sem  * HYPOL_SCALE;
        ttl_hyp_mean_d  = ttl_hyp_mean  * HYPOL_SCALE;
        ttl_hyp_sem_d   = ttl_hyp_sem   * HYPOL_SCALE;

        % Depol (positive): shaded SEM + line
        draw_shaded_line(ax, positions, ctrl_dep_mean, ctrl_dep_sem, ...
            COL_CTRL, COL_CTRL_FILL, ALPHA_CTRL, '-', 1.5);
        draw_shaded_line(ax, positions, ttl_dep_mean, ttl_dep_sem, ...
            COL_TTL, COL_TTL_FILL, ALPHA_TTL, '-', 1.2);

        % Hyperpol (negative, ×2 scaled): shaded SEM + line
        draw_shaded_line(ax, positions, ctrl_hyp_mean_d, ctrl_hyp_sem_d, ...
            COL_CTRL, COL_CTRL_FILL, ALPHA_CTRL, '-', 1.5);
        draw_shaded_line(ax, positions, ttl_hyp_mean_d, ttl_hyp_sem_d, ...
            COL_TTL, COL_TTL_FILL, ALPHA_TTL, '-', 1.2);

        % Markers on mean lines (all filled)
        plot(ax, positions, ctrl_dep_mean, 'o', 'Color', COL_CTRL, ...
            'MarkerFaceColor', COL_CTRL, 'MarkerSize', 3.5, 'LineWidth', 0.5);
        plot(ax, positions, ctrl_hyp_mean_d, 'o', 'Color', COL_CTRL, ...
            'MarkerFaceColor', COL_CTRL, 'MarkerSize', 3.5, 'LineWidth', 0.5);
        plot(ax, positions, ttl_dep_mean, 'o', 'Color', COL_TTL, ...
            'MarkerFaceColor', COL_TTL, 'MarkerSize', 3.5, 'LineWidth', 0.5);
        plot(ax, positions, ttl_hyp_mean_d, 'o', 'Color', COL_TTL, ...
            'MarkerFaceColor', COL_TTL, 'MarkerSize', 3.5, 'LineWidth', 0.5);

        % Zero line
        plot(ax, [-5.5 5.5], [0 0], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);

        % FWHM lines on depol curves
        if ~isnan(fwhm_ctrl_mean)
            tick_h = hm_ctrl * 0.12;
            plot(ax, [fwhm_ctrl_L fwhm_ctrl_R], [hm_ctrl hm_ctrl], '--', ...
                'Color', COL_CTRL, 'LineWidth', 0.8);
            plot(ax, [fwhm_ctrl_L fwhm_ctrl_L], [hm_ctrl-tick_h hm_ctrl+tick_h], '-', ...
                'Color', COL_CTRL, 'LineWidth', 0.7);
            plot(ax, [fwhm_ctrl_R fwhm_ctrl_R], [hm_ctrl-tick_h hm_ctrl+tick_h], '-', ...
                'Color', COL_CTRL, 'LineWidth', 0.7);
        end
        if ~isnan(fwhm_ttl_mean)
            tick_h = hm_ttl * 0.12;
            plot(ax, [fwhm_ttl_L fwhm_ttl_R], [hm_ttl hm_ttl], '--', ...
                'Color', COL_TTL, 'LineWidth', 0.8);
            plot(ax, [fwhm_ttl_L fwhm_ttl_L], [hm_ttl-tick_h hm_ttl+tick_h], '-', ...
                'Color', COL_TTL, 'LineWidth', 0.7);
            plot(ax, [fwhm_ttl_R fwhm_ttl_R], [hm_ttl-tick_h hm_ttl+tick_h], '-', ...
                'Color', COL_TTL, 'LineWidth', 0.7);
        end

        % Asterisks — depol above, hyperpol below (using display-scaled coords)
        dep_ymax   = max([ctrl_dep_mean + ctrl_dep_sem, ttl_dep_mean + ttl_dep_sem], [], 'omitnan');
        hyp_ymin_d = min([ctrl_hyp_mean_d - ctrl_hyp_sem_d, ttl_hyp_mean_d - ttl_hyp_sem_d], [], 'omitnan');

        for i = 1:numel(dep_sm_stats)
            if dep_sm_stats(i).sig
                center = dep_sm_stats(i).center;
                ci = center + 6;
                y_star = max([ctrl_dep_mean(ci) + ctrl_dep_sem(ci), ...
                              ttl_dep_mean(ci) + ttl_dep_sem(ci)]) + 0.04;
                star_str = sig_stars(dep_sm_stats(i).p);
                text(ax, center, y_star, star_str, ...
                    'FontSize', FONT_STAT+2, 'Color', STAT_COLOR, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            end
        end
        for i = 1:numel(hyp_sm_stats)
            if hyp_sm_stats(i).sig
                center = hyp_sm_stats(i).center;
                ci = center + 6;
                y_star = min([ctrl_hyp_mean_d(ci) - ctrl_hyp_sem_d(ci), ...
                              ttl_hyp_mean_d(ci) - ttl_hyp_sem_d(ci)]) - 0.04;
                star_str = sig_stars(hyp_sm_stats(i).p);
                text(ax, center, y_star, star_str, ...
                    'FontSize', FONT_STAT+2, 'Color', STAT_COLOR, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            end
        end

        % Axes formatting
        xlim(ax, [-5.5 5.5]);
        set(ax, 'XTick', -5:5, 'FontSize', FONT_AX, 'TickDir', 'out', 'Box', 'off');
        xtl = arrayfun(@(x) sprintf('%+d', x), -5:5, 'uni', false);
        xtl{6} = ac.center_label;
        set(ax, 'XTickLabel', xtl);
        xlabel(ax, 'Position (bar widths from center)', 'FontSize', FONT_LABEL);
        ylabel(ax, 'Normalized response', 'FontSize', FONT_LABEL);

        % y-limits with padding (using display-scaled hyperpol)
        all_vals_d = [ctrl_dep_mean + ctrl_dep_sem; ttl_dep_mean + ttl_dep_sem; ...
                      ctrl_hyp_mean_d - ctrl_hyp_sem_d; ttl_hyp_mean_d - ttl_hyp_sem_d];
        y_lo = min(all_vals_d(~isnan(all_vals_d)));
        y_hi = max(all_vals_d(~isnan(all_vals_d)));
        pad  = 0.20 * (y_hi - y_lo);
        ylim(ax, [y_lo - pad, y_hi + pad]);

        % Custom y-tick labels: negative ticks show true (un-scaled) values
        ytk = get(ax, 'YTick');
        ytl_new = cell(size(ytk));
        for ti = 1:numel(ytk)
            if ytk(ti) < 0
                ytl_new{ti} = sprintf('%.1f', ytk(ti) / HYPOL_SCALE);
            else
                ytl_new{ti} = sprintf('%.1f', ytk(ti));
            end
        end
        set(ax, 'YTickLabel', ytl_new);

        % Hyperpol scale annotation
        neg_ticks = ytk(ytk < 0);
        if ~isempty(neg_ticks)
            text(ax, 5.3, mean(neg_ticks), sprintf('\\times%d', HYPOL_SCALE), ...
                'FontSize', 7, 'Color', [0.4 0.4 0.4], 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
        end

        % FWHM text annotation (upper left)
        fwhm_str = sprintf('FWHM: ctrl %.1f, TTL %.1f bw', ...
            fwhm_ctrl_mean, fwhm_ttl_mean);
        if ~isnan(p_fwhm)
            if p_fwhm < 0.001
                fwhm_str = sprintf('%s  p<0.001', fwhm_str);
            else
                fwhm_str = sprintf('%s  p=%.3f', fwhm_str, p_fwhm);
            end
        end
        text(ax, -5.2, y_hi + pad*0.7, fwhm_str, ...
            'FontSize', 5.5, 'FontName', FONT_NAME, ...
            'VerticalAlignment', 'top');

        % Legend
        h1 = plot(ax, NaN, NaN, '-o', 'Color', COL_CTRL, ...
            'MarkerFaceColor', COL_CTRL, 'MarkerSize', 3, 'LineWidth', 1.2);
        h2 = plot(ax, NaN, NaN, '-o', 'Color', COL_TTL, ...
            'MarkerFaceColor', COL_TTL, 'MarkerSize', 3, 'LineWidth', 1.0);
        leg = legend(ax, [h1 h2], ...
            {sprintf('ctrl (n=%d)', nc), ...
             sprintf('%s (n=%d)', TTL_TEX, nt)}, ...
            'FontSize', FONT_LEGEND, 'Location', 'northeast', 'Box', 'off', ...
            'Interpreter', 'tex');
        leg.ItemTokenSize = [10 10];

        % Title
        title_str = sprintf('%s (%s) — %s', label, type_str, ac.name);
        title(ax, title_str, 'FontSize', FONT_TITLE, 'FontWeight', 'bold');

        % Distal / Proximal labels
        yl = ylim(ax);
        text(ax, -5, yl(2) + 0.02*(yl(2)-yl(1)), 'distal', ...
            'FontSize', 6, 'FontAngle', 'italic', 'Color', [0.3 0.3 0.3], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        text(ax, 5, yl(2) + 0.02*(yl(2)-yl(1)), 'proximal', ...
            'FontSize', 6, 'FontAngle', 'italic', 'Color', [0.3 0.3 0.3], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

        hold(ax, 'off');

        % Export
        export_panel(fig, out_dir, sprintf('fig_ds_panel_%s_gruntman_%s_%s', label, type_str, ac.name));
        close(fig);
    end
end

fclose(fid);
fprintf('\nSaved stats: %s\n', stats_file);
fprintf('=== Done ===\n');


%% ========================= Local Functions ===============================

function [dep_mat, hyp_mat, n_dep_zeroed, n_hyp_zeroed] = extract_position_profiles( ...
    results_sub, trace_field, stim_on, resp_end, dep_pct, hyp_pct, dep_sd_mult, hyp_sd_mult)
% Extract per-cell raw amplitude profiles with Gruntman-style SD thresholding.
%   dep_mat: n_cells × 11 (99.5th pctile in [stim_on, resp_end])
%   hyp_mat: n_cells × 11 (0.5th pctile in [stim_on, end])
%   n_dep_zeroed: count of positions zeroed by depol SD threshold
%   n_hyp_zeroed: count of positions zeroed by hyperpol SD threshold
%
%   SD threshold: compute pooled baseline SD from samples [1:stim_on-1]
%   across all 11 positions for each cell. Zero positions where:
%     depol < dep_sd_mult × baseline_SD
%     |hyperpol| < hyp_sd_mult × baseline_SD

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

        % Compute pooled baseline SD across all 11 positions
        baseline_samples = [];
        for pos = 1:11
            row = traces(pos, :);
            if all(isnan(row)), continue; end
            bl = row(1:stim_on-1);
            bl = bl(~isnan(bl));
            baseline_samples = [baseline_samples, bl]; %#ok<AGROW>
        end
        if numel(baseline_samples) < 100
            continue;  % not enough baseline data
        end
        baseline_sd = std(baseline_samples);

        % Extract amplitudes at each position
        for pos = 1:11
            row = traces(pos, :);
            if all(isnan(row)), continue; end

            % Depolarization: 99.5th pctile in [stim_on, resp_end]
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

            % Hyperpolarization: 0.5th pctile in [stim_on, end]
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
% Normalize per-cell profiles.
%   'depol':   threshold negatives to 0, divide by max → [0, 1]
%   'hyperpol': keep values negative, normalize by abs(min) → [-1, 0]
    [n, p] = size(raw_mat);
    norm_mat = NaN(n, p);

    for k = 1:n
        profile = raw_mat(k, :);
        valid = ~isnan(profile);
        if sum(valid) < 3, continue; end  % need >=3 valid positions

        if strcmp(polarity, 'depol')
            % Threshold negatives to 0, normalize by peak
            profile(profile < 0) = 0;
            pk = max(profile);
            if pk > 0
                norm_mat(k, :) = profile / pk;
            end
        else  % hyperpol
            % Keep values negative; normalize by abs(min) so most negative -> -1.0
            % Zero out any positive values
            profile(profile > 0) = 0;
            mn = min(profile);
            if mn < 0
                norm_mat(k, :) = profile / abs(mn);  % divides by |min|, so min -> -1
            end
        end
    end
end


function [m, s, n] = mean_sem(data)
% Column-wise mean, SEM, and count ignoring NaN.
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


function draw_shaded_line(ax, x, y_mean, y_sem, line_col, fill_col, alpha, style, lw)
% Draw a line with shaded SEM band.
    valid = ~isnan(y_mean) & ~isnan(y_sem);
    if sum(valid) < 2, return; end

    xv = x(valid);
    yv = y_mean(valid);
    sv = y_sem(valid);

    % SEM shading
    fill_x = [xv, fliplr(xv)];
    fill_y = [yv + sv, fliplr(yv - sv)];
    fill(ax, fill_x, fill_y, fill_col, 'FaceAlpha', alpha, 'EdgeColor', 'none');

    % Mean line
    plot(ax, xv, yv, style, 'Color', line_col, 'LineWidth', lw);
end


function stats = sliding_mean_wilcoxon(ctrl_norm, ttl_norm)
% 3-position sliding-mean Wilcoxon rank-sum.
%   9 centers: positions -4 to +4 (indices 2:10)
%   For each cell: nanmean of 3 adjacent normalized values.
%   Wilcoxon rank-sum ctrl vs TTL on these sliding means.
%   BH FDR q=0.05 across 9 tests.

    centers = -4:4;  % 9 centers
    n_centers = numel(centers);
    stats = struct('center', num2cell(centers), ...
        'n_ctrl', 0, 'n_ttl', 0, ...
        'mean_ctrl', NaN, 'mean_ttl', NaN, ...
        'p', NaN, 'sig', false);

    raw_p = NaN(1, n_centers);
    for i = 1:n_centers
        ci = centers(i) + 6;  % center index in 1:11
        cols = (ci-1):(ci+1);  % 3 adjacent columns

        % Per-cell sliding mean (ctrl)
        ctrl_vals = nanmean(ctrl_norm(:, cols), 2);  %#ok<NANMEAN>
        ctrl_vals = ctrl_vals(~isnan(ctrl_vals));

        % Per-cell sliding mean (TTL)
        ttl_vals = nanmean(ttl_norm(:, cols), 2);  %#ok<NANMEAN>
        ttl_vals = ttl_vals(~isnan(ttl_vals));

        stats(i).n_ctrl = numel(ctrl_vals);
        stats(i).n_ttl  = numel(ttl_vals);
        if ~isempty(ctrl_vals), stats(i).mean_ctrl = mean(ctrl_vals); end
        if ~isempty(ttl_vals),  stats(i).mean_ttl  = mean(ttl_vals); end

        if numel(ctrl_vals) >= 2 && numel(ttl_vals) >= 2
            raw_p(i) = ranksum(ctrl_vals, ttl_vals);
            stats(i).p = raw_p(i);
        end
    end

    % Benjamini-Hochberg FDR correction
    valid_idx = find(~isnan(raw_p));
    if ~isempty(valid_idx)
        [sorted_p, sort_order] = sort(raw_p(valid_idx));
        m = numel(sorted_p);
        q = 0.05;
        threshold = (1:m)' / m * q;
        sig_mask = sorted_p(:) <= threshold;
        last_sig = find(sig_mask, 1, 'last');
        if ~isempty(last_sig)
            sig_indices = sort_order(1:last_sig);
            for j = sig_indices(:)'
                stats(valid_idx(j)).sig = true;
            end
        end
    end
end


function str = sig_stars(p)
% Convert p-value to asterisk string.
    if p < 0.001,     str = '***';
    elseif p < 0.01,  str = '**';
    else,              str = '*';
    end
end


function export_panel(fig, out_dir, name)
% Export panel as timestamped PNG only.
    ts = datestr(now, 'yyyymmdd_HHMM');
    png_ts = fullfile(out_dir, sprintf('%s_%s.png', name, ts));
    exportgraphics(fig, png_ts, 'Resolution', 300);
    fprintf('  Saved: %s\n', png_ts);
end


function write_descriptive_stats(fid, metric_label, positions, ...
    ctrl_mean, ctrl_sem, ctrl_n, ttl_mean, ttl_sem, ttl_n)
% Write per-position descriptive stats (no p-values).
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
% Write sliding-mean Wilcoxon stats table.
    fprintf(fid, '%s:\n', metric_label);
    fprintf(fid, '%-7s  %5s  %5s  %8s  %8s  %8s  %s\n', ...
        'Center', 'nCtrl', 'nTTL', 'MnCtrl', 'MnTTL', 'p', 'Sig');
    fprintf(fid, '%s\n', repmat('-', 1, 55));
    for i = 1:numel(sm_stats)
        sig_str = '';
        if sm_stats(i).sig
            sig_str = sig_stars(sm_stats(i).p);
        end
        fprintf(fid, '%+3d      %5d  %5d  %8.3f  %8.3f  %8.4f  %s\n', ...
            sm_stats(i).center, sm_stats(i).n_ctrl, sm_stats(i).n_ttl, ...
            sm_stats(i).mean_ctrl, sm_stats(i).mean_ttl, ...
            sm_stats(i).p, sig_str);
    end
    fprintf(fid, '\n');
end


function [fwhm, left_x, right_x, half_max] = compute_profile_fwhm(positions, profile)
% Compute FWHM of a 1D profile by pchip interpolation and half-max walk.
%   positions: 1×N position vector (e.g. -5:5)
%   profile:   1×N amplitude vector
%   Returns: fwhm (bar widths), left_x, right_x, half_max
    fwhm = NaN;  left_x = NaN;  right_x = NaN;  half_max = NaN;
    valid = ~isnan(positions) & ~isnan(profile);
    if sum(valid) < 3, return; end
    pos_v = positions(valid);
    pro_v = profile(valid);
    [pk, ~] = max(pro_v);
    if pk <= 0, return; end
    half_max = pk / 2;

    % Interpolate to fine grid
    x_fine = linspace(pos_v(1), pos_v(end), 1000);
    y_fine = interp1(pos_v, pro_v, x_fine, 'pchip');
    [~, pk_fine] = max(y_fine);

    % Walk left from peak
    left_x = x_fine(1);
    for ki = pk_fine:-1:2
        if y_fine(ki-1) <= half_max && y_fine(ki) > half_max
            frac = (half_max - y_fine(ki-1)) / (y_fine(ki) - y_fine(ki-1));
            left_x = x_fine(ki-1) + frac * (x_fine(ki) - x_fine(ki-1));
            break;
        end
    end

    % Walk right from peak
    right_x = x_fine(end);
    for ki = pk_fine:length(y_fine)-1
        if y_fine(ki) > half_max && y_fine(ki+1) <= half_max
            frac = (half_max - y_fine(ki+1)) / (y_fine(ki) - y_fine(ki+1));
            right_x = x_fine(ki+1) + frac * (x_fine(ki) - x_fine(ki+1));
            break;
        end
    end

    fwhm = right_x - left_x;
    if fwhm <= 0 || fwhm > (pos_v(end) - pos_v(1))
        fwhm = NaN;  left_x = NaN;  right_x = NaN;
    end
end

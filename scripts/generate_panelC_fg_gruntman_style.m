% generate_panelC_fg_gruntman_style.m
%
% Two deliverables:
%   1) Standalone Panel C  – DSI using Gruntman convention:
%      DSI = (PDmax − NDmax) / PDmax
%      with PDmax / NDmax = 99.5th percentile of sweep responses
%      (current pipeline uses 98th pctile & (PD−ND)/(PD+ND))
%      Uses ALL 48 cells: 25 late (bar flash protocol) + 23 early (pre-bar-flash)
%
%   2) F/G panels in Gruntman visual style (connected lines, mean±SEM,
%      depol up / hyperpol down, black=ctrl red=TTL) for M2 and M5,
%      using raw mV per-position amplitudes (not normalized).
%      Uses only 25 late cells (which have bar flash data).
%
% Outputs  → /Users/reiserm/Documents/ttl_1DRF/manuscript_figures/
%   fig_ds_panelC_gruntman_DSI_*.pdf/.png
%   fig_ds_panel_F_gruntman_raw_ON_{M2,M5}_*.pdf/.png
%   fig_ds_panel_G_gruntman_raw_OFF_{M2,M5}_*.pdf/.png
%   fig_ds_panelC_fg_gruntman_stats.txt

%% ====== paths ==========================================================
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root      = '/Users/reiserm/Documents/ttl_1DRF';
out_dir        = fullfile(data_root, 'manuscript_figures');
late_res_file  = fullfile(data_root, 'population_results', 'batch_results.mat');
early_res_file = fullfile(data_root, 'pre-bar-flash', 'population_results', 'batch_results_pre_bf.mat');

stamp = datestr(now, 'yyyymmdd_HHMM');

% Load late dataset (25 cells — bar flash protocol)
S_late  = load(late_res_file, 'results');
results_late = S_late.results;
n_late = numel(results_late);

% Load early dataset (23 cells — pre-bar-flash, sweeps only)
S_early = load(early_res_file, 'results');
results_early = S_early.results;
n_early = numel(results_early);

n_all = n_late + n_early;
fprintf('Loaded %d late + %d early = %d total cells\n', n_late, n_early, n_all);

%% ====== Part 1: Recompute DSI with Gruntman convention =================
% Reload each cell's bar sweeps with 99.5th percentile, compute
% DSI_gruntman = (PD − ND) / PD
% Process both late and early datasets.

plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
sweep_opts_new = struct('baseline_range', [1000 9000], ...
                        'stim_trim_end',  7000, ...
                        'percentile',     99.5);   % ← Gruntman

lut_file = fullfile(fileparts(which('batch_analyze_1DRF')), 'bar_lut.mat');
L = load(lut_file);

% Preallocate for all cells
dsi_gruntman = nan(n_all, 1);
dsi_current  = nan(n_all, 1);
all_groups   = cell(n_all, 1);
all_labels   = cell(n_all, 1);   % 'late' or 'early'

%% --- Process late cells (25) ---
fprintf('\n--- Reprocessing LATE cells (n=%d) with 99.5th pctile ---\n', n_late);
for ci = 1:n_late
    folder = fullfile(data_root, results_late(ci).folder);
    all_groups{ci} = results_late(ci).group;
    all_labels{ci} = 'late';
    try
        orig_dir = pwd;
        [~, ~, Log, ~, ~] = load_protocol2_data(folder);
        cd(orig_dir);
        f_data = Log.ADC.Volts(1, :);
        v_data = Log.ADC.Volts(2, :) * 10;
        ce = load(fullfile(folder, 'currentExp.mat'), 'pattern_order', 'func_order');
        bar_data = parse_bar_data(f_data, v_data);

        max_v_new = compute_bar_sweep_responses(bar_data, plot_order, sweep_opts_new);

        [lut_directions, ~, ~, ~] = verify_lut_directions( ...
            L.Tbl, ce.pattern_order, ce.func_order, plot_order);
        lut_dirs_ordered = lut_directions(plot_order);
        [~, sort_idx] = sort(lut_dirs_ordered);
        max_v_sorted = max_v_new(sort_idx);
        max_v_polar  = [max_v_sorted; max_v_sorted(1)];

        [d_aligned, ~] = find_PD_and_order_idx(max_v_polar, 0);

        pd_val = d_aligned(5, 2);
        nd_val = d_aligned(13, 2);

        if pd_val > 0
            dsi_gruntman(ci) = (pd_val - nd_val) / pd_val;
        end
        dsi_current(ci) = results_late(ci).dsi_pdnd;

        fprintf('  [L%02d] %s %-12s DSI_old=%.3f DSI_g=%.3f\n', ...
            ci, results_late(ci).date_str, results_late(ci).group, ...
            dsi_current(ci), dsi_gruntman(ci));
    catch ME
        fprintf('  [L%02d] %s: ERROR %s\n', ci, results_late(ci).date_str, ME.message);
    end
end

%% --- Process early cells (23) ---
fprintf('\n--- Reprocessing EARLY cells (n=%d) with 99.5th pctile ---\n', n_early);
for ci = 1:n_early
    gi = n_late + ci;   % global index
    folder = results_early(ci).folder;   % early stores full path
    all_groups{gi} = results_early(ci).group;
    all_labels{gi} = 'early';
    try
        orig_dir = pwd;
        [date_str_e, ~, Log, ~, ~] = load_protocol2_data(folder);
        cd(orig_dir);
        f_data = Log.ADC.Volts(1, :);
        v_data = Log.ADC.Volts(2, :) * 10;
        ce = load(fullfile(folder, 'currentExp.mat'), 'pattern_order', 'func_order');

        % Early protocol uses parse_bar_data_pre_bf (2-speed, no bar flashes)
        bar_data = parse_bar_data_pre_bf(f_data, v_data);

        % Apply dark-bar polarity correction for pre-Oct-15 OFF cells
        is_off = ~results_early(ci).is_on;
        [bar_data, ~] = correct_off_polarity_swap( ...
            bar_data, ce.pattern_order, ce.func_order, date_str_e, is_off);

        max_v_new = compute_bar_sweep_responses(bar_data, plot_order, sweep_opts_new);

        [lut_directions, ~, ~, ~] = verify_lut_directions( ...
            L.Tbl, ce.pattern_order, ce.func_order, plot_order);
        lut_dirs_ordered = lut_directions(plot_order);
        [~, sort_idx] = sort(lut_dirs_ordered);
        max_v_sorted = max_v_new(sort_idx);
        max_v_polar  = [max_v_sorted; max_v_sorted(1)];

        [d_aligned, ~] = find_PD_and_order_idx(max_v_polar, 0);

        pd_val = d_aligned(5, 2);
        nd_val = d_aligned(13, 2);

        if pd_val > 0
            dsi_gruntman(gi) = (pd_val - nd_val) / pd_val;
        end
        dsi_current(gi) = results_early(ci).dsi_pdnd;

        fprintf('  [E%02d] %s %-12s DSI_old=%.3f DSI_g=%.3f\n', ...
            ci, date_str_e, results_early(ci).group, ...
            dsi_current(gi), dsi_gruntman(gi));
    catch ME
        fprintf('  [E%02d] %s: ERROR %s\n', ci, results_early(ci).folder, ME.message);
    end
end

%% Build combined struct for Panel C plotting (all 48 cells)
groups = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
group_labels = {'T4 ctrl', 'T4 tutl^{-}', 'T5 ctrl', 'T5 tutl^{-}'};
group_colors = {[0 0 0], [0.85 0 0], [0 0 0], [0.85 0 0]};

combined = struct('group', {}, 'dsi_gruntman', {}, 'dsi_pdnd', {}, 'batch', {});
for ci = 1:n_all
    s = struct();
    s.group        = all_groups{ci};
    s.dsi_gruntman = dsi_gruntman(ci);
    s.dsi_pdnd     = dsi_current(ci);
    s.batch        = all_labels{ci};
    combined(ci)   = s;
end

%% ---- Plot Panel C: Gruntman DSI vs current DSI (all 48 cells) ----
fig_c = figure('Units', 'centimeters', 'Position', [2 2 10 12], ...
    'Color', 'w', 'Renderer', 'painters');

% Top: Gruntman DSI
ax1 = axes(fig_c, 'Position', [0.18 0.55 0.75 0.38]);
draw_dsi_boxplot(ax1, combined, 'dsi_gruntman', ...
    'DSI = (PD-ND)/PD  [99.5th pctile]', groups, group_labels, group_colors);

% Bottom: Current DSI for comparison
ax2 = axes(fig_c, 'Position', [0.18 0.08 0.75 0.38]);
draw_dsi_boxplot(ax2, combined, 'dsi_pdnd', ...
    'DSI = (PD-ND)/(PD+ND)  [98th pctile]', groups, group_labels, group_colors);

% Labels
annotation(fig_c, 'textbox', [0.01 0.93 0.08 0.06], 'String', 'C', ...
    'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica', ...
    'EdgeColor', 'none');
annotation(fig_c, 'textbox', [0.01 0.46 0.12 0.06], 'String', 'C (old)', ...
    'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Helvetica', ...
    'EdgeColor', 'none');

fname_c = sprintf('fig_ds_panelC_gruntman_DSI_%s', stamp);
exportgraphics(fig_c, fullfile(out_dir, [fname_c '.png']), 'Resolution', 300);
fprintf('\nSaved Panel C: %s\n', fname_c);
close(fig_c);


%% ====== Part 2: F/G panels — Gruntman visual style, raw mV ============
% (Uses only the 25 late cells which have bar flash data)
% Connected lines, mean±SEM, depol up / hyperpol down, per M2 and M5

COL_CTRL = [0 0 0];
COL_TTL  = [0.85 0 0];
COL_CTRL_FILL = [0.65 0.65 0.65];
COL_TTL_FILL  = [1.0 0.7 0.7];

align_configs = struct( ...
    'name',  {'M2', 'M5'}, ...
    'field', {'pd_flash_peak_aligned', 'pd_flash_m5_aligned'}, ...
    'center_label', {'0 (M2 peak)', '0 (M5 centroid)'});

cell_types = struct( ...
    'name',     {'ON', 'OFF'}, ...
    'panel',    {'F', 'G'}, ...
    'ctrl_grp', {'on_control', 'on_control'}, ...
    'ttl_grp',  {'on_ttl', 'on_ttl'});
cell_types(1).ctrl_grp = 'on_control';
cell_types(1).ttl_grp  = 'on_ttl';
cell_types(2).ctrl_grp = 'off_control';
cell_types(2).ttl_grp  = 'off_ttl';

% Amplitude extraction params
DEPOL_WIN  = [5001, 6551];
HYPOL_WIN  = [5001, 10001];   % to end of extracted trace (10001 samples)
REJECT_THR = 0.5;             % mV rejection threshold
HYPOL_SCALE = 2;              % display magnification for hyperpol

% Open stats file
stats_file = fullfile(out_dir, sprintf('fig_ds_panelC_fg_gruntman_stats_%s.txt', stamp));
fid = fopen(stats_file, 'w');
fprintf(fid, 'Panel C (Gruntman DSI, ALL 48 cells) + F/G (raw mV, 25 late cells only)\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Late cells: %d,  Early cells: %d,  Total: %d\n\n', n_late, n_early, n_all);

% --- Write Panel C stats ---
fprintf(fid, '========== Panel C: Gruntman DSI = (PD-ND)/PD [99.5th pctile] ==========\n');
fprintf(fid, 'All 48 cells (25 late + 23 early)\n\n');
fprintf(fid, '%-14s  %3s  %5s %5s  %8s %8s   %8s %8s\n', ...
    'Group', 'n', 'nLat', 'nEar', 'Mean_g', 'SEM_g', 'Mn_old', 'SEM_old');
fprintf(fid, '%s\n', repmat('-', 1, 78));

for gi = 1:4
    mask = strcmp({combined.group}, groups{gi});
    mask_late  = mask & strcmp({combined.batch}, 'late');
    mask_early = mask & strcmp({combined.batch}, 'early');
    vals_new = [combined(mask).dsi_gruntman]; vals_new = vals_new(~isnan(vals_new));
    vals_old = [combined(mask).dsi_pdnd];     vals_old = vals_old(~isnan(vals_old));
    fprintf(fid, '%-14s  %3d  %5d %5d  %8.3f %8.3f   %8.3f %8.3f\n', ...
        groups{gi}, numel(vals_new), sum(mask_late), sum(mask_early), ...
        mean(vals_new), std(vals_new)/sqrt(numel(vals_new)), ...
        mean(vals_old), std(vals_old)/sqrt(numel(vals_old)));
end

% Wilcoxon rank-sum: ctrl vs TTL within ON and OFF
for ct = 1:2
    if ct == 1
        ctrl_mask = strcmp({combined.group}, 'on_control');
        ttl_mask  = strcmp({combined.group}, 'on_ttl');
        label = 'T4 (ON)';
    else
        ctrl_mask = strcmp({combined.group}, 'off_control');
        ttl_mask  = strcmp({combined.group}, 'off_ttl');
        label = 'T5 (OFF)';
    end
    v_c = [combined(ctrl_mask).dsi_gruntman]; v_c = v_c(~isnan(v_c));
    v_t = [combined(ttl_mask).dsi_gruntman];  v_t = v_t(~isnan(v_t));
    if numel(v_c) >= 2 && numel(v_t) >= 2
        p = ranksum(v_c, v_t);
    else
        p = NaN;
    end
    fprintf(fid, '\n%s ctrl vs TTL: p = %.4f  (n_ctrl=%d, n_ttl=%d, Wilcoxon rank-sum)\n', ...
        label, p, numel(v_c), numel(v_t));
end
fprintf(fid, '\n');

%% Generate F/G panels (late cells only — they have bar flash data)
results = results_late;   % F/G panels use late cells
for ai = 1:numel(align_configs)
    ac = align_configs(ai);
    for ct = 1:numel(cell_types)
        ctype = cell_types(ct);

        ctrl_mask = strcmp({results.group}, ctype.ctrl_grp);
        ttl_mask  = strcmp({results.group}, ctype.ttl_grp);
        ctrl_idx  = find(ctrl_mask);
        ttl_idx   = find(ttl_mask);
        n_ctrl = numel(ctrl_idx);
        n_ttl  = numel(ttl_idx);

        fprintf(fid, '========== Panel %s (%s): %s raw mV ==========\n', ...
            ctype.panel, ac.name, ctype.name);
        fprintf(fid, 'ctrl n=%d, tutl- n=%d\n\n', n_ctrl, n_ttl);

        % Extract per-position amplitudes
        positions = -5:5;
        n_pos = 11;

        [dep_ctrl, hyp_ctrl] = extract_group_amps(results, ctrl_idx, ac.field, ...
            DEPOL_WIN, HYPOL_WIN, REJECT_THR);
        [dep_ttl, hyp_ttl]   = extract_group_amps(results, ttl_idx, ac.field, ...
            DEPOL_WIN, HYPOL_WIN, REJECT_THR);

        % Compute mean ± SEM per position
        [dep_ctrl_mn, dep_ctrl_se, dep_ctrl_n] = position_stats(dep_ctrl);
        [dep_ttl_mn,  dep_ttl_se,  dep_ttl_n]  = position_stats(dep_ttl);
        [hyp_ctrl_mn, hyp_ctrl_se, hyp_ctrl_n] = position_stats(hyp_ctrl);
        [hyp_ttl_mn,  hyp_ttl_se,  hyp_ttl_n]  = position_stats(hyp_ttl);

        % Wilcoxon rank-sum at each position
        dep_p = nan(n_pos, 1);
        hyp_p = nan(n_pos, 1);
        for pi = 1:n_pos
            dc = dep_ctrl(:, pi); dc = dc(~isnan(dc));
            dt = dep_ttl(:, pi);  dt = dt(~isnan(dt));
            if numel(dc) >= 2 && numel(dt) >= 2
                dep_p(pi) = ranksum(dc, dt);
            end
            hc = hyp_ctrl(:, pi); hc = hc(~isnan(hc));
            ht = hyp_ttl(:, pi);  ht = ht(~isnan(ht));
            if numel(hc) >= 2 && numel(ht) >= 2
                hyp_p(pi) = ranksum(hc, ht);
            end
        end

        % Write stats
        fprintf(fid, 'Depolarization (mV):\n');
        fprintf(fid, 'Pos  nCtrl nTTL  MnCtrl  SEctrl   MnTTL   SEttl       p  Sig\n');
        fprintf(fid, '%s\n', repmat('-', 1, 72));
        for pi = 1:n_pos
            sig_str = '';
            if dep_p(pi) < 0.01, sig_str = '**';
            elseif dep_p(pi) < 0.05, sig_str = '*'; end
            fprintf(fid, '%+3d  %4d  %4d  %7.2f %7.2f  %7.2f %7.2f  %7.4f  %s\n', ...
                positions(pi), dep_ctrl_n(pi), dep_ttl_n(pi), ...
                dep_ctrl_mn(pi), dep_ctrl_se(pi), ...
                dep_ttl_mn(pi), dep_ttl_se(pi), dep_p(pi), sig_str);
        end

        fprintf(fid, '\nHyperpolarization (mV):\n');
        fprintf(fid, 'Pos  nCtrl nTTL  MnCtrl  SEctrl   MnTTL   SEttl       p  Sig\n');
        fprintf(fid, '%s\n', repmat('-', 1, 72));
        for pi = 1:n_pos
            sig_str = '';
            if hyp_p(pi) < 0.01, sig_str = '**';
            elseif hyp_p(pi) < 0.05, sig_str = '*'; end
            fprintf(fid, '%+3d  %4d  %4d  %7.2f %7.2f  %7.2f %7.2f  %7.4f  %s\n', ...
                positions(pi), hyp_ctrl_n(pi), hyp_ttl_n(pi), ...
                hyp_ctrl_mn(pi), hyp_ctrl_se(pi), ...
                hyp_ttl_mn(pi), hyp_ttl_se(pi), hyp_p(pi), sig_str);
        end
        fprintf(fid, '\n');

        % --- FWHM of depolarization per cell (for statistics) ---
        fwhm_ctrl_cells = nan(n_ctrl, 1);
        for ci_f = 1:n_ctrl
            fwhm_ctrl_cells(ci_f) = compute_profile_fwhm(positions, dep_ctrl(ci_f, :));
        end
        fwhm_ttl_cells = nan(n_ttl, 1);
        for ci_f = 1:n_ttl
            fwhm_ttl_cells(ci_f) = compute_profile_fwhm(positions, dep_ttl(ci_f, :));
        end
        fwhm_ctrl_v = fwhm_ctrl_cells(~isnan(fwhm_ctrl_cells));
        fwhm_ttl_v  = fwhm_ttl_cells(~isnan(fwhm_ttl_cells));

        % FWHM of mean curves (for drawing on plot)
        [fwhm_ctrl_mean, fwhm_ctrl_L, fwhm_ctrl_R, hm_ctrl] = ...
            compute_profile_fwhm(positions, dep_ctrl_mn);
        [fwhm_ttl_mean, fwhm_ttl_L, fwhm_ttl_R, hm_ttl] = ...
            compute_profile_fwhm(positions, dep_ttl_mn);

        % Wilcoxon rank-sum on per-cell FWHM
        if numel(fwhm_ctrl_v) >= 2 && numel(fwhm_ttl_v) >= 2
            p_fwhm = ranksum(fwhm_ctrl_v, fwhm_ttl_v);
        else
            p_fwhm = NaN;
        end

        fprintf(fid, 'Depol FWHM (bar widths):\n');
        fprintf(fid, '  ctrl: %.2f +/- %.2f (n=%d)\n', ...
            mean(fwhm_ctrl_v), std(fwhm_ctrl_v)/sqrt(numel(fwhm_ctrl_v)), numel(fwhm_ctrl_v));
        fprintf(fid, '  TTL:  %.2f +/- %.2f (n=%d)\n', ...
            mean(fwhm_ttl_v), std(fwhm_ttl_v)/sqrt(numel(fwhm_ttl_v)), numel(fwhm_ttl_v));
        fprintf(fid, '  Mean curve FWHM: ctrl=%.2f, TTL=%.2f bw\n', ...
            fwhm_ctrl_mean, fwhm_ttl_mean);
        fprintf(fid, '  ctrl vs TTL: p = %.4f (Wilcoxon rank-sum)\n\n', p_fwhm);

        % ---- Plot ----
        fig = figure('Units', 'centimeters', 'Position', [2 2 8 7], ...
            'Color', 'w', 'Renderer', 'painters');
        ax = axes(fig, 'Position', [0.15 0.18 0.78 0.72]);
        hold(ax, 'on');

        % Scale hyperpol for display (×2 magnification)
        hyp_ctrl_mn_d = hyp_ctrl_mn * HYPOL_SCALE;
        hyp_ctrl_se_d = hyp_ctrl_se * HYPOL_SCALE;
        hyp_ttl_mn_d  = hyp_ttl_mn  * HYPOL_SCALE;
        hyp_ttl_se_d  = hyp_ttl_se  * HYPOL_SCALE;

        % Depolarization (positive): shaded SEM + line
        draw_shaded_line(ax, positions, dep_ctrl_mn, dep_ctrl_se, ...
            COL_CTRL, COL_CTRL_FILL);
        draw_shaded_line(ax, positions, dep_ttl_mn, dep_ttl_se, ...
            COL_TTL, COL_TTL_FILL);

        % Hyperpolarization (negative, ×2 scaled): shaded SEM + line
        draw_shaded_line(ax, positions, hyp_ctrl_mn_d, hyp_ctrl_se_d, ...
            COL_CTRL, COL_CTRL_FILL);
        draw_shaded_line(ax, positions, hyp_ttl_mn_d, hyp_ttl_se_d, ...
            COL_TTL, COL_TTL_FILL);

        % Zero line
        plot(ax, [-5.5 5.5], [0 0], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);

        % --- FWHM lines on depol curves ---
        if ~isnan(fwhm_ctrl_mean)
            tick_h = hm_ctrl * 0.12;   % proportional tick height
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

        % Asterisks — depol above, hyperpol below (using display coords)
        dep_ymax  = max([dep_ctrl_mn + dep_ctrl_se, dep_ttl_mn + dep_ttl_se], [], 'omitnan');
        hyp_ymin_d = min([hyp_ctrl_mn_d - hyp_ctrl_se_d, hyp_ttl_mn_d - hyp_ttl_se_d], [], 'omitnan');
        ast_dep_y = dep_ymax + 1.5;
        ast_hyp_y = hyp_ymin_d - 1.5;

        for pi = 1:n_pos
            if dep_p(pi) < 0.05
                text(ax, positions(pi), ast_dep_y, '*', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12, ...
                    'FontWeight', 'bold', 'Color', [0.2 0.2 0.8]);
            end
            if hyp_p(pi) < 0.05
                text(ax, positions(pi), ast_hyp_y, '*', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12, ...
                    'FontWeight', 'bold', 'Color', [0.2 0.2 0.8]);
            end
        end

        % Axes formatting
        xlim(ax, [-5.5 5.5]);
        set(ax, 'XTick', -5:5);
        xtl = arrayfun(@(x) sprintf('%+d', x), -5:5, 'uni', false);
        xtl{6} = ac.center_label;
        set(ax, 'XTickLabel', xtl, 'FontSize', 7, 'FontName', 'Helvetica');
        xlabel(ax, 'Position (bar widths from center)', 'FontSize', 8);
        ylabel(ax, 'Amplitude (mV)', 'FontSize', 8);

        % y-limits with padding (using display-scaled hyperpol)
        all_vals_d = [dep_ctrl_mn + dep_ctrl_se; dep_ttl_mn + dep_ttl_se; ...
                    hyp_ctrl_mn_d - hyp_ctrl_se_d; hyp_ttl_mn_d - hyp_ttl_se_d];
        y_lo = min(all_vals_d(~isnan(all_vals_d)));
        y_hi = max(all_vals_d(~isnan(all_vals_d)));
        pad  = 0.20 * (y_hi - y_lo);
        ylim(ax, [y_lo - pad, y_hi + pad]);

        % Custom y-tick labels: negative ticks show true (un-scaled) mV values
        ytk = get(ax, 'YTick');
        ytl_new = cell(size(ytk));
        for ti = 1:numel(ytk)
            if ytk(ti) < 0
                ytl_new{ti} = sprintf('%.1f', ytk(ti) / HYPOL_SCALE);
            else
                ytl_new{ti} = sprintf('%.0f', ytk(ti));
            end
        end
        set(ax, 'YTickLabel', ytl_new);

        % Hyperpol scale annotation (near negative axis)
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
            'FontSize', 5.5, 'FontName', 'Helvetica', ...
            'VerticalAlignment', 'top');

        % Legend
        h1 = plot(ax, NaN, NaN, '-', 'Color', COL_CTRL, 'LineWidth', 1.5);
        h2 = plot(ax, NaN, NaN, '-', 'Color', COL_TTL,  'LineWidth', 1.5);
        legend(ax, [h1 h2], {sprintf('ctrl (n=%d)', n_ctrl), ...
            sprintf('tutl^- (n=%d)', n_ttl)}, ...
            'FontSize', 6, 'Location', 'northeast', 'Box', 'off');

        title_str = sprintf('%s (%s) — %s', ctype.panel, ctype.name, ac.name);
        title(ax, title_str, 'FontSize', 9, 'FontName', 'Helvetica');

        % n-count along bottom
        for pi = 1:n_pos
            nc = dep_ctrl_n(pi); nt = dep_ttl_n(pi);
            if nc > 0 || nt > 0
                text(ax, positions(pi), y_lo - pad*0.6, ...
                    sprintf('%d/%d', nc, nt), ...
                    'HorizontalAlignment', 'center', 'FontSize', 5, ...
                    'Color', [0.4 0.4 0.4]);
            end
        end

        % Distal / Proximal labels
        yl = ylim(ax);
        text(ax, -5, yl(2) + 0.02*(yl(2)-yl(1)), 'distal', ...
            'FontSize', 6, 'FontAngle', 'italic', 'Color', [0.3 0.3 0.3], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        text(ax, 5, yl(2) + 0.02*(yl(2)-yl(1)), 'proximal', ...
            'FontSize', 6, 'FontAngle', 'italic', 'Color', [0.3 0.3 0.3], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

        hold(ax, 'off');
        box(ax, 'off');

        fname_fg = sprintf('fig_ds_panel_%s_gruntman_raw_%s_%s_%s', ...
            ctype.panel, ctype.name, ac.name, stamp);
        exportgraphics(fig, fullfile(out_dir, [fname_fg '.png']), 'Resolution', 300);
        fprintf('Saved: %s\n', fname_fg);
        close(fig);
    end
end

fclose(fid);
fprintf('\nStats: %s\n', stats_file);
fprintf('=== Done ===\n');


%% ====================== LOCAL FUNCTIONS =================================

function draw_dsi_boxplot(ax, combined, field, title_str, groups, group_labels, group_colors)
% Draw a 4-group boxplot with individual dots and Wilcoxon brackets

    hold(ax, 'on');
    x_positions = [1 2 3.5 4.5];   % gap between ON and OFF

    for gi = 1:4
        mask = strcmp({combined.group}, groups{gi});
        vals = [combined(mask).(field)];
        vals = vals(~isnan(vals));
        x = x_positions(gi);

        % Jittered dots
        jitter = 0.15 * (rand(size(vals)) - 0.5);
        scatter(ax, x + jitter, vals, 18, group_colors{gi}, 'filled', ...
            'MarkerFaceAlpha', 0.6);

        % Mean ± SEM bar
        mn = mean(vals);
        se = std(vals) / sqrt(numel(vals));
        plot(ax, [x-0.2, x+0.2], [mn mn], '-', 'Color', group_colors{gi}, ...
            'LineWidth', 1.5);
        plot(ax, [x x], [mn-se, mn+se], '-', 'Color', group_colors{gi}, ...
            'LineWidth', 1);
    end

    % Wilcoxon brackets: ON ctrl vs TTL, OFF ctrl vs TTL
    bracket_pairs = {[1 2], [3 4]};
    for bi = 1:2
        idx = bracket_pairs{bi};
        g1_mask = strcmp({combined.group}, groups{idx(1)});
        g2_mask = strcmp({combined.group}, groups{idx(2)});
        v1 = [combined(g1_mask).(field)]; v1 = v1(~isnan(v1));
        v2 = [combined(g2_mask).(field)]; v2 = v2(~isnan(v2));
        if numel(v1) >= 2 && numel(v2) >= 2
            p = ranksum(v1, v2);
            y_bracket = max([v1, v2]) + 0.08;
            x1 = x_positions(idx(1)); x2 = x_positions(idx(2));
            plot(ax, [x1 x1 x2 x2], ...
                [y_bracket y_bracket+0.03 y_bracket+0.03 y_bracket], ...
                'k-', 'LineWidth', 0.8);
            text(ax, (x1+x2)/2, y_bracket+0.05, sprintf('p=%.3f', p), ...
                'HorizontalAlignment', 'center', 'FontSize', 6);
        end
    end

    set(ax, 'XTick', x_positions, 'XTickLabel', group_labels, 'FontSize', 7);
    ylabel(ax, title_str, 'FontSize', 7);
    xlim(ax, [0.3 5.2]);
    box(ax, 'off');
    hold(ax, 'off');
end


function [dep_amps, hyp_amps] = extract_group_amps(results, idx, field, ...
    depol_win, hypol_win, reject_thr)
% Extract per-cell depol/hyperpol amplitudes at each of 11 positions
% Returns n_cells × 11 matrices (NaN where rejected or missing)

    n = numel(idx);
    dep_amps = nan(n, 11);
    hyp_amps = nan(n, 11);

    for ci = 1:n
        traces = results(idx(ci)).(field);   % 11 × N_samples
        if isempty(traces), continue; end

        for pi = 1:11
            tr = traces(pi, :);
            if all(isnan(tr)), continue; end

            % Depolarization: 99.5th percentile in depol window
            d_win = tr(depol_win(1):min(depol_win(2), length(tr)));
            dep_val = prctile(d_win, 99.5);

            % Hyperpolarization: 0.5th percentile in hypol window
            h_win = tr(hypol_win(1):min(hypol_win(2), length(tr)));
            hyp_val = prctile(h_win, 0.5);

            % Rejection: mean response too small
            mean_resp = mean(tr(depol_win(1):end), 'omitnan');
            if abs(mean_resp) < reject_thr
                dep_val = NaN;
                hyp_val = NaN;
            end

            dep_amps(ci, pi) = max(dep_val, 0);    % depol is positive
            hyp_amps(ci, pi) = min(hyp_val, 0);    % hyperpol is negative
        end
    end
end


function [mn, se, n_valid] = position_stats(amps)
% Compute mean, SEM, n at each of 11 positions (ignoring NaNs)
    n_pos = size(amps, 2);
    mn = nan(1, n_pos);
    se = nan(1, n_pos);
    n_valid = zeros(1, n_pos);
    for pi = 1:n_pos
        vals = amps(:, pi);
        vals = vals(~isnan(vals));
        n_valid(pi) = numel(vals);
        if numel(vals) >= 1
            mn(pi) = mean(vals);
            se(pi) = std(vals) / sqrt(numel(vals));
        end
    end
end


function draw_shaded_line(ax, x, mn, se, line_col, fill_col)
% Draw mean line with shaded SEM band
    valid = ~isnan(mn) & ~isnan(se);
    if sum(valid) < 2, return; end
    xv = x(valid);
    mv = mn(valid);
    sv = se(valid);

    % Shaded region
    fill_x = [xv, fliplr(xv)];
    fill_y = [mv + sv, fliplr(mv - sv)];
    fill(ax, fill_x, fill_y, fill_col, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.35);

    % Mean line
    plot(ax, xv, mv, '-o', 'Color', line_col, 'LineWidth', 1.5, ...
        'MarkerSize', 4, 'MarkerFaceColor', line_col);
end


function [fwhm, left_x, right_x, half_max] = compute_profile_fwhm(positions, profile)
% Compute full-width at half-maximum of a 1-D spatial profile.
%   [fwhm, left_x, right_x, half_max] = compute_profile_fwhm(positions, profile)
%
%   positions : 1×N vector of spatial positions (e.g. -5:5)
%   profile   : 1×N vector of amplitudes (expected non-negative for depol)
%
%   Returns:
%     fwhm     — width in position units (bar widths) at half-max
%     left_x   — interpolated left half-max crossing position
%     right_x  — interpolated right half-max crossing position
%     half_max — the half-maximum value used

    fwhm = NaN;  left_x = NaN;  right_x = NaN;  half_max = NaN;

    % Remove NaN entries
    valid = ~isnan(positions) & ~isnan(profile);
    if sum(valid) < 3, return; end

    pos_v = positions(valid);
    pro_v = profile(valid);

    % Peak and half-max
    [pk, pk_idx] = max(pro_v);
    if pk <= 0, return; end
    half_max = pk / 2;

    % Interpolate to fine grid for smooth crossing detection
    x_fine = linspace(pos_v(1), pos_v(end), 1000);
    y_fine = interp1(pos_v, pro_v, x_fine, 'pchip');

    % Find fine-grid index of peak
    [~, pk_fine] = max(y_fine);

    % Walk left from peak to find half-max crossing
    left_x = x_fine(1);   % default: clamp at left edge
    for ki = pk_fine:-1:2
        if y_fine(ki-1) <= half_max && y_fine(ki) > half_max
            % Linear interpolation between ki-1 and ki
            frac = (half_max - y_fine(ki-1)) / (y_fine(ki) - y_fine(ki-1));
            left_x = x_fine(ki-1) + frac * (x_fine(ki) - x_fine(ki-1));
            break;
        end
    end

    % Walk right from peak to find half-max crossing
    right_x = x_fine(end);  % default: clamp at right edge
    for ki = pk_fine:length(y_fine)-1
        if y_fine(ki) > half_max && y_fine(ki+1) <= half_max
            frac = (half_max - y_fine(ki+1)) / (y_fine(ki) - y_fine(ki+1));
            right_x = x_fine(ki+1) + frac * (x_fine(ki) - x_fine(ki+1));
            break;
        end
    end

    fwhm = right_x - left_x;

    % Sanity: FWHM should be positive and reasonable
    if fwhm <= 0 || fwhm > (pos_v(end) - pos_v(1))
        fwhm = NaN;  left_x = NaN;  right_x = NaN;
    end
end

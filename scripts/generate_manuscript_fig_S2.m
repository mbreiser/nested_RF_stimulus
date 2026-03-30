% GENERATE_MANUSCRIPT_FIG_S2  Supplemental Figure S2: Slow vs Fast flash RF comparison.
%
%   Two-row figure comparing slow (80ms) and fast (14ms) bar flash 1D RF
%   traces and amplitude-by-position plots, both aligned with the same M6
%   center per cell.
%
%   Toggle ALIGN_MODE to generate two versions:
%     'fast_m6' — both rows aligned to fast-flash M6 center
%     'comb_m6' — both rows aligned to combined (avg of normalized slow+fast) M6 center
%
%   Layout (each half = one speed):
%     PD axis:    [T4 traces 1x11] [T5 traces 1x11]  |  [T4 amp] [T5 amp]
%     Ortho axis: [T4 traces 1x11] [T5 traces 1x11]  |  [T4 amp] [T5 amp]
%
%   Usage:
%     run('scripts/generate_manuscript_fig_S2.m')

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
out_dir   = fullfile(data_root, 'manuscript_figures');
if ~isfolder(out_dir), mkdir(out_dir); end

%% ========================= ALIGNMENT MODE ================================
% Change this to generate different versions:
ALIGN_MODE = 'fast_m6';   % 'fast_m6' or 'comb_m6'

%% Load batch results
res_file = fullfile(data_root, 'population_results', 'batch_results.mat');
fprintf('Loading: %s\n', res_file);
S = load(res_file, 'results');
results = S.results;
fprintf('Loaded %d cells.  Alignment: %s\n', numel(results), ALIGN_MODE);

%% ========================= Constants =====================================
FONT_NAME  = 'Helvetica';
FONT_AX    = 6;
FONT_LABEL = 7;
FONT_TITLE = 8;
FONT_STAT  = 5.5;

% T4 (ON) colours
COL_T4 = struct('ctrl_line', [0 0 0],           'ttl_line', [1 0 0], ...
                'ctrl_fill', [0.80 0.80 0.80],   'ttl_fill', [1 0.70 0.70], ...
                'alpha', 0.35, 'stim_line', [0.2 0.7 0.2]);

% T5 (OFF) colours
COL_T5 = struct('ctrl_line', [0.4 0.4 0.4],     'ttl_line', [0.8 0.2 0.2], ...
                'ctrl_fill', [0.70 0.70 0.70],   'ttl_fill', [0.90 0.60 0.60], ...
                'alpha', 0.35, 'stim_line', [0.2 0.7 0.2]);

STAT_COLOR = [0.15 0.15 0.6];
ALPHA_AMP  = 0.20;
LINE_W_TRACE = 0.5;
LINE_W_AMP   = 0.5;
MARKER_SZ    = 2;

% Slow timing
SLOW.stim_onset  = 5001;
SLOW.stim_offset = 5801;
SLOW.dep_window  = [5001, 5001 + 2000 - 1];
SLOW.hyp_window  = [5001 + 1000, Inf];
SLOW.trace_start = 3701;

% Fast timing
FAST.stim_onset  = 2501;
FAST.stim_offset = 2641;   % 2501 + 140
FAST.dep_window  = [2501, 2501 + 2000 - 1];
FAST.hyp_window  = [2501 + 1000, Inf];
FAST.trace_start = 1501;

DEP_PCTILE    = 99.9;
HYP_PCTILE    = 0.1;
REJECT_THRESH = 0.5;

Y_LIM    = [-5 25];
YLIM_HYP = [-5  0];
STAT_METHOD = 'mean_sem';
N_POS       = 11;
TTL_TEX = '{\ittutl}^{-}';

%% ========================= Trace field names =============================
switch ALIGN_MODE
    case 'fast_m6'
        slow_pd_field    = 'pd_flash_slow_m6fast_aligned';
        slow_ortho_field = 'ortho_flash_slow_m6fast_aligned';
        fast_pd_field    = 'pd_flash_fast_m6fast_aligned';
        fast_ortho_field = 'ortho_flash_fast_m6fast_aligned';
        align_label      = 'fast-M6';
    case 'comb_m6'
        slow_pd_field    = 'pd_flash_slow_m6comb_aligned';
        slow_ortho_field = 'ortho_flash_slow_m6comb_aligned';
        fast_pd_field    = 'pd_flash_fast_m6comb_aligned';
        fast_ortho_field = 'ortho_flash_fast_m6comb_aligned';
        align_label      = 'combined-M6';
    otherwise
        error('Unknown ALIGN_MODE: %s', ALIGN_MODE);
end

%% Group masks
on_ctrl  = [results.is_on] & ~[results.is_ttl];
on_ttl   = [results.is_on] &  [results.is_ttl];
off_ctrl = ~[results.is_on] & ~[results.is_ttl];
off_ttl  = ~[results.is_on] &  [results.is_ttl];

n_on_c  = sum(on_ctrl);   n_on_t  = sum(on_ttl);
n_off_c = sum(off_ctrl);  n_off_t = sum(off_ttl);
fprintf('T4 (ON):  ctrl=%d  tutl-=%d\n', n_on_c, n_on_t);
fprintf('T5 (OFF): ctrl=%d  tutl-=%d\n', n_off_c, n_off_t);

%% ========================= Trace statistics ==============================
% Slow PD
[spd_on_c, spd_on_t]   = compute_all_stats({results(on_ctrl).(slow_pd_field)},  {results(on_ttl).(slow_pd_field)},  N_POS, STAT_METHOD);
[spd_off_c, spd_off_t] = compute_all_stats({results(off_ctrl).(slow_pd_field)}, {results(off_ttl).(slow_pd_field)}, N_POS, STAT_METHOD);
% Slow Ortho
[so_on_c, so_on_t]   = compute_all_stats({results(on_ctrl).(slow_ortho_field)},  {results(on_ttl).(slow_ortho_field)},  N_POS, STAT_METHOD);
[so_off_c, so_off_t] = compute_all_stats({results(off_ctrl).(slow_ortho_field)}, {results(off_ttl).(slow_ortho_field)}, N_POS, STAT_METHOD);
% Fast PD
[fpd_on_c, fpd_on_t]   = compute_all_stats({results(on_ctrl).(fast_pd_field)},  {results(on_ttl).(fast_pd_field)},  N_POS, STAT_METHOD);
[fpd_off_c, fpd_off_t] = compute_all_stats({results(off_ctrl).(fast_pd_field)}, {results(off_ttl).(fast_pd_field)}, N_POS, STAT_METHOD);
% Fast Ortho
[fo_on_c, fo_on_t]   = compute_all_stats({results(on_ctrl).(fast_ortho_field)},  {results(on_ttl).(fast_ortho_field)},  N_POS, STAT_METHOD);
[fo_off_c, fo_off_t] = compute_all_stats({results(off_ctrl).(fast_ortho_field)}, {results(off_ttl).(fast_ortho_field)}, N_POS, STAT_METHOD);

%% ========================= Amplitude matrices ============================
positions = -5:5;

% Slow PD
[spd_on_c_dep, spd_on_c_hyp]   = extract_robust_amps(results(on_ctrl),  slow_pd_field, SLOW.dep_window, SLOW.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[spd_on_t_dep, spd_on_t_hyp]   = extract_robust_amps(results(on_ttl),   slow_pd_field, SLOW.dep_window, SLOW.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[spd_off_c_dep, spd_off_c_hyp] = extract_robust_amps(results(off_ctrl), slow_pd_field, SLOW.dep_window, SLOW.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[spd_off_t_dep, spd_off_t_hyp] = extract_robust_amps(results(off_ttl),  slow_pd_field, SLOW.dep_window, SLOW.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
% Slow Ortho
[so_on_c_dep, so_on_c_hyp]   = extract_robust_amps(results(on_ctrl),  slow_ortho_field, SLOW.dep_window, SLOW.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[so_on_t_dep, so_on_t_hyp]   = extract_robust_amps(results(on_ttl),   slow_ortho_field, SLOW.dep_window, SLOW.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[so_off_c_dep, so_off_c_hyp] = extract_robust_amps(results(off_ctrl), slow_ortho_field, SLOW.dep_window, SLOW.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[so_off_t_dep, so_off_t_hyp] = extract_robust_amps(results(off_ttl),  slow_ortho_field, SLOW.dep_window, SLOW.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
% Fast PD
[fpd_on_c_dep, fpd_on_c_hyp]   = extract_robust_amps(results(on_ctrl),  fast_pd_field, FAST.dep_window, FAST.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[fpd_on_t_dep, fpd_on_t_hyp]   = extract_robust_amps(results(on_ttl),   fast_pd_field, FAST.dep_window, FAST.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[fpd_off_c_dep, fpd_off_c_hyp] = extract_robust_amps(results(off_ctrl), fast_pd_field, FAST.dep_window, FAST.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[fpd_off_t_dep, fpd_off_t_hyp] = extract_robust_amps(results(off_ttl),  fast_pd_field, FAST.dep_window, FAST.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
% Fast Ortho
[fo_on_c_dep, fo_on_c_hyp]   = extract_robust_amps(results(on_ctrl),  fast_ortho_field, FAST.dep_window, FAST.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[fo_on_t_dep, fo_on_t_hyp]   = extract_robust_amps(results(on_ttl),   fast_ortho_field, FAST.dep_window, FAST.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[fo_off_c_dep, fo_off_c_hyp] = extract_robust_amps(results(off_ctrl), fast_ortho_field, FAST.dep_window, FAST.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[fo_off_t_dep, fo_off_t_hyp] = extract_robust_amps(results(off_ttl),  fast_ortho_field, FAST.dep_window, FAST.hyp_window, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

%% ========================= Pooled rank-sum ===============================
[pool_centers] = -4:4;
% Slow PD
[spd_on_dep_pp,~]  = compute_pooled_ranksum(spd_on_c_dep,  spd_on_t_dep);
[spd_on_hyp_pp,~]  = compute_pooled_ranksum(spd_on_c_hyp,  spd_on_t_hyp);
[spd_off_dep_pp,~] = compute_pooled_ranksum(spd_off_c_dep, spd_off_t_dep);
[spd_off_hyp_pp,~] = compute_pooled_ranksum(spd_off_c_hyp, spd_off_t_hyp);
% Slow Ortho
[so_on_dep_pp,~]  = compute_pooled_ranksum(so_on_c_dep,  so_on_t_dep);
[so_on_hyp_pp,~]  = compute_pooled_ranksum(so_on_c_hyp,  so_on_t_hyp);
[so_off_dep_pp,~] = compute_pooled_ranksum(so_off_c_dep, so_off_t_dep);
[so_off_hyp_pp,~] = compute_pooled_ranksum(so_off_c_hyp, so_off_t_hyp);
% Fast PD
[fpd_on_dep_pp,~]  = compute_pooled_ranksum(fpd_on_c_dep,  fpd_on_t_dep);
[fpd_on_hyp_pp,~]  = compute_pooled_ranksum(fpd_on_c_hyp,  fpd_on_t_hyp);
[fpd_off_dep_pp,~] = compute_pooled_ranksum(fpd_off_c_dep, fpd_off_t_dep);
[fpd_off_hyp_pp,~] = compute_pooled_ranksum(fpd_off_c_hyp, fpd_off_t_hyp);
% Fast Ortho
[fo_on_dep_pp,~]  = compute_pooled_ranksum(fo_on_c_dep,  fo_on_t_dep);
[fo_on_hyp_pp,~]  = compute_pooled_ranksum(fo_on_c_hyp,  fo_on_t_hyp);
[fo_off_dep_pp,~] = compute_pooled_ranksum(fo_off_c_dep, fo_off_t_dep);
[fo_off_hyp_pp,~] = compute_pooled_ranksum(fo_off_c_hyp, fo_off_t_hyp);

%% ========================= Figure ========================================
FIG_W = 18;  FIG_H = 26;

fig = figure('Units', 'centimeters', 'Position', [2 2 FIG_W FIG_H], ...
    'PaperUnits', 'centimeters', 'PaperSize', [FIG_W FIG_H], ...
    'PaperPosition', [0 0 FIG_W FIG_H], 'Color', 'w');
set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultTextFontName', FONT_NAME);

%% ========================= Layout geometry ===============================
% Left: trace tiles
TRACE_L  = 0.06;
TRACE_R  = 0.53;
TILE_GAP = 0.001;
TILE_W   = (TRACE_R - TRACE_L - (N_POS-1)*TILE_GAP) / N_POS;

% Right: amplitude panels
AMP_T4_L = 0.58;
AMP_T5_L = 0.73;
AMP_W    = 0.13;

% Vertical layout — 2 halves
HALF_H  = 0.44;    % each half occupies 44% of figure
HALF_GAP = 0.06;   % gap between halves for speed label
SLOW_BASE = 0.52;  % bottom of slow half
FAST_BASE = 0.02;  % bottom of fast half

ROW_H       = 0.08;
ROW_GAP     = 0.01;
SECTION_GAP = 0.03;

%% ========================= Draw function =================================
% This helper draws one full EFGH-like block (PD + Ortho, T4 + T5)
draw_half(fig, SLOW_BASE, HALF_H, ...
    TRACE_L, TRACE_R, TILE_W, TILE_GAP, AMP_T4_L, AMP_T5_L, AMP_W, ...
    ROW_H, ROW_GAP, SECTION_GAP, ...
    spd_on_c, spd_on_t, spd_off_c, spd_off_t, ...   % PD trace stats
    so_on_c, so_on_t, so_off_c, so_off_t, ...        % Ortho trace stats
    spd_on_c_dep, spd_on_t_dep, spd_off_c_dep, spd_off_t_dep, ...  % PD dep amps
    spd_on_c_hyp, spd_on_t_hyp, spd_off_c_hyp, spd_off_t_hyp, ...  % PD hyp amps
    so_on_c_dep, so_on_t_dep, so_off_c_dep, so_off_t_dep, ...       % Ortho dep amps
    so_on_c_hyp, so_on_t_hyp, so_off_c_hyp, so_off_t_hyp, ...      % Ortho hyp amps
    spd_on_dep_pp, spd_on_hyp_pp, spd_off_dep_pp, spd_off_hyp_pp, ...
    so_on_dep_pp, so_on_hyp_pp, so_off_dep_pp, so_off_hyp_pp, ...
    pool_centers, positions, N_POS, Y_LIM, YLIM_HYP, ...
    COL_T4, COL_T5, LINE_W_TRACE, LINE_W_AMP, MARKER_SZ, ALPHA_AMP, ...
    STAT_COLOR, FONT_AX, FONT_LABEL, FONT_TITLE, FONT_STAT, ...
    SLOW.stim_onset, SLOW.stim_offset, SLOW.trace_start, ...
    TTL_TEX, 'S', 200);  % 'S' prefix for tags, 200ms scale bar

draw_half(fig, FAST_BASE, HALF_H, ...
    TRACE_L, TRACE_R, TILE_W, TILE_GAP, AMP_T4_L, AMP_T5_L, AMP_W, ...
    ROW_H, ROW_GAP, SECTION_GAP, ...
    fpd_on_c, fpd_on_t, fpd_off_c, fpd_off_t, ...
    fo_on_c, fo_on_t, fo_off_c, fo_off_t, ...
    fpd_on_c_dep, fpd_on_t_dep, fpd_off_c_dep, fpd_off_t_dep, ...
    fpd_on_c_hyp, fpd_on_t_hyp, fpd_off_c_hyp, fpd_off_t_hyp, ...
    fo_on_c_dep, fo_on_t_dep, fo_off_c_dep, fo_off_t_dep, ...
    fo_on_c_hyp, fo_on_t_hyp, fo_off_c_hyp, fo_off_t_hyp, ...
    fpd_on_dep_pp, fpd_on_hyp_pp, fpd_off_dep_pp, fpd_off_hyp_pp, ...
    fo_on_dep_pp, fo_on_hyp_pp, fo_off_dep_pp, fo_off_hyp_pp, ...
    pool_centers, positions, N_POS, Y_LIM, YLIM_HYP, ...
    COL_T4, COL_T5, LINE_W_TRACE, LINE_W_AMP, MARKER_SZ, ALPHA_AMP, ...
    STAT_COLOR, FONT_AX, FONT_LABEL, FONT_TITLE, FONT_STAT, ...
    FAST.stim_onset, FAST.stim_offset, FAST.trace_start, ...
    TTL_TEX, 'F', 100);  % 'F' prefix for tags, 100ms scale bar

%% ========================= Panel labels ==================================
lbl_fsize = 12;
% Slow half panels
slow_pd_top  = SLOW_BASE + HALF_H - 0.01;
slow_ort_top = SLOW_BASE + HALF_H/2 - 0.01;
fast_pd_top  = FAST_BASE + HALF_H - 0.01;
fast_ort_top = FAST_BASE + HALF_H/2 - 0.01;

labels = {'A','B','C','D','E','F','G','H'};
y_vals = [slow_pd_top, slow_pd_top, slow_ort_top, slow_ort_top, ...
          fast_pd_top, fast_pd_top, fast_ort_top, fast_ort_top];
x_vals = [0.00, AMP_T4_L-0.04, 0.00, AMP_T4_L-0.04, ...
          0.00, AMP_T4_L-0.04, 0.00, AMP_T4_L-0.04];

for li = 1:8
    annotation(fig, 'textbox', [x_vals(li), y_vals(li), 0.04, 0.04], ...
        'String', labels{li}, 'FontSize', lbl_fsize, 'FontWeight', 'bold', ...
        'FontName', FONT_NAME, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

%% ========================= Speed headers =================================
annotation(fig, 'textbox', [TRACE_L, SLOW_BASE + HALF_H + 0.002, 0.50, 0.025], ...
    'String', sprintf('80 ms flash (slow)  —  aligned: %s', align_label), ...
    'FontSize', FONT_TITLE + 1, 'FontWeight', 'bold', 'FontName', FONT_NAME, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

annotation(fig, 'textbox', [TRACE_L, FAST_BASE + HALF_H + 0.002, 0.50, 0.025], ...
    'String', sprintf('14 ms flash (fast)  —  aligned: %s', align_label), ...
    'FontSize', FONT_TITLE + 1, 'FontWeight', 'bold', 'FontName', FONT_NAME, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

%% ========================= Export ========================================
ts = datestr(now, 'yyyymmdd_HHMM');
pdf_file = fullfile(out_dir, sprintf('fig_S2_flash_%s_%s.pdf', ALIGN_MODE, ts));
png_file = fullfile(out_dir, sprintf('fig_S2_flash_%s_%s.png', ALIGN_MODE, ts));

exportgraphics(fig, pdf_file, 'ContentType', 'vector');
exportgraphics(fig, png_file, 'Resolution', 300);
fprintf('Saved: %s\n', pdf_file);
fprintf('Saved: %s\n', png_file);
fprintf('=== Done (%s alignment) ===\n', ALIGN_MODE);


%% =========================================================================
%%                          LOCAL FUNCTIONS
%% =========================================================================

function draw_half(fig, base_y, half_h, ...
    trace_l, trace_r, tile_w, tile_gap, amp_t4_l, amp_t5_l, amp_w, ...
    row_h, row_gap, section_gap, ...
    pd_on_c_st, pd_on_t_st, pd_off_c_st, pd_off_t_st, ...
    ort_on_c_st, ort_on_t_st, ort_off_c_st, ort_off_t_st, ...
    pd_on_c_dep, pd_on_t_dep, pd_off_c_dep, pd_off_t_dep, ...
    pd_on_c_hyp, pd_on_t_hyp, pd_off_c_hyp, pd_off_t_hyp, ...
    ort_on_c_dep, ort_on_t_dep, ort_off_c_dep, ort_off_t_dep, ...
    ort_on_c_hyp, ort_on_t_hyp, ort_off_c_hyp, ort_off_t_hyp, ...
    pd_on_dep_pp, pd_on_hyp_pp, pd_off_dep_pp, pd_off_hyp_pp, ...
    ort_on_dep_pp, ort_on_hyp_pp, ort_off_dep_pp, ort_off_hyp_pp, ...
    pool_centers, positions, n_pos, y_lim, ylim_hyp, ...
    col_t4, col_t5, lw_trace, lw_amp, marker_sz, alpha_amp, ...
    stat_color, font_ax, font_label, font_title, font_stat, ...
    stim_onset, stim_offset, trace_start, ...
    ttl_tex, tag_prefix, scalebar_ms)

    % Vertical positions within this half
    pd_t4_y  = base_y + half_h - row_h;
    pd_t5_y  = pd_t4_y - row_h - row_gap;
    ort_t4_y = pd_t5_y - section_gap - row_h;
    ort_t5_y = ort_t4_y - row_h - row_gap;

    pos_labels = arrayfun(@(x) sprintf('%+d', x), -5:5, 'UniformOutput', false);
    pos_labels{6} = '0';

    % --- Trace rows ---
    draw_trace_row(fig, trace_l, pd_t4_y, tile_w, row_h, tile_gap, ...
        pd_on_c_st, pd_on_t_st, n_pos, y_lim, col_t4, lw_trace, ...
        pos_labels, true, font_ax, font_label, stim_onset, stim_offset, trace_start, [tag_prefix '_PD_T4']);
    add_row_label(fig, trace_l - 0.04, pd_t4_y + row_h/2, 'T4', font_title);
    add_n_label(fig, trace_r + 0.005, pd_t4_y + row_h*0.7, pd_on_c_st(6).n, pd_on_t_st(6).n, col_t4, font_ax);

    draw_trace_row(fig, trace_l, pd_t5_y, tile_w, row_h, tile_gap, ...
        pd_off_c_st, pd_off_t_st, n_pos, y_lim, col_t5, lw_trace, ...
        {}, false, font_ax, font_label, stim_onset, stim_offset, trace_start, [tag_prefix '_PD_T5']);
    add_row_label(fig, trace_l - 0.04, pd_t5_y + row_h/2, 'T5', font_title);
    add_n_label(fig, trace_r + 0.005, pd_t5_y + row_h*0.7, pd_off_c_st(6).n, pd_off_t_st(6).n, col_t5, font_ax);

    draw_trace_row(fig, trace_l, ort_t4_y, tile_w, row_h, tile_gap, ...
        ort_on_c_st, ort_on_t_st, n_pos, y_lim, col_t4, lw_trace, ...
        pos_labels, true, font_ax, font_label, stim_onset, stim_offset, trace_start, [tag_prefix '_O_T4']);
    add_row_label(fig, trace_l - 0.04, ort_t4_y + row_h/2, 'T4', font_title);
    add_n_label(fig, trace_r + 0.005, ort_t4_y + row_h*0.7, ort_on_c_st(6).n, ort_on_t_st(6).n, col_t4, font_ax);

    draw_trace_row(fig, trace_l, ort_t5_y, tile_w, row_h, tile_gap, ...
        ort_off_c_st, ort_off_t_st, n_pos, y_lim, col_t5, lw_trace, ...
        {}, false, font_ax, font_label, stim_onset, stim_offset, trace_start, [tag_prefix '_O_T5']);
    add_row_label(fig, trace_l - 0.04, ort_t5_y + row_h/2, 'T5', font_title);
    add_n_label(fig, trace_r + 0.005, ort_t5_y + row_h*0.7, ort_off_c_st(6).n, ort_off_t_st(6).n, col_t5, font_ax);

    % --- Section headers ---
    annotation(fig, 'textbox', [trace_l, pd_t4_y + row_h + 0.003, 0.25, 0.02], ...
        'String', 'PD axis', 'FontSize', font_title, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    annotation(fig, 'textbox', [trace_l, ort_t4_y + row_h + 0.003, 0.25, 0.02], ...
        'String', 'Orthogonal axis', 'FontSize', font_title, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

    % --- Amplitude panels ---
    % PD section vertical positions
    pd_top = pd_t4_y + row_h;  pd_bot = pd_t5_y;
    pd_h = pd_top - pd_bot;
    pd_dep_h = pd_h * 0.64;  pd_hyp_h = pd_h * 0.26;
    pd_dep_y = pd_top - pd_dep_h;
    pd_hyp_y = pd_bot;

    ort_top = ort_t4_y + row_h;  ort_bot = ort_t5_y;
    ort_h = ort_top - ort_bot;
    ort_dep_h = ort_h * 0.64;  ort_hyp_h = ort_h * 0.26;
    ort_dep_y = ort_top - ort_dep_h;
    ort_hyp_y = ort_bot;

    % PD Dep T4
    ax = axes(fig, 'Position', [amp_t4_l, pd_dep_y, amp_w, pd_dep_h]);
    st = draw_amp_line(ax, positions, pd_on_c_dep, pd_on_t_dep, col_t4.ctrl_line, col_t4.ttl_line, lw_amp, marker_sz, font_ax, col_t4.ctrl_fill, col_t4.ttl_fill, alpha_amp);
    draw_pooled_asterisks(ax, pool_centers, pd_on_dep_pp, st, stat_color, font_stat);
    format_amp(ax, y_lim, font_ax, false, false); ylabel(ax, 'mV', 'FontSize', font_label);
    add_fwhm_bars(ax, positions, st, col_t4.ctrl_line, col_t4.ttl_line, y_lim, font_stat, ttl_tex, pd_on_c_dep, pd_on_t_dep);

    % PD Dep T5
    ax = axes(fig, 'Position', [amp_t5_l, pd_dep_y, amp_w, pd_dep_h]);
    st = draw_amp_line(ax, positions, pd_off_c_dep, pd_off_t_dep, col_t5.ctrl_line, col_t5.ttl_line, lw_amp, marker_sz, font_ax, col_t5.ctrl_fill, col_t5.ttl_fill, alpha_amp);
    draw_pooled_asterisks(ax, pool_centers, pd_off_dep_pp, st, stat_color, font_stat);
    format_amp(ax, y_lim, font_ax, false, true);
    add_fwhm_bars(ax, positions, st, col_t5.ctrl_line, col_t5.ttl_line, y_lim, font_stat, ttl_tex, pd_off_c_dep, pd_off_t_dep);

    % PD Hyp T4
    ax = axes(fig, 'Position', [amp_t4_l, pd_hyp_y, amp_w, pd_hyp_h]);
    st = draw_amp_line(ax, positions, pd_on_c_hyp, pd_on_t_hyp, col_t4.ctrl_line, col_t4.ttl_line, lw_amp, marker_sz, font_ax, col_t4.ctrl_fill, col_t4.ttl_fill, alpha_amp);
    draw_pooled_asterisks(ax, pool_centers, pd_on_hyp_pp, st, stat_color, font_stat);
    format_amp(ax, ylim_hyp, font_ax, false, false); ylabel(ax, 'mV', 'FontSize', font_label);

    % PD Hyp T5
    ax = axes(fig, 'Position', [amp_t5_l, pd_hyp_y, amp_w, pd_hyp_h]);
    st = draw_amp_line(ax, positions, pd_off_c_hyp, pd_off_t_hyp, col_t5.ctrl_line, col_t5.ttl_line, lw_amp, marker_sz, font_ax, col_t5.ctrl_fill, col_t5.ttl_fill, alpha_amp);
    draw_pooled_asterisks(ax, pool_centers, pd_off_hyp_pp, st, stat_color, font_stat);
    format_amp(ax, ylim_hyp, font_ax, false, true);

    % Ortho Dep T4
    ax = axes(fig, 'Position', [amp_t4_l, ort_dep_y, amp_w, ort_dep_h]);
    st = draw_amp_line(ax, positions, ort_on_c_dep, ort_on_t_dep, col_t4.ctrl_line, col_t4.ttl_line, lw_amp, marker_sz, font_ax, col_t4.ctrl_fill, col_t4.ttl_fill, alpha_amp);
    draw_pooled_asterisks(ax, pool_centers, ort_on_dep_pp, st, stat_color, font_stat);
    format_amp(ax, y_lim, font_ax, false, false); ylabel(ax, 'mV', 'FontSize', font_label);
    add_fwhm_bars(ax, positions, st, col_t4.ctrl_line, col_t4.ttl_line, y_lim, font_stat, ttl_tex, ort_on_c_dep, ort_on_t_dep);

    % Ortho Dep T5
    ax = axes(fig, 'Position', [amp_t5_l, ort_dep_y, amp_w, ort_dep_h]);
    st = draw_amp_line(ax, positions, ort_off_c_dep, ort_off_t_dep, col_t5.ctrl_line, col_t5.ttl_line, lw_amp, marker_sz, font_ax, col_t5.ctrl_fill, col_t5.ttl_fill, alpha_amp);
    draw_pooled_asterisks(ax, pool_centers, ort_off_dep_pp, st, stat_color, font_stat);
    format_amp(ax, y_lim, font_ax, false, true);
    add_fwhm_bars(ax, positions, st, col_t5.ctrl_line, col_t5.ttl_line, y_lim, font_stat, ttl_tex, ort_off_c_dep, ort_off_t_dep);

    % Ortho Hyp T4
    ax = axes(fig, 'Position', [amp_t4_l, ort_hyp_y, amp_w, ort_hyp_h]);
    st = draw_amp_line(ax, positions, ort_on_c_hyp, ort_on_t_hyp, col_t4.ctrl_line, col_t4.ttl_line, lw_amp, marker_sz, font_ax, col_t4.ctrl_fill, col_t4.ttl_fill, alpha_amp);
    draw_pooled_asterisks(ax, pool_centers, ort_on_hyp_pp, st, stat_color, font_stat);
    format_amp(ax, ylim_hyp, font_ax, true, false); ylabel(ax, 'mV', 'FontSize', font_label);

    % Ortho Hyp T5
    ax = axes(fig, 'Position', [amp_t5_l, ort_hyp_y, amp_w, ort_hyp_h]);
    st = draw_amp_line(ax, positions, ort_off_c_hyp, ort_off_t_hyp, col_t5.ctrl_line, col_t5.ttl_line, lw_amp, marker_sz, font_ax, col_t5.ctrl_fill, col_t5.ttl_fill, alpha_amp);
    draw_pooled_asterisks(ax, pool_centers, ort_off_hyp_pp, st, stat_color, font_stat);
    format_amp(ax, ylim_hyp, font_ax, true, true);

    % Column headers over amplitude panels (only for top section of this half)
    annotation(fig, 'textbox', [amp_t4_l, pd_dep_y + pd_dep_h + 0.003, amp_w, 0.02], ...
        'String', 'T4 (ON)', 'FontSize', font_title, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    annotation(fig, 'textbox', [amp_t5_l, pd_dep_y + pd_dep_h + 0.003, amp_w, 0.02], ...
        'String', 'T5 (OFF)', 'FontSize', font_title, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    % Scale bar on bottom-left tile (Ortho T5, position 1)
    ax_sb = findobj(fig, 'Type', 'axes', 'Tag', sprintf('%s_O_T5_1', tag_prefix));
    if ~isempty(ax_sb)
        add_horiz_scale_bar(ax_sb(1), y_lim, font_ax, trace_start, scalebar_ms);
    end
end


% ==================== Trace panel helpers =================================

function [ctrl_stats, ttl_stats] = compute_all_stats(traces_ctrl, traces_ttl, n_pos, stat_method)
    ctrl_stats = struct('center', cell(1, n_pos), 'spread', cell(1, n_pos), 'n', num2cell(zeros(1, n_pos)));
    ttl_stats  = struct('center', cell(1, n_pos), 'spread', cell(1, n_pos), 'n', num2cell(zeros(1, n_pos)));
    for pi = 1:n_pos
        [ctrl_stats(pi).center, ctrl_stats(pi).spread, ctrl_stats(pi).n] = stack_and_stat(traces_ctrl, pi, stat_method);
        [ttl_stats(pi).center,  ttl_stats(pi).spread,  ttl_stats(pi).n]  = stack_and_stat(traces_ttl,  pi, stat_method);
    end
end

function [center, spread, n] = stack_and_stat(traces_cell, pos_idx, stat_method)
    center = []; spread = []; n = 0;
    if isempty(traces_cell), return; end
    min_len = Inf;
    for k = 1:numel(traces_cell)
        mat = traces_cell{k};
        if ~isempty(mat) && pos_idx <= size(mat, 1)
            min_len = min(min_len, size(mat, 2));
        end
    end
    if isinf(min_len), return; end
    all_tr = [];
    for k = 1:numel(traces_cell)
        mat = traces_cell{k};
        if ~isempty(mat) && pos_idx <= size(mat, 1)
            row = mat(pos_idx, 1:min_len);
            if ~all(isnan(row)), all_tr = [all_tr; row]; end %#ok<AGROW>
        end
    end
    if isempty(all_tr), return; end
    n = size(all_tr, 1);
    if strcmpi(stat_method, 'mean_sem')
        center = mean(all_tr, 1, 'omitnan');
        spread = std(all_tr, 0, 1, 'omitnan') ./ sqrt(n);
    else
        center = median(all_tr, 1, 'omitnan');
        spread = mad(all_tr, 1, 1);
    end
end


function draw_trace_row(fig, x_left, y_bottom, tile_w, tile_h, tile_gap, ...
    ctrl_stats, ttl_stats, n_pos, y_lim, col, line_w, ...
    pos_labels, show_titles, font_size, label_font, stim_onset, stim_offset, trace_start, row_tag)
    for pi = 1:n_pos
        x = x_left + (pi-1) * (tile_w + tile_gap);
        ax = axes(fig, 'Position', [x, y_bottom, tile_w, tile_h]); %#ok<LAXES>
        if ~isempty(row_tag), ax.Tag = sprintf('%s_%d', row_tag, pi); end
        hold(ax, 'on');
        ctr_c = ctrl_stats(pi).center;  spr_c = ctrl_stats(pi).spread;
        ctr_t = ttl_stats(pi).center;   spr_t = ttl_stats(pi).spread;
        x_len = 0;
        if ~isempty(ctr_c), x_len = max(x_len, numel(ctr_c)); end
        if ~isempty(ctr_t), x_len = max(x_len, numel(ctr_t)); end
        if x_len == 0
            ylim(ax, y_lim); set(ax, 'XTick', [], 'YTick', []); box(ax, 'off'); continue;
        end
        fill(ax, [stim_onset, stim_offset, stim_offset, stim_onset], ...
            [-2, -2, -1, -1], col.stim_line, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        xv = 1:x_len;
        if ~isempty(ctr_c) && ctrl_stats(pi).n >= 2
            fill(ax, [xv, fliplr(xv)], [ctr_c + spr_c, fliplr(ctr_c - spr_c)], ...
                col.ctrl_fill, 'FaceAlpha', col.alpha, 'EdgeColor', 'none');
            plot(ax, xv, ctr_c, '-', 'Color', col.ctrl_line, 'LineWidth', line_w);
        end
        if ~isempty(ctr_t) && ttl_stats(pi).n >= 2
            fill(ax, [xv, fliplr(xv)], [ctr_t + spr_t, fliplr(ctr_t - spr_t)], ...
                col.ttl_fill, 'FaceAlpha', col.alpha, 'EdgeColor', 'none');
            plot(ax, xv, ctr_t, '-', 'Color', col.ttl_line, 'LineWidth', line_w);
        end
        ylim(ax, y_lim); xlim(ax, [trace_start, x_len]);
        set(ax, 'XTick', []);
        if pi == 1
            ylabel(ax, '\DeltamV', 'FontSize', label_font);
            set(ax, 'FontSize', font_size, 'TickDir', 'out', 'TickLength', [0.02 0.02]);
        else
            set(ax, 'YTick', []);
        end
        if show_titles && ~isempty(pos_labels)
            title(ax, pos_labels{pi}, 'FontSize', font_size, 'FontWeight', 'normal');
        end
        box(ax, 'off'); ax.LineWidth = 0.5;
    end
end


function add_row_label(fig, x, y, label_str, font_size)
    annotation(fig, 'textbox', [x, y - 0.04, 0.04, 0.08], ...
        'String', label_str, 'FontSize', font_size + 1, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Rotation', 90);
end


function add_n_label(fig, x, y, n_ctrl, n_ttl, col, font_size)
    str = sprintf('n=%d/%d', n_ctrl, n_ttl);
    annotation(fig, 'textbox', [x, y - 0.02, 0.05, 0.04], ...
        'String', str, 'FontSize', font_size, ...
        'FontName', 'Helvetica', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'Color', col.ctrl_line * 0.6 + 0.4);
end


function add_horiz_scale_bar(ax, y_lim, font_size, trace_start, ms)
    v_range = diff(y_lim);
    bar_t = ms * 10;  % ms to samples at 10kHz
    x_left   = trace_start + 0.05 * (xlim(ax) * [0;1] - trace_start);
    y_bottom = y_lim(1) + 0.08 * v_range;
    plot(ax, [x_left, x_left + bar_t], [y_bottom, y_bottom], 'k-', ...
        'LineWidth', 0.8, 'Clipping', 'off');
    text(ax, x_left, y_bottom - 0.06*v_range, sprintf('%d ms', ms), ...
        'HorizontalAlignment', 'left', 'FontSize', font_size, 'Clipping', 'off');
end


% ==================== Amplitude panel helpers =============================

function [dep_mat, hyp_mat] = extract_robust_amps(results_sub, trace_field, ...
    dep_win, hyp_win, dep_pct, hyp_pct, reject_thresh)
    n = numel(results_sub);
    dep_mat = NaN(n, 11); hyp_mat = NaN(n, 11);
    for k = 1:n
        if ~isfield(results_sub(k), trace_field), continue; end
        traces = results_sub(k).(trace_field);
        if isempty(traces), continue; end
        n_samp = size(traces, 2);
        for pos = 1:11
            row = traces(pos, :);
            if all(isnan(row)), continue; end
            de = min(dep_win(2), n_samp);
            if ~isinf(de)
                dep_mat(k, pos) = prctile(row(dep_win(1):de), dep_pct);
            else
                dep_mat(k, pos) = prctile(row(dep_win(1):end), dep_pct);
            end
            he = hyp_win(2); if isinf(he), he = n_samp; end
            he = min(he, n_samp);
            hval = prctile(row(hyp_win(1):he), hyp_pct);
            hyp_mat(k, pos) = min(hval, 0);
        end
    end
end


function stats = draw_amp_line(ax, positions, ctrl_data, ttl_data, ...
    col_c, col_t, line_w, marker_sz, font_ax, col_cf, col_tf, alpha_fill)
    hold(ax, 'on');
    n_pos = numel(positions);
    stats = struct('pos', num2cell(positions), 'n_ctrl', 0, 'n_ttl', 0, ...
        'mean_ctrl', NaN, 'mean_ttl', NaN, 'sem_ctrl', NaN, 'sem_ttl', NaN);
    ctrl_mean = NaN(1, n_pos); ctrl_sem = NaN(1, n_pos);
    ttl_mean  = NaN(1, n_pos); ttl_sem  = NaN(1, n_pos);
    for i = 1:n_pos
        cv = ctrl_data(:, i); cv = cv(~isnan(cv));
        tv = ttl_data(:, i);  tv = tv(~isnan(tv));
        stats(i).n_ctrl = numel(cv); stats(i).n_ttl = numel(tv);
        if numel(cv) >= 2
            ctrl_mean(i) = mean(cv); ctrl_sem(i) = std(cv)/sqrt(numel(cv));
        elseif numel(cv) == 1
            ctrl_mean(i) = cv; ctrl_sem(i) = 0;
        end
        stats(i).mean_ctrl = ctrl_mean(i); stats(i).sem_ctrl = ctrl_sem(i);
        if numel(tv) >= 2
            ttl_mean(i) = mean(tv); ttl_sem(i) = std(tv)/sqrt(numel(tv));
        elseif numel(tv) == 1
            ttl_mean(i) = tv; ttl_sem(i) = 0;
        end
        stats(i).mean_ttl = ttl_mean(i); stats(i).sem_ttl = ttl_sem(i);
    end
    valid_c = ~isnan(ctrl_mean); valid_t = ~isnan(ttl_mean);
    if any(valid_c)
        xc = positions(valid_c);
        fill(ax, [xc, fliplr(xc)], ...
            [ctrl_mean(valid_c)+ctrl_sem(valid_c), fliplr(ctrl_mean(valid_c)-ctrl_sem(valid_c))], ...
            col_cf, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none');
    end
    if any(valid_t)
        xt = positions(valid_t);
        fill(ax, [xt, fliplr(xt)], ...
            [ttl_mean(valid_t)+ttl_sem(valid_t), fliplr(ttl_mean(valid_t)-ttl_sem(valid_t))], ...
            col_tf, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none');
    end
    if any(valid_c)
        plot(ax, positions(valid_c), ctrl_mean(valid_c), '-o', ...
            'Color', col_c, 'MarkerFaceColor', col_c, 'MarkerSize', marker_sz, 'LineWidth', line_w);
    end
    if any(valid_t)
        plot(ax, positions(valid_t), ttl_mean(valid_t), '-o', ...
            'Color', col_t, 'MarkerFaceColor', col_t, 'MarkerSize', marker_sz, 'LineWidth', line_w);
    end
end


function format_amp(ax, y_lim, font_ax, show_xlabel, hide_yticklabels)
    xline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    xlim(ax, [-5.5 5.5]); ylim(ax, y_lim);
    set(ax, 'FontSize', font_ax, 'TickDir', 'out', 'LineWidth', 0.4, 'XTick', -4:2:4);
    box(ax, 'off');
    if show_xlabel
        xlabel(ax, 'Position', 'FontSize', font_ax);
    else
        set(ax, 'XTickLabel', []);
    end
    if hide_yticklabels, set(ax, 'YTickLabel', []); end
end


function [pool_pvals, pool_centers] = compute_pooled_ranksum(ctrl_mat, ttl_mat)
    pool_centers = -4:4;
    pool_pvals = NaN(1, 9);
    for p = 1:9
        cols = p:p+2;
        c = ctrl_mat(:, cols); c = c(:); c = c(~isnan(c));
        t = ttl_mat(:, cols);  t = t(:); t = t(~isnan(t));
        if numel(c) >= 2 && numel(t) >= 2
            pool_pvals(p) = ranksum(c, t);
        end
    end
end


function draw_pooled_asterisks(ax, pool_centers, pool_pvals, stats_11, col_stat, font_stat)
    pos11 = -5:5;
    yl = ylim(ax);
    is_hyp = yl(2) <= 0;
    for p = 1:numel(pool_centers)
        if isnan(pool_pvals(p)) || pool_pvals(p) >= 0.05, continue; end
        x = pool_centers(p);
        idx = find(pos11 == x, 1);
        if isempty(idx), continue; end
        y_c = stats_11(idx).mean_ctrl;  sem_c = stats_11(idx).sem_ctrl;
        y_t = stats_11(idx).mean_ttl;   sem_t = stats_11(idx).sem_ttl;
        if isnan(y_c) && isnan(y_t), continue; end
        if is_hyp
            y_bot = min([y_c - sem_c, y_t - sem_t], [], 'omitnan');
            y_star = y_bot - 0.08 * diff(yl); va = 'top';
        else
            y_top = max([y_c + sem_c, y_t + sem_t], [], 'omitnan');
            y_star = y_top + 0.06 * diff(yl);
            y_star = min(y_star, yl(2) - 0.02 * diff(yl)); va = 'bottom';
        end
        if pool_pvals(p) < 0.001, str = '***';
        elseif pool_pvals(p) < 0.01, str = '**';
        else, str = '*'; end
        text(ax, x, y_star, str, 'FontSize', font_stat, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', va, ...
            'Color', col_stat, 'FontWeight', 'bold', 'Clipping', 'off');
    end
end


function [fw, left_x, right_x] = compute_fwhm_positions(positions, mean_vals)
    fw = NaN; left_x = NaN; right_x = NaN;
    valid = ~isnan(mean_vals);
    if sum(valid) < 3, return; end
    pos_v = positions(valid); val_v = mean_vals(valid);
    [pk, pk_idx] = max(val_v);
    if pk <= 0, return; end
    half_max = pk / 2;
    for j = pk_idx:-1:2
        if val_v(j-1) <= half_max
            frac = (half_max - val_v(j-1)) / (val_v(j) - val_v(j-1));
            left_x = pos_v(j-1) + frac * (pos_v(j) - pos_v(j-1)); break;
        end
    end
    if isnan(left_x), left_x = pos_v(1); end
    for j = pk_idx:numel(val_v)-1
        if val_v(j+1) <= half_max
            frac = (half_max - val_v(j+1)) / (val_v(j) - val_v(j+1));
            right_x = pos_v(j+1) - frac * (pos_v(j+1) - pos_v(j)); break;
        end
    end
    if isnan(right_x), right_x = pos_v(end); end
    fw = right_x - left_x;
end


function fwhm_vec = compute_percell_fwhm(positions, amp_mat)
    n = size(amp_mat, 1); fwhm_vec = NaN(n, 1);
    for k = 1:n
        [fw, ~, ~] = compute_fwhm_positions(positions, amp_mat(k, :));
        fwhm_vec(k) = fw;
    end
end


function add_fwhm_bars(ax, positions, stats, col_c, col_t, y_lim, font_stat, ttl_tex, ...
    ctrl_dep_mat, ttl_dep_mat)
    ctrl_mn = [stats.mean_ctrl]; ttl_mn = [stats.mean_ttl];
    [fw_c, lx_c, rx_c] = compute_fwhm_positions(positions, ctrl_mn);
    [fw_t, lx_t, rx_t] = compute_fwhm_positions(positions, ttl_mn);
    fwhm_y_base = y_lim(1) + 0.5; fwhm_y_gap = 1.0;
    if ~isnan(fw_c)
        y_c = fwhm_y_base;
        plot(ax, [lx_c rx_c], [y_c y_c], '-', 'Color', col_c, 'LineWidth', 2.5);
        plot(ax, [lx_c lx_c], y_c+[-0.3 0.3], '-', 'Color', col_c, 'LineWidth', 1.0);
        plot(ax, [rx_c rx_c], y_c+[-0.3 0.3], '-', 'Color', col_c, 'LineWidth', 1.0);
    end
    if ~isnan(fw_t)
        y_t = fwhm_y_base + fwhm_y_gap;
        plot(ax, [lx_t rx_t], [y_t y_t], '-', 'Color', col_t, 'LineWidth', 2.5);
        plot(ax, [lx_t lx_t], y_t+[-0.3 0.3], '-', 'Color', col_t, 'LineWidth', 1.0);
        plot(ax, [rx_t rx_t], y_t+[-0.3 0.3], '-', 'Color', col_t, 'LineWidth', 1.0);
    end
    fwhm_parts = {};
    if ~isnan(fw_c), fwhm_parts{end+1} = sprintf('c=%.1f', fw_c); end
    if ~isnan(fw_t), fwhm_parts{end+1} = sprintf('%s=%.1f', ttl_tex, fw_t); end
    cv_fw = compute_percell_fwhm(positions, ctrl_dep_mat);
    tv_fw = compute_percell_fwhm(positions, ttl_dep_mat);
    cv_fw = cv_fw(~isnan(cv_fw)); tv_fw = tv_fw(~isnan(tv_fw));
    if numel(cv_fw) >= 2 && numel(tv_fw) >= 2
        fwhm_parts{end+1} = sprintf('p=%.3f', ranksum(cv_fw, tv_fw));
    end
    if ~isempty(fwhm_parts)
        text(ax, -5.3, fwhm_y_base + fwhm_y_gap + 1.5, ...
            ['FW: ' strjoin(fwhm_parts, ', ')], ...
            'FontSize', font_stat, 'Interpreter', 'tex', 'VerticalAlignment', 'bottom');
    end
end

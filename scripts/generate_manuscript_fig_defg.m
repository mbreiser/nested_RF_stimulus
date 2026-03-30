% GENERATE_MANUSCRIPT_FIG_DEFG  DS figure bottom half: Panels E-H.
%
%   Left side:   1D RF trace panels (1x11 tiles, mean +/- SEM)
%   Right side:  Gruntman-style amplitude-by-position (connected-line mean
%                with SEM bands, FWHM bars, pooled rank-sum asterisks)
%
%   Layout:
%     PD axis section (Panel E):
%       Left:  [T4 PD traces (1x11)]     Right:  [T4 Dep] [T5 Dep]  <- taller
%              [T5 PD traces (1x11)]              [T4 Hyp] [T5 Hyp]  <- squished
%
%     Orthogonal axis section (Panel G):
%       Left:  [T4 Ortho traces (1x11)]  Right:  [T4 Dep] [T5 Dep]
%              [T5 Ortho traces (1x11)]           [T4 Hyp] [T5 Hyp]
%
%   Traces:   thin lines (0.5), SEM bands, y capped at [-5, 25].
%   Amplitude: SEM bands, FWHM bars, pooled rank-sum asterisks on all panels.
%   Colors:   T4 ctrl=black, tutl-=red;  T5 ctrl=gray, tutl-=burgundy.
%   Scale bar: 200 ms horizontal only (bottom-left tile of G).
%
%   Usage:
%     run('scripts/generate_manuscript_fig_defg.m')

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
out_dir   = fullfile(data_root, 'manuscript_figures');
if ~isfolder(out_dir), mkdir(out_dir); end

%% Load batch results
res_file = fullfile(data_root, 'population_results', 'batch_results.mat');
fprintf('Loading: %s\n', res_file);
S = load(res_file, 'results');
results = S.results;
fprintf('Loaded %d cells.\n', numel(results));

%% ========================= Constants =====================================
FONT_NAME  = 'Helvetica';
FONT_AX    = 6;
FONT_LABEL = 7;
FONT_TITLE = 8;
FONT_STAT  = 5.5;

% T4 (ON) colours — same as Panel A in ABC
COL_T4 = struct('ctrl_line', [0 0 0],           'ttl_line', [1 0 0], ...
                'ctrl_fill', [0.80 0.80 0.80],   'ttl_fill', [1 0.70 0.70], ...
                'alpha', 0.35, 'stim_line', [0.2 0.7 0.2]);

% T5 (OFF) colours — same as Panel B in ABC
COL_T5 = struct('ctrl_line', [0.4 0.4 0.4],     'ttl_line', [0.8 0.2 0.2], ...
                'ctrl_fill', [0.70 0.70 0.70],   'ttl_fill', [0.90 0.60 0.60], ...
                'alpha', 0.35, 'stim_line', [0.2 0.7 0.2]);

STAT_COLOR = [0.15 0.15 0.6];   % dark blue asterisks
STAT_COLOR_HYP = [0.15 0.15 0.6]; % same for hyp panels
ALPHA_AMP  = 0.20;              % SEM band alpha on amplitude plots

LINE_W_TRACE = 0.5;   % thin trace lines
LINE_W_AMP   = 0.5;   % thin amplitude lines
MARKER_SZ    = 2;

% Timing / amplitude constants
STIM_ONSET  = 5001;
STIM_OFFSET = 5801;
DEP_WINDOW  = [STIM_ONSET, STIM_ONSET + 2000 - 1];
HYP_WINDOW  = [STIM_ONSET + 1000, Inf];
DEP_PCTILE  = 99.9;
HYP_PCTILE  = 0.1;
REJECT_THRESH = 0.5;

TRACE_START = 3701;  % crop x-axis to ~130 ms pre-stimulus (20 ms trimmed for tile gaps)

Y_LIM      = [-5 25];    % shared for traces + amp depol
YLIM_HYP   = [-5  0];    % amp hyperpol

STAT_METHOD = 'mean_sem';
N_POS       = 11;

TTL_TEX = '{\ittutl}^{-}';

%% Group masks
on_ctrl  = [results.is_on] & ~[results.is_ttl];
on_ttl   = [results.is_on] &  [results.is_ttl];
off_ctrl = ~[results.is_on] & ~[results.is_ttl];
off_ttl  = ~[results.is_on] &  [results.is_ttl];

n_on_c  = sum(on_ctrl);   n_on_t  = sum(on_ttl);
n_off_c = sum(off_ctrl);  n_off_t = sum(off_ttl);
fprintf('T4 (ON):  ctrl=%d  tutl-=%d\n', n_on_c, n_on_t);
fprintf('T5 (OFF): ctrl=%d  tutl-=%d\n', n_off_c, n_off_t);

%% ========================= Trace statistics (left panels) ================
% PD-axis
pd_on_ctrl_tr  = {results(on_ctrl).pd_flash_m6_aligned};
pd_on_ttl_tr   = {results(on_ttl).pd_flash_m6_aligned};
pd_off_ctrl_tr = {results(off_ctrl).pd_flash_m6_aligned};
pd_off_ttl_tr  = {results(off_ttl).pd_flash_m6_aligned};

[pd_on_c_st, pd_on_t_st]   = compute_all_stats(pd_on_ctrl_tr,  pd_on_ttl_tr,  N_POS, STAT_METHOD);
[pd_off_c_st, pd_off_t_st] = compute_all_stats(pd_off_ctrl_tr, pd_off_ttl_tr, N_POS, STAT_METHOD);

% Ortho-axis
ort_on_ctrl_tr  = {results(on_ctrl).ortho_flash_m6_aligned};
ort_on_ttl_tr   = {results(on_ttl).ortho_flash_m6_aligned};
ort_off_ctrl_tr = {results(off_ctrl).ortho_flash_m6_aligned};
ort_off_ttl_tr  = {results(off_ttl).ortho_flash_m6_aligned};

[ort_on_c_st, ort_on_t_st]   = compute_all_stats(ort_on_ctrl_tr,  ort_on_ttl_tr,  N_POS, STAT_METHOD);
[ort_off_c_st, ort_off_t_st] = compute_all_stats(ort_off_ctrl_tr, ort_off_ttl_tr, N_POS, STAT_METHOD);

%% ========================= Amplitude matrices (right panels) =============
positions = -5:5;

% PD axis
[pd_on_c_dep, pd_on_c_hyp]   = extract_robust_amps(results(on_ctrl),  'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[pd_on_t_dep, pd_on_t_hyp]   = extract_robust_amps(results(on_ttl),   'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[pd_off_c_dep, pd_off_c_hyp] = extract_robust_amps(results(off_ctrl), 'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[pd_off_t_dep, pd_off_t_hyp] = extract_robust_amps(results(off_ttl),  'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

% Ortho axis
[ort_on_c_dep, ort_on_c_hyp]   = extract_robust_amps(results(on_ctrl),  'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[ort_on_t_dep, ort_on_t_hyp]   = extract_robust_amps(results(on_ttl),   'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[ort_off_c_dep, ort_off_c_hyp] = extract_robust_amps(results(off_ctrl), 'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[ort_off_t_dep, ort_off_t_hyp] = extract_robust_amps(results(off_ttl),  'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

%% ========================= Pooled rank-sum ===============================
[pd_on_dep_pp, pool_centers]  = compute_pooled_ranksum(pd_on_c_dep,  pd_on_t_dep);
[pd_on_hyp_pp, ~]             = compute_pooled_ranksum(pd_on_c_hyp,  pd_on_t_hyp);
[pd_off_dep_pp, ~]            = compute_pooled_ranksum(pd_off_c_dep, pd_off_t_dep);
[pd_off_hyp_pp, ~]            = compute_pooled_ranksum(pd_off_c_hyp, pd_off_t_hyp);

[ort_on_dep_pp, ~]  = compute_pooled_ranksum(ort_on_c_dep,  ort_on_t_dep);
[ort_on_hyp_pp, ~]  = compute_pooled_ranksum(ort_on_c_hyp,  ort_on_t_hyp);
[ort_off_dep_pp, ~] = compute_pooled_ranksum(ort_off_c_dep, ort_off_t_dep);
[ort_off_hyp_pp, ~] = compute_pooled_ranksum(ort_off_c_hyp, ort_off_t_hyp);

%% ========================= Figure ========================================
FIG_W = 18;  FIG_H = 13;

fig = figure('Units', 'centimeters', 'Position', [2 2 FIG_W FIG_H], ...
    'PaperUnits', 'centimeters', 'PaperSize', [FIG_W FIG_H], ...
    'PaperPosition', [0 0 FIG_W FIG_H], 'Color', 'w');
set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultTextFontName', FONT_NAME);

%% ========================= Layout geometry ===============================
% --- Left: trace tiles ---
TRACE_L  = 0.06;
TRACE_R  = 0.53;
TILE_GAP = 0.001;
TILE_W   = (TRACE_R - TRACE_L - (N_POS-1)*TILE_GAP) / N_POS;

% --- Right: amplitude panels (two columns, ~30% narrower) ---
AMP_T4_L = 0.58;     % left edge of T4 column
AMP_T5_L = 0.73;     % left edge of T5 column
AMP_W    = 0.13;     % width of each amplitude column

% --- Vertical: 4 trace rows ---
ROW_H       = 0.17;
ROW_GAP     = 0.02;   % gap between T4 and T5 rows within a section
SECTION_GAP = 0.06;   % gap between PD and Ortho sections

E_T4_y = 0.72;                                 % PD section, T4 row
E_T5_y = E_T4_y - ROW_H - ROW_GAP;            % PD section, T5 row
G_T4_y = E_T5_y - SECTION_GAP - ROW_H;        % Ortho section, T4 row
G_T5_y = G_T4_y - ROW_H - ROW_GAP;            % Ortho section, T5 row

% --- Right: amplitude panel vertical positions (aligned with sections) ---
% PD section — aligned to E_T4 + E_T5
PD_top = E_T4_y + ROW_H;  PD_bot = E_T5_y;
PD_h   = PD_top - PD_bot;
PD_DEP_H = PD_h * 0.64;   PD_HYP_H = PD_h * 0.26;
PD_DEP_Y = PD_top - PD_DEP_H;
PD_HYP_Y = PD_bot;

% Ortho section — aligned to G_T4 + G_T5
ORT_top = G_T4_y + ROW_H;  ORT_bot = G_T5_y;
ORT_h   = ORT_top - ORT_bot;
ORT_DEP_H = ORT_h * 0.64;  ORT_HYP_H = ORT_h * 0.26;
ORT_DEP_Y = ORT_top - ORT_DEP_H;
ORT_HYP_Y = ORT_bot;

%% ========================= Position labels ===============================
pos_labels = arrayfun(@(x) sprintf('%+d', x), -5:5, 'UniformOutput', false);
pos_labels{6} = '0';

%% ========================= Draw TRACE panels (left) ======================

% --- Panel E: PD-axis ---
draw_trace_row(fig, TRACE_L, E_T4_y, TILE_W, ROW_H, TILE_GAP, ...
    pd_on_c_st, pd_on_t_st, N_POS, Y_LIM, COL_T4, LINE_W_TRACE, ...
    pos_labels, true, FONT_AX, FONT_LABEL, STIM_ONSET, STIM_OFFSET, TRACE_START, 'E_T4');
add_row_label(fig, TRACE_L - 0.04, E_T4_y + ROW_H/2, 'T4', FONT_TITLE);
add_n_label(fig, TRACE_R + 0.005, E_T4_y + ROW_H*0.7, ...
    pd_on_c_st(6).n, pd_on_t_st(6).n, COL_T4, FONT_AX);

draw_trace_row(fig, TRACE_L, E_T5_y, TILE_W, ROW_H, TILE_GAP, ...
    pd_off_c_st, pd_off_t_st, N_POS, Y_LIM, COL_T5, LINE_W_TRACE, ...
    {}, false, FONT_AX, FONT_LABEL, STIM_ONSET, STIM_OFFSET, TRACE_START, 'E_T5');
add_row_label(fig, TRACE_L - 0.04, E_T5_y + ROW_H/2, 'T5', FONT_TITLE);
add_n_label(fig, TRACE_R + 0.005, E_T5_y + ROW_H*0.7, ...
    pd_off_c_st(6).n, pd_off_t_st(6).n, COL_T5, FONT_AX);

% --- Panel G: Ortho-axis ---
draw_trace_row(fig, TRACE_L, G_T4_y, TILE_W, ROW_H, TILE_GAP, ...
    ort_on_c_st, ort_on_t_st, N_POS, Y_LIM, COL_T4, LINE_W_TRACE, ...
    pos_labels, true, FONT_AX, FONT_LABEL, STIM_ONSET, STIM_OFFSET, TRACE_START, 'G_T4');
add_row_label(fig, TRACE_L - 0.04, G_T4_y + ROW_H/2, 'T4', FONT_TITLE);
add_n_label(fig, TRACE_R + 0.005, G_T4_y + ROW_H*0.7, ...
    ort_on_c_st(6).n, ort_on_t_st(6).n, COL_T4, FONT_AX);

draw_trace_row(fig, TRACE_L, G_T5_y, TILE_W, ROW_H, TILE_GAP, ...
    ort_off_c_st, ort_off_t_st, N_POS, Y_LIM, COL_T5, LINE_W_TRACE, ...
    {}, false, FONT_AX, FONT_LABEL, STIM_ONSET, STIM_OFFSET, TRACE_START, 'G_T5');
add_row_label(fig, TRACE_L - 0.04, G_T5_y + ROW_H/2, 'T5', FONT_TITLE);
add_n_label(fig, TRACE_R + 0.005, G_T5_y + ROW_H*0.7, ...
    ort_off_c_st(6).n, ort_off_t_st(6).n, COL_T5, FONT_AX);

%% ========================= Draw AMPLITUDE panels (right) =================

% ===== PD axis: Depolarization =====
% T4 (ON)
ax = axes(fig, 'Position', [AMP_T4_L, PD_DEP_Y, AMP_W, PD_DEP_H]);
dep_st_pd_t4 = draw_amp_line(ax, positions, pd_on_c_dep, pd_on_t_dep, ...
    COL_T4.ctrl_line, COL_T4.ttl_line, LINE_W_AMP, MARKER_SZ, FONT_AX, ...
    COL_T4.ctrl_fill, COL_T4.ttl_fill, ALPHA_AMP);
draw_pooled_asterisks(ax, pool_centers, pd_on_dep_pp, dep_st_pd_t4, STAT_COLOR, FONT_STAT);
format_amp(ax, Y_LIM, FONT_AX, false, false);
ylabel(ax, 'mV', 'FontSize', FONT_LABEL);
add_fwhm_bars(ax, positions, dep_st_pd_t4, COL_T4.ctrl_line, COL_T4.ttl_line, ...
    Y_LIM, FONT_STAT, TTL_TEX, pd_on_c_dep, pd_on_t_dep);

% T5 (OFF)
ax = axes(fig, 'Position', [AMP_T5_L, PD_DEP_Y, AMP_W, PD_DEP_H]);
dep_st_pd_t5 = draw_amp_line(ax, positions, pd_off_c_dep, pd_off_t_dep, ...
    COL_T5.ctrl_line, COL_T5.ttl_line, LINE_W_AMP, MARKER_SZ, FONT_AX, ...
    COL_T5.ctrl_fill, COL_T5.ttl_fill, ALPHA_AMP);
draw_pooled_asterisks(ax, pool_centers, pd_off_dep_pp, dep_st_pd_t5, STAT_COLOR, FONT_STAT);
format_amp(ax, Y_LIM, FONT_AX, false, true);
add_fwhm_bars(ax, positions, dep_st_pd_t5, COL_T5.ctrl_line, COL_T5.ttl_line, ...
    Y_LIM, FONT_STAT, TTL_TEX, pd_off_c_dep, pd_off_t_dep);

% ===== PD axis: Hyperpolarization =====
% T4 (ON)
ax = axes(fig, 'Position', [AMP_T4_L, PD_HYP_Y, AMP_W, PD_HYP_H]);
hyp_st_pd_t4 = draw_amp_line(ax, positions, pd_on_c_hyp, pd_on_t_hyp, ...
    COL_T4.ctrl_line, COL_T4.ttl_line, LINE_W_AMP, MARKER_SZ, FONT_AX, ...
    COL_T4.ctrl_fill, COL_T4.ttl_fill, ALPHA_AMP);
draw_pooled_asterisks(ax, pool_centers, pd_on_hyp_pp, hyp_st_pd_t4, STAT_COLOR, FONT_STAT);
format_amp(ax, YLIM_HYP, FONT_AX, false, false);
ylabel(ax, 'mV', 'FontSize', FONT_LABEL);

% T5 (OFF)
ax = axes(fig, 'Position', [AMP_T5_L, PD_HYP_Y, AMP_W, PD_HYP_H]);
hyp_st_pd_t5 = draw_amp_line(ax, positions, pd_off_c_hyp, pd_off_t_hyp, ...
    COL_T5.ctrl_line, COL_T5.ttl_line, LINE_W_AMP, MARKER_SZ, FONT_AX, ...
    COL_T5.ctrl_fill, COL_T5.ttl_fill, ALPHA_AMP);
draw_pooled_asterisks(ax, pool_centers, pd_off_hyp_pp, hyp_st_pd_t5, STAT_COLOR, FONT_STAT);
format_amp(ax, YLIM_HYP, FONT_AX, false, true);

% ===== Ortho axis: Depolarization =====
% T4 (ON)
ax = axes(fig, 'Position', [AMP_T4_L, ORT_DEP_Y, AMP_W, ORT_DEP_H]);
dep_st_ort_t4 = draw_amp_line(ax, positions, ort_on_c_dep, ort_on_t_dep, ...
    COL_T4.ctrl_line, COL_T4.ttl_line, LINE_W_AMP, MARKER_SZ, FONT_AX, ...
    COL_T4.ctrl_fill, COL_T4.ttl_fill, ALPHA_AMP);
draw_pooled_asterisks(ax, pool_centers, ort_on_dep_pp, dep_st_ort_t4, STAT_COLOR, FONT_STAT);
format_amp(ax, Y_LIM, FONT_AX, false, false);
ylabel(ax, 'mV', 'FontSize', FONT_LABEL);
add_fwhm_bars(ax, positions, dep_st_ort_t4, COL_T4.ctrl_line, COL_T4.ttl_line, ...
    Y_LIM, FONT_STAT, TTL_TEX, ort_on_c_dep, ort_on_t_dep);

% T5 (OFF)
ax = axes(fig, 'Position', [AMP_T5_L, ORT_DEP_Y, AMP_W, ORT_DEP_H]);
dep_st_ort_t5 = draw_amp_line(ax, positions, ort_off_c_dep, ort_off_t_dep, ...
    COL_T5.ctrl_line, COL_T5.ttl_line, LINE_W_AMP, MARKER_SZ, FONT_AX, ...
    COL_T5.ctrl_fill, COL_T5.ttl_fill, ALPHA_AMP);
draw_pooled_asterisks(ax, pool_centers, ort_off_dep_pp, dep_st_ort_t5, STAT_COLOR, FONT_STAT);
format_amp(ax, Y_LIM, FONT_AX, false, true);
add_fwhm_bars(ax, positions, dep_st_ort_t5, COL_T5.ctrl_line, COL_T5.ttl_line, ...
    Y_LIM, FONT_STAT, TTL_TEX, ort_off_c_dep, ort_off_t_dep);

% ===== Ortho axis: Hyperpolarization =====
% T4 (ON)
ax = axes(fig, 'Position', [AMP_T4_L, ORT_HYP_Y, AMP_W, ORT_HYP_H]);
hyp_st_ort_t4 = draw_amp_line(ax, positions, ort_on_c_hyp, ort_on_t_hyp, ...
    COL_T4.ctrl_line, COL_T4.ttl_line, LINE_W_AMP, MARKER_SZ, FONT_AX, ...
    COL_T4.ctrl_fill, COL_T4.ttl_fill, ALPHA_AMP);
draw_pooled_asterisks(ax, pool_centers, ort_on_hyp_pp, hyp_st_ort_t4, STAT_COLOR, FONT_STAT);
format_amp(ax, YLIM_HYP, FONT_AX, true, false);
ylabel(ax, 'mV', 'FontSize', FONT_LABEL);

% T5 (OFF)
ax = axes(fig, 'Position', [AMP_T5_L, ORT_HYP_Y, AMP_W, ORT_HYP_H]);
hyp_st_ort_t5 = draw_amp_line(ax, positions, ort_off_c_hyp, ort_off_t_hyp, ...
    COL_T5.ctrl_line, COL_T5.ttl_line, LINE_W_AMP, MARKER_SZ, FONT_AX, ...
    COL_T5.ctrl_fill, COL_T5.ttl_fill, ALPHA_AMP);
draw_pooled_asterisks(ax, pool_centers, ort_off_hyp_pp, hyp_st_ort_t5, STAT_COLOR, FONT_STAT);
format_amp(ax, YLIM_HYP, FONT_AX, true, true);

%% ========================= Diagnostic: p-value summary ===================
fprintf('\nAsterisk summary (p<0.05 count):\n');
fprintf('  PD T4 dep:  %d sig\n', sum(pd_on_dep_pp < 0.05));
fprintf('  PD T5 dep:  %d sig\n', sum(pd_off_dep_pp < 0.05));
fprintf('  PD T4 hyp:  %d sig\n', sum(pd_on_hyp_pp < 0.05));
fprintf('  PD T5 hyp:  %d sig\n', sum(pd_off_hyp_pp < 0.05));
fprintf('  Ort T4 dep: %d sig\n', sum(ort_on_dep_pp < 0.05));
fprintf('  Ort T5 dep: %d sig\n', sum(ort_off_dep_pp < 0.05));
fprintf('  Ort T4 hyp: %d sig\n', sum(ort_on_hyp_pp < 0.05));
fprintf('  Ort T5 hyp: %d sig\n', sum(ort_off_hyp_pp < 0.05));

%% ========================= Labels & annotations ==========================
lbl_fsize = 12;

% Panel E (PD traces)
annotation(fig, 'textbox', [0.00, E_T4_y + ROW_H - 0.02, 0.04, 0.06], ...
    'String', 'E', 'FontSize', lbl_fsize, 'FontWeight', 'bold', ...
    'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Panel F (PD amplitude)
annotation(fig, 'textbox', [AMP_T4_L - 0.04, PD_DEP_Y + PD_DEP_H - 0.02, 0.04, 0.06], ...
    'String', 'F', 'FontSize', lbl_fsize, 'FontWeight', 'bold', ...
    'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Panel G (Ortho traces)
annotation(fig, 'textbox', [0.00, G_T4_y + ROW_H - 0.02, 0.04, 0.06], ...
    'String', 'G', 'FontSize', lbl_fsize, 'FontWeight', 'bold', ...
    'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Panel H (Ortho amplitude)
annotation(fig, 'textbox', [AMP_T4_L - 0.04, ORT_DEP_Y + ORT_DEP_H - 0.02, 0.04, 0.06], ...
    'String', 'H', 'FontSize', lbl_fsize, 'FontWeight', 'bold', ...
    'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Section headers (above trace area)
annotation(fig, 'textbox', [TRACE_L, E_T4_y + ROW_H + 0.005, 0.30, 0.03], ...
    'String', 'PD axis', 'FontSize', FONT_TITLE, ...
    'FontWeight', 'bold', 'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

annotation(fig, 'textbox', [TRACE_L, G_T4_y + ROW_H + 0.005, 0.30, 0.03], ...
    'String', 'Orthogonal axis', 'FontSize', FONT_TITLE, ...
    'FontWeight', 'bold', 'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

% Column headers over amplitude panels
annotation(fig, 'textbox', [AMP_T4_L, PD_DEP_Y + PD_DEP_H + 0.005, AMP_W, 0.03], ...
    'String', 'T4 (ON)', 'FontSize', FONT_TITLE, ...
    'FontWeight', 'bold', 'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
annotation(fig, 'textbox', [AMP_T5_L, PD_DEP_Y + PD_DEP_H + 0.005, AMP_W, 0.03], ...
    'String', 'T5 (OFF)', 'FontSize', FONT_TITLE, ...
    'FontWeight', 'bold', 'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

%% ========================= Scale bar =====================================
% 200 ms horizontal bar only, on bottom-left trace tile (G, T5, position 1)
ax_sb = findobj(fig, 'Type', 'axes', 'Tag', 'G_T5_1');
if ~isempty(ax_sb)
    add_horiz_scale_bar(ax_sb(1), Y_LIM, FONT_AX, TRACE_START);
end

%% ========================= Export ========================================
ts = datestr(now, 'yyyymmdd_HHMM');
pdf_file = fullfile(out_dir, sprintf('fig_ds_panels_EFGH_%s.pdf', ts));
png_file = fullfile(out_dir, sprintf('fig_ds_panels_EFGH_%s.png', ts));

exportgraphics(fig, pdf_file, 'ContentType', 'vector');
exportgraphics(fig, png_file, 'Resolution', 300);
fprintf('Saved: %s\n', pdf_file);
fprintf('Saved: %s\n', png_file);
fprintf('=== Done ===\n');


%% =========================================================================
%%                          LOCAL FUNCTIONS
%% =========================================================================

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
        if ~isempty(row_tag)
            ax.Tag = sprintf('%s_%d', row_tag, pi);
        end
        hold(ax, 'on');

        ctr_c = ctrl_stats(pi).center;  spr_c = ctrl_stats(pi).spread;
        ctr_t = ttl_stats(pi).center;   spr_t = ttl_stats(pi).spread;

        x_len = 0;
        if ~isempty(ctr_c), x_len = max(x_len, numel(ctr_c)); end
        if ~isempty(ctr_t), x_len = max(x_len, numel(ctr_t)); end
        if x_len == 0
            ylim(ax, y_lim); set(ax, 'XTick', [], 'YTick', []); box(ax, 'off'); continue;
        end

        % Green stimulus rectangle (below traces, from -2 to -1 mV)
        fill(ax, [stim_onset, stim_offset, stim_offset, stim_onset], ...
            [-2, -2, -1, -1], col.stim_line, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

        % Ctrl: SEM band + mean line
        xv = 1:x_len;
        if ~isempty(ctr_c) && ctrl_stats(pi).n >= 2
            fill(ax, [xv, fliplr(xv)], [ctr_c + spr_c, fliplr(ctr_c - spr_c)], ...
                col.ctrl_fill, 'FaceAlpha', col.alpha, 'EdgeColor', 'none');
            plot(ax, xv, ctr_c, '-', 'Color', col.ctrl_line, 'LineWidth', line_w);
        end

        % TTL: SEM band + mean line
        if ~isempty(ctr_t) && ttl_stats(pi).n >= 2
            fill(ax, [xv, fliplr(xv)], [ctr_t + spr_t, fliplr(ctr_t - spr_t)], ...
                col.ttl_fill, 'FaceAlpha', col.alpha, 'EdgeColor', 'none');
            plot(ax, xv, ctr_t, '-', 'Color', col.ttl_line, 'LineWidth', line_w);
        end

        ylim(ax, y_lim);
        xlim(ax, [trace_start, x_len]);
        set(ax, 'XTick', []);
        if pi == 1
            ylabel(ax, '\DeltamV', 'FontSize', label_font);
            set(ax, 'FontSize', font_size, 'TickDir', 'out', 'TickLength', [0.02 0.02]);
        else
            set(ax, 'YTick', []);   % remove all ticks for non-leftmost tiles
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
        'Color', col.ctrl_line * 0.6 + 0.4);  % muted version of ctrl color
end


function add_horiz_scale_bar(ax, y_lim, font_size, trace_start)
% Horizontal-only 200 ms scale bar (2000 samples at 10 kHz).
    v_range = diff(y_lim);
    bar_t = 2000;  % 200 ms

    x_left   = trace_start + 0.05 * (xlim(ax) * [0;1] - trace_start);
    y_bottom = y_lim(1) + 0.08 * v_range;

    plot(ax, [x_left, x_left + bar_t], [y_bottom, y_bottom], 'k-', ...
        'LineWidth', 0.8, 'Clipping', 'off');
    text(ax, x_left, y_bottom - 0.06*v_range, '200 ms', ...
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
% Draw Gruntman-style connected-line amplitude with SEM bands.
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

    % SEM bands
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

    % Connected mean lines
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
    is_hyp = yl(2) <= 0;   % hyperpolarization panel (all values ≤ 0)
    for p = 1:numel(pool_centers)
        if isnan(pool_pvals(p)) || pool_pvals(p) >= 0.05, continue; end
        x = pool_centers(p);
        idx = find(pos11 == x, 1);
        if isempty(idx), continue; end
        y_c = stats_11(idx).mean_ctrl;  sem_c = stats_11(idx).sem_ctrl;
        y_t = stats_11(idx).mean_ttl;   sem_t = stats_11(idx).sem_ttl;
        if isnan(y_c) && isnan(y_t), continue; end
        if is_hyp
            % Place asterisks below the data (more negative)
            y_bot = min([y_c - sem_c, y_t - sem_t], [], 'omitnan');
            y_star = y_bot - 0.08 * diff(yl);
            va = 'top';
        else
            % Place asterisks above the data
            y_top = max([y_c + sem_c, y_t + sem_t], [], 'omitnan');
            y_star = y_top + 0.06 * diff(yl);
            % Cap so it doesn't go above ylim
            y_star = min(y_star, yl(2) - 0.02 * diff(yl));
            va = 'bottom';
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

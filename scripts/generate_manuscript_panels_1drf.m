% GENERATE_MANUSCRIPT_PANELS_1DRF  Manuscript-quality 1D RF panels.
%
%   Generates standalone panels for the bottom half of the DS figure:
%     Panel D/E (PD):    M6-aligned 1D RF traces (PD axis)
%     Panel D/E (Ortho): M6-aligned 1D RF traces (Orthogonal axis)
%     Panel F/G:         Per-position amplitudes, PD (Gruntman connected-line style)
%     Panel F2/G2:       Per-position amplitudes, Ortho (Gruntman connected-line style)
%   All panels: Combined T4 (ON) + T5 (OFF), mean +/- SEM,
%               per-position Wilcoxon rank-sum (uncorrected).
%
%   Uses the 25-cell late dataset (batch_results.mat).
%   Alignment: M6 (68%-area centroid) for both PD and orthogonal axes.
%
%   Outputs (in manuscript_figures/):
%     fig_ds_panel_DE_aligned_1drf_ON_OFF_<ts>.png
%     fig_ds_panel_DE_aligned_1drf_ON_OFF_ortho_<ts>.png
%     fig_ds_panel_F_amplitudes_ON_<ts>.png
%     fig_ds_panel_G_amplitudes_OFF_<ts>.png
%     fig_ds_panel_F2_amplitudes_ON_ortho_<ts>.png
%     fig_ds_panel_G2_amplitudes_OFF_ortho_<ts>.png
%     fig_ds_1drf_stats.txt  — per-position stats summary
%
%   Usage:
%     run('scripts/generate_manuscript_panels_1drf.m')

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
on_mask  = [results.is_on];
off_mask = ~on_mask;
ctrl_mask = ~[results.is_ttl];
ttl_mask  = [results.is_ttl];

n_on_ctrl  = sum(on_mask & ctrl_mask);
n_on_ttl   = sum(on_mask & ttl_mask);
n_off_ctrl = sum(off_mask & ctrl_mask);
n_off_ttl  = sum(off_mask & ttl_mask);
fprintf('  ON ctrl=%d, ON tutl-=%d, OFF ctrl=%d, OFF tutl-=%d\n', ...
    n_on_ctrl, n_on_ttl, n_off_ctrl, n_off_ttl);

%% Manuscript defaults
FONT_NAME    = 'Helvetica';
FONT_AX      = 7;     % axis tick labels
FONT_LABEL   = 8;     % axis labels
FONT_TITLE   = 9;     % panel titles
FONT_STAT    = 6;     % stat annotations

COL_CTRL     = [0 0 0];           % black
COL_CTRL_F   = [0.75 0.75 0.75];  % light gray fill
COL_TTL      = [0.85 0 0];        % dark red
COL_TTL_F    = [1 0.70 0.70];     % light red fill
ALPHA_FILL   = 0.20;              % light SEM shading (Gruntman style)

COL_STIM     = [0.2 0.7 0.2];     % green stim lines
COL_RESP_END = [0.3 0.3 0.3];     % gray response end line
STIM_ON      = 5001;
STIM_OFF     = 5801;
RESP_END     = 6551;
X_START      = 4001;              % 100ms pre-stim (crop long baseline)

STAT_COLOR   = [0.15 0.15 0.6];   % dark blue for asterisks

% Shared y-axis range for all depol panels (DE traces + F/G depol)
YLIM_SHARED  = [-5, 25];

% Label helpers
TTL_TEX = '{\ittutl}^{-}';

%% ========================================================================
%  Panels D/E: M6-aligned 1D RF traces (PD axis only)
%  Combined figure: top row = T4 (ON), bottom row = T5 (OFF)
%  ~150mm wide × ~70mm tall (Gruntman style: thin lines, wide aspect)
% =========================================================================

fprintf('\n=== Generating combined ON+OFF aligned 1D RF panel ===\n');

% Open stats file
stats_file = fullfile(out_dir, 'fig_ds_1drf_stats.txt');
fid = fopen(stats_file, 'w');
fprintf(fid, 'Per-position statistics for aligned 1D RF panels\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

% --- Extract PD traces for both ON and OFF ---
on_ctrl_sel  = on_mask & ctrl_mask;
on_ttl_sel   = on_mask & ttl_mask;
off_ctrl_sel = off_mask & ctrl_mask;
off_ttl_sel  = off_mask & ttl_mask;

on_pd_ctrl  = {results(on_ctrl_sel).pd_flash_m6_aligned};
on_pd_ttl   = {results(on_ttl_sel).pd_flash_m6_aligned};
off_pd_ctrl = {results(off_ctrl_sel).pd_flash_m6_aligned};
off_pd_ttl  = {results(off_ttl_sel).pd_flash_m6_aligned};

fprintf('ON (T4):  ctrl n=%d, tutl- n=%d\n', n_on_ctrl, n_on_ttl);
fprintf('OFF (T5): ctrl n=%d, tutl- n=%d\n', n_off_ctrl, n_off_ttl);

% Fixed y-limits for consistent range across all panels
y_lim = YLIM_SHARED;

% Figure: 15 cm wide × 7 cm tall (narrower → tiles closer together)
fig = figure('Units', 'centimeters', 'Position', [1 2 15 7], ...
    'Color', 'w', 'PaperUnits', 'centimeters', ...
    'PaperSize', [15 7], 'PaperPosition', [0 0 15 7]);
set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultTextFontName', FONT_NAME);

% 2 rows × 11 columns: top=ON PD, bottom=OFF PD
% 'tight' spacing → small gaps between tiles (~100ms equiv)
tl = tiledlayout(fig, 2, 11, 'TileSpacing', 'tight', 'Padding', 'compact');

% Position labels (top row only)
pos_labels = cell(1, 11);
for p = 1:11
    if p == 1,      pos_labels{p} = 'Distal';
    elseif p == 6,  pos_labels{p} = '0 (M6)';
    elseif p == 11, pos_labels{p} = 'Proximal';
    else,            pos_labels{p} = num2str(p - 6);
    end
end

% Per-position Wilcoxon tests
on_pd_stats  = compute_position_stats(on_pd_ctrl, on_pd_ttl, STIM_ON, RESP_END);
off_pd_stats = compute_position_stats(off_pd_ctrl, off_pd_ttl, STIM_ON, RESP_END);

% --- Top row: T4 (ON) PD traces ---
for p = 1:11
    ax = nexttile(tl, p);
    hold(ax, 'on');
    draw_1drf_tile(ax, on_pd_ctrl, on_pd_ttl, p, y_lim, ...
        COL_CTRL, COL_CTRL_F, COL_TTL, COL_TTL_F, ALPHA_FILL, ...
        STIM_ON, STIM_OFF, RESP_END, COL_STIM, COL_RESP_END, X_START);

    % Asterisk
    if on_pd_stats(p).sig
        draw_position_asterisk(ax, on_pd_stats(p).p, y_lim, STAT_COLOR, FONT_STAT);
    end

    % Axis styling: remove all lines except leftmost y-axis
    set(ax, 'FontSize', FONT_AX, 'XTick', [], 'LineWidth', 0.4);
    ax.XAxis.Visible = 'off';
    if p == 1
        ylabel(ax, '\DeltamV', 'FontSize', FONT_LABEL);
        ax.YAxis.Visible = 'on';
        ax.TickDir = 'out';
    else
        ax.YAxis.Visible = 'off';
    end
    box(ax, 'off');

    % Position labels on top row
    title(ax, pos_labels{p}, 'FontSize', FONT_STAT, 'FontWeight', 'normal');
    if p == 6
        title(ax, pos_labels{p}, 'FontSize', FONT_STAT, 'FontWeight', 'bold');
    end
end

% Row label: T4 (ON)
annotation(fig, 'textbox', [0.0, 0.55, 0.04, 0.20], ...
    'String', 'T4 (ON)', 'FontSize', FONT_LABEL, 'FontWeight', 'bold', ...
    'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'Rotation', 90);

% --- Bottom row: T5 (OFF) PD traces ---
for p = 1:11
    ax = nexttile(tl, 11 + p);
    hold(ax, 'on');
    draw_1drf_tile(ax, off_pd_ctrl, off_pd_ttl, p, y_lim, ...
        COL_CTRL, COL_CTRL_F, COL_TTL, COL_TTL_F, ALPHA_FILL, ...
        STIM_ON, STIM_OFF, RESP_END, COL_STIM, COL_RESP_END, X_START);

    % Asterisk
    if off_pd_stats(p).sig
        draw_position_asterisk(ax, off_pd_stats(p).p, y_lim, STAT_COLOR, FONT_STAT);
    end

    % Axis styling: remove all lines except leftmost y-axis
    set(ax, 'FontSize', FONT_AX, 'XTick', [], 'LineWidth', 0.4);
    ax.XAxis.Visible = 'off';
    if p == 1
        ylabel(ax, '\DeltamV', 'FontSize', FONT_LABEL);
        ax.YAxis.Visible = 'on';
        ax.TickDir = 'out';

        % Time scale bar: 100 ms = 1000 samples, bottom-left
        yl = ylim(ax);
        xl = xlim(ax);
        bar_y = yl(1) + 0.08 * diff(yl);
        bar_x0 = xl(1) + 0.02 * diff(xl);
        bar_x1 = bar_x0 + 1000;  % 100 ms at 10 kHz
        plot(ax, [bar_x0 bar_x1], [bar_y bar_y], '-k', 'LineWidth', 1.0);
        text(ax, (bar_x0 + bar_x1)/2, bar_y - 0.06*diff(yl), '100 ms', ...
            'FontSize', FONT_STAT, 'FontName', FONT_NAME, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    else
        ax.YAxis.Visible = 'off';
    end
    box(ax, 'off');
end

% Row label: T5 (OFF)
annotation(fig, 'textbox', [0.0, 0.08, 0.04, 0.20], ...
    'String', 'T5 (OFF)', 'FontSize', FONT_LABEL, 'FontWeight', 'bold', ...
    'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'Rotation', 90);

% Legend in top-right corner (ON n-counts)
ax_leg = nexttile(tl, 11);
yl = ylim(ax_leg);
xl = xlim(ax_leg);
text(ax_leg, xl(2), yl(2), sprintf('ctrl (n=%d)', n_on_ctrl), ...
    'FontSize', FONT_STAT, 'Color', COL_CTRL, 'FontName', FONT_NAME, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
text(ax_leg, xl(2), yl(2) - 0.18*diff(yl), sprintf('%s (n=%d)', TTL_TEX, n_on_ttl), ...
    'FontSize', FONT_STAT, 'Color', COL_TTL, 'FontName', FONT_NAME, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'Interpreter', 'tex');

% OFF n-counts in bottom-right corner
ax_leg2 = nexttile(tl, 22);
yl2 = ylim(ax_leg2);
xl2 = xlim(ax_leg2);
text(ax_leg2, xl2(2), yl2(2), sprintf('ctrl (n=%d)', n_off_ctrl), ...
    'FontSize', FONT_STAT, 'Color', COL_CTRL, 'FontName', FONT_NAME, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
text(ax_leg2, xl2(2), yl2(2) - 0.18*diff(yl2), sprintf('%s (n=%d)', TTL_TEX, n_off_ttl), ...
    'FontSize', FONT_STAT, 'Color', COL_TTL, 'FontName', FONT_NAME, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'Interpreter', 'tex');

% Panel title
title(tl, 'M6-aligned 1D RF (PD axis)', ...
    'FontSize', FONT_TITLE, 'FontWeight', 'bold', 'FontName', FONT_NAME);

% Export
export_panel(fig, out_dir, 'fig_ds_panel_DE_aligned_1drf_ON_OFF');
close(fig);

% Write stats
fprintf(fid, '=== Panel D/E: M6-aligned 1D RF (PD axis) ===\n\n');
fprintf(fid, 'ON (T4): ctrl n=%d, tutl- n=%d\n', n_on_ctrl, n_on_ttl);
write_pos_stats(fid, 'ON PD axis', on_pd_stats, pos_labels);
fprintf(fid, 'OFF (T5): ctrl n=%d, tutl- n=%d\n', n_off_ctrl, n_off_ttl);
write_pos_stats(fid, 'OFF PD axis', off_pd_stats, pos_labels);

%% ========================================================================
%  Panel D/E (ortho): M6-aligned 1D RF — orthogonal (OD) axis
%  Same layout as PD panel above, using ortho_flash_m6_aligned traces
% =========================================================================

fprintf('\n=== Generating M6-aligned trace panels (OD axis) ===\n');

% --- Extract ortho traces for both ON and OFF ---
on_od_ctrl  = {results(on_ctrl_sel).ortho_flash_m6_aligned};
on_od_ttl   = {results(on_ttl_sel).ortho_flash_m6_aligned};
off_od_ctrl = {results(off_ctrl_sel).ortho_flash_m6_aligned};
off_od_ttl  = {results(off_ttl_sel).ortho_flash_m6_aligned};

% Fixed y-limits — same as PD panel for visual consistency
y_lim_od = YLIM_SHARED;

% Figure: same dimensions as PD panel
fig_od = figure('Units', 'centimeters', 'Position', [1 2 15 7], ...
    'Color', 'w', 'PaperUnits', 'centimeters', ...
    'PaperSize', [15 7], 'PaperPosition', [0 0 15 7]);
set(fig_od, 'DefaultAxesFontName', FONT_NAME, 'DefaultTextFontName', FONT_NAME);

tl_od = tiledlayout(fig_od, 2, 11, 'TileSpacing', 'tight', 'Padding', 'compact');

% Position labels for ortho axis (numbered, no Proximal/Distal)
od_labels = cell(1, 11);
for p = 1:11
    if p == 6,  od_labels{p} = '0 (M6)';
    else,        od_labels{p} = num2str(p - 6);
    end
end

% Per-position Wilcoxon tests
on_od_stats  = compute_position_stats(on_od_ctrl, on_od_ttl, STIM_ON, RESP_END);
off_od_stats = compute_position_stats(off_od_ctrl, off_od_ttl, STIM_ON, RESP_END);

% --- Top row: T4 (ON) ortho traces ---
for p = 1:11
    ax = nexttile(tl_od, p);
    hold(ax, 'on');
    draw_1drf_tile(ax, on_od_ctrl, on_od_ttl, p, y_lim_od, ...
        COL_CTRL, COL_CTRL_F, COL_TTL, COL_TTL_F, ALPHA_FILL, ...
        STIM_ON, STIM_OFF, RESP_END, COL_STIM, COL_RESP_END, X_START);

    % Asterisk
    if on_od_stats(p).sig
        draw_position_asterisk(ax, on_od_stats(p).p, y_lim_od, STAT_COLOR, FONT_STAT);
    end

    % Axis styling
    set(ax, 'FontSize', FONT_AX, 'XTick', [], 'LineWidth', 0.4);
    ax.XAxis.Visible = 'off';
    if p == 1
        ylabel(ax, '\DeltamV', 'FontSize', FONT_LABEL);
        ax.YAxis.Visible = 'on';
        ax.TickDir = 'out';
    else
        ax.YAxis.Visible = 'off';
    end
    box(ax, 'off');

    % Position labels on top row
    title(ax, od_labels{p}, 'FontSize', FONT_STAT, 'FontWeight', 'normal');
    if p == 6
        title(ax, od_labels{p}, 'FontSize', FONT_STAT, 'FontWeight', 'bold');
    end
end

% Row label: T4 (ON)
annotation(fig_od, 'textbox', [0.0, 0.55, 0.04, 0.20], ...
    'String', 'T4 (ON)', 'FontSize', FONT_LABEL, 'FontWeight', 'bold', ...
    'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'Rotation', 90);

% --- Bottom row: T5 (OFF) ortho traces ---
for p = 1:11
    ax = nexttile(tl_od, 11 + p);
    hold(ax, 'on');
    draw_1drf_tile(ax, off_od_ctrl, off_od_ttl, p, y_lim_od, ...
        COL_CTRL, COL_CTRL_F, COL_TTL, COL_TTL_F, ALPHA_FILL, ...
        STIM_ON, STIM_OFF, RESP_END, COL_STIM, COL_RESP_END, X_START);

    % Asterisk
    if off_od_stats(p).sig
        draw_position_asterisk(ax, off_od_stats(p).p, y_lim_od, STAT_COLOR, FONT_STAT);
    end

    % Axis styling
    set(ax, 'FontSize', FONT_AX, 'XTick', [], 'LineWidth', 0.4);
    ax.XAxis.Visible = 'off';
    if p == 1
        ylabel(ax, '\DeltamV', 'FontSize', FONT_LABEL);
        ax.YAxis.Visible = 'on';
        ax.TickDir = 'out';

        % Time scale bar: 100 ms = 1000 samples, bottom-left
        yl = ylim(ax);
        xl = xlim(ax);
        bar_y = yl(1) + 0.08 * diff(yl);
        bar_x0 = xl(1) + 0.02 * diff(xl);
        bar_x1 = bar_x0 + 1000;  % 100 ms at 10 kHz
        plot(ax, [bar_x0 bar_x1], [bar_y bar_y], '-k', 'LineWidth', 1.0);
        text(ax, (bar_x0 + bar_x1)/2, bar_y - 0.06*diff(yl), '100 ms', ...
            'FontSize', FONT_STAT, 'FontName', FONT_NAME, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    else
        ax.YAxis.Visible = 'off';
    end
    box(ax, 'off');
end

% Row label: T5 (OFF)
annotation(fig_od, 'textbox', [0.0, 0.08, 0.04, 0.20], ...
    'String', 'T5 (OFF)', 'FontSize', FONT_LABEL, 'FontWeight', 'bold', ...
    'FontName', FONT_NAME, 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'Rotation', 90);

% Legend in top-right corner (ON n-counts)
ax_leg_od = nexttile(tl_od, 11);
yl_od = ylim(ax_leg_od);
xl_od = xlim(ax_leg_od);
text(ax_leg_od, xl_od(2), yl_od(2), sprintf('ctrl (n=%d)', n_on_ctrl), ...
    'FontSize', FONT_STAT, 'Color', COL_CTRL, 'FontName', FONT_NAME, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
text(ax_leg_od, xl_od(2), yl_od(2) - 0.18*diff(yl_od), sprintf('%s (n=%d)', TTL_TEX, n_on_ttl), ...
    'FontSize', FONT_STAT, 'Color', COL_TTL, 'FontName', FONT_NAME, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'Interpreter', 'tex');

% OFF n-counts in bottom-right corner
ax_leg2_od = nexttile(tl_od, 22);
yl2_od = ylim(ax_leg2_od);
xl2_od = xlim(ax_leg2_od);
text(ax_leg2_od, xl2_od(2), yl2_od(2), sprintf('ctrl (n=%d)', n_off_ctrl), ...
    'FontSize', FONT_STAT, 'Color', COL_CTRL, 'FontName', FONT_NAME, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
text(ax_leg2_od, xl2_od(2), yl2_od(2) - 0.18*diff(yl2_od), sprintf('%s (n=%d)', TTL_TEX, n_off_ttl), ...
    'FontSize', FONT_STAT, 'Color', COL_TTL, 'FontName', FONT_NAME, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'Interpreter', 'tex');

% Panel title
title(tl_od, 'M6-aligned 1D RF (Orthogonal axis)', ...
    'FontSize', FONT_TITLE, 'FontWeight', 'bold', 'FontName', FONT_NAME);

% Export
export_panel(fig_od, out_dir, 'fig_ds_panel_DE_aligned_1drf_ON_OFF_ortho');
close(fig_od);

% Write ortho stats
fprintf(fid, '\n=== Panel D/E: M6-aligned 1D RF (Orthogonal axis) ===\n\n');
fprintf(fid, 'ON (T4): ctrl n=%d, tutl- n=%d\n', n_on_ctrl, n_on_ttl);
write_pos_stats(fid, 'ON Ortho axis', on_od_stats, od_labels);
fprintf(fid, 'OFF (T5): ctrl n=%d, tutl- n=%d\n', n_off_ctrl, n_off_ttl);
write_pos_stats(fid, 'OFF Ortho axis', off_od_stats, od_labels);

%% ========================================================================
%  Panels F/G: Position-by-position amplitudes
%  Depol (99.5th pctile) + Hyperpol (0.5th pctile) at each M6-aligned pos
%  ~9cm wide × 5cm tall per panel
% =========================================================================

fprintf('\n=== Generating position-by-position amplitude panels ===\n');

DEP_WINDOW   = [STIM_ON, STIM_ON + 2000 - 1];   % 200 ms post-stim onset
HYP_WINDOW   = [STIM_ON + 1000, Inf];            % 100 ms post-stim onset to end
DEP_PCTILE   = 99.9;
HYP_PCTILE   = 0.1;
REJECT_THRESH = 0.5;  % mV

for panel_idx = 1:2
    if panel_idx == 1
        label = 'F'; type_str = 'ON'; type_long = 'T4 (ON)';
        mask = on_mask;
    else
        label = 'G'; type_str = 'OFF'; type_long = 'T5 (OFF)';
        mask = off_mask;
    end

    ctrl_sel = mask & ctrl_mask;
    ttl_sel  = mask & ttl_mask;
    nc = sum(ctrl_sel);
    nt = sum(ttl_sel);

    fprintf('Panel %s: %s — ctrl n=%d, tutl- n=%d\n', label, type_long, nc, nt);

    % Extract per-cell robust amplitude metrics
    [ctrl_dep, ctrl_hyp] = extract_robust_amps(results(ctrl_sel), ...
        'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
    [ttl_dep, ttl_hyp]   = extract_robust_amps(results(ttl_sel), ...
        'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

    % Per-cell FWHM (from each cell's 11-position amplitude profile)
    ctrl_fwhm = compute_percell_fwhm(-5:5, ctrl_dep);
    ttl_fwhm  = compute_percell_fwhm(-5:5, ttl_dep);

    % Figure — taller depol, shorter hyperpol (manual 3:1 axis positioning)
    fig = figure('Units', 'centimeters', 'Position', [2 2 9 5.5], ...
        'Color', 'w', 'PaperUnits', 'centimeters', ...
        'PaperSize', [9 5.5], 'PaperPosition', [0 0 9 5.5]);
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultTextFontName', FONT_NAME);

    % Manual axes: depol gets ~50% height, hyperpol ~18%, clear gap for title
    ax1 = axes(fig, 'Position', [0.14 0.36 0.82 0.48]);  % top: depol
    ax2 = axes(fig, 'Position', [0.14 0.10 0.82 0.18]);  % bottom: hyperpol

    positions = -5:5;  % 11 positions relative to M6 center

    % --- Top: Depolarization ---
    hold(ax1, 'on');
    dep_stats = draw_amplitude_panel(ax1, positions, ctrl_dep, ttl_dep, ...
        COL_CTRL, COL_TTL, STAT_COLOR, FONT_AX, FONT_STAT, ...
        COL_CTRL_F, COL_TTL_F, ALPHA_FILL);
    xline(ax1, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    ylabel(ax1, sprintf('Depol. (%.1f%% pctile, mV)', DEP_PCTILE), 'FontSize', FONT_LABEL);
    set(ax1, 'FontSize', FONT_AX, 'TickDir', 'out', 'LineWidth', 0.4);
    xlim(ax1, [-5.5 5.5]);
    ylim(ax1, YLIM_SHARED);
    box(ax1, 'off');
    set(ax1, 'XTickLabel', []);

    % FWHM annotation — bars projected onto x-axis (bottom of depol plot)
    ctrl_mn = [dep_stats.mean_ctrl];
    ttl_mn  = [dep_stats.mean_ttl];
    [fw_c, lx_c, rx_c] = compute_fwhm_positions(positions, ctrl_mn);
    [fw_t, lx_t, rx_t] = compute_fwhm_positions(positions, ttl_mn);
    fwhm_y_base = YLIM_SHARED(1) + 0.5;  % just above y-axis floor
    fwhm_y_gap  = 1.0;                    % vertical offset between ctrl and TTL bars
    if ~isnan(fw_c)
        y_c = fwhm_y_base;
        plot(ax1, [lx_c rx_c], [y_c y_c], '-', 'Color', COL_CTRL, 'LineWidth', 2.5);
        plot(ax1, [lx_c lx_c], y_c + [-0.3 0.3], '-', 'Color', COL_CTRL, 'LineWidth', 1.0);
        plot(ax1, [rx_c rx_c], y_c + [-0.3 0.3], '-', 'Color', COL_CTRL, 'LineWidth', 1.0);
    end
    if ~isnan(fw_t)
        y_t = fwhm_y_base + fwhm_y_gap;
        plot(ax1, [lx_t rx_t], [y_t y_t], '-', 'Color', COL_TTL, 'LineWidth', 2.5);
        plot(ax1, [lx_t lx_t], y_t + [-0.3 0.3], '-', 'Color', COL_TTL, 'LineWidth', 1.0);
        plot(ax1, [rx_t rx_t], y_t + [-0.3 0.3], '-', 'Color', COL_TTL, 'LineWidth', 1.0);
    end
    % FWHM text with per-cell stats (positioned near x-axis bars)
    fwhm_str = 'FWHM: ';
    fwhm_parts = {};
    if ~isnan(fw_c), fwhm_parts{end+1} = sprintf('ctrl=%.1f', fw_c); end
    if ~isnan(fw_t), fwhm_parts{end+1} = sprintf('%s=%.1f', TTL_TEX, fw_t); end
    % Rank-sum on per-cell FWHM
    cv_fw = ctrl_fwhm(~isnan(ctrl_fwhm));
    tv_fw = ttl_fwhm(~isnan(ttl_fwhm));
    if numel(cv_fw) >= 2 && numel(tv_fw) >= 2
        p_fw = ranksum(cv_fw, tv_fw);
        fwhm_parts{end+1} = sprintf('p=%.3f', p_fw);
    end
    if ~isempty(fwhm_parts)
        text(ax1, -5.3, fwhm_y_base + fwhm_y_gap + 1.5, [fwhm_str strjoin(fwhm_parts, ', ')], ...
            'FontSize', FONT_STAT, 'Interpreter', 'tex', 'VerticalAlignment', 'bottom');
    end

    % Legend
    h1 = plot(ax1, NaN, NaN, '-o', 'Color', COL_CTRL, 'MarkerFaceColor', COL_CTRL, 'MarkerSize', 3, 'LineWidth', 0.8);
    h2 = plot(ax1, NaN, NaN, '-o', 'Color', COL_TTL, 'MarkerFaceColor', COL_TTL, 'MarkerSize', 3, 'LineWidth', 0.8);
    legend(ax1, [h1 h2], {sprintf('ctrl (n=%d)', nc), sprintf('%s (n=%d)', TTL_TEX, nt)}, ...
        'FontSize', FONT_STAT, 'Location', 'northeast', 'Box', 'off', 'Interpreter', 'tex');

    % --- Bottom: Hyperpolarization (shown as negative) ---
    hold(ax2, 'on');
    hyp_stats = draw_amplitude_panel(ax2, positions, ctrl_hyp, ttl_hyp, ...
        COL_CTRL, COL_TTL, STAT_COLOR, FONT_AX, FONT_STAT, ...
        COL_CTRL_F, COL_TTL_F, ALPHA_FILL);
    xline(ax2, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    ylabel(ax2, sprintf('Hyperpol. (%.1f%% pctile, mV)', HYP_PCTILE), 'FontSize', FONT_LABEL);
    xlabel(ax2, 'Position relative to M6 center', 'FontSize', FONT_LABEL);
    set(ax2, 'FontSize', FONT_AX, 'TickDir', 'out', 'LineWidth', 0.4);
    xlim(ax2, [-5.5 5.5]);
    box(ax2, 'off');

    % Panel title (using annotation above axes)
    annotation(fig, 'textbox', [0.05 0.86 0.90 0.13], ...
        'String', sprintf('%s — Per-Position Amplitudes (M6-aligned)', type_long), ...
        'FontSize', FONT_TITLE, 'FontWeight', 'bold', 'FontName', FONT_NAME, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');

    % Export
    export_panel(fig, out_dir, sprintf('fig_ds_panel_%s_amplitudes_%s', label, type_str));
    close(fig);

    % Write stats
    fprintf(fid, '=== Panel %s: %s Per-Position Amplitudes ===\n', label, type_long);
    fprintf(fid, 'ctrl n=%d, tutl- n=%d\n\n', nc, nt);
    write_amp_stats(fid, 'Depolarization', positions, ctrl_dep, ttl_dep, dep_stats);
    write_amp_stats(fid, 'Hyperpolarization', positions, ctrl_hyp, ttl_hyp, hyp_stats);
end

%% ========================================================================
%  Panels F2/G2: Position-by-position amplitudes — ORTHOGONAL axis
%  Same Gruntman style as F/G but using ortho_flash_m6_aligned traces
% =========================================================================

fprintf('\n=== Generating position-by-position amplitude panels (Ortho) ===\n');

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

    fprintf('Panel %s: %s — ctrl n=%d, tutl- n=%d\n', label, type_long, nc, nt);

    % Extract per-cell robust amplitude metrics (ortho axis)
    [ctrl_dep, ctrl_hyp] = extract_robust_amps(results(ctrl_sel), ...
        'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
    [ttl_dep, ttl_hyp]   = extract_robust_amps(results(ttl_sel), ...
        'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

    % Per-cell FWHM (from each cell's 11-position amplitude profile)
    ctrl_fwhm = compute_percell_fwhm(-5:5, ctrl_dep);
    ttl_fwhm  = compute_percell_fwhm(-5:5, ttl_dep);

    % Figure — taller depol, shorter hyperpol (manual 3:1 axis positioning)
    fig = figure('Units', 'centimeters', 'Position', [2 2 9 5.5], ...
        'Color', 'w', 'PaperUnits', 'centimeters', ...
        'PaperSize', [9 5.5], 'PaperPosition', [0 0 9 5.5]);
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultTextFontName', FONT_NAME);

    % Manual axes: depol gets ~48% height, hyperpol ~18%, clear gap for title
    ax1 = axes(fig, 'Position', [0.14 0.36 0.82 0.48]);  % top: depol
    ax2 = axes(fig, 'Position', [0.14 0.10 0.82 0.18]);  % bottom: hyperpol

    positions = -5:5;  % 11 positions relative to M6 center

    % --- Top: Depolarization ---
    hold(ax1, 'on');
    dep_stats = draw_amplitude_panel(ax1, positions, ctrl_dep, ttl_dep, ...
        COL_CTRL, COL_TTL, STAT_COLOR, FONT_AX, FONT_STAT, ...
        COL_CTRL_F, COL_TTL_F, ALPHA_FILL);
    xline(ax1, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    ylabel(ax1, sprintf('Depol. (%.1f%% pctile, mV)', DEP_PCTILE), 'FontSize', FONT_LABEL);
    set(ax1, 'FontSize', FONT_AX, 'TickDir', 'out', 'LineWidth', 0.4);
    xlim(ax1, [-5.5 5.5]);
    ylim(ax1, YLIM_SHARED);
    box(ax1, 'off');
    set(ax1, 'XTickLabel', []);

    % FWHM annotation — bars projected onto x-axis (ortho)
    ctrl_mn = [dep_stats.mean_ctrl];
    ttl_mn  = [dep_stats.mean_ttl];
    [fw_c, lx_c, rx_c] = compute_fwhm_positions(positions, ctrl_mn);
    [fw_t, lx_t, rx_t] = compute_fwhm_positions(positions, ttl_mn);
    fwhm_y_base = YLIM_SHARED(1) + 0.5;  % just above y-axis floor
    fwhm_y_gap  = 1.0;                    % vertical offset between ctrl and TTL bars
    if ~isnan(fw_c)
        y_c = fwhm_y_base;
        plot(ax1, [lx_c rx_c], [y_c y_c], '-', 'Color', COL_CTRL, 'LineWidth', 2.5);
        plot(ax1, [lx_c lx_c], y_c + [-0.3 0.3], '-', 'Color', COL_CTRL, 'LineWidth', 1.0);
        plot(ax1, [rx_c rx_c], y_c + [-0.3 0.3], '-', 'Color', COL_CTRL, 'LineWidth', 1.0);
    end
    if ~isnan(fw_t)
        y_t = fwhm_y_base + fwhm_y_gap;
        plot(ax1, [lx_t rx_t], [y_t y_t], '-', 'Color', COL_TTL, 'LineWidth', 2.5);
        plot(ax1, [lx_t lx_t], y_t + [-0.3 0.3], '-', 'Color', COL_TTL, 'LineWidth', 1.0);
        plot(ax1, [rx_t rx_t], y_t + [-0.3 0.3], '-', 'Color', COL_TTL, 'LineWidth', 1.0);
    end
    % FWHM text with per-cell stats (positioned near x-axis bars)
    fwhm_str = 'FWHM: ';
    fwhm_parts = {};
    if ~isnan(fw_c), fwhm_parts{end+1} = sprintf('ctrl=%.1f', fw_c); end
    if ~isnan(fw_t), fwhm_parts{end+1} = sprintf('%s=%.1f', TTL_TEX, fw_t); end
    % Rank-sum on per-cell FWHM
    cv_fw = ctrl_fwhm(~isnan(ctrl_fwhm));
    tv_fw = ttl_fwhm(~isnan(ttl_fwhm));
    if numel(cv_fw) >= 2 && numel(tv_fw) >= 2
        p_fw = ranksum(cv_fw, tv_fw);
        fwhm_parts{end+1} = sprintf('p=%.3f', p_fw);
    end
    if ~isempty(fwhm_parts)
        text(ax1, -5.3, fwhm_y_base + fwhm_y_gap + 1.5, [fwhm_str strjoin(fwhm_parts, ', ')], ...
            'FontSize', FONT_STAT, 'Interpreter', 'tex', 'VerticalAlignment', 'bottom');
    end

    % Legend
    h1 = plot(ax1, NaN, NaN, '-o', 'Color', COL_CTRL, 'MarkerFaceColor', COL_CTRL, 'MarkerSize', 3, 'LineWidth', 0.8);
    h2 = plot(ax1, NaN, NaN, '-o', 'Color', COL_TTL, 'MarkerFaceColor', COL_TTL, 'MarkerSize', 3, 'LineWidth', 0.8);
    legend(ax1, [h1 h2], {sprintf('ctrl (n=%d)', nc), sprintf('%s (n=%d)', TTL_TEX, nt)}, ...
        'FontSize', FONT_STAT, 'Location', 'northeast', 'Box', 'off', 'Interpreter', 'tex');

    % --- Bottom: Hyperpolarization ---
    hold(ax2, 'on');
    hyp_stats = draw_amplitude_panel(ax2, positions, ctrl_hyp, ttl_hyp, ...
        COL_CTRL, COL_TTL, STAT_COLOR, FONT_AX, FONT_STAT, ...
        COL_CTRL_F, COL_TTL_F, ALPHA_FILL);
    xline(ax2, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    ylabel(ax2, sprintf('Hyperpol. (%.1f%% pctile, mV)', HYP_PCTILE), 'FontSize', FONT_LABEL);
    xlabel(ax2, 'Position relative to M6 center (ortho)', 'FontSize', FONT_LABEL);
    set(ax2, 'FontSize', FONT_AX, 'TickDir', 'out', 'LineWidth', 0.4);
    xlim(ax2, [-5.5 5.5]);
    box(ax2, 'off');

    % Panel title (using annotation above axes)
    annotation(fig, 'textbox', [0.05 0.86 0.90 0.13], ...
        'String', sprintf('%s — Per-Position Amplitudes, Ortho (M6-aligned)', type_long), ...
        'FontSize', FONT_TITLE, 'FontWeight', 'bold', 'FontName', FONT_NAME, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');

    % Export
    export_panel(fig, out_dir, sprintf('fig_ds_panel_%s_amplitudes_%s_ortho', label, type_str));
    close(fig);

    % Write stats
    fprintf(fid, '=== Panel %s: %s Per-Position Amplitudes (Ortho) ===\n', label, type_long);
    fprintf(fid, 'ctrl n=%d, tutl- n=%d\n\n', nc, nt);
    write_amp_stats(fid, 'Depolarization (ortho)', positions, ctrl_dep, ttl_dep, dep_stats);
    write_amp_stats(fid, 'Hyperpolarization (ortho)', positions, ctrl_hyp, ttl_hyp, hyp_stats);
end

%% ========================================================================
%  Panels F/G POOLED-STATS: Same 11-position amplitude plots as F/G,
%  but asterisks computed from 3-position pooled rank-sum instead of
%  per-position rank-sum. 9 pools centered at positions -4 to +4.
% =========================================================================

fprintf('\n=== Generating amplitude panels with pooled stats (PD) ===\n');

for panel_idx = 1:2
    if panel_idx == 1
        label = 'F'; type_str = 'ON'; type_long = 'T4 (ON)';
        mask = on_mask;
    else
        label = 'G'; type_str = 'OFF'; type_long = 'T5 (OFF)';
        mask = off_mask;
    end

    ctrl_sel = mask & ctrl_mask;
    ttl_sel  = mask & ttl_mask;
    nc = sum(ctrl_sel);
    nt = sum(ttl_sel);

    fprintf('Pooled-stats Panel %s: %s — ctrl n=%d, tutl- n=%d\n', label, type_long, nc, nt);

    % Extract per-cell robust amplitude metrics (same 11-position data as F/G)
    [ctrl_dep, ctrl_hyp] = extract_robust_amps(results(ctrl_sel), ...
        'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
    [ttl_dep, ttl_hyp]   = extract_robust_amps(results(ttl_sel), ...
        'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

    % Compute pooled p-values: 3-position sliding window rank-sum
    [dep_pool_p, pool_centers] = compute_pooled_ranksum(ctrl_dep, ttl_dep);
    [hyp_pool_p, ~]            = compute_pooled_ranksum(ctrl_hyp, ttl_hyp);

    % Per-cell FWHM (from each cell's 11-position amplitude profile)
    ctrl_fwhm = compute_percell_fwhm(-5:5, ctrl_dep);
    ttl_fwhm  = compute_percell_fwhm(-5:5, ttl_dep);

    % Figure — identical layout to F/G
    fig = figure('Units', 'centimeters', 'Position', [2 2 9 5.5], ...
        'Color', 'w', 'PaperUnits', 'centimeters', ...
        'PaperSize', [9 5.5], 'PaperPosition', [0 0 9 5.5]);
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultTextFontName', FONT_NAME);

    ax1 = axes(fig, 'Position', [0.14 0.36 0.82 0.48]);  % top: depol
    ax2 = axes(fig, 'Position', [0.14 0.10 0.82 0.18]);  % bottom: hyperpol

    positions = -5:5;  % 11 positions relative to M6 center

    % --- Top: Depolarization (original 11-pos data, pooled-stats asterisks) ---
    hold(ax1, 'on');
    dep_stats = draw_amplitude_panel_no_stars(ax1, positions, ctrl_dep, ttl_dep, ...
        COL_CTRL, COL_TTL, FONT_AX, COL_CTRL_F, COL_TTL_F, ALPHA_FILL);
    draw_pooled_asterisks(ax1, pool_centers, dep_pool_p, dep_stats, STAT_COLOR, FONT_STAT);
    xline(ax1, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    ylabel(ax1, sprintf('Depol. (%.1f%% pctile, mV)', DEP_PCTILE), 'FontSize', FONT_LABEL);
    set(ax1, 'FontSize', FONT_AX, 'TickDir', 'out', 'LineWidth', 0.4);
    xlim(ax1, [-5.5 5.5]);
    ylim(ax1, YLIM_SHARED);
    box(ax1, 'off');
    set(ax1, 'XTickLabel', []);

    % FWHM annotation — bars projected onto x-axis (same as F/G)
    ctrl_mn = [dep_stats.mean_ctrl];
    ttl_mn  = [dep_stats.mean_ttl];
    [fw_c, lx_c, rx_c] = compute_fwhm_positions(positions, ctrl_mn);
    [fw_t, lx_t, rx_t] = compute_fwhm_positions(positions, ttl_mn);
    fwhm_y_base = YLIM_SHARED(1) + 0.5;
    fwhm_y_gap  = 1.0;
    if ~isnan(fw_c)
        y_c = fwhm_y_base;
        plot(ax1, [lx_c rx_c], [y_c y_c], '-', 'Color', COL_CTRL, 'LineWidth', 2.5);
        plot(ax1, [lx_c lx_c], y_c + [-0.3 0.3], '-', 'Color', COL_CTRL, 'LineWidth', 1.0);
        plot(ax1, [rx_c rx_c], y_c + [-0.3 0.3], '-', 'Color', COL_CTRL, 'LineWidth', 1.0);
    end
    if ~isnan(fw_t)
        y_t = fwhm_y_base + fwhm_y_gap;
        plot(ax1, [lx_t rx_t], [y_t y_t], '-', 'Color', COL_TTL, 'LineWidth', 2.5);
        plot(ax1, [lx_t lx_t], y_t + [-0.3 0.3], '-', 'Color', COL_TTL, 'LineWidth', 1.0);
        plot(ax1, [rx_t rx_t], y_t + [-0.3 0.3], '-', 'Color', COL_TTL, 'LineWidth', 1.0);
    end
    fwhm_str = 'FWHM: ';
    fwhm_parts = {};
    if ~isnan(fw_c), fwhm_parts{end+1} = sprintf('ctrl=%.1f', fw_c); end
    if ~isnan(fw_t), fwhm_parts{end+1} = sprintf('%s=%.1f', TTL_TEX, fw_t); end
    cv_fw = ctrl_fwhm(~isnan(ctrl_fwhm));
    tv_fw = ttl_fwhm(~isnan(ttl_fwhm));
    if numel(cv_fw) >= 2 && numel(tv_fw) >= 2
        p_fw = ranksum(cv_fw, tv_fw);
        fwhm_parts{end+1} = sprintf('p=%.3f', p_fw);
    end
    if ~isempty(fwhm_parts)
        text(ax1, -5.3, fwhm_y_base + fwhm_y_gap + 1.5, [fwhm_str strjoin(fwhm_parts, ', ')], ...
            'FontSize', FONT_STAT, 'Interpreter', 'tex', 'VerticalAlignment', 'bottom');
    end

    % Legend
    h1 = plot(ax1, NaN, NaN, '-o', 'Color', COL_CTRL, 'MarkerFaceColor', COL_CTRL, 'MarkerSize', 3, 'LineWidth', 0.8);
    h2 = plot(ax1, NaN, NaN, '-o', 'Color', COL_TTL, 'MarkerFaceColor', COL_TTL, 'MarkerSize', 3, 'LineWidth', 0.8);
    legend(ax1, [h1 h2], {sprintf('ctrl (n=%d)', nc), sprintf('%s (n=%d)', TTL_TEX, nt)}, ...
        'FontSize', FONT_STAT, 'Location', 'northeast', 'Box', 'off', 'Interpreter', 'tex');

    % --- Bottom: Hyperpolarization (original 11-pos data, pooled-stats asterisks) ---
    hold(ax2, 'on');
    hyp_stats = draw_amplitude_panel_no_stars(ax2, positions, ctrl_hyp, ttl_hyp, ...
        COL_CTRL, COL_TTL, FONT_AX, COL_CTRL_F, COL_TTL_F, ALPHA_FILL);
    draw_pooled_asterisks(ax2, pool_centers, hyp_pool_p, hyp_stats, STAT_COLOR, FONT_STAT);
    xline(ax2, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    ylabel(ax2, sprintf('Hyperpol. (%.1f%% pctile, mV)', HYP_PCTILE), 'FontSize', FONT_LABEL);
    xlabel(ax2, 'Position relative to M6 center', 'FontSize', FONT_LABEL);
    set(ax2, 'FontSize', FONT_AX, 'TickDir', 'out', 'LineWidth', 0.4);
    xlim(ax2, [-5.5 5.5]);
    box(ax2, 'off');

    % Panel title
    annotation(fig, 'textbox', [0.05 0.86 0.90 0.13], ...
        'String', sprintf('%s — Amplitudes, pooled stats (M6-aligned)', type_long), ...
        'FontSize', FONT_TITLE, 'FontWeight', 'bold', 'FontName', FONT_NAME, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');

    % Export
    export_panel(fig, out_dir, sprintf('fig_ds_panel_%s_amplitudes_%s_pooled', label, type_str));
    close(fig);

    % Write stats
    fprintf(fid, '=== Panel %s POOLED-STATS: %s Amplitudes with 3-pos pooled rank-sum ===\n', label, type_long);
    fprintf(fid, 'ctrl n=%d, tutl- n=%d\n', nc, nt);
    fprintf(fid, 'Pooled stat: for each pool center, rank-sum on concatenated values from 3 consecutive positions\n\n');
    write_pooled_stats(fid, 'Depolarization (pooled stats)', pool_centers, dep_pool_p);
    write_pooled_stats(fid, 'Hyperpolarization (pooled stats)', pool_centers, hyp_pool_p);
    fprintf(fid, '\n');
end

%% ========================================================================
%  Panels F2/G2 POOLED-STATS: Same 11-pos ortho amplitude plots,
%  but asterisks from 3-position pooled rank-sum.
% =========================================================================

fprintf('\n=== Generating amplitude panels with pooled stats (Ortho) ===\n');

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

    fprintf('Pooled-stats Panel %s: %s — ctrl n=%d, tutl- n=%d\n', label, type_long, nc, nt);

    % Extract per-cell robust amplitude metrics (ortho)
    [ctrl_dep, ctrl_hyp] = extract_robust_amps(results(ctrl_sel), ...
        'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
    [ttl_dep, ttl_hyp]   = extract_robust_amps(results(ttl_sel), ...
        'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

    % Compute pooled p-values: 3-position sliding window rank-sum
    [dep_pool_p, pool_centers] = compute_pooled_ranksum(ctrl_dep, ttl_dep);
    [hyp_pool_p, ~]            = compute_pooled_ranksum(ctrl_hyp, ttl_hyp);

    % Per-cell FWHM
    ctrl_fwhm = compute_percell_fwhm(-5:5, ctrl_dep);
    ttl_fwhm  = compute_percell_fwhm(-5:5, ttl_dep);

    % Figure — identical layout to F2/G2
    fig = figure('Units', 'centimeters', 'Position', [2 2 9 5.5], ...
        'Color', 'w', 'PaperUnits', 'centimeters', ...
        'PaperSize', [9 5.5], 'PaperPosition', [0 0 9 5.5]);
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultTextFontName', FONT_NAME);

    ax1 = axes(fig, 'Position', [0.14 0.36 0.82 0.48]);
    ax2 = axes(fig, 'Position', [0.14 0.10 0.82 0.18]);

    positions = -5:5;

    % --- Top: Depolarization (original 11-pos data, pooled-stats asterisks) ---
    hold(ax1, 'on');
    dep_stats = draw_amplitude_panel_no_stars(ax1, positions, ctrl_dep, ttl_dep, ...
        COL_CTRL, COL_TTL, FONT_AX, COL_CTRL_F, COL_TTL_F, ALPHA_FILL);
    draw_pooled_asterisks(ax1, pool_centers, dep_pool_p, dep_stats, STAT_COLOR, FONT_STAT);
    xline(ax1, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    ylabel(ax1, sprintf('Depol. (%.1f%% pctile, mV)', DEP_PCTILE), 'FontSize', FONT_LABEL);
    set(ax1, 'FontSize', FONT_AX, 'TickDir', 'out', 'LineWidth', 0.4);
    xlim(ax1, [-5.5 5.5]);
    ylim(ax1, YLIM_SHARED);
    box(ax1, 'off');
    set(ax1, 'XTickLabel', []);

    % FWHM annotation — bars projected onto x-axis (ortho)
    ctrl_mn = [dep_stats.mean_ctrl];
    ttl_mn  = [dep_stats.mean_ttl];
    [fw_c, lx_c, rx_c] = compute_fwhm_positions(positions, ctrl_mn);
    [fw_t, lx_t, rx_t] = compute_fwhm_positions(positions, ttl_mn);
    fwhm_y_base = YLIM_SHARED(1) + 0.5;
    fwhm_y_gap  = 1.0;
    if ~isnan(fw_c)
        y_c = fwhm_y_base;
        plot(ax1, [lx_c rx_c], [y_c y_c], '-', 'Color', COL_CTRL, 'LineWidth', 2.5);
        plot(ax1, [lx_c lx_c], y_c + [-0.3 0.3], '-', 'Color', COL_CTRL, 'LineWidth', 1.0);
        plot(ax1, [rx_c rx_c], y_c + [-0.3 0.3], '-', 'Color', COL_CTRL, 'LineWidth', 1.0);
    end
    if ~isnan(fw_t)
        y_t = fwhm_y_base + fwhm_y_gap;
        plot(ax1, [lx_t rx_t], [y_t y_t], '-', 'Color', COL_TTL, 'LineWidth', 2.5);
        plot(ax1, [lx_t lx_t], y_t + [-0.3 0.3], '-', 'Color', COL_TTL, 'LineWidth', 1.0);
        plot(ax1, [rx_t rx_t], y_t + [-0.3 0.3], '-', 'Color', COL_TTL, 'LineWidth', 1.0);
    end
    fwhm_str = 'FWHM: ';
    fwhm_parts = {};
    if ~isnan(fw_c), fwhm_parts{end+1} = sprintf('ctrl=%.1f', fw_c); end
    if ~isnan(fw_t), fwhm_parts{end+1} = sprintf('%s=%.1f', TTL_TEX, fw_t); end
    cv_fw = ctrl_fwhm(~isnan(ctrl_fwhm));
    tv_fw = ttl_fwhm(~isnan(ttl_fwhm));
    if numel(cv_fw) >= 2 && numel(tv_fw) >= 2
        p_fw = ranksum(cv_fw, tv_fw);
        fwhm_parts{end+1} = sprintf('p=%.3f', p_fw);
    end
    if ~isempty(fwhm_parts)
        text(ax1, -5.3, fwhm_y_base + fwhm_y_gap + 1.5, [fwhm_str strjoin(fwhm_parts, ', ')], ...
            'FontSize', FONT_STAT, 'Interpreter', 'tex', 'VerticalAlignment', 'bottom');
    end

    % Legend
    h1 = plot(ax1, NaN, NaN, '-o', 'Color', COL_CTRL, 'MarkerFaceColor', COL_CTRL, 'MarkerSize', 3, 'LineWidth', 0.8);
    h2 = plot(ax1, NaN, NaN, '-o', 'Color', COL_TTL, 'MarkerFaceColor', COL_TTL, 'MarkerSize', 3, 'LineWidth', 0.8);
    legend(ax1, [h1 h2], {sprintf('ctrl (n=%d)', nc), sprintf('%s (n=%d)', TTL_TEX, nt)}, ...
        'FontSize', FONT_STAT, 'Location', 'northeast', 'Box', 'off', 'Interpreter', 'tex');

    % --- Bottom: Hyperpolarization (original 11-pos data, pooled-stats asterisks) ---
    hold(ax2, 'on');
    hyp_stats = draw_amplitude_panel_no_stars(ax2, positions, ctrl_hyp, ttl_hyp, ...
        COL_CTRL, COL_TTL, FONT_AX, COL_CTRL_F, COL_TTL_F, ALPHA_FILL);
    draw_pooled_asterisks(ax2, pool_centers, hyp_pool_p, hyp_stats, STAT_COLOR, FONT_STAT);
    xline(ax2, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    ylabel(ax2, sprintf('Hyperpol. (%.1f%% pctile, mV)', HYP_PCTILE), 'FontSize', FONT_LABEL);
    xlabel(ax2, 'Position relative to M6 center (ortho)', 'FontSize', FONT_LABEL);
    set(ax2, 'FontSize', FONT_AX, 'TickDir', 'out', 'LineWidth', 0.4);
    xlim(ax2, [-5.5 5.5]);
    box(ax2, 'off');

    % Panel title
    annotation(fig, 'textbox', [0.05 0.86 0.90 0.13], ...
        'String', sprintf('%s — Amplitudes, pooled stats, Ortho (M6-aligned)', type_long), ...
        'FontSize', FONT_TITLE, 'FontWeight', 'bold', 'FontName', FONT_NAME, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');

    % Export
    export_panel(fig, out_dir, sprintf('fig_ds_panel_%s_amplitudes_%s_ortho_pooled', label, type_str));
    close(fig);

    % Write stats
    fprintf(fid, '=== Panel %s POOLED-STATS: %s Amplitudes, Ortho, with 3-pos pooled rank-sum ===\n', label, type_long);
    fprintf(fid, 'ctrl n=%d, tutl- n=%d\n', nc, nt);
    fprintf(fid, 'Pooled stat: for each pool center, rank-sum on concatenated values from 3 consecutive positions\n\n');
    write_pooled_stats(fid, 'Depolarization, ortho (pooled stats)', pool_centers, dep_pool_p);
    write_pooled_stats(fid, 'Hyperpolarization, ortho (pooled stats)', pool_centers, hyp_pool_p);
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('\nSaved: %s\n', stats_file);
fprintf('=== Done ===\n');


%% ========================= Local Functions ===============================

function y_lim = compute_shared_ylim_pd(on_ctrl, on_ttl, off_ctrl, off_ttl, x_start)
% Compute shared y-limits across ON + OFF PD trace groups (no ortho).
% Only considers samples from x_start onward (visible window).
    all_vals = [];
    for traces = {on_ctrl, on_ttl, off_ctrl, off_ttl}
        for c = 1:numel(traces{1})
            mat = traces{1}{c};
            if ~isempty(mat)
                vis = mat(:, x_start:end);  % visible window only
                mn = mean(vis(:), 'omitnan');
                sd = std(vis(:), 'omitnan');
                all_vals = [all_vals; mn - 2*sd; mn + 2*sd]; %#ok<AGROW>
            end
        end
    end
    if isempty(all_vals)
        y_lim = [-15 35];
    else
        y_lim = [min(all_vals), max(all_vals)];
        pad = 0.10 * diff(y_lim);
        y_lim = [y_lim(1) - pad, y_lim(2) + pad];
    end
end


function draw_1drf_tile(ax, traces_ctrl, traces_ttl, pos, y_lim, ...
    col_c, col_cf, col_t, col_tf, alpha_f, ...
    stim_on, stim_off, resp_end, col_stim, col_resp, x_start)
% Draw a single tile of the 1D RF panel.
% x_start: first sample to display (crops pre-stimulus baseline).

    % Compute stats
    [ctr_c, spr_c, n_c] = compute_group_stats(traces_ctrl, pos);
    [ctr_t, spr_t, n_t] = compute_group_stats(traces_ttl, pos);

    % Stim timing lines (Gruntman style: subtle)
    if ~isempty(ctr_c) || ~isempty(ctr_t)
        xline(ax, stim_on, '-', 'Color', col_stim, 'LineWidth', 0.3, 'Alpha', 0.6);
        xline(ax, stim_off, '-', 'Color', col_stim, 'LineWidth', 0.3, 'Alpha', 0.6);
        xline(ax, resp_end, '-', 'Color', col_resp, 'LineWidth', 0.2, 'Alpha', 0.4);
    end

    % Control traces (Gruntman style: thin mean line)
    if ~isempty(ctr_c) && n_c >= 2
        x = 1:numel(ctr_c);
        fill_x = [x, fliplr(x)];
        fill_y = [ctr_c + spr_c, fliplr(ctr_c - spr_c)];
        fill(ax, fill_x, fill_y, col_cf, 'FaceAlpha', alpha_f, 'EdgeColor', 'none');
        plot(ax, x, ctr_c, '-', 'Color', col_c, 'LineWidth', 0.6);
    end

    % TTL traces (Gruntman style: thin mean line)
    if ~isempty(ctr_t) && n_t >= 2
        x = 1:numel(ctr_t);
        fill_x = [x, fliplr(x)];
        fill_y = [ctr_t + spr_t, fliplr(ctr_t - spr_t)];
        fill(ax, fill_x, fill_y, col_tf, 'FaceAlpha', alpha_f, 'EdgeColor', 'none');
        plot(ax, x, ctr_t, '-', 'Color', col_t, 'LineWidth', 0.6);
    end

    ylim(ax, y_lim);

    % Crop x-axis to visible window
    trace_len = max([numel(ctr_c), numel(ctr_t), 1]);
    xlim(ax, [x_start, trace_len]);
end


function [ctr, spr, n] = compute_group_stats(traces, pos)
% Mean ± SEM for a given position across all cells.
% Handles variable trace lengths by truncating to minimum common length.
    ctr = []; spr = []; n = 0;
    n_cells = numel(traces);
    if n_cells == 0, return; end

    % First pass: collect valid rows and find min length
    valid_rows = {};
    for c = 1:n_cells
        mat = traces{c};
        if isempty(mat), continue; end
        if pos > size(mat, 1), continue; end
        row = mat(pos, :);
        if all(isnan(row)), continue; end
        valid_rows{end+1} = row; %#ok<AGROW>
    end
    n = numel(valid_rows);
    if n < 1, return; end

    min_len = min(cellfun(@numel, valid_rows));
    rows = NaN(n, min_len);
    for r = 1:n
        rows(r, :) = valid_rows{r}(1:min_len);
    end
    ctr = mean(rows, 1, 'omitnan');
    spr = std(rows, 0, 1, 'omitnan') ./ sqrt(n);
end


function stats = compute_position_stats(traces_ctrl, traces_ttl, stim_on, resp_end)
% Per-position Wilcoxon rank-sum on 99.5th pctile peak amplitude.
    n_pos = 11;
    stats = struct('pos', num2cell(1:n_pos), 'n_ctrl', 0, 'n_ttl', 0, ...
        'mean_ctrl', NaN, 'mean_ttl', NaN, 'p', NaN, 'sig', false);

    for p = 1:n_pos
        ctrl_vals = extract_peak_at_pos(traces_ctrl, p, stim_on, resp_end);
        ttl_vals  = extract_peak_at_pos(traces_ttl, p, stim_on, resp_end);
        stats(p).n_ctrl = numel(ctrl_vals);
        stats(p).n_ttl  = numel(ttl_vals);
        stats(p).mean_ctrl = mean(ctrl_vals);
        stats(p).mean_ttl  = mean(ttl_vals);

        if numel(ctrl_vals) >= 2 && numel(ttl_vals) >= 2
            stats(p).p = ranksum(ctrl_vals, ttl_vals);
            stats(p).sig = stats(p).p < 0.05;
        end
    end
end


function vals = extract_peak_at_pos(traces, pos, stim_on, resp_end)
% Extract 99.5th percentile peak for each cell at a given position.
    vals = [];
    for c = 1:numel(traces)
        mat = traces{c};
        if isempty(mat) || pos > size(mat, 1), continue; end
        row = mat(pos, :);
        if all(isnan(row)), continue; end
        win_end = min(resp_end, numel(row));
        if stim_on > numel(row), continue; end
        win = row(stim_on:win_end);
        win = win(~isnan(win));
        if isempty(win), continue; end
        vals(end+1) = prctile(win, 99.5); %#ok<AGROW>
    end
end


function draw_position_asterisk(ax, p_val, y_lim, col, font_size)
% Draw asterisk at top of tile for significant position.
    if p_val < 0.001
        str = '***';
    elseif p_val < 0.01
        str = '**';
    else
        str = '*';
    end
    xl = xlim(ax);
    x_mid = mean(xl);
    y_top = y_lim(2) - 0.02 * diff(y_lim);
    text(ax, x_mid, y_top, str, 'FontSize', font_size, 'FontWeight', 'bold', ...
        'Color', col, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end


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

            % Depolarization
            de = min(dep_win(2), n_samp);
            if ~isinf(de)
                dep_mat(k, pos) = prctile(row(dep_win(1):de), dep_pct);
            else
                dep_mat(k, pos) = prctile(row(dep_win(1):end), dep_pct);
            end

            % Hyperpolarization — clamp positive values to zero
            % (no hyperpolarization detected → report zero, not a positive artifact)
            he = hyp_win(2);
            if isinf(he), he = n_samp; end
            he = min(he, n_samp);
            hval = prctile(row(hyp_win(1):he), hyp_pct);
            hyp_mat(k, pos) = min(hval, 0);
        end
    end
end


function stats = draw_amplitude_panel(ax, positions, ctrl_data, ttl_data, ...
    col_c, col_t, col_stat, font_ax, font_stat, col_cf, col_tf, alpha_fill)
% Draw Gruntman-style connected-line amplitude panel.
%   Connected points (mean) with SEM shading, Wilcoxon asterisks.
%   col_cf, col_tf: fill colors for SEM bands.
%   alpha_fill: fill transparency.
    n_pos = numel(positions);
    stats = struct('pos', num2cell(positions), 'n_ctrl', 0, 'n_ttl', 0, ...
        'mean_ctrl', NaN, 'mean_ttl', NaN, 'sem_ctrl', NaN, 'sem_ttl', NaN, ...
        'p', NaN, 'sig', false);

    % Pre-compute means and SEMs
    ctrl_mean = NaN(1, n_pos);
    ctrl_sem  = NaN(1, n_pos);
    ttl_mean  = NaN(1, n_pos);
    ttl_sem   = NaN(1, n_pos);

    for i = 1:n_pos
        cv = ctrl_data(:, i); cv = cv(~isnan(cv));
        tv = ttl_data(:, i);  tv = tv(~isnan(tv));
        stats(i).n_ctrl = numel(cv);
        stats(i).n_ttl  = numel(tv);

        if numel(cv) >= 2
            ctrl_mean(i) = mean(cv);
            ctrl_sem(i)  = std(cv) / sqrt(numel(cv));
        elseif numel(cv) == 1
            ctrl_mean(i) = cv;
            ctrl_sem(i)  = 0;
        end
        stats(i).mean_ctrl = ctrl_mean(i);
        stats(i).sem_ctrl  = ctrl_sem(i);

        if numel(tv) >= 2
            ttl_mean(i) = mean(tv);
            ttl_sem(i)  = std(tv) / sqrt(numel(tv));
        elseif numel(tv) == 1
            ttl_mean(i) = tv;
            ttl_sem(i)  = 0;
        end
        stats(i).mean_ttl = ttl_mean(i);
        stats(i).sem_ttl  = ttl_sem(i);

        % Wilcoxon rank-sum
        if numel(cv) >= 2 && numel(tv) >= 2
            p = ranksum(cv, tv);
            stats(i).p = p;
            stats(i).sig = p < 0.05;
        end
    end

    % Draw SEM bands (ctrl first, then TTL on top)
    valid_c = ~isnan(ctrl_mean);
    if any(valid_c)
        xc = positions(valid_c);
        fill_x = [xc, fliplr(xc)];
        fill_y = [ctrl_mean(valid_c) + ctrl_sem(valid_c), ...
                  fliplr(ctrl_mean(valid_c) - ctrl_sem(valid_c))];
        fill(ax, fill_x, fill_y, col_cf, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none');
    end

    valid_t = ~isnan(ttl_mean);
    if any(valid_t)
        xt = positions(valid_t);
        fill_x = [xt, fliplr(xt)];
        fill_y = [ttl_mean(valid_t) + ttl_sem(valid_t), ...
                  fliplr(ttl_mean(valid_t) - ttl_sem(valid_t))];
        fill(ax, fill_x, fill_y, col_tf, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none');
    end

    % Draw connected mean lines with markers
    if any(valid_c)
        plot(ax, positions(valid_c), ctrl_mean(valid_c), '-o', ...
            'Color', col_c, 'MarkerFaceColor', col_c, 'MarkerSize', 3, ...
            'LineWidth', 0.8);
    end
    if any(valid_t)
        plot(ax, positions(valid_t), ttl_mean(valid_t), '-o', ...
            'Color', col_t, 'MarkerFaceColor', col_t, 'MarkerSize', 3, ...
            'LineWidth', 0.8);
    end

    % Asterisks
    for i = 1:n_pos
        if stats(i).sig
            x = positions(i);
            y_top = max([ctrl_mean(i) + ctrl_sem(i), ttl_mean(i) + ttl_sem(i)]);
            if isnan(y_top), continue; end
            yl = ylim(ax);
            y_star = y_top + 0.06 * diff(yl);
            if stats(i).p < 0.001,    str = '***';
            elseif stats(i).p < 0.01, str = '**';
            else,                       str = '*';
            end
            text(ax, x, y_star, str, 'FontSize', font_stat, ...
                'HorizontalAlignment', 'center', 'Color', col_stat, ...
                'FontWeight', 'bold');
        end
    end
end


function export_panel(fig, out_dir, name)
% Export panel as timestamped PNG only.
    ts = datestr(now, 'yyyymmdd_HHMM');
    png_ts = fullfile(out_dir, sprintf('%s_%s.png', name, ts));
    exportgraphics(fig, png_ts, 'Resolution', 300);
    fprintf('  Saved: %s\n', png_ts);
end


function write_pos_stats(fid, axis_label, stats, labels)
% Write per-position stats table to file.
    fprintf(fid, '%s:\n', axis_label);
    fprintf(fid, '%-8s  %5s  %5s  %8s  %8s  %8s  %s\n', ...
        'Pos', 'nCtrl', 'nTTL', 'MnCtrl', 'MnTTL', 'p', 'Sig');
    fprintf(fid, '%s\n', repmat('-', 1, 60));
    for i = 1:numel(stats)
        sig_str = '';
        if stats(i).sig, sig_str = '*'; end
        if ~isnan(stats(i).p) && stats(i).p < 0.01, sig_str = '**'; end
        if ~isnan(stats(i).p) && stats(i).p < 0.001, sig_str = '***'; end
        fprintf(fid, '%-8s  %5d  %5d  %8.2f  %8.2f  %8.4f  %s\n', ...
            labels{i}, stats(i).n_ctrl, stats(i).n_ttl, ...
            stats(i).mean_ctrl, stats(i).mean_ttl, stats(i).p, sig_str);
    end
    fprintf(fid, '\n');
end


function write_amp_stats(fid, metric_label, positions, ctrl_data, ttl_data, stats)
% Write amplitude stats table to file.
    fprintf(fid, '%s:\n', metric_label);
    fprintf(fid, '%-5s  %5s  %5s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %s\n', ...
        'Pos', 'nCtrl', 'nTTL', 'MnCtrl', 'MnTTL', 'SEMctrl', 'SEMttl', ...
        'MdCtrl', 'MdTTL', 'p', 'Sig');
    fprintf(fid, '%s\n', repmat('-', 1, 100));
    for i = 1:numel(positions)
        cv = ctrl_data(:, i); cv = cv(~isnan(cv));
        tv = ttl_data(:, i);  tv = tv(~isnan(tv));

        mn_c = mean(cv); mn_t = mean(tv);
        se_c = std(cv)/sqrt(max(1,numel(cv)));
        se_t = std(tv)/sqrt(max(1,numel(tv)));
        md_c = median(cv); md_t = median(tv);

        sig_str = '';
        if stats(i).sig, sig_str = '*'; end
        if ~isnan(stats(i).p) && stats(i).p < 0.01, sig_str = '**'; end
        if ~isnan(stats(i).p) && stats(i).p < 0.001, sig_str = '***'; end

        fprintf(fid, '%+3d    %5d  %5d  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f  %8.4f  %s\n', ...
            positions(i), stats(i).n_ctrl, stats(i).n_ttl, ...
            mn_c, mn_t, se_c, se_t, md_c, md_t, stats(i).p, sig_str);
    end
    fprintf(fid, '\n');
end


function [fw, left_x, right_x] = compute_fwhm_positions(positions, mean_vals)
% Compute FWHM of a mean amplitude curve using linear interpolation.
%   Returns FWHM in position units, and the left/right half-max x-coords.
%   Works on the depol (positive) amplitude curves.
    fw = NaN; left_x = NaN; right_x = NaN;
    valid = ~isnan(mean_vals);
    if sum(valid) < 3, return; end

    pos_v = positions(valid);
    val_v = mean_vals(valid);

    [pk, pk_idx] = max(val_v);
    if pk <= 0, return; end
    half_max = pk / 2;

    % Walk left from peak
    for j = pk_idx:-1:2
        if val_v(j-1) <= half_max
            % Linear interpolation between j-1 and j
            frac = (half_max - val_v(j-1)) / (val_v(j) - val_v(j-1));
            left_x = pos_v(j-1) + frac * (pos_v(j) - pos_v(j-1));
            break;
        end
    end
    if isnan(left_x), left_x = pos_v(1); end  % clamp to edge

    % Walk right from peak
    for j = pk_idx:numel(val_v)-1
        if val_v(j+1) <= half_max
            frac = (half_max - val_v(j+1)) / (val_v(j) - val_v(j+1));
            right_x = pos_v(j+1) - frac * (pos_v(j+1) - pos_v(j));
            break;
        end
    end
    if isnan(right_x), right_x = pos_v(end); end  % clamp to edge

    fw = right_x - left_x;
end


function fwhm_vec = compute_percell_fwhm(positions, amp_mat)
% Compute FWHM for each cell's amplitude profile.
%   positions: 1×N vector of position centers
%   amp_mat:   n_cells × N matrix of amplitude values
%   Returns:   n_cells × 1 vector of FWHM values (NaN if not computable)
    n = size(amp_mat, 1);
    fwhm_vec = NaN(n, 1);
    for k = 1:n
        row = amp_mat(k, :);
        [fw, ~, ~] = compute_fwhm_positions(positions, row);
        fwhm_vec(k) = fw;
    end
end


function stats = draw_amplitude_panel_no_stars(ax, positions, ctrl_data, ttl_data, ...
    col_c, col_t, font_ax, col_cf, col_tf, alpha_fill)
% Draw Gruntman-style connected-line amplitude panel WITHOUT asterisks.
%   Same as draw_amplitude_panel but omits per-position Wilcoxon asterisks.
%   Returns stats struct with mean/SEM for each position (used for asterisk placement).
    n_pos = numel(positions);
    stats = struct('pos', num2cell(positions), 'n_ctrl', 0, 'n_ttl', 0, ...
        'mean_ctrl', NaN, 'mean_ttl', NaN, 'sem_ctrl', NaN, 'sem_ttl', NaN);

    ctrl_mean = NaN(1, n_pos);  ctrl_sem = NaN(1, n_pos);
    ttl_mean  = NaN(1, n_pos);  ttl_sem  = NaN(1, n_pos);

    for i = 1:n_pos
        cv = ctrl_data(:, i); cv = cv(~isnan(cv));
        tv = ttl_data(:, i);  tv = tv(~isnan(tv));
        stats(i).n_ctrl = numel(cv);
        stats(i).n_ttl  = numel(tv);

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

    % Draw SEM bands
    valid_c = ~isnan(ctrl_mean);
    if any(valid_c)
        xc = positions(valid_c);
        fill(ax, [xc, fliplr(xc)], ...
            [ctrl_mean(valid_c)+ctrl_sem(valid_c), fliplr(ctrl_mean(valid_c)-ctrl_sem(valid_c))], ...
            col_cf, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none');
    end
    valid_t = ~isnan(ttl_mean);
    if any(valid_t)
        xt = positions(valid_t);
        fill(ax, [xt, fliplr(xt)], ...
            [ttl_mean(valid_t)+ttl_sem(valid_t), fliplr(ttl_mean(valid_t)-ttl_sem(valid_t))], ...
            col_tf, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none');
    end

    % Draw connected mean lines with markers
    if any(valid_c)
        plot(ax, positions(valid_c), ctrl_mean(valid_c), '-o', ...
            'Color', col_c, 'MarkerFaceColor', col_c, 'MarkerSize', 3, 'LineWidth', 0.8);
    end
    if any(valid_t)
        plot(ax, positions(valid_t), ttl_mean(valid_t), '-o', ...
            'Color', col_t, 'MarkerFaceColor', col_t, 'MarkerSize', 3, 'LineWidth', 0.8);
    end
end


function [pool_pvals, pool_centers] = compute_pooled_ranksum(ctrl_mat, ttl_mat)
% Compute rank-sum p-values using 3-position sliding-window pooling.
%   For each pool center (9 pools from 11 positions), concatenate amplitude
%   values from 3 consecutive positions across all cells, then rank-sum.
%   ctrl_mat, ttl_mat: n_cells × 11 matrices.
%   Returns: 1×9 p-value vector and 1×9 pool center positions.
    pool_centers = -4:4;  % 9 pools
    n_pools = numel(pool_centers);
    pool_pvals = NaN(1, n_pools);

    for p = 1:n_pools
        cols = p : p+2;  % 3 consecutive columns (1-indexed into 11-col matrix)

        % Concatenate values from 3 positions for each group
        c_vals = ctrl_mat(:, cols);
        c_vals = c_vals(:);
        c_vals = c_vals(~isnan(c_vals));

        t_vals = ttl_mat(:, cols);
        t_vals = t_vals(:);
        t_vals = t_vals(~isnan(t_vals));

        if numel(c_vals) >= 2 && numel(t_vals) >= 2
            pool_pvals(p) = ranksum(c_vals, t_vals);
        end
    end
end


function draw_pooled_asterisks(ax, pool_centers, pool_pvals, stats_11pos, col_stat, font_stat)
% Draw asterisks at pool center positions based on pooled p-values.
%   pool_centers: 1×9 vector of pool center positions (-4 to +4)
%   pool_pvals:   1×9 vector of p-values
%   stats_11pos:  1×11 struct array with mean_ctrl/sem_ctrl/mean_ttl/sem_ttl
%                 (used to determine asterisk y-position from the 11-pos means)
    positions_11 = -5:5;

    for p = 1:numel(pool_centers)
        if isnan(pool_pvals(p)) || pool_pvals(p) >= 0.05
            continue;
        end

        x = pool_centers(p);

        % Find the corresponding index in the 11-pos stats
        idx = find(positions_11 == x, 1);
        if isempty(idx), continue; end

        % y-position: above the higher of the two group means+SEM at this position
        y_c = stats_11pos(idx).mean_ctrl + stats_11pos(idx).sem_ctrl;
        y_t = stats_11pos(idx).mean_ttl  + stats_11pos(idx).sem_ttl;
        y_top = max([y_c, y_t]);
        if isnan(y_top), continue; end

        yl = ylim(ax);
        y_star = y_top + 0.06 * diff(yl);

        if pool_pvals(p) < 0.001,    str = '***';
        elseif pool_pvals(p) < 0.01, str = '**';
        else,                          str = '*';
        end
        text(ax, x, y_star, str, 'FontSize', font_stat, ...
            'HorizontalAlignment', 'center', 'Color', col_stat, 'FontWeight', 'bold');
    end
end


function write_pooled_stats(fid, label, pool_centers, pool_pvals)
% Write pooled-stats table to file.
    fprintf(fid, '%s:\n', label);
    fprintf(fid, '%-8s  %8s  %s\n', 'PoolCtr', 'p', 'Sig');
    fprintf(fid, '%s\n', repmat('-', 1, 30));
    for i = 1:numel(pool_centers)
        sig_str = '';
        if ~isnan(pool_pvals(i))
            if pool_pvals(i) < 0.001, sig_str = '***';
            elseif pool_pvals(i) < 0.01, sig_str = '**';
            elseif pool_pvals(i) < 0.05, sig_str = '*';
            end
        end
        fprintf(fid, '%+3d       %8.4f  %s\n', pool_centers(i), pool_pvals(i), sig_str);
    end
    fprintf(fid, '\n');
end


function [pooled_mat, pool_centers] = pool_3position(amp_mat)
% Pool amplitude values over sliding window of 3 consecutive positions.
%   amp_mat:      n_cells × 11 matrix (positions -5 to +5)
%   pooled_mat:   n_cells × 9 matrix (pool centers -4 to +4)
%   pool_centers: 1×9 vector [-4, -3, ..., +4]
%
%   Each pool(k) = mean of positions (k), (k+1), (k+2) for each cell.
%   NaN-aware: if all 3 positions are NaN, pool is NaN.
    pool_centers = -4:4;  % 9 pools
    n_cells = size(amp_mat, 1);
    n_pools = numel(pool_centers);
    pooled_mat = NaN(n_cells, n_pools);

    for p = 1:n_pools
        % Columns in the 11-position matrix: positions are -5:5 → cols 1:11
        % Pool center at pool_centers(p) corresponds to original columns
        %   (pool_centers(p)+6-1), (pool_centers(p)+6), (pool_centers(p)+6+1)
        col_start = p;  % = pool_centers(p) - (-5) = pool_centers(p) + 5, but 1-indexed → p
        cols = col_start : col_start + 2;  % 3 consecutive columns
        for c = 1:n_cells
            vals = amp_mat(c, cols);
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                pooled_mat(c, p) = mean(vals);
            end
        end
    end
end

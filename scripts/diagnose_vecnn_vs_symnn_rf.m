% DIAGNOSE_VECNN_VS_SYMNN_RF  Phase 0B: Compare VecNN vs SymNN PD for 4 ambiguous cells.
%
%   For 4 cells at DS grid boundaries where VecNN and SymNN disagree by 22.5°,
%   this script:
%     1. Per-cell diagnostic: side-by-side 1x11 bar flash traces under both PDs
%     2. Population-level propagation: re-run M2-aligned flash analysis twice
%        (VecNN PD vs substituted SymNN PD for affected cells), compare stats
%
%   The 4 affected cells (from ds_alignment_exploration.m):
%     #1  2025_10_23_10_31  ON TTL    VecNN=292.5° SymNN=270.0° (Δ=-22.5°)
%     #6  2025_10_28_13_54  ON TTL    VecNN=225.0° SymNN=202.5° (Δ=-22.5°)
%     #17 2025_11_05_15_54  OFF TTL   VecNN=135.0° SymNN=112.5° (Δ=-22.5°)
%     #20 2025_11_10_11_43  OFF ctrl  VecNN=67.5°  SymNN=90.0°  (Δ=+22.5°)
%
%   Outputs:
%     diag_vecnn_vs_symnn_per_cell.png  — 4×2 figure, per-cell bar flash traces
%     diag_vecnn_vs_symnn_pop_on.png    — population ON: VecNN vs SymNN
%     diag_vecnn_vs_symnn_pop_off.png   — population OFF: VecNN vs SymNN
%     diag_vecnn_vs_symnn_stats.txt     — per-position p-value comparison
%
%   Usage:
%     run('scripts/diagnose_vecnn_vs_symnn_rf.m')

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root   = '/Users/reiserm/Documents/ttl_1DRF';
preview_dir = fullfile(data_root, 'figure_previews');
if ~isfolder(preview_dir), mkdir(preview_dir); end

plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

% Load LUT
lut_path = fullfile(fileparts(mfilename('fullpath')), ...
    '..', 'src', 'analysis', 'protocol2', 'bar_lut.mat');
S_lut = load(lut_path, 'Tbl');
Tbl = S_lut.Tbl;

% Load batch results (late dataset only — has bar flash data)
S_batch = load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
results = S_batch.results;
n_cells = numel(results);
fprintf('Loaded %d cells from batch_results.mat\n', n_cells);

%% Identify the 4 affected cells
% Match by date_str prefix (batch results use date_str from load_protocol2_data)
affected_dates = {'2025_10_23_10_31', '2025_10_28_13_54', ...
                  '2025_11_05_15_54', '2025_11_10_11_43'};
affected_symnn_pd = [270.0, 202.5, 112.5, 90.0];  % SymNN PD in degrees

% Find indices in results
affected_idx = NaN(4, 1);
for a = 1:4
    for k = 1:n_cells
        if contains(results(k).date_str, affected_dates{a})
            affected_idx(a) = k;
            break;
        end
    end
end

fprintf('\nAffected cells:\n');
for a = 1:4
    k = affected_idx(a);
    if isnan(k)
        fprintf('  #%d: %s — NOT FOUND\n', a, affected_dates{a});
    else
        fprintf('  #%d: %s [%s] VecNN=%.1f° SymNN=%.1f°\n', ...
            a, results(k).date_str, results(k).group, ...
            results(k).pd_direction, affected_symnn_pd(a));
    end
end

%% Part 1: Per-cell diagnostic — side-by-side bar flash traces
fprintf('\n=== Part 1: Per-cell VecNN vs SymNN bar flash comparison ===\n');

sweep_opts.baseline_range = [1000 9000];
sweep_opts.stim_trim_end  = 7000;
sweep_opts.percentile     = 98;
on_threshold = 129;

fig_cells = figure('Position', [50 50 1200 800], 'Color', 'w');
t_cells = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for a = 1:4
    k = affected_idx(a);
    if isnan(k), continue; end

    exp_folder = results(k).folder;
    ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
        'pattern_order', 'func_order', 'metadata');

    orig_dir = pwd;
    cleanup_dir = onCleanup(@() cd(orig_dir));
    [~, ~, Log, ~, ~] = load_protocol2_data(exp_folder);
    f_data = Log.ADC.Volts(1, :);
    v_data = Log.ADC.Volts(2, :) * 10;

    is_on = ce.metadata.Frame > on_threshold;

    % Parse bar flash data
    bar_flash = parse_bar_flash_data(f_data, v_data, 0.5);

    % Get LUT info for VecNN PD
    pd_vecnn = results(k).pd_direction;
    [pd_info_vecnn] = find_pd_from_lut(Tbl, ce.pattern_order, ce.func_order, pd_vecnn);

    % Get LUT info for SymNN PD
    pd_symnn = affected_symnn_pd(a);
    [pd_info_symnn] = find_pd_from_lut(Tbl, ce.pattern_order, ce.func_order, pd_symnn);

    % Extract PD bar flash traces for both
    traces_vecnn = extract_pd_flash(bar_flash, pd_info_vecnn, is_on);
    traces_symnn = extract_pd_flash(bar_flash, pd_info_symnn, is_on);

    % Plot VecNN traces (left column)
    ax1 = nexttile(t_cells);
    plot_flash_panel(ax1, traces_vecnn, ...
        sprintf('#%d %s VecNN PD=%.1f°', k, results(k).date_str, pd_vecnn));

    % Plot SymNN traces (right column)
    ax2 = nexttile(t_cells);
    plot_flash_panel(ax2, traces_symnn, ...
        sprintf('#%d %s SymNN PD=%.1f°', k, results(k).date_str, pd_symnn));
end

title(t_cells, 'VecNN (left) vs SymNN (right) — PD Bar Flash Traces', 'FontSize', 12);

fname_cells = 'diag_vecnn_vs_symnn_per_cell.png';
exportgraphics(fig_cells, fullfile(preview_dir, fname_cells), 'Resolution', 150);
fprintf('Saved: %s\n', fname_cells);
close(fig_cells);

%% Part 2: Population-level propagation
%   Run M2-aligned flash pipeline twice:
%     A) All cells with original VecNN PD
%     B) 4 affected cells get SymNN PD; rest unchanged
%   Then compare per-position ctrl vs TTL rank-sum p-values.

fprintf('\n=== Part 2: Population-level propagation ===\n');

% Re-extract flash traces for all 25 cells with both PD assignments
% We need to reload data for affected cells only (others use existing results)

% First, compute SymNN PDs for the 4 affected cells
symnn_pds = NaN(n_cells, 1);  % NaN means "use VecNN" (no change)
for a = 1:4
    k = affected_idx(a);
    if ~isnan(k)
        symnn_pds(k) = affected_symnn_pd(a);
    end
end

% Collect M2-aligned PD flash traces for both scenarios
% Scenario A: VecNN (already in results — pd_flash_peak_aligned)
% Scenario B: SymNN for 4 cells (need to re-extract and re-align)

% For scenario B, we only need to re-extract the 4 affected cells
pd_flash_b = cell(n_cells, 1);  % will hold 11×N matrices per cell
peak_pos_b = NaN(n_cells, 1);

for k = 1:n_cells
    if isnan(symnn_pds(k))
        % Use existing VecNN result
        pd_flash_b{k} = results(k).pd_flash_peak_aligned;
        peak_pos_b(k) = results(k).peak_pos;
    else
        % Re-extract with SymNN PD
        fprintf('  Re-extracting cell %d (%s) with SymNN PD=%.1f° ...\n', ...
            k, results(k).date_str, symnn_pds(k));

        exp_folder = results(k).folder;
        ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
            'pattern_order', 'func_order', 'metadata');

        orig_dir = pwd;
        cleanup_dir = onCleanup(@() cd(orig_dir));
        [~, ~, Log, ~, ~] = load_protocol2_data(exp_folder);
        f_data = Log.ADC.Volts(1, :);
        v_data = Log.ADC.Volts(2, :) * 10;

        is_on = ce.metadata.Frame > on_threshold;

        bar_flash = parse_bar_flash_data(f_data, v_data, 0.5);

        pd_info = find_pd_from_lut(Tbl, ce.pattern_order, ce.func_order, symnn_pds(k));

        % Extract PD flash traces (same as batch pipeline)
        [flash_traces, ~] = extract_flash_traces_for_cell(bar_flash, pd_info, is_on);

        % Compute peak position (M2)
        A = compute_peak_amplitudes(flash_traces);
        [~, pk] = max(A);
        peak_pos_b(k) = pk;

        % Reindex to peak
        pd_flash_b{k} = reindex_to_peak(flash_traces, pk);

        fprintf('    Peak pos: VecNN=%d, SymNN=%d\n', results(k).peak_pos, pk);
    end
end

% Now compare population flash traces: Scenario A (VecNN) vs B (SymNN)
% Group by ON/OFF and ctrl/TTL

for on_off = ["ON", "OFF"]
    if on_off == "ON"
        ctrl_mask = [results.is_on]' & ~[results.is_ttl]';
        ttl_mask  = [results.is_on]' &  [results.is_ttl]';
    else
        ctrl_mask = ~[results.is_on]' & ~[results.is_ttl]';
        ttl_mask  = ~[results.is_on]' &  [results.is_ttl]';
    end

    n_ctrl = sum(ctrl_mask);
    n_ttl  = sum(ttl_mask);
    fprintf('\n%s cells: ctrl=%d, TTL=%d\n', on_off, n_ctrl, n_ttl);

    % Compute per-position peak amplitudes for both scenarios
    fig_pop = figure('Position', [50 50 1400 600], 'Color', 'w');
    t_pop = tiledlayout(2, 11, 'TileSpacing', 'compact', 'Padding', 'compact');

    stats_a = cell(11, 1);
    stats_b = cell(11, 1);

    for pos = 1:11
        % Scenario A (VecNN)
        ctrl_vals_a = NaN(n_ctrl, 1);
        ttl_vals_a  = NaN(n_ttl, 1);
        ci = 0; ti = 0;
        for k = 1:n_cells
            aligned = results(k).pd_flash_peak_aligned;
            if isempty(aligned), continue; end
            val = prctile(mean(aligned(pos, :), 2, 'omitnan'), 99.5);
            if ctrl_mask(k), ci = ci + 1; ctrl_vals_a(ci) = val; end
            if ttl_mask(k),  ti = ti + 1; ttl_vals_a(ti)  = val; end
        end
        ctrl_vals_a = ctrl_vals_a(1:ci);
        ttl_vals_a  = ttl_vals_a(1:ti);

        % Scenario B (SymNN for 4 cells)
        ctrl_vals_b = NaN(n_ctrl, 1);
        ttl_vals_b  = NaN(n_ttl, 1);
        ci = 0; ti = 0;
        for k = 1:n_cells
            aligned = pd_flash_b{k};
            if isempty(aligned), continue; end
            val = prctile(mean(aligned(pos, :), 2, 'omitnan'), 99.5);
            if ctrl_mask(k), ci = ci + 1; ctrl_vals_b(ci) = val; end
            if ttl_mask(k),  ti = ti + 1; ttl_vals_b(ti)  = val; end
        end
        ctrl_vals_b = ctrl_vals_b(1:ci);
        ttl_vals_b  = ttl_vals_b(1:ti);

        % Stats
        if numel(ctrl_vals_a) >= 2 && numel(ttl_vals_a) >= 2
            p_a = ranksum(ctrl_vals_a, ttl_vals_a);
        else
            p_a = NaN;
        end
        if numel(ctrl_vals_b) >= 2 && numel(ttl_vals_b) >= 2
            p_b = ranksum(ctrl_vals_b, ttl_vals_b);
        else
            p_b = NaN;
        end

        stats_a{pos} = struct('ctrl', ctrl_vals_a, 'ttl', ttl_vals_a, 'p', p_a);
        stats_b{pos} = struct('ctrl', ctrl_vals_b, 'ttl', ttl_vals_b, 'p', p_b);

        % Plot Scenario A (top row)
        ax1 = nexttile(t_pop, pos);
        hold(ax1, 'on');
        if ~isempty(ctrl_vals_a)
            scatter(ax1, ones(size(ctrl_vals_a))*0.8 + 0.1*randn(size(ctrl_vals_a)), ...
                ctrl_vals_a, 15, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
        end
        if ~isempty(ttl_vals_a)
            scatter(ax1, ones(size(ttl_vals_a))*1.2 + 0.1*randn(size(ttl_vals_a)), ...
                ttl_vals_a, 15, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
        end
        title(ax1, sprintf('pos %d', pos - 6), 'FontSize', 7);
        set(ax1, 'XTick', []);
        if pos == 1, ylabel(ax1, 'VecNN (A)', 'FontSize', 8); end
        if ~isnan(p_a)
            text(ax1, 1, max(ax1.YLim)*0.9, sprintf('p=%.3f', p_a), ...
                'FontSize', 6, 'HorizontalAlignment', 'center');
        end

        % Plot Scenario B (bottom row)
        ax2 = nexttile(t_pop, 11 + pos);
        hold(ax2, 'on');
        if ~isempty(ctrl_vals_b)
            scatter(ax2, ones(size(ctrl_vals_b))*0.8 + 0.1*randn(size(ctrl_vals_b)), ...
                ctrl_vals_b, 15, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
        end
        if ~isempty(ttl_vals_b)
            scatter(ax2, ones(size(ttl_vals_b))*1.2 + 0.1*randn(size(ttl_vals_b)), ...
                ttl_vals_b, 15, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
        end
        set(ax2, 'XTick', []);
        if pos == 1, ylabel(ax2, 'SymNN (B)', 'FontSize', 8); end
        if ~isnan(p_b)
            text(ax2, 1, max(ax2.YLim)*0.9, sprintf('p=%.3f', p_b), ...
                'FontSize', 6, 'HorizontalAlignment', 'center');
        end
    end

    title(t_pop, sprintf('%s — Per-Position Peak Amplitude: VecNN (top) vs SymNN (bottom)', ...
        on_off), 'FontSize', 11);

    fname_pop = sprintf('diag_vecnn_vs_symnn_pop_%s.png', lower(char(on_off)));
    exportgraphics(fig_pop, fullfile(preview_dir, fname_pop), 'Resolution', 150);
    fprintf('Saved: %s\n', fname_pop);
    close(fig_pop);
end

%% Part 3: Save stats comparison table
fid = fopen(fullfile(preview_dir, 'diag_vecnn_vs_symnn_stats.txt'), 'w');
fprintf(fid, 'VecNN vs SymNN — Per-Position Rank-Sum P-Values\n');
fprintf(fid, '================================================\n\n');

fprintf(fid, 'Affected cells:\n');
for a = 1:4
    k = affected_idx(a);
    if isnan(k), continue; end
    fprintf(fid, '  %s [%s] VecNN=%.1f° SymNN=%.1f° M2: VecNN=%d SymNN=%d\n', ...
        results(k).date_str, results(k).group, ...
        results(k).pd_direction, affected_symnn_pd(a), ...
        results(k).peak_pos, peak_pos_b(k));
end

fprintf(fid, '\nConclusion: ');
if all(peak_pos_b(~isnan(symnn_pds)) == [results(~isnan(symnn_pds')).peak_pos]')
    fprintf(fid, 'M2 peak positions identical → SymNN has no effect on RF alignment.\n');
else
    fprintf(fid, 'Some M2 peak positions differ → SymNN affects RF alignment.\n');
end

fclose(fid);
fprintf('\nSaved: diag_vecnn_vs_symnn_stats.txt\n');

fprintf('\n=== Phase 0B Complete ===\n');


%% ========================= Helper Functions ============================

function traces = extract_pd_flash(bar_flash, pd_info, is_on)
% EXTRACT_PD_FLASH  Extract 1x11 PD bar flash traces using pd_info mapping.
%   Returns 11×N matrix (rows = positions, cols = repetitions).

    pd_orient = pd_info.pd_orientation;
    pos_order = pd_info.pos_order;

    % Find the PD orientation index
    % bar_flash is organized by orientation then position
    % Each orientation has 11 positions × N reps
    n_reps = size(bar_flash, 2) - 1;  % last col is mean

    traces = [];
    for pos = 1:11
        flash_row = (pd_orient - 1) * 11 + pos_order(pos);
        if flash_row > size(bar_flash, 1), continue; end
        mean_trace = bar_flash{flash_row, n_reps + 1};
        if isempty(mean_trace), continue; end

        % Baseline-subtract
        bl = mean(mean_trace(1:min(5000, numel(mean_trace))));
        traces(pos, :) = mean_trace - bl; %#ok<AGROW>
    end
end


function plot_flash_panel(ax, traces, title_str)
% PLOT_FLASH_PANEL  Plot 1x11 bar flash traces in a compact panel.

    hold(ax, 'on');
    if isempty(traces)
        title(ax, [title_str ' — NO DATA'], 'FontSize', 7);
        return;
    end

    n_pos = size(traces, 1);
    n_samp = size(traces, 2);
    t_ms = (0:n_samp-1) * 0.1;  % ms

    cmap = parula(n_pos);
    y_offset = 0;
    offsets = zeros(n_pos, 1);

    % Stack traces with offsets for visibility
    max_amp = max(abs(traces(:)));
    if max_amp == 0, max_amp = 1; end
    spacing = max_amp * 0.6;

    for pos = 1:n_pos
        offsets(pos) = -(pos - 6) * spacing;
        plot(ax, t_ms, traces(pos, :) + offsets(pos), ...
            'Color', cmap(pos, :), 'LineWidth', 0.8);
    end

    % Stim timing
    xline(ax, 500.1, 'g-', 'LineWidth', 0.5);  % onset = sample 5001 = 500.1ms
    xline(ax, 580.1, 'g-', 'LineWidth', 0.5);  % offset = sample 5801

    % Labels
    yticks_pos = offsets;
    set(ax, 'YTick', yticks_pos, 'YTickLabel', arrayfun(@(x) sprintf('%d', x-6), 1:n_pos, 'uni', 0));
    xlabel(ax, 'Time (ms)', 'FontSize', 6);
    title(ax, title_str, 'FontSize', 7);
    set(ax, 'FontSize', 6);
end


function [flash_traces, bl_subtracted] = extract_flash_traces_for_cell(bar_flash, pd_info, is_on)
% EXTRACT_FLASH_TRACES_FOR_CELL  Extract 11-position PD bar flash traces.
%   Matches the batch_analyze_1DRF extract_flash_traces logic.

    pd_orient = pd_info.pd_orientation;
    pos_order = pd_info.pos_order;
    n_reps = size(bar_flash, 2) - 1;

    flash_traces = [];
    bl_subtracted = [];

    for pos = 1:11
        flash_row = (pd_orient - 1) * 11 + pos_order(pos);
        if flash_row > size(bar_flash, 1), continue; end

        % Average across reps
        traces_all = [];
        for rep = 1:n_reps
            tr = bar_flash{flash_row, rep};
            if isempty(tr), continue; end
            traces_all = [traces_all; tr(:)']; %#ok<AGROW>
        end

        if isempty(traces_all), continue; end
        mean_trace = mean(traces_all, 1);

        % Baseline-subtract (samples 1:5000)
        bl = mean(mean_trace(1:min(5000, numel(mean_trace))));
        bl_sub = mean_trace - bl;

        if isempty(flash_traces)
            flash_traces = NaN(11, numel(mean_trace));
            bl_subtracted = NaN(11, numel(bl_sub));
        end

        flash_traces(pos, 1:numel(mean_trace)) = mean_trace;
        bl_subtracted(pos, 1:numel(bl_sub)) = bl_sub;
    end

    flash_traces = bl_subtracted;  % return baseline-subtracted
end


function A = compute_peak_amplitudes(flash_traces)
% COMPUTE_PEAK_AMPLITUDES  99.5th percentile in response window per position.

    A = NaN(size(flash_traces, 1), 1);
    for pos = 1:size(flash_traces, 1)
        tr = flash_traces(pos, :);
        if all(isnan(tr)), continue; end
        resp_win = tr(5001:min(6551, numel(tr)));
        bl = mean(tr(1:min(5000, numel(tr))), 'omitnan');
        A(pos) = prctile(resp_win, 99.5) - bl;
    end

    % Threshold negatives to 0
    A(A < 0) = 0;
end


function aligned = reindex_to_peak(traces, peak_pos)
% REINDEX_TO_PEAK  Shift 11×N traces so peak_pos maps to row 6 (center).

    aligned = NaN(size(traces));
    shift = 6 - peak_pos;

    for pos = 1:11
        src = pos - shift;
        if src >= 1 && src <= 11
            aligned(pos, :) = traces(src, :);
        end
    end
end

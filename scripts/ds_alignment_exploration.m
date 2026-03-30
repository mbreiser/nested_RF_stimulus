% DS_ALIGNMENT_EXPLORATION  Compare 4 PD extraction/alignment methods.
%
%   Loads early (pre-bar-flash, 23 cells) and late (main 1DRF, 25 cells)
%   batch results, then computes 4 PD methods per cell:
%     1. VecNN     — vector sum, snap to 22.5 deg grid (current method)
%     2. SymNN     — symmetry axis from 16 candidates, snap to 22.5 deg
%     3. VecInterp — vector sum (continuous), 2.5 deg interpolated grid
%     4. SymInterp — brute-force symmetry axis at 2.5 deg resolution
%
%   Outputs (~12 PNGs + disagreement table):
%     Individual cell polar plots with 4 PD arrows
%     Population polar tuning: 4 methods x ON/OFF (ctrl vs TTL)
%     FWHM comparison: 4 methods x 4 groups
%     Symmetry score comparison: 4 methods x 4 groups
%     Disagreement table (text file)
%
%   Usage:
%     run('scripts/ds_alignment_exploration.m')

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root   = '/Users/reiserm/Documents/ttl_1DRF';
preview_dir = fullfile(data_root, 'figure_previews');
if ~isfolder(preview_dir), mkdir(preview_dir); end

%% Load both batch results
fprintf('Loading batch results...\n');

S_late = load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
late = S_late.results;
fprintf('  Late batch (main): %d cells\n', numel(late));

S_early = load(fullfile(data_root, 'pre-bar-flash', 'population_results', ...
    'batch_results_pre_bf.mat'), 'results');
early = S_early.results;
fprintf('  Early batch (pre-BF): %d cells\n', numel(early));

%% Build combined cell list with common fields
all_cells = build_combined_cells(early, late);
n_cells = numel(all_cells);
fprintf('  Combined: %d cells\n\n', n_cells);

%% Compute 4 PD methods for each cell
fprintf('=== Computing 4 PD methods per cell ===\n\n');

method_names = {'vec_nn', 'sym_nn', 'vec_interp', 'sym_interp'};
method_labels = {'VecNN', 'SymNN', 'VecInterp', 'SymInterp'};

for k = 1:n_cells
    d = all_cells(k).max_v_aligned;
    angles_16    = d(:, 1);
    responses_16 = d(:, 2);

    all_cells(k).methods = compute_pd_four_methods(angles_16, responses_16);
end

fprintf('Done. All %d cells processed.\n\n', n_cells);

%% Compute per-cell max disagreement (needed for flagged cell selection)
max_delta = compute_max_delta(all_cells, method_names);

%% Output A: Individual cell polar plots for ALL flagged cells (MaxΔ > 22.5°)
fprintf('=== Generating individual cell polar plots (flagged cells) ===\n');

flagged_idx = find(max_delta > 22.5);
fprintf('  %d flagged cells with MaxDelta > 22.5 deg\n', numel(flagged_idx));

for ei = 1:numel(flagged_idx)
    k = flagged_idx(ei);
    c = all_cells(k);
    fig = plot_individual_4pd(c, max_delta(k));
    fname = sprintf('ds_explore_flagged_%s.png', c.date_str);
    exportgraphics(fig, fullfile(preview_dir, fname), 'Resolution', 150);
    close(fig);
    fprintf('  Saved: %s  (MaxDelta=%.1f°)\n', fname, max_delta(k));
end

%% Output B: Disagreement table
fprintf('\n=== Method disagreement table ===\n\n');

table_path = fullfile(preview_dir, 'ds_explore_disagreement.txt');
write_disagreement_table(all_cells, method_names, method_labels, table_path, max_delta);
fprintf('  Saved: ds_explore_disagreement.txt\n');

%% Output C: Population polar plots (4 methods x ON/OFF)
fprintf('\n=== Generating population polar plots ===\n');

polar_opts = struct();
polar_opts.stat_method = 'mean_sem';

for mi = 1:4
    mname = method_names{mi};
    mlabel = method_labels{mi};

    for on_off = ["ON", "OFF"]
        if on_off == "ON"
            mask = [all_cells.is_on];
        else
            mask = ~[all_cells.is_on];
        end

        ctrl_mask = mask & ~[all_cells.is_ttl];
        ttl_mask  = mask &  [all_cells.is_ttl];

        if ~any(ctrl_mask) && ~any(ttl_mask)
            continue;
        end

        aligned_ctrl = extract_aligned(all_cells(ctrl_mask), mname);
        aligned_ttl  = extract_aligned(all_cells(ttl_mask), mname);

        ptitle = sprintf('%s — %s PD-Aligned Polar (ctrl n=%d, TTL n=%d)', ...
            mlabel, on_off, sum(ctrl_mask), sum(ttl_mask));

        fig = plot_polar_population(aligned_ctrl, aligned_ttl, ptitle, polar_opts);

        fname = sprintf('ds_explore_polar_%s_%s.png', mlabel, lower(char(on_off)));
        exportgraphics(fig, fullfile(preview_dir, fname), 'Resolution', 150);
        close(fig);
        fprintf('  Saved: %s\n', fname);
    end
end

%% Output D: FWHM and symmetry comparison
fprintf('\n=== Generating FWHM and symmetry comparison plots ===\n');

fig_fwhm = plot_metric_comparison(all_cells, method_names, method_labels, ...
    'fwhm_deg', 'FWHM (deg)', 'Direction Tuning FWHM — Method Comparison');
exportgraphics(fig_fwhm, fullfile(preview_dir, 'ds_explore_fwhm_comparison.png'), ...
    'Resolution', 150);
close(fig_fwhm);
fprintf('  Saved: ds_explore_fwhm_comparison.png\n');

fig_sym = plot_metric_comparison(all_cells, method_names, method_labels, ...
    'sym_score', 'Symmetry Score (0=perfect)', ...
    'Tuning Curve Symmetry — Method Comparison');
exportgraphics(fig_sym, fullfile(preview_dir, 'ds_explore_symmetry_comparison.png'), ...
    'Resolution', 150);
close(fig_sym);
fprintf('  Saved: ds_explore_symmetry_comparison.png\n');

%% Output E: DSI (PD-ND) comparison — 4x1 layout (wider panels)
fprintf('\n=== Generating per-method DSI comparison ===\n');

fig_dsi = plot_metric_comparison(all_cells, method_names, method_labels, ...
    'dsi_pdnd', 'DSI (PD−ND)', 'Direction Selectivity (PD−ND) — Method Comparison');
exportgraphics(fig_dsi, fullfile(preview_dir, 'ds_explore_dsi_comparison.png'), ...
    'Resolution', 150);
close(fig_dsi);
fprintf('  Saved: ds_explore_dsi_comparison.png\n');

%% Output F: Aspect ratio comparison — 2x2 layout (squarer panels)
fprintf('\n=== Generating per-method aspect ratio comparison ===\n');

fig_ar = plot_metric_comparison_2x2(all_cells, method_names, method_labels, ...
    'aspect_ratio', 'Aspect Ratio (PD / mean ortho)', ...
    'Tuning Aspect Ratio — Method Comparison');
exportgraphics(fig_ar, fullfile(preview_dir, 'ds_explore_aspect_ratio.png'), ...
    'Resolution', 150);
close(fig_ar);
fprintf('  Saved: ds_explore_aspect_ratio.png\n');

%% Output G: Combined 2x2 polar summary — all 4 groups × 4 methods
fprintf('\n=== Generating combined 2x2 polar summary ===\n');

fig_polar_summ = plot_polar_4method_summary(all_cells, method_names, method_labels);
exportgraphics(fig_polar_summ, fullfile(preview_dir, 'ds_explore_polar_summary_2x2.png'), ...
    'Resolution', 150);
close(fig_polar_summ);
fprintf('  Saved: ds_explore_polar_summary_2x2.png\n');

%% Output H: Population alignment quality — vector strength and symmetry of mean curves
fprintf('\n=== Computing population alignment quality (mean curve metrics) ===\n');

fig_qual = plot_alignment_quality(all_cells, method_names, method_labels);
exportgraphics(fig_qual, fullfile(preview_dir, 'ds_explore_alignment_quality.png'), ...
    'Resolution', 150);
close(fig_qual);
fprintf('  Saved: ds_explore_alignment_quality.png\n');

%% Output I: Publication-quality SymInterp polar + statistics
fprintf('\n=== Generating publication-quality SymInterp figures ===\n');

fig_polar = plot_polar_summary_single(all_cells, 'sym_interp');
exportgraphics(fig_polar, fullfile(preview_dir, 'ds_summary_polar_syminterp.png'), ...
    'Resolution', 300);
close(fig_polar);
fprintf('  Saved: ds_summary_polar_syminterp.png (300 DPI)\n');

fig_stats = plot_dsi_ar_summary(all_cells, 'sym_interp');
exportgraphics(fig_stats, fullfile(preview_dir, 'ds_summary_stats_syminterp.png'), ...
    'Resolution', 300);
close(fig_stats);
fprintf('  Saved: ds_summary_stats_syminterp.png (300 DPI)\n');

fprintf('\n=== Done. All outputs in: %s ===\n', preview_dir);


%% ========================= Helper Functions ============================

function cells = build_combined_cells(early, late)
% BUILD_COMBINED_CELLS  Merge early and late results with common fields.

    cells = struct([]);
    for k = 1:numel(early)
        s.date_str      = early(k).date_str;
        s.is_on         = early(k).is_on;
        s.is_ttl        = early(k).is_ttl;
        s.group         = early(k).group;
        s.max_v_aligned = early(k).max_v_aligned;
        s.batch         = 'early';
        s.methods       = struct();
        if isempty(cells)
            cells = s;
        else
            cells(end + 1) = s; %#ok<AGROW>
        end
    end
    for k = 1:numel(late)
        s.date_str      = late(k).date_str;
        s.is_on         = late(k).is_on;
        s.is_ttl        = late(k).is_ttl;
        s.group         = late(k).group;
        s.max_v_aligned = late(k).max_v_aligned;
        s.batch         = 'late';
        s.methods       = struct();
        cells(end + 1) = s; %#ok<AGROW>
    end

end


function max_delta = compute_max_delta(cells, method_names)
% COMPUTE_MAX_DELTA  Max pairwise angular difference across 4 methods.

    n = numel(cells);
    n_methods = numel(method_names);
    max_delta = NaN(n, 1);

    for k = 1:n
        pds = NaN(1, n_methods);
        for mi = 1:n_methods
            pds(mi) = cells(k).methods.(method_names{mi}).pd_deg;
        end
        max_d = 0;
        for i = 1:n_methods
            for j = i+1:n_methods
                d = abs(pds(i) - pds(j));
                d = min(d, 360 - d);
                max_d = max(max_d, d);
            end
        end
        max_delta(k) = max_d;
    end

end


function fig = plot_individual_4pd(cell_data, max_delta_val)
% PLOT_INDIVIDUAL_4PD  Polar plot with 4 colored PD direction arrows.

    d = cell_data.max_v_aligned;
    angles = d(:, 1);
    responses = d(:, 2);
    m = cell_data.methods;

    % Close the curve for plotting
    angles_closed = [angles; angles(1) + 2*pi];
    resp_closed   = [responses; responses(1)];

    fig = figure('Position', [100 200 550 500]);
    ax = polaraxes;
    hold(ax, 'on');

    % Plot raw tuning curve
    polarplot(ax, angles_closed, resp_closed, '-', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    polarscatter(ax, angles, responses, 40, [0.5 0.5 0.5], 'filled');

    % PD arrows — 4 methods
    pd_colors = [0.12 0.47 0.71;    % blue  — VecNN
                 0.84 0.15 0.16;    % red   — SymNN
                 0.17 0.63 0.17;    % green — VecInterp
                 0.58 0.40 0.74];   % purple — SymInterp
    method_fields = {'vec_nn', 'sym_nn', 'vec_interp', 'sym_interp'};
    method_labels = {'VecNN', 'SymNN', 'VecInterp', 'SymInterp'};

    rmax = max(responses) * 1.15;
    legend_entries = gobjects(4, 1);

    for mi = 1:4
        pd_rad = m.(method_fields{mi}).pd_angle;
        pd_deg = m.(method_fields{mi}).pd_deg;
        h = polarplot(ax, [pd_rad pd_rad], [0 rmax], '-', ...
            'Color', pd_colors(mi, :), 'LineWidth', 2.5);
        legend_entries(mi) = h;
        method_labels{mi} = sprintf('%s (%.1f°)', method_labels{mi}, pd_deg);
    end

    legend(ax, legend_entries, method_labels, ...
        'Location', 'southoutside', 'NumColumns', 2, 'FontSize', 9);

    if nargin >= 2 && ~isempty(max_delta_val)
        title(ax, sprintf('%s [%s] MaxΔ=%.1f°', cell_data.date_str, ...
            strrep(cell_data.group, '_', ' '), max_delta_val), 'FontSize', 12);
    else
        title(ax, sprintf('%s [%s]', cell_data.date_str, ...
            strrep(cell_data.group, '_', ' ')), 'FontSize', 12);
    end

end


function aligned_cell_array = extract_aligned(cells, method_name)
% EXTRACT_ALIGNED  Get aligned tuning data for a specific method.

    n = numel(cells);
    aligned_cell_array = cell(1, n);
    for k = 1:n
        aligned_cell_array{k} = cells(k).methods.(method_name).aligned;
    end

end


function write_disagreement_table(cells, method_names, method_labels, filepath, max_delta)
% WRITE_DISAGREEMENT_TABLE  Console + text file of PD angles and disagreements.

    n = numel(cells);
    n_methods = numel(method_names);

    % Collect PD angles
    pd_angles = NaN(n, n_methods);
    for k = 1:n
        for mi = 1:n_methods
            pd_angles(k, mi) = cells(k).methods.(method_names{mi}).pd_deg;
        end
    end

    % Write to file and console
    fid = fopen(filepath, 'w');

    header = sprintf('%-4s  %-20s  %-12s', '#', 'Date', 'Group');
    for mi = 1:n_methods
        header = [header, sprintf('  %10s', method_labels{mi})]; %#ok<AGROW>
    end
    header = [header, sprintf('  %8s  %4s', 'MaxDelta', 'Flag')];

    separator = repmat('-', 1, numel(header) + 5);

    for dest = {1, fid}  % console then file
        d = dest{1};
        fprintf(d, '\nDS ALIGNMENT EXPLORATION — PD DISAGREEMENT TABLE\n');
        fprintf(d, '(All PD angles in aligned-curve coordinates, degrees)\n\n');
        fprintf(d, '%s\n', header);
        fprintf(d, '%s\n', separator);

        for k = 1:n
            flag = '';
            if max_delta(k) > 22.5
                flag = '*';
            end
            line = sprintf('%-4d  %-20s  %-12s', k, cells(k).date_str, ...
                strrep(cells(k).group, '_', ' '));
            for mi = 1:n_methods
                line = [line, sprintf('  %10.1f', pd_angles(k, mi))]; %#ok<AGROW>
            end
            line = [line, sprintf('  %8.1f  %4s', max_delta(k), flag)]; %#ok<AGROW>
            fprintf(d, '%s\n', line);
        end

        fprintf(d, '%s\n', separator);
        n_flagged = sum(max_delta > 22.5);
        fprintf(d, '\nFlagged (MaxDelta > 22.5 deg): %d / %d cells\n', n_flagged, n);
        fprintf(d, 'Mean MaxDelta: %.1f deg\n', mean(max_delta));

        % Per-group summary
        groups = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
        group_labels = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
        fprintf(d, '\nMean MaxDelta by group:\n');
        for g = 1:numel(groups)
            mask = strcmp({cells.group}, groups{g});
            if any(mask)
                fprintf(d, '  %-12s: %.1f deg (n=%d)\n', ...
                    group_labels{g}, mean(max_delta(mask)), sum(mask));
            end
        end
        fprintf(d, '\n');
    end

    fclose(fid);

end


function fig = plot_metric_comparison(cells, method_names, method_labels, ...
    metric_field, ylabel_str, title_str)
% PLOT_METRIC_COMPARISON  1x4 tiled box+dot plot: one panel per method.

    groups = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
    group_labels = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
    colors = [0 0 0; 1 0 0; 0.4 0.4 0.4; 0.8 0.2 0.2];

    n_methods = numel(method_names);
    n_groups = numel(groups);

    fig = figure('Name', title_str, 'Position', [50 200 1200 350]);
    t = tiledlayout(1, n_methods, 'TileSpacing', 'compact', 'Padding', 'compact');

    for mi = 1:n_methods
        nexttile; hold on;

        all_vals = [];
        all_grp_idx = [];
        group_data = cell(1, n_groups);

        for g = 1:n_groups
            mask = strcmp({cells.group}, groups{g});
            grp_idx = find(mask);
            vals = NaN(1, numel(grp_idx));
            for ki = 1:numel(grp_idx)
                vals(ki) = cells(grp_idx(ki)).methods.(method_names{mi}).(metric_field);
            end
            vals = vals(~isnan(vals));
            group_data{g} = vals;
            ng = numel(vals);
            all_vals = [all_vals, vals]; %#ok<AGROW>
            all_grp_idx = [all_grp_idx, repmat(g, 1, ng)]; %#ok<AGROW>
        end

        if isempty(all_vals)
            title(method_labels{mi});
            continue;
        end

        boxchart(all_grp_idx(:), all_vals(:), ...
            'BoxFaceColor', [0.9 0.9 0.9], 'MarkerStyle', 'none', ...
            'BoxWidth', 0.5);

        for g = 1:n_groups
            idx = all_grp_idx == g;
            vals = all_vals(idx);
            x = g + 0.15 * (rand(size(vals)) - 0.5);
            scatter(x, vals, 50, colors(g, :), 'filled', ...
                'MarkerFaceAlpha', 0.8);
        end

        set(gca, 'XTick', 1:n_groups, 'XTickLabel', group_labels);
        xlim([0.5, n_groups + 0.5]);
        box off;
        set(gca, 'TickDir', 'out', 'FontSize', 10);
        title(method_labels{mi}, 'FontSize', 12);

        if mi == 1
            ylabel(ylabel_str);
        end

        % Wilcoxon rank-sum: ON ctrl vs TTL, OFF ctrl vs TTL
        yl = ylim;
        y_stat = yl(2) + 0.05 * diff(yl);
        dy = 0.08 * diff(yl);

        if numel(group_data{1}) >= 2 && numel(group_data{2}) >= 2
            p = ranksum(group_data{1}, group_data{2});
            draw_bracket_local(1, 2, y_stat, p);
        end
        if numel(group_data{3}) >= 2 && numel(group_data{4}) >= 2
            p = ranksum(group_data{3}, group_data{4});
            draw_bracket_local(3, 4, y_stat + dy, p);
        end
        ylim([yl(1), y_stat + 2.5*dy]);
    end

    title(t, title_str, 'FontSize', 14);

end


function draw_bracket_local(x1, x2, y, p)
% DRAW_BRACKET_LOCAL  Bracket with p-value annotation.

    line([x1 x1 x2 x2], [y - 0.02*abs(y), y, y, y - 0.02*abs(y)], ...
        'Color', 'k', 'LineWidth', 1);
    if p < 0.001
        p_str = 'p < 0.001';
    elseif p < 0.01
        p_str = sprintf('p = %.3f', p);
    else
        p_str = sprintf('p = %.2f', p);
    end
    text(mean([x1 x2]), y, p_str, 'FontSize', 8, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

end


function fig = plot_metric_comparison_2x2(cells, method_names, method_labels, ...
    metric_field, ylabel_str, title_str)
% PLOT_METRIC_COMPARISON_2X2  2x2 tiled box+dot plot: one panel per method.
%   Same content as plot_metric_comparison but in a 2x2 layout for squarer
%   panels (better for metrics like aspect ratio).

    groups = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
    group_labels = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
    colors = [0 0 0; 1 0 0; 0.4 0.4 0.4; 0.8 0.2 0.2];

    n_methods = numel(method_names);
    n_groups = numel(groups);

    fig = figure('Name', title_str, 'Position', [50 100 900 700]);
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for mi = 1:n_methods
        nexttile; hold on;

        all_vals = [];
        all_grp_idx = [];
        group_data = cell(1, n_groups);

        for g = 1:n_groups
            mask = strcmp({cells.group}, groups{g});
            grp_idx = find(mask);
            vals = NaN(1, numel(grp_idx));
            for ki = 1:numel(grp_idx)
                vals(ki) = cells(grp_idx(ki)).methods.(method_names{mi}).(metric_field);
            end
            vals = vals(~isnan(vals));
            group_data{g} = vals;
            ng = numel(vals);
            all_vals = [all_vals, vals]; %#ok<AGROW>
            all_grp_idx = [all_grp_idx, repmat(g, 1, ng)]; %#ok<AGROW>
        end

        if isempty(all_vals)
            title(method_labels{mi});
            continue;
        end

        boxchart(all_grp_idx(:), all_vals(:), ...
            'BoxFaceColor', [0.9 0.9 0.9], 'MarkerStyle', 'none', ...
            'BoxWidth', 0.5);

        for g = 1:n_groups
            idx = all_grp_idx == g;
            vals = all_vals(idx);
            x = g + 0.15 * (rand(size(vals)) - 0.5);
            scatter(x, vals, 50, colors(g, :), 'filled', ...
                'MarkerFaceAlpha', 0.8);
        end

        set(gca, 'XTick', 1:n_groups, 'XTickLabel', group_labels);
        xlim([0.5, n_groups + 0.5]);
        box off;
        set(gca, 'TickDir', 'out', 'FontSize', 10);
        title(method_labels{mi}, 'FontSize', 12);
        ylabel(ylabel_str);

        % Wilcoxon rank-sum: ON ctrl vs TTL, OFF ctrl vs TTL
        yl = ylim;
        y_stat = yl(2) + 0.05 * diff(yl);
        dy = 0.08 * diff(yl);

        if numel(group_data{1}) >= 2 && numel(group_data{2}) >= 2
            p = ranksum(group_data{1}, group_data{2});
            draw_bracket_local(1, 2, y_stat, p);
        end
        if numel(group_data{3}) >= 2 && numel(group_data{4}) >= 2
            p = ranksum(group_data{3}, group_data{4});
            draw_bracket_local(3, 4, y_stat + dy, p);
        end
        ylim([yl(1), y_stat + 2.5*dy]);
    end

    title(t, title_str, 'FontSize', 14);

end


function fig = plot_alignment_quality(all_cells, method_names, method_labels)
% PLOT_ALIGNMENT_QUALITY  Compare population vector strength and symmetry.
%   For each method x group, stacks all aligned curves, computes the
%   population mean, then measures vector strength (R) and symmetry.
%   A method that aligns cells better produces a more concentrated
%   (higher R) and more symmetric population mean tuning curve.

    groups     = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
    grp_labels = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
    colors     = [0 0 0; 0.85 0 0; 0.15 0.35 0.70; 0.90 0.50 0.10];

    n_methods = numel(method_names);
    n_groups  = numel(groups);

    vs_mat  = NaN(n_methods, n_groups);
    sym_mat = NaN(n_methods, n_groups);

    for mi = 1:n_methods
        for g = 1:n_groups
            mask    = strcmp({all_cells.group}, groups{g});
            cells_g = all_cells(mask);
            ng      = numel(cells_g);
            if ng < 2, continue; end

            % Stack aligned curves
            theta = [];
            vals  = [];
            for k = 1:ng
                a = cells_g(k).methods.(method_names{mi}).aligned;
                if isempty(theta), theta = a(:, 1); end
                vals = [vals, a(:, 2)]; %#ok<AGROW>
            end

            mean_resp = mean(vals, 2, 'omitnan');

            % Vector strength of population mean curve
            vs_mat(mi, g) = compute_vector_strength(theta, mean_resp);

            % Symmetry of population mean curve
            sym_mat(mi, g) = compute_symmetry_of_mean(theta, mean_resp);
        end
    end

    % --- Figure ---
    fig = figure('Name', 'Alignment Quality Comparison', ...
        'Position', [80 80 950 700], 'Color', 'w');

    bw = 0.18;  % bar width
    x  = 1:n_methods;

    % Top panel: Vector strength
    subplot(2, 1, 1); hold on;
    for g = 1:n_groups
        bx = x + (g - (n_groups + 1) / 2) * bw;
        bar(bx, vs_mat(:, g), bw, ...
            'FaceColor', colors(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    end
    % Value labels on bars
    for mi = 1:n_methods
        for g = 1:n_groups
            bx = mi + (g - (n_groups + 1) / 2) * bw;
            v = vs_mat(mi, g);
            if ~isnan(v)
                text(bx, v + 0.01, sprintf('%.3f', v), ...
                    'FontSize', 7, 'HorizontalAlignment', 'center', ...
                    'Color', colors(g, :));
            end
        end
    end
    set(gca, 'XTick', x, 'XTickLabel', method_labels, 'FontSize', 11);
    ylabel('Vector strength R (0–1)', 'FontSize', 11);
    title('Population Vector Strength — higher = tighter alignment', 'FontSize', 13);
    legend(grp_labels, 'Location', 'northeast', 'FontSize', 9, 'Box', 'on');
    box off; set(gca, 'TickDir', 'out');

    % Bottom panel: Symmetry
    subplot(2, 1, 2); hold on;
    for g = 1:n_groups
        bx = x + (g - (n_groups + 1) / 2) * bw;
        bar(bx, sym_mat(:, g), bw, ...
            'FaceColor', colors(g, :), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    end
    % Value labels on bars
    for mi = 1:n_methods
        for g = 1:n_groups
            bx = mi + (g - (n_groups + 1) / 2) * bw;
            v = sym_mat(mi, g);
            if ~isnan(v)
                text(bx, v + 0.005, sprintf('%.3f', v), ...
                    'FontSize', 7, 'HorizontalAlignment', 'center', ...
                    'Color', colors(g, :));
            end
        end
    end
    set(gca, 'XTick', x, 'XTickLabel', method_labels, 'FontSize', 11);
    ylabel('Symmetry score (0 = perfect)', 'FontSize', 11);
    title('Population Symmetry — lower = more symmetric alignment', 'FontSize', 13);
    legend(grp_labels, 'Location', 'northeast', 'FontSize', 9, 'Box', 'on');
    box off; set(gca, 'TickDir', 'out');

    sgtitle({'Alignment Quality: Vector Strength and Symmetry of Mean Tuning Curves', ...
             '(computed on population-averaged aligned data)'}, ...
        'FontSize', 14, 'FontWeight', 'bold');

    % --- Print summary table to console ---
    fprintf('\nAlignment Quality — Vector strength of population mean curves (R):\n');
    fprintf('%-12s', '');
    for g = 1:n_groups, fprintf('  %10s', grp_labels{g}); end
    fprintf('\n');
    for mi = 1:n_methods
        fprintf('%-12s', method_labels{mi});
        for g = 1:n_groups
            fprintf('  %10.4f', vs_mat(mi, g));
        end
        fprintf('\n');
    end

    fprintf('\nAlignment Quality — Symmetry of population mean curves:\n');
    fprintf('%-12s', '');
    for g = 1:n_groups, fprintf('  %10s', grp_labels{g}); end
    fprintf('\n');
    for mi = 1:n_methods
        fprintf('%-12s', method_labels{mi});
        for g = 1:n_groups
            fprintf('  %10.4f', sym_mat(mi, g));
        end
        fprintf('\n');
    end
    fprintf('\n');

end


function R = compute_vector_strength(theta, mean_resp)
% COMPUTE_VECTOR_STRENGTH  Mean resultant length of a directional response.
%   R = |sum(resp .* exp(i*theta))| / sum(resp)
%   Ranges [0,1]: 0 = uniform, 1 = all response at one direction.
%   Standard circular statistic for concentration (complement of circular
%   variance: CV = 1 - R). Higher R = tighter/better-aligned peak.

    resp = max(mean_resp(:), 0);  % threshold negatives to zero
    total = sum(resp);
    if total == 0
        R = NaN;
        return;
    end
    R = abs(sum(resp .* exp(1i * theta(:)))) / total;

end


function score = compute_symmetry_of_mean(theta, mean_resp)
% Symmetry score of a mean curve: fold around PD at pi/2.
%   Returns 0 for perfectly symmetric, higher for more asymmetric.

    aw = [theta(end) - 2*pi; theta; theta(1) + 2*pi];
    rw = [mean_resp(end);     mean_resp; mean_resp(1)];

    % Mirror each angle across pi/2: mirror = pi - angle
    mirror_angles = mod(pi - theta, 2*pi);
    mirror_resp   = interp1(aw, rw, mirror_angles, 'linear');

    ssd = sum((mean_resp - mirror_resp).^2);
    ss  = sum(mean_resp.^2);

    if ss > 0
        score = ssd / ss;
    else
        score = NaN;
    end

end

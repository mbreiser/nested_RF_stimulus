% RUN_PRE_BAR_FLASH_ANALYSIS  Run batch pipeline + generate population plots
%   for the pre-bar-flash dataset (summer 2025, 2-speed protocol).
%
%   Produces:
%     - batch_results_pre_bf.mat  (per-cell metrics)
%     - pre_bf_polar_on.png       (ON polar tuning: ctrl vs TTL)
%     - pre_bf_polar_off.png      (OFF polar tuning: ctrl vs TTL)
%     - pre_bf_dsi_comparison.png (DSI + FWHM by group)
%     - pre_bf_aspect_ratio.png   (T4 vs T5 x treatment)
%     - pre_bf_metrics_table.txt  (per-cell metrics text file)
%
%   Usage:
%     run('scripts/run_pre_bar_flash_analysis.m')

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root   = '/Users/reiserm/Documents/ttl_1DRF/pre-bar-flash';
preview_dir = fullfile(data_root, 'figure_previews');
results_dir = fullfile(data_root, 'population_results');

if ~isfolder(preview_dir), mkdir(preview_dir); end
if ~isfolder(results_dir), mkdir(results_dir); end

%% Run batch pipeline (or load saved results)
results_path = fullfile(results_dir, 'batch_results_pre_bf.mat');

run_fresh = true;  % Set to false to load saved results

if run_fresh || ~isfile(results_path)
    fprintf('Running batch pipeline...\n');
    opts = struct();
    opts.save_results = true;
    opts.save_dir = results_dir;
    results = batch_analyze_pre_bar_flash(data_root, opts);
else
    fprintf('Loading saved results from: %s\n', results_path);
    S = load(results_path, 'results');
    results = S.results;
    fprintf('Loaded %d cells.\n', numel(results));
end

if isempty(results)
    error('No results — cannot generate plots.');
end

fprintf('\n=== Generating population plots ===\n\n');

%% 1. Polar tuning plots — ON and OFF cells
for on_off_label = ["ON", "OFF"]
    if on_off_label == "ON"
        mask = [results.is_on];
    else
        mask = ~[results.is_on];
    end

    ctrl_mask = mask & ~[results.is_ttl];
    ttl_mask  = mask &  [results.is_ttl];

    if ~any(mask)
        fprintf('  No %s cells — skipping polar plot.\n', on_off_label);
        continue;
    end

    aligned_ctrl = {results(ctrl_mask).max_v_aligned};
    aligned_ttl  = {results(ttl_mask).max_v_aligned};

    polar_opts = struct();
    polar_opts.stat_method = 'mean_sem';
    polar_title = sprintf('Pre-BF %s Cells — PD-Aligned Polar Tuning', on_off_label);
    fig = plot_polar_population(aligned_ctrl, aligned_ttl, polar_title, polar_opts);

    fname = sprintf('pre_bf_polar_%s.png', lower(char(on_off_label)));
    exportgraphics(fig, fullfile(preview_dir, fname), 'Resolution', 150);
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

%% 2. DSI comparison (3-panel: DSI vector, DSI PD-ND, FWHM)
fig_dsi = plot_dsi_comparison(results);
exportgraphics(fig_dsi, fullfile(preview_dir, 'pre_bf_dsi_comparison.png'), 'Resolution', 150);
close(fig_dsi);
fprintf('  Saved: pre_bf_dsi_comparison.png\n');

%% 3. Aspect ratio comparison (T4 vs T5 x ctrl vs TTL)
fig_ar = plot_aspect_ratio_comparison(results);
exportgraphics(fig_ar, fullfile(preview_dir, 'pre_bf_aspect_ratio.png'), 'Resolution', 150);
close(fig_ar);
fprintf('  Saved: pre_bf_aspect_ratio.png\n');

%% 4. Save per-cell metrics table
table_path = fullfile(preview_dir, 'pre_bf_metrics_table.txt');
write_metrics_table(results, table_path);
fprintf('  Saved: pre_bf_metrics_table.txt\n');

fprintf('\n=== Done. All outputs in: %s ===\n', preview_dir);


%% ========================= Helper Functions ============================

function write_metrics_table(results, filepath)
% WRITE_METRICS_TABLE  Write per-cell metrics to a text file.

    fid = fopen(filepath, 'w');
    if fid == -1
        warning('Could not open %s for writing.', filepath);
        return;
    end

    fprintf(fid, '================================================================================\n');
    fprintf(fid, '  PRE-BAR-FLASH DATASET — PER-CELL METRICS\n');
    fprintf(fid, '  Generated: %s\n', datestr(now)); %#ok<TNOW1,DATST>
    fprintf(fid, '================================================================================\n\n');

    % Group-by-group output
    for grp = ["on_control", "on_ttl", "off_control", "off_ttl"]
        grp_mask = strcmp({results.group}, grp);
        grp_idx  = find(grp_mask);
        grp_label = upper(strrep(char(grp), '_', ' '));

        fprintf(fid, '--- %s (n=%d) ---\n', grp_label, numel(grp_idx));
        fprintf(fid, '%-4s  %-20s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %-10s\n', ...
            '#', 'Date_str', 'PD', 'DSI_v', 'DSI_pn', 'FWHM', 'AR', 'SymR', 'CV', 'Strain');
        fprintf(fid, '%s\n', repmat('-', 1, 100));

        for k = 1:numel(grp_idx)
            i = grp_idx(k);
            r = results(i);
            fprintf(fid, '%-4d  %-20s  %6.1f  %6.3f  %6.3f  %6.1f  %6.2f  %6.3f  %6.3f  %-10s\n', ...
                i, r.date_str, r.pd_direction, ...
                r.dsi_vector, r.dsi_pdnd, r.dir_tuning_fwhm, ...
                r.tuning_aspect_ratio, r.sym_ratio, r.dir_tuning_cv, ...
                r.strain);
        end
        fprintf(fid, '\n');
    end

    % Summary statistics
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '  SUMMARY STATISTICS BY GROUP (mean +/- std)\n');
    fprintf(fid, '================================================================================\n\n');

    fprintf(fid, '%-14s  %4s  %12s  %12s  %12s  %12s\n', ...
        'Group', 'n', 'DSI_vector', 'DSI_pdnd', 'FWHM', 'Aspect_Ratio');
    fprintf(fid, '%s\n', repmat('-', 1, 75));

    for grp = ["on_control", "on_ttl", "off_control", "off_ttl"]
        grp_mask = strcmp({results.group}, grp);
        if ~any(grp_mask), continue; end

        r_grp = results(grp_mask);
        n = numel(r_grp);

        dsi_v  = [r_grp.dsi_vector];
        dsi_pn = [r_grp.dsi_pdnd];
        fwhm   = [r_grp.dir_tuning_fwhm];
        ar     = [r_grp.tuning_aspect_ratio];

        fprintf(fid, '%-14s  %4d  %5.3f+/-%.3f  %5.3f+/-%.3f  %5.1f+/-%.1f  %5.2f+/-%.2f\n', ...
            upper(strrep(char(grp), '_', ' ')), n, ...
            mean(dsi_v), std(dsi_v), ...
            mean(dsi_pn), std(dsi_pn), ...
            mean(fwhm), std(fwhm), ...
            mean(ar, 'omitnan'), std(ar, 'omitnan'));
    end

    fprintf(fid, '\n');
    fclose(fid);

end

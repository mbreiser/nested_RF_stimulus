% GENERATE_FIGURE_PREVIEWS  Master script for all population figure PNGs.
%
%   Runs the updated batch pipeline on all 25 experiments, then generates
%   preview PNGs for each population figure panel. PNGs are stored in
%   <data_root>/figure_previews/ for use with the figure layout composer.
%
%   FIGURES GENERATED:
%     Population polar tuning:
%       pop_polar_ON_mean_sem.png
%       pop_polar_OFF_mean_sem.png
%
%     Population bar flash (unaligned, ND-to-PD ordering):
%       pop_pd_flash_ON_mean_sem.png
%       pop_pd_flash_OFF_mean_sem.png
%       pop_ortho_flash_ON_mean_sem.png
%       pop_ortho_flash_OFF_mean_sem.png
%
%     Population bar flash (M2-aligned, with n>=2 filter):
%       pop_pd_flash_ON_M2_aligned.png
%       pop_pd_flash_OFF_M2_aligned.png
%       pop_ortho_flash_ON_M2_aligned.png
%       pop_ortho_flash_OFF_M2_aligned.png
%
%     Population bar flash (M5-aligned, with n>=2 filter):
%       pop_pd_flash_ON_M5_aligned.png
%       pop_pd_flash_OFF_M5_aligned.png
%       pop_ortho_flash_ON_M5_aligned.png
%       pop_ortho_flash_OFF_M5_aligned.png
%
%     Summary figures:
%       pop_baseline_voltage.png
%       pop_rf_width.png
%       pop_dsi_comparison.png
%       pop_ortho_width_comparison.png
%       pop_timing_ON_M2.png, pop_timing_OFF_M2.png
%       pop_timing_ON_M5.png, pop_timing_OFF_M5.png
%
%   See also BATCH_ANALYZE_1DRF, GENERATE_EXAMPLE_CELL_PANELS

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root   = '/Users/reiserm/Documents/ttl_1DRF';
preview_dir = fullfile(data_root, 'figure_previews');
if ~isfolder(preview_dir), mkdir(preview_dir); end

export_res = 150;  % dpi for preview PNGs

%% Step 1: Run batch pipeline
fprintf('=== Running batch pipeline ===\n');
batch_opts = struct();
batch_opts.stat_method = 'mean_sem';
batch_opts.save_figs   = false;  % we'll save PNGs manually below
batch_opts.save_dir    = fullfile(data_root, 'population_results');

results = batch_analyze_1DRF(data_root, batch_opts);

% Save results for reuse
if ~isfolder(batch_opts.save_dir), mkdir(batch_opts.save_dir); end
save(fullfile(batch_opts.save_dir, 'batch_results.mat'), 'results');
fprintf('Results saved.\n');

%% Step 2: Generate population polar + unaligned flash PNGs
fprintf('\n=== Generating population plots ===\n');

for on_off_label = ["ON", "OFF"]
    if on_off_label == "ON"
        type_mask = [results.is_on];
    else
        type_mask = ~[results.is_on];
    end
    ctrl_mask = type_mask & ~[results.is_ttl];
    ttl_mask  = type_mask &  [results.is_ttl];
    label = lower(char(on_off_label));

    if ~any(type_mask)
        fprintf('  No %s cells — skipping.\n', on_off_label);
        continue;
    end

    % Polar tuning
    aligned_ctrl = {results(ctrl_mask).max_v_aligned};
    aligned_ttl  = {results(ttl_mask).max_v_aligned};
    polar_opts.stat_method = 'mean_sem';
    fig_polar = plot_polar_population(aligned_ctrl, aligned_ttl, ...
        sprintf('%s Cells — PD-Aligned Polar Tuning', on_off_label), polar_opts);
    exportgraphics(fig_polar, fullfile(preview_dir, ...
        sprintf('pop_polar_%s_mean_sem.png', label)), 'Resolution', export_res);
    close(fig_polar);

    % Unaligned PD-ND flash
    flash_opts.y_limits    = [-15 35];
    flash_opts.plot_type   = 'pd_nd';
    flash_opts.stat_method = 'mean_sem';
    fig_pd = plot_flash_1x11_population(...
        {results(ctrl_mask).pd_flash_bl}, {results(ttl_mask).pd_flash_bl}, ...
        sprintf('%s Cells — PD-ND Bar Flash', on_off_label), flash_opts);
    exportgraphics(fig_pd, fullfile(preview_dir, ...
        sprintf('pop_pd_flash_%s_mean_sem.png', label)), 'Resolution', export_res);
    close(fig_pd);

    % Unaligned orthogonal flash
    flash_opts.plot_type = 'orthogonal';
    fig_ortho = plot_flash_1x11_population(...
        {results(ctrl_mask).ortho_flash_bl}, {results(ttl_mask).ortho_flash_bl}, ...
        sprintf('%s Cells — Orthogonal Bar Flash', on_off_label), flash_opts);
    exportgraphics(fig_ortho, fullfile(preview_dir, ...
        sprintf('pop_ortho_flash_%s_mean_sem.png', label)), 'Resolution', export_res);
    close(fig_ortho);

    % --- Common aligned plot options ---
    aligned_flash_opts = struct();
    aligned_flash_opts.y_limits          = [-15 35];
    aligned_flash_opts.stat_method       = 'mean_sem';
    aligned_flash_opts.show_n_per_pos    = true;
    aligned_flash_opts.require_both_n2   = true;
    aligned_flash_opts.stim_onset_sample  = 5001;
    aligned_flash_opts.stim_offset_sample = 5801;
    aligned_flash_opts.resp_end_sample    = 6551;
    aligned_flash_opts.show_ordinal_ranks = true;
    aligned_flash_opts.col_labels        = arrayfun(@(x) sprintf('%+d', x), ...
        -5:5, 'UniformOutput', false);

    % === M2-aligned PD-ND flash ===
    m2_opts = aligned_flash_opts;
    m2_opts.plot_type     = 'pd_nd';
    m2_opts.col_labels{6} = '0 (M2 peak)';

    fig_pd_m2 = plot_flash_1x11_population(...
        {results(ctrl_mask).pd_flash_peak_aligned}, ...
        {results(ttl_mask).pd_flash_peak_aligned}, ...
        sprintf('%s Cells — M2-Aligned PD-ND', on_off_label), m2_opts);
    exportgraphics(fig_pd_m2, fullfile(preview_dir, ...
        sprintf('pop_pd_flash_%s_M2_aligned.png', label)), 'Resolution', export_res);
    close(fig_pd_m2);

    % === M2-aligned orthogonal flash ===
    m2_opts.plot_type = 'orthogonal';
    fig_ortho_m2 = plot_flash_1x11_population(...
        {results(ctrl_mask).ortho_flash_peak_aligned}, ...
        {results(ttl_mask).ortho_flash_peak_aligned}, ...
        sprintf('%s Cells — M2-Aligned Orthogonal', on_off_label), m2_opts);
    exportgraphics(fig_ortho_m2, fullfile(preview_dir, ...
        sprintf('pop_ortho_flash_%s_M2_aligned.png', label)), 'Resolution', export_res);
    close(fig_ortho_m2);

    % === M5-aligned PD-ND flash ===
    m5_opts = aligned_flash_opts;
    m5_opts.plot_type     = 'pd_nd';
    m5_opts.col_labels{6} = '0 (M5 centroid)';

    fig_pd_m5 = plot_flash_1x11_population(...
        {results(ctrl_mask).pd_flash_m5_aligned}, ...
        {results(ttl_mask).pd_flash_m5_aligned}, ...
        sprintf('%s Cells — M5-Aligned PD-ND', on_off_label), m5_opts);
    exportgraphics(fig_pd_m5, fullfile(preview_dir, ...
        sprintf('pop_pd_flash_%s_M5_aligned.png', label)), 'Resolution', export_res);
    close(fig_pd_m5);

    % === M5-aligned orthogonal flash ===
    m5_opts.plot_type = 'orthogonal';
    fig_ortho_m5 = plot_flash_1x11_population(...
        {results(ctrl_mask).ortho_flash_m5_aligned}, ...
        {results(ttl_mask).ortho_flash_m5_aligned}, ...
        sprintf('%s Cells — M5-Aligned Orthogonal', on_off_label), m5_opts);
    exportgraphics(fig_ortho_m5, fullfile(preview_dir, ...
        sprintf('pop_ortho_flash_%s_M5_aligned.png', label)), 'Resolution', export_res);
    close(fig_ortho_m5);
end

%% Step 3: Baseline voltage comparison
fprintf('\n=== Generating voltage comparison ===\n');
fig_v = plot_baseline_voltage_comparison(results);
exportgraphics(fig_v, fullfile(preview_dir, 'pop_baseline_voltage.png'), ...
    'Resolution', export_res);
close(fig_v);

%% Step 4: RF width comparison
fprintf('\n=== Generating RF width comparison ===\n');
fig_w = plot_rf_width_comparison(results);
exportgraphics(fig_w, fullfile(preview_dir, 'pop_rf_width.png'), ...
    'Resolution', export_res);
close(fig_w);

%% Step 5: DSI comparison
fprintf('\n=== Generating DSI comparison ===\n');
fig_dsi = plot_dsi_comparison(results);
exportgraphics(fig_dsi, fullfile(preview_dir, 'pop_dsi_comparison.png'), ...
    'Resolution', export_res);
close(fig_dsi);

%% Step 5b: Tuning aspect ratio
fprintf('\n=== Generating tuning aspect ratio ===\n');
fig_ar = plot_aspect_ratio_comparison(results);
exportgraphics(fig_ar, fullfile(preview_dir, 'pop_tuning_aspect_ratio.png'), ...
    'Resolution', export_res);
close(fig_ar);

%% Step 6: PD vs Orthogonal RF width
fprintf('\n=== Generating PD vs orthogonal RF width ===\n');
fig_ow = plot_ortho_width_comparison(results);
exportgraphics(fig_ow, fullfile(preview_dir, 'pop_ortho_width_comparison.png'), ...
    'Resolution', export_res);
close(fig_ow);

%% Step 7: Timing metrics by position (M2 and M5)
fprintf('\n=== Generating timing figures ===\n');
for cell_type = ["ON", "OFF"]
    % M2-based timing
    t_opts = struct();
    t_opts.cell_type = char(cell_type);
    t_opts.temporal_field = 'temporal_metrics';
    fig_t2 = plot_timing_by_position(results, t_opts);
    exportgraphics(fig_t2, fullfile(preview_dir, ...
        sprintf('pop_timing_%s_M2.png', lower(char(cell_type)))), ...
        'Resolution', export_res);
    close(fig_t2);

    % M5-based timing
    t_opts.temporal_field = 'temporal_metrics_m5';
    fig_t5 = plot_timing_by_position(results, t_opts);
    exportgraphics(fig_t5, fullfile(preview_dir, ...
        sprintf('pop_timing_%s_M5.png', lower(char(cell_type)))), ...
        'Resolution', export_res);
    close(fig_t5);
end

%% Summary
png_files = dir(fullfile(preview_dir, '*.png'));
fprintf('\n=== Done ===\n');
fprintf('Generated %d PNGs in: %s\n', numel(png_files), preview_dir);
for k = 1:numel(png_files)
    fprintf('  %s\n', png_files(k).name);
end

% RUN_T4_VS_T5_COMPARISON  Generate T4 vs T5 tuning comparison PNGs.
%
%   Loads batch_results.mat and generates 3 preview PNGs:
%     pop_polar_ctrl_t4_vs_t5.png  — Polar overlay: T4-ctrl vs T5-ctrl
%     pop_polar_ttl_t4_vs_t5.png   — Polar overlay: T4-TTL vs T5-TTL
%     pop_tuning_t4_vs_t5.png      — FWHM + DSI: T4 vs T5 box plots
%
%   Colors: T4 (ON) = blue, T5 (OFF) = orange
%
%   USAGE:
%     run('scripts/run_t4_vs_t5_comparison.m')
%
%   REQUIRES: batch_results.mat in population_results/
%
%   See also PLOT_POLAR_POPULATION, PLOT_TUNING_T4_VS_T5

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));

data_root   = '/Users/reiserm/Documents/ttl_1DRF';
preview_dir = fullfile(data_root, 'figure_previews');
results_file = fullfile(data_root, 'population_results', 'batch_results.mat');
if ~isfolder(preview_dir), mkdir(preview_dir); end

export_res = 150;

fprintf('Loading batch results...\n');
S = load(results_file, 'results');
results = S.results;
fprintf('Loaded %d cells.\n', numel(results));

%% Define groups
on_ctrl_mask  = [results.is_on] & ~[results.is_ttl];
off_ctrl_mask = ~[results.is_on] & ~[results.is_ttl];
on_ttl_mask   = [results.is_on] &  [results.is_ttl];
off_ttl_mask  = ~[results.is_on] &  [results.is_ttl];

% Colors: T4 = blue, T5 = orange
col_t4_line = [0.12 0.47 0.71];       % blue
col_t4_fill = [0.55 0.73 0.90];       % light blue
col_t5_line = [0.85 0.37 0.01];       % orange
col_t5_fill = [0.95 0.65 0.40];       % light orange

%% Polar plot 1: Control — T4 vs T5
fprintf('\n=== Polar: ctrl T4 vs T5 ===\n');

polar_opts = struct();
polar_opts.stat_method      = 'mean_sem';
polar_opts.group1_label     = 'T4 (ON)';
polar_opts.group2_label     = 'T5 (OFF)';
polar_opts.group1_line_color = col_t4_line;
polar_opts.group1_fill_color = col_t4_fill;
polar_opts.group2_line_color = col_t5_line;
polar_opts.group2_fill_color = col_t5_fill;

aligned_on_ctrl  = {results(on_ctrl_mask).max_v_aligned};
aligned_off_ctrl = {results(off_ctrl_mask).max_v_aligned};

fig_pc = plot_polar_population(aligned_on_ctrl, aligned_off_ctrl, ...
    'Control — T4 vs T5 Polar Tuning', polar_opts);
exportgraphics(fig_pc, fullfile(preview_dir, 'pop_polar_ctrl_t4_vs_t5.png'), ...
    'Resolution', export_res);
close(fig_pc);
fprintf('  Saved: pop_polar_ctrl_t4_vs_t5.png\n');

%% Polar plot 2: TTL — T4 vs T5
fprintf('\n=== Polar: TTL T4 vs T5 ===\n');

aligned_on_ttl  = {results(on_ttl_mask).max_v_aligned};
aligned_off_ttl = {results(off_ttl_mask).max_v_aligned};

fig_pt = plot_polar_population(aligned_on_ttl, aligned_off_ttl, ...
    'TTL — T4 vs T5 Polar Tuning', polar_opts);
exportgraphics(fig_pt, fullfile(preview_dir, 'pop_polar_ttl_t4_vs_t5.png'), ...
    'Resolution', export_res);
close(fig_pt);
fprintf('  Saved: pop_polar_ttl_t4_vs_t5.png\n');

%% Tuning metrics: T4 vs T5
fprintf('\n=== Tuning metrics: T4 vs T5 ===\n');

fig_tm = plot_tuning_t4_vs_t5(results);
exportgraphics(fig_tm, fullfile(preview_dir, 'pop_tuning_t4_vs_t5.png'), ...
    'Resolution', export_res);
close(fig_tm);
fprintf('  Saved: pop_tuning_t4_vs_t5.png\n');

%% Summary
fprintf('\n=== Done ===\n');
fprintf('Generated 3 T4 vs T5 comparison PNGs in: %s\n', preview_dir);

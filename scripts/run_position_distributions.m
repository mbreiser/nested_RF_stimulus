% RUN_POSITION_DISTRIBUTIONS  Generate per-position amplitude and histogram PNGs.
%
%   Loads batch_results.mat and generates 8 preview PNGs:
%     pop_amp_dist_{on,off}_{M2,M5}.png   — depol/hyperpol per position
%     pop_resp_hist_{on,off}_{M2,M5}.png  — full response histograms
%
%   USAGE:
%     run('scripts/run_position_distributions.m')
%
%   REQUIRES: batch_results.mat in population_results/
%
%   See also PLOT_AMPLITUDE_BY_POSITION, PLOT_RESPONSE_HISTOGRAMS_BY_POSITION

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

%% Generate amplitude distribution plots
fprintf('\n=== Per-position amplitude distributions ===\n');

configs = {
    'ON',  'pd_flash_peak_aligned',  'M2';
    'ON',  'pd_flash_m5_aligned',    'M5';
    'OFF', 'pd_flash_peak_aligned',  'M2';
    'OFF', 'pd_flash_m5_aligned',    'M5';
};

for c = 1:size(configs, 1)
    cell_type  = configs{c, 1};
    trace_field = configs{c, 2};
    align_label = configs{c, 3};

    amp_opts = struct();
    amp_opts.cell_type       = cell_type;
    amp_opts.trace_field     = trace_field;
    amp_opts.dep_window      = [5001, 6551];
    amp_opts.hyp_window      = [5001, Inf];
    amp_opts.dep_pctile      = 99.5;
    amp_opts.hyp_pctile      = 0.5;
    amp_opts.reject_threshold = 0.5;
    amp_opts.alignment_label = sprintf('0 (%s)', align_label);

    fig = plot_amplitude_by_position(results, amp_opts);
    fname = sprintf('pop_amp_dist_%s_%s.png', lower(cell_type), align_label);
    exportgraphics(fig, fullfile(preview_dir, fname), 'Resolution', export_res);
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

%% Generate response histogram plots
fprintf('\n=== Per-position response histograms ===\n');

for c = 1:size(configs, 1)
    cell_type  = configs{c, 1};
    trace_field = configs{c, 2};
    align_label = configs{c, 3};

    hist_opts = struct();
    hist_opts.cell_type    = cell_type;
    hist_opts.trace_field  = trace_field;
    hist_opts.stim_onset   = 5001;
    hist_opts.n_bins       = 20;
    hist_opts.reject_threshold = 0;

    fig = plot_response_histograms_by_position(results, hist_opts);
    fname = sprintf('pop_resp_hist_%s_%s.png', lower(cell_type), align_label);
    exportgraphics(fig, fullfile(preview_dir, fname), 'Resolution', export_res);
    close(fig);
    fprintf('  Saved: %s\n', fname);
end

%% Summary
fprintf('\n=== Done ===\n');
png_files = dir(fullfile(preview_dir, 'pop_amp_dist_*.png'));
png_files = [png_files; dir(fullfile(preview_dir, 'pop_resp_hist_*.png'))];
fprintf('Generated %d position distribution PNGs:\n', numel(png_files));
for k = 1:numel(png_files)
    fprintf('  %s\n', png_files(k).name);
end

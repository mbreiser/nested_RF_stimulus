% run_batch_m6_update.m — Re-run batch pipeline to add M6 fields
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
opts = struct();
opts.save_figs = true;   % triggers .mat save
opts.stat_method = 'mean_sem';
opts.skip_plots = true;  % just recompute results, skip figure generation

batch_analyze_1DRF(data_root, opts);
fprintf('\nDone. Checking M6 fields...\n');

S = load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
r1 = S.results(1);
fprintf('Fields present: pd_flash_m6_aligned=%d, ortho_flash_m6_aligned=%d\n', ...
    isfield(r1, 'pd_flash_m6_aligned'), isfield(r1, 'ortho_flash_m6_aligned'));
fprintf('centroid_m6=%.2f, centroid_m6_rounded=%d\n', r1.centroid_m6, r1.centroid_m6_rounded);
fprintf('ortho_centroid_m6=%.2f, ortho_centroid_m6_rounded=%d\n', r1.ortho_centroid_m6, r1.ortho_centroid_m6_rounded);
fprintf('m6_bump_width=%d, ortho_m6_bump_width=%d\n', r1.m6_bump_width, r1.ortho_m6_bump_width);

% Re-run batch pipeline to add fast flash fields
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
opts = struct();
opts.save_figs = true;  % needed to trigger .mat save (also generates plots)
opts.save_dir = fullfile(data_root, 'population_results');

fprintf('Starting batch pipeline re-run (no figures)...\n');
tic;
batch_analyze_1DRF(data_root, opts);
fprintf('\nBatch pipeline complete in %.1f s.\n', toc);

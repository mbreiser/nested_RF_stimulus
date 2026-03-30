%% Test analyze_single_experiment_mr
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

exp_folder = '/Users/reiserm/Documents/ttl_1DRF/2025_11_10_10_17';
opts = struct();
opts.visualize_everything = true;
opts.save_figs = true;  % just view, don't save PDFs

analyze_single_experiment_mr(exp_folder, opts);

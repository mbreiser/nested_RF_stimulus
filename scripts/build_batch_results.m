% BUILD_BATCH_RESULTS  One-time setup: build the two batch .mat files
% required by the manuscript figure pipeline.
%
% Edit DATA_ROOT and CIRCSTAT_PATH below to match your local setup, then
% run this script once. It produces:
%
%   <DATA_ROOT>/population_results/batch_results.mat
%       (late-batch dataset; built by batch_analyze_1DRF.m from the
%        per-experiment folders directly under DATA_ROOT)
%
%   <DATA_ROOT>/pre-bar-flash/population_results/batch_results_pre_bf.mat
%       (early-batch dataset; built by batch_analyze_pre_bar_flash.m
%        from the {control,ttl}/{ON,OFF}/<exp_folders> tree under
%        DATA_ROOT/pre-bar-flash/)
%
% After both .mat files exist, run scripts/generate_manuscript_fig_main.m
% and scripts/generate_manuscript_fig_supp.m to produce the figures.

DATA_ROOT     = '/Users/reiserm/Documents/ttl_1DRF';                                            % <-- edit for your setup
CIRCSTAT_PATH = '/Users/reiserm/Library/CloudStorage/Dropbox-HHMI/Michael Reiser/Matlab_work/CircStat2012a';  % <-- edit for your setup

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));
addpath(CIRCSTAT_PATH);

%% Late-batch (full 1DRF protocol)
fprintf('\n========== Building batch_results.mat (late dataset) ==========\n');
late_opts = struct('save_figs', false);
results = batch_analyze_1DRF(DATA_ROOT, late_opts); %#ok<NASGU>
late_dir = fullfile(DATA_ROOT, 'population_results');
if ~isfolder(late_dir), mkdir(late_dir); end
late_out = fullfile(late_dir, 'batch_results.mat');
save(late_out, 'results', '-v7.3');
fprintf('Wrote %s\n', late_out);
clear results;

%% Early-batch (pre-bar-flash, summer 2025 2-speed protocol)
fprintf('\n========== Building batch_results_pre_bf.mat (early dataset) ==========\n');
early_root = fullfile(DATA_ROOT, 'pre-bar-flash');
early_opts = struct('save_results', true, ...
    'save_dir', fullfile(early_root, 'population_results'));
batch_analyze_pre_bar_flash(early_root, early_opts);
fprintf('Wrote %s\n', fullfile(early_opts.save_dir, 'batch_results_pre_bf.mat'));

fprintf('\n========== Done ==========\n');
fprintf('You can now run:\n');
fprintf('  run(''scripts/generate_manuscript_fig_main.m'')\n');
fprintf('  run(''scripts/generate_manuscript_fig_supp.m'')\n');

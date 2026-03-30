% TEST_SINGLE_CELL_MR  Quick test of Enhancement G on one cell.
%
%   Tests analyze_single_experiment_mr with the default test experiment
%   (2025_11_10_10_17) and verifies that cell_info is returned correctly.

addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

exp_folder = '/Users/reiserm/Documents/ttl_1DRF/2025_11_10_10_17';

opts = struct();
opts.save_figs = true;
opts.save_dir  = fullfile(exp_folder, 'analysis_output_mr');
opts.visualize_everything = true;

fprintf('=== Testing analyze_single_experiment_mr with Enhancement G ===\n');
fprintf('Experiment: %s\n\n', exp_folder);

cell_info = analyze_single_experiment_mr(exp_folder, opts);

fprintf('\n=== Returned cell_info struct ===\n');
disp(cell_info);

fprintf('\n=== Field verification ===\n');
required_fields = {'exp_folder', 'date_str', 'strain', 'frame', 'on_off', ...
    'pd_direction', 'pd_orientation', 'centroid_idx', 'centroid_deg', 'peak_amplitudes'};
all_ok = true;
for i = 1:numel(required_fields)
    if isfield(cell_info, required_fields{i})
        fprintf('  %-18s  OK\n', required_fields{i});
    else
        fprintf('  %-18s  MISSING!\n', required_fields{i});
        all_ok = false;
    end
end

if all_ok
    fprintf('\nAll fields present. Test PASSED.\n');
else
    fprintf('\nSome fields missing. Test FAILED.\n');
end

close all;
fprintf('\nDone.\n');

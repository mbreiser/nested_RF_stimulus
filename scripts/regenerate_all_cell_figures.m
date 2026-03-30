% REGENERATE_ALL_CELL_FIGURES  Re-run analyze_single_experiment_mr on all 25 cells.
%
%   Regenerates all individual cell figures (5 PDFs per cell) with the latest
%   code. Each cell's figures are saved to <exp_folder>/analysis_output_mr/.
%
%   See also ANALYZE_SINGLE_EXPERIMENT_MR

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';

%% All 25 cells
all_cells = {
    '2025_10_29_09_21',   'on_control'
    '2025_11_10_10_17',   'on_control'
    '2025_11_11_12_04',   'on_control'
    '2025_11_13_11_03',   'on_control'
    '2025_11_13_12_21',   'on_control'
    '2025_10_23_10_31',   'on_ttl'
    '2025_10_23_13_53',   'on_ttl'
    '2025_10_28_13_54',   'on_ttl'
    '2025_11_03_11_54',   'on_ttl'
    '2025_11_17_09_44',   'on_ttl'
    '2025_10_24_12_49',   'off_control'
    '2025_10_24_17_50',   'off_control'
    '2025_10_29_11_06',   'off_control'
    '2025_10_31_10_16',   'off_control'
    '2025_11_04_10_01',   'off_control'
    '2025_11_04_11_56',   'off_control'
    '2025_11_10_11_43',   'off_control'
    '2025_10_23_15_03',   'off_ttl'
    '2025_10_28_17_37',   'off_ttl'
    '2025_10_30_10_47',   'off_ttl'
    '2025_11_03_15_50',   'off_ttl'
    '2025_11_05_10_10',   'off_ttl'
    '2025_11_05_15_54',   'off_ttl'
    '2025_11_05_18_33',   'off_ttl'
    '2025_11_17_14_24',   'off_ttl'
};

n_total = size(all_cells, 1);
t_start = tic;

fprintf('=== Regenerating figures for %d cells ===\n\n', n_total);

errors = {};
for i = 1:n_total
    folder_name = all_cells{i, 1};
    group_label = all_cells{i, 2};
    exp_folder  = fullfile(data_root, folder_name);

    fprintf('[%d/%d] %s (%s)... ', i, n_total, folder_name, group_label);

    opts = struct();
    opts.save_figs = true;
    opts.save_dir  = fullfile(exp_folder, 'analysis_output_mr');
    opts.visualize_everything = true;

    try
        analyze_single_experiment_mr(exp_folder, opts);
        fprintf('OK (%.0f s)\n', toc(t_start) / i * 1);
    catch ME
        fprintf('ERROR: %s\n', ME.message);
        errors{end+1} = sprintf('%s: %s', folder_name, ME.message); %#ok<SAGROW>
    end

    close all;
end

elapsed = toc(t_start);
fprintf('\n=== All %d cells processed in %.0f s (%.1f s/cell) ===\n', ...
    n_total, elapsed, elapsed / n_total);

if ~isempty(errors)
    fprintf('\nErrors:\n');
    for k = 1:numel(errors)
        fprintf('  %s\n', errors{k});
    end
end

fprintf('\nDone.\n');

% RUN_OFF_CELLS_MR  Batch run analyze_single_experiment_mr on 15 OFF cells.
%
%   Loops over all 15 OFF cells (7 control + 8 TTL), runs the enhanced
%   analysis with visualize_everything=true, collects cell_info structs,
%   and saves a centering summary.
%
%   See also ANALYZE_SINGLE_EXPERIMENT_MR, RUN_ON_CELLS_MR

%% Setup paths
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
results_dir = fullfile(data_root, 'population_results');
if ~isfolder(results_dir), mkdir(results_dir); end

%% Define 15 OFF cells (7 control + 8 TTL)
% Extracted from batch_summary.txt
off_cells = {
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

n_cells = size(off_cells, 1);

%% Run analysis on each cell
cell_infos = cell(n_cells, 1);
t_start = tic;

for i = 1:n_cells
    folder_name = off_cells{i, 1};
    group_label = off_cells{i, 2};
    exp_folder  = fullfile(data_root, folder_name);

    fprintf('\n========================================================\n');
    fprintf('  [%d/%d] %s (%s)\n', i, n_cells, folder_name, group_label);
    fprintf('========================================================\n');

    opts = struct();
    opts.save_figs = true;
    opts.save_dir  = fullfile(exp_folder, 'analysis_output_mr');
    opts.visualize_everything = true;

    try
        cell_infos{i} = analyze_single_experiment_mr(exp_folder, opts);
    catch ME
        fprintf('  *** ERROR: %s\n', ME.message);
        cell_infos{i} = struct('exp_folder', exp_folder, ...
            'date_str', folder_name, 'error', ME.message);
    end

    close all;
end

elapsed = toc(t_start);
fprintf('\n\nAll %d OFF cells processed in %.1f s\n', n_cells, elapsed);

%% Print centering summary table
header = sprintf('\n%-28s  %-12s  %7s  %8s  %10s  %10s  %8s  %6s', ...
    'Folder', 'Group', 'PD_dir', 'Peak@', 'Centroid', 'Offset', 'Bump', 'Delta');
separator = repmat('-', 1, numel(header));

fprintf('\n%s\n%s\n', header, separator);

summary_lines = cell(n_cells + 3, 1);
summary_lines{1} = 'OFF Cells Centering Summary';
summary_lines{2} = header;
summary_lines{3} = separator;

for i = 1:n_cells
    ci = cell_infos{i};
    if isfield(ci, 'error')
        line = sprintf('%-28s  %-12s  ERROR: %s', ...
            off_cells{i, 1}, off_cells{i, 2}, ci.error);
    else
        line = sprintf('%-28s  %-12s  %7.1f  %5d    %9.2f  %+9.1f°  %4d-%-2d  %5.2f', ...
            off_cells{i, 1}, off_cells{i, 2}, ...
            ci.pd_direction, ci.peak_pos, ci.centroid_idx, ci.centroid_deg, ...
            ci.bump_range(1), ci.bump_range(2), ci.centroid_peak_delta);
    end
    fprintf('%s\n', line);
    summary_lines{i + 3} = line;
end

%% Save summary
txt_path = fullfile(results_dir, 'off_cells_centering_summary.txt');
fid = fopen(txt_path, 'w');
for i = 1:numel(summary_lines)
    fprintf(fid, '%s\n', summary_lines{i});
end
fclose(fid);
fprintf('\nSummary saved to: %s\n', txt_path);

mat_path = fullfile(results_dir, 'off_cells_centering_summary.mat');
save(mat_path, 'cell_infos', 'off_cells');
fprintf('Data saved to:    %s\n', mat_path);

fprintf('\nDone.\n');

% RUN_ON_CELLS_MR  Batch run analyze_single_experiment_mr on 10 ON cells.
%
%   Loops over all 10 ON cells (5 control + 5 TTL), runs the enhanced
%   analysis with visualize_everything=true, collects cell_info structs,
%   and prints/saves a centering summary table.
%
%   OUTPUT:
%     - 5 PDFs per cell in <exp_folder>/analysis_output_mr/
%     - Centering summary table printed to console
%     - Summary saved to population_results/on_cells_centering_summary.txt
%     - Summary saved to population_results/on_cells_centering_summary.mat
%
%   See also ANALYZE_SINGLE_EXPERIMENT_MR

%% Setup paths
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
results_dir = fullfile(data_root, 'population_results');
if ~isfolder(results_dir), mkdir(results_dir); end

%% Define 10 ON cells (5 control + 5 TTL)
% Extracted from batch_summary.txt
on_cells = {
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
};

n_cells = size(on_cells, 1);

%% Run analysis on each cell
cell_infos = cell(n_cells, 1);
t_start = tic;

for i = 1:n_cells
    folder_name = on_cells{i, 1};
    group_label = on_cells{i, 2};
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

    % Close figures between cells to manage memory
    close all;
end

elapsed = toc(t_start);
fprintf('\n\nAll %d ON cells processed in %.1f s\n', n_cells, elapsed);

%% Print centering summary table
header = sprintf('\n%-28s  %-12s  %7s  %8s  %10s  %10s  %8s  %6s', ...
    'Folder', 'Group', 'PD_dir', 'Peak@', 'Centroid', 'Offset', 'Bump', 'Delta');
separator = repmat('-', 1, numel(header));

fprintf('\n%s\n%s\n', header, separator);

summary_lines = cell(n_cells + 3, 1);
summary_lines{1} = 'ON Cells Centering Summary';
summary_lines{2} = header;
summary_lines{3} = separator;

for i = 1:n_cells
    ci = cell_infos{i};
    if isfield(ci, 'error')
        line = sprintf('%-28s  %-12s  ERROR: %s', ...
            on_cells{i, 1}, on_cells{i, 2}, ci.error);
    else
        line = sprintf('%-28s  %-12s  %7.1f  %5d    %9.2f  %+9.1f°  %4d-%-2d  %5.2f', ...
            on_cells{i, 1}, on_cells{i, 2}, ...
            ci.pd_direction, ci.peak_pos, ci.centroid_idx, ci.centroid_deg, ...
            ci.bump_range(1), ci.bump_range(2), ci.centroid_peak_delta);
    end
    fprintf('%s\n', line);
    summary_lines{i + 3} = line;
end

%% Save summary
% Text file
txt_path = fullfile(results_dir, 'on_cells_centering_summary.txt');
fid = fopen(txt_path, 'w');
for i = 1:numel(summary_lines)
    fprintf(fid, '%s\n', summary_lines{i});
end
fclose(fid);
fprintf('\nSummary saved to: %s\n', txt_path);

% MAT file
mat_path = fullfile(results_dir, 'on_cells_centering_summary.mat');
save(mat_path, 'cell_infos', 'on_cells');
fprintf('Data saved to:    %s\n', mat_path);

fprintf('\nDone.\n');

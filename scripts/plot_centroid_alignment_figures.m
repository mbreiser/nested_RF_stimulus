% PLOT_CENTROID_ALIGNMENT_FIGURES  Generate centroid-aligned population figures.
%
%   Loads saved batch results and centering data for ON and OFF cells,
%   applies nearest-neighbor centroid shifts onto a 15-column canvas, and
%   generates per-cell-type:
%
%     Figure Type 1 (individual cell rows):
%       - Control cells: original and centroid-shifted
%       - TTL cells: original and centroid-shifted
%
%     Figure Type 2 (population average):
%       - Original alignment: mean ± SD
%       - Centroid-shifted: mean ± SD with per-position n-counts
%
%   OUTPUT:
%     12 PDFs in population_results/centroid_aligned/
%     (6 for ON cells + 6 for OFF cells)
%
%   See also SHIFT_TRACES_TO_CANVAS, PLOT_FLASH_1X11_INDIVIDUAL_ROWS,
%            PLOT_FLASH_1X11_POPULATION

%% Setup paths
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));

data_root   = '/Users/reiserm/Documents/ttl_1DRF';
results_dir = fullfile(data_root, 'population_results');
save_dir    = fullfile(results_dir, 'centroid_aligned');
if ~isfolder(save_dir), mkdir(save_dir); end

%% Load batch results (shared across ON/OFF)
S1 = load(fullfile(results_dir, 'batch_results.mat'), 'results');
results = S1.results;
batch_folders = {results.folder};

%% Degree label templates
deg_labels_11 = make_degree_labels(11, 6);
deg_labels_15 = make_degree_labels(15, 8);

export_opts = {'ContentType', 'image', 'Resolution', 300};

%% Process ON and OFF cell types
cell_types = {'on', 'off'};

for ct = 1:numel(cell_types)
    cell_type = cell_types{ct};
    upper_type = upper(cell_type);

    fprintf('\n============================================================\n');
    fprintf('  Processing %s cells\n', upper_type);
    fprintf('============================================================\n');

    %% Load centering data for this cell type
    summary_file = fullfile(results_dir, ...
        sprintf('%s_cells_centering_summary.mat', cell_type));
    if ~isfile(summary_file)
        fprintf('  WARNING: %s not found — skipping %s cells\n', ...
            summary_file, upper_type);
        continue;
    end

    S2 = load(summary_file);
    cell_infos = S2.cell_infos;
    % The variable name is 'on_cells' or 'off_cells' in the .mat
    if isfield(S2, 'on_cells')
        cell_list = S2.on_cells;
    elseif isfield(S2, 'off_cells')
        cell_list = S2.off_cells;
    else
        error('Cannot find cell list variable in %s', summary_file);
    end

    n_cells = size(cell_list, 1);

    %% Match cells and apply shifts
    folder_list    = cell(n_cells, 1);
    group_list     = cell(n_cells, 1);
    centroid_list  = zeros(n_cells, 1);
    shift_list     = zeros(n_cells, 1);
    orig_traces    = cell(n_cells, 1);
    shifted_traces = cell(n_cells, 1);
    valid          = true(n_cells, 1);

    for i = 1:n_cells
        folder = cell_list{i, 1};
        group  = cell_list{i, 2};

        batch_idx = find(strcmp(batch_folders, folder));
        if isempty(batch_idx)
            fprintf('  WARNING: %s not found in batch_results\n', folder);
            valid(i) = false;
            continue;
        end

        ci = cell_infos{i};
        if isfield(ci, 'error')
            fprintf('  WARNING: %s had error\n', folder);
            valid(i) = false;
            continue;
        end

        folder_list{i}   = folder;
        group_list{i}    = group;
        centroid_list(i)  = ci.centroid_idx;
        orig_traces{i}    = results(batch_idx).pd_flash_bl;

        [shifted_traces{i}, shift_list(i)] = ...
            shift_traces_to_canvas(orig_traces{i}, ci.centroid_idx);

        fprintf('  [%2d] %-25s  %-14s  c=%.2f  s=%+d\n', ...
            i, folder, group, ci.centroid_idx, shift_list(i));
    end

    % Remove invalid
    folder_list    = folder_list(valid);
    group_list     = group_list(valid);
    centroid_list  = centroid_list(valid);
    shift_list     = shift_list(valid);
    orig_traces    = orig_traces(valid);
    shifted_traces = shifted_traces(valid);

    %% Split into control and TTL
    ctrl_mask = contains(group_list, 'control');
    ttl_mask  = contains(group_list, 'ttl');

    ctrl_orig    = orig_traces(ctrl_mask);
    ctrl_shifted = shifted_traces(ctrl_mask);
    ctrl_folders = folder_list(ctrl_mask);
    ctrl_centroids = centroid_list(ctrl_mask);
    ctrl_shifts  = shift_list(ctrl_mask);

    ttl_orig     = orig_traces(ttl_mask);
    ttl_shifted  = shifted_traces(ttl_mask);
    ttl_folders  = folder_list(ttl_mask);
    ttl_centroids = centroid_list(ttl_mask);
    ttl_shifts   = shift_list(ttl_mask);

    fprintf('\n  %s Control: %d cells, %s TTL: %d cells\n', ...
        upper_type, sum(ctrl_mask), upper_type, sum(ttl_mask));

    %% Build row labels
    ctrl_labels = build_row_labels(ctrl_folders, ctrl_centroids, ctrl_shifts);
    ttl_labels  = build_row_labels(ttl_folders, ttl_centroids, ttl_shifts);

    %% ==================================================================
    %  FIGURE TYPE 1: Individual Cell Rows
    %  ==================================================================
    fprintf('\n  --- Individual cell row figures ---\n');

    opts_orig = struct('col_labels', {deg_labels_11}, 'center_col', 6);
    opts_shift = struct('col_labels', {deg_labels_15}, 'center_col', 8);

    % Control original
    fname = sprintf('%s_ctrl_individual_original.pdf', cell_type);
    fig = plot_flash_1x11_individual_rows(ctrl_orig, ctrl_labels, ...
        sprintf('%s Control — Original Alignment (PD-ND)', upper_type), opts_orig);
    exportgraphics(fig, fullfile(save_dir, fname), export_opts{:});
    fprintf('    Saved: %s\n', fname);
    close(fig);

    % Control shifted
    fname = sprintf('%s_ctrl_individual_shifted.pdf', cell_type);
    fig = plot_flash_1x11_individual_rows(ctrl_shifted, ctrl_labels, ...
        sprintf('%s Control — Centroid-Shifted Alignment', upper_type), opts_shift);
    exportgraphics(fig, fullfile(save_dir, fname), export_opts{:});
    fprintf('    Saved: %s\n', fname);
    close(fig);

    % TTL original
    fname = sprintf('%s_ttl_individual_original.pdf', cell_type);
    fig = plot_flash_1x11_individual_rows(ttl_orig, ttl_labels, ...
        sprintf('%s TTL — Original Alignment (PD-ND)', upper_type), opts_orig);
    exportgraphics(fig, fullfile(save_dir, fname), export_opts{:});
    fprintf('    Saved: %s\n', fname);
    close(fig);

    % TTL shifted
    fname = sprintf('%s_ttl_individual_shifted.pdf', cell_type);
    fig = plot_flash_1x11_individual_rows(ttl_shifted, ttl_labels, ...
        sprintf('%s TTL — Centroid-Shifted Alignment', upper_type), opts_shift);
    exportgraphics(fig, fullfile(save_dir, fname), export_opts{:});
    fprintf('    Saved: %s\n', fname);
    close(fig);

    %% ==================================================================
    %  FIGURE TYPE 2: Population Average — mean ± SD
    %  ==================================================================
    fprintf('\n  --- Population comparison figures ---\n');

    pop_base = struct('y_limits', [-15 35], 'plot_type', 'pd_nd', ...
        'stat_method', 'mean_sd', 'show_n_per_pos', true);

    % Original
    pop_opts = pop_base;
    pop_opts.col_labels = deg_labels_11;
    fname = sprintf('%s_population_original_mean_sd.pdf', cell_type);
    fig = plot_flash_1x11_population(ctrl_orig, ttl_orig, ...
        sprintf('%s Cells — Original PD-ND (mean \\pm SD)', upper_type), pop_opts);
    exportgraphics(fig, fullfile(save_dir, fname), export_opts{:});
    fprintf('    Saved: %s\n', fname);
    close(fig);

    % Centroid-shifted
    pop_opts = pop_base;
    pop_opts.col_labels = deg_labels_15;
    fname = sprintf('%s_population_centroid_aligned_mean_sd.pdf', cell_type);
    fig = plot_flash_1x11_population(ctrl_shifted, ttl_shifted, ...
        sprintf('%s Cells — Centroid-Aligned PD-ND (mean \\pm SD)', upper_type), pop_opts);
    exportgraphics(fig, fullfile(save_dir, fname), export_opts{:});
    fprintf('    Saved: %s\n', fname);
    close(fig);

end  % cell_types loop

%% Summary
fprintf('\n=== All figures saved to: %s ===\n', save_dir);
fprintf('Done.\n');


%% ========================= Local Functions ============================

function labels = make_degree_labels(n_cols, center_col)
% MAKE_DEGREE_LABELS  Generate degree label strings for column headers.

    labels = cell(1, n_cols);
    for p = 1:n_cols
        deg_val = (p - center_col) * 2.5;
        if p == 1
            labels{p} = sprintf('%.1f° (Lead)', deg_val);
        elseif p == center_col
            labels{p} = '0° (Center)';
        elseif p == n_cols
            labels{p} = sprintf('+%.1f° (Trail)', deg_val);
        elseif deg_val > 0
            labels{p} = sprintf('+%.1f°', deg_val);
        else
            labels{p} = sprintf('%.1f°', deg_val);
        end
    end

end


function labels = build_row_labels(folders, centroids, shifts)
% BUILD_ROW_LABELS  Format row labels with folder name, centroid, shift.

    labels = cell(numel(folders), 1);
    for i = 1:numel(folders)
        labels{i} = sprintf('%s (c=%.1f, s=%+d)', ...
            strrep(folders{i}, '_', '-'), centroids(i), shifts(i));
    end

end

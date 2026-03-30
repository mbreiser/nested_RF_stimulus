% GENERATE_EXAMPLE_CELL_PANELS  Generate sweep + polar PNGs for example cells.
%
%   Produces preview PNGs of the bar sweep polar timeseries figure for one
%   representative cell from each of the 4 groups (ON ctrl, ON TTL, OFF ctrl,
%   OFF TTL). These PNGs are intended for use in the figure layout composer.
%
%   USAGE:
%     1. Fill in the 4 experiment folder names below
%     2. Run this script
%     3. PNGs saved to <data_root>/figure_previews/
%
%   See also ANALYZE_SINGLE_EXPERIMENT_MR, GENERATE_FIGURE_PREVIEWS

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root   = '/Users/reiserm/Documents/ttl_1DRF';
preview_dir = fullfile(data_root, 'figure_previews');
if ~isfolder(preview_dir), mkdir(preview_dir); end

%% === FILL IN THESE FOLDER NAMES ===
% Select one representative experiment from each group.
% These should be cells with clear direction selectivity.
example_cells = struct();
example_cells.on_control  = '';  % e.g. '2025_11_10_10_17'
example_cells.on_ttl      = '';  % e.g. '2025_11_12_09_45'
example_cells.off_control = '';  % e.g. '2025_11_14_11_30'
example_cells.off_ttl     = '';  % e.g. '2025_11_15_10_00'

%% Generate panels
group_fields = fieldnames(example_cells);
panel_labels = {'a', 'b', 'c', 'd'};

for g = 1:numel(group_fields)
    grp = group_fields{g};
    folder_name = example_cells.(grp);

    if isempty(folder_name)
        fprintf('Skipping %s — no folder specified.\n', grp);
        continue;
    end

    exp_folder = fullfile(data_root, folder_name);
    if ~isfolder(exp_folder)
        fprintf('WARNING: folder not found: %s\n', exp_folder);
        continue;
    end

    fprintf('Processing %s: %s\n', grp, folder_name);

    % Run single-experiment analysis
    opts = struct();
    opts.save_figs = false;
    opts.visualize_everything = true;
    cell_info = analyze_single_experiment_mr(exp_folder, opts);

    % Find the polar/sweep figure (Figure 1)
    fig_sweep = findobj('Type', 'figure', '-regexp', 'Name', 'Slow Bar Sweep');
    if ~isempty(fig_sweep)
        fig_sweep = fig_sweep(1);
        fname = sprintf('fig_example_%s_%s_sweep_polar_%s.png', ...
            panel_labels{g}, strrep(grp, '_', ''), folder_name);
        exportgraphics(fig_sweep, fullfile(preview_dir, fname), 'Resolution', 150);
        fprintf('  Saved: %s\n', fname);
    end

    close all;
end

fprintf('\nDone. PNGs saved to: %s\n', preview_dir);

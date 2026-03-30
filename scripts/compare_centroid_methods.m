% COMPARE_CENTROID_METHODS  Compare old vs new RF centroid on all 25 cells.
%
%   Loads the old centering summaries (from previous runs), then re-runs
%   the updated compute_rf_centroid on a single test cell to verify, and
%   finally runs all 25 cells to produce a comparison table.
%
%   The comparison table shows:
%     - Old centroid (all-position center-of-mass, max over full trace)
%     - New centroid (bump-based, 99.5th pctile in response window)
%     - Peak position, bump width, centroid-peak delta
%     - Whether the shift amount changed

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
results_dir = fullfile(data_root, 'population_results');

%% Load old centering summaries
old_on  = load(fullfile(results_dir, 'on_cells_centering_summary.mat'));
old_off = load(fullfile(results_dir, 'off_cells_centering_summary.mat'));

%% Quick single-cell test first
fprintf('=== Single cell test (2025_11_10_10_17) ===\n');
test_folder = fullfile(data_root, '2025_11_10_10_17');
opts = struct();
opts.save_figs = false;
opts.visualize_everything = false;
ci_test = analyze_single_experiment_mr(test_folder, opts);
close all;

fprintf('\n  New centroid: %.2f (%.1f°)\n', ci_test.centroid_idx, ci_test.centroid_deg);
fprintf('  Peak position: %d\n', ci_test.peak_pos);
fprintf('  Bump: %d-%d (%d wide)\n', ci_test.bump_range(1), ci_test.bump_range(2), ci_test.bump_width);
fprintf('  Centroid-peak delta: %.2f\n', ci_test.centroid_peak_delta);

%% Run all 25 cells with new algorithm (no figures)
all_cells = [old_on.on_cells; old_off.off_cells];
old_infos = [old_on.cell_infos; old_off.cell_infos];
n_total = size(all_cells, 1);

new_infos = cell(n_total, 1);
fprintf('\n=== Running new algorithm on all %d cells ===\n', n_total);

for i = 1:n_total
    folder_name = all_cells{i, 1};
    exp_folder = fullfile(data_root, folder_name);

    opts = struct();
    opts.save_figs = false;
    opts.visualize_everything = false;

    try
        new_infos{i} = analyze_single_experiment_mr(exp_folder, opts);
    catch ME
        fprintf('  ERROR on %s: %s\n', folder_name, ME.message);
        new_infos{i} = struct('centroid_idx', NaN, 'centroid_deg', NaN, ...
            'peak_pos', NaN, 'bump_width', NaN, 'bump_range', [NaN NaN], ...
            'centroid_peak_delta', NaN);
    end
    close all;
end

%% Print comparison table
fprintf('\n\n=== CENTROID COMPARISON: Old (all-pos, max) vs New (bump, q99.5) ===\n');
fprintf('%-28s  %-12s  %8s  %8s  %6s  %5s  %8s  %6s  %8s\n', ...
    'Folder', 'Group', 'Old_ctr', 'New_ctr', 'Delta', 'Peak', 'Bump', 'Width', 'Ctr-Pk');
fprintf('%s\n', repmat('-', 1, 110));

shift_changed = 0;
for i = 1:n_total
    oi = old_infos{i};
    ni = new_infos{i};

    old_ctr = oi.centroid_idx;
    new_ctr = ni.centroid_idx;
    old_shift = 6 - round(old_ctr);
    new_shift = 6 - round(new_ctr);
    changed = old_shift ~= new_shift;
    if changed, shift_changed = shift_changed + 1; end

    marker = '';
    if changed, marker = ' ***'; end

    fprintf('%-28s  %-12s  %8.2f  %8.2f  %+5.2f  %5d  %4d-%-2d  %5d  %7.2f%s\n', ...
        all_cells{i, 1}, all_cells{i, 2}, ...
        old_ctr, new_ctr, new_ctr - old_ctr, ...
        ni.peak_pos, ni.bump_range(1), ni.bump_range(2), ni.bump_width, ...
        ni.centroid_peak_delta, marker);
end

fprintf('\n%d of %d cells had shift amount change (marked ***)\n', shift_changed, n_total);

%% Summary statistics by group
groups = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
fprintf('\n=== Bump width by group ===\n');
for g = 1:numel(groups)
    mask = strcmp(all_cells(:, 2), groups{g});
    widths = cellfun(@(x) x.bump_width, new_infos(mask));
    fprintf('  %-12s: mean=%.1f SD=%.1f range=[%d, %d] (n=%d)\n', ...
        groups{g}, mean(widths), std(widths), min(widths), max(widths), sum(mask));
end

fprintf('\nDone.\n');

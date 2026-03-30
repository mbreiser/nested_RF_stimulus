% RUN_COMBINED_BATCH_COMPARISON  Compare early (pre-bar-flash) and late (main)
%   datasets, plus combined. Generates DSI, aspect ratio, and polar tuning
%   plots for 3 groupings: all combined, early only, late only.
%
%   Produces 12 PNGs in the combined figure_previews directory:
%     Polar plots (6):  comb_polar_{on,off}_{combined,early,late}.png
%     DSI (3):          comb_dsi_{combined,early,late}.png
%     Aspect ratio (3): comb_aspect_ratio_{combined,early,late}.png
%
%   Usage:
%     run('scripts/run_combined_batch_comparison.m')

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root   = '/Users/reiserm/Documents/ttl_1DRF';
preview_dir = fullfile(data_root, 'figure_previews');
if ~isfolder(preview_dir), mkdir(preview_dir); end

%% Load both batch results
fprintf('Loading batch results...\n');

S_late = load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
late = S_late.results;
fprintf('  Late batch (main): %d cells\n', numel(late));

S_early = load(fullfile(data_root, 'pre-bar-flash', 'population_results', ...
    'batch_results_pre_bf.mat'), 'results');
early = S_early.results;
fprintf('  Early batch (pre-BF): %d cells\n', numel(early));

%% Harmonize struct fields for concatenation
%   Both have the key shared fields. We create a minimal common struct
%   with just the fields needed for DSI, aspect ratio, and polar plots.

combined = harmonize_and_merge(early, late);
early_h  = harmonize_results(early, 'early');
late_h   = harmonize_results(late, 'late');

fprintf('  Combined: %d cells\n\n', numel(combined));

%% Define the 3 dataset groupings
datasets = struct( ...
    'label',   {'Combined (n=%d)', 'Early batch (n=%d)', 'Late batch (n=%d)'}, ...
    'tag',     {'combined', 'early', 'late'}, ...
    'results', {combined, early_h, late_h} ...
);

%% Generate all plots
fprintf('=== Generating comparison plots ===\n\n');

for d = 1:numel(datasets)
    ds = datasets(d);
    n_cells = numel(ds.results);
    ds_label = sprintf(ds.label, n_cells);
    fprintf('--- %s ---\n', ds_label);

    % --- Polar tuning plots (ON and OFF) ---
    for on_off = ["ON", "OFF"]
        if on_off == "ON"
            mask = [ds.results.is_on];
        else
            mask = ~[ds.results.is_on];
        end

        ctrl_mask = mask & ~[ds.results.is_ttl];
        ttl_mask  = mask &  [ds.results.is_ttl];

        n_ctrl = sum(ctrl_mask);
        n_ttl  = sum(ttl_mask);

        if n_ctrl == 0 && n_ttl == 0
            fprintf('  No %s cells — skipping polar.\n', on_off);
            continue;
        end

        aligned_ctrl = {ds.results(ctrl_mask).max_v_aligned};
        aligned_ttl  = {ds.results(ttl_mask).max_v_aligned};

        polar_opts = struct();
        polar_opts.stat_method = 'mean_sem';
        polar_title = sprintf('%s %s — PD-Aligned Polar (ctrl n=%d, TTL n=%d)', ...
            ds_label, on_off, n_ctrl, n_ttl);

        fig = plot_polar_population(aligned_ctrl, aligned_ttl, polar_title, polar_opts);

        fname = sprintf('comb_polar_%s_%s.png', lower(char(on_off)), ds.tag);
        exportgraphics(fig, fullfile(preview_dir, fname), 'Resolution', 150);
        close(fig);
        fprintf('  Saved: %s\n', fname);
    end

    % --- DSI comparison ---
    fig_dsi = plot_dsi_comparison(ds.results);
    % Update title
    ax_children = findall(fig_dsi, 'Type', 'tiledlayout');
    if ~isempty(ax_children)
        title(ax_children, sprintf('Direction Selectivity — %s', ds_label), 'FontSize', 14);
    end

    fname_dsi = sprintf('comb_dsi_%s.png', ds.tag);
    exportgraphics(fig_dsi, fullfile(preview_dir, fname_dsi), 'Resolution', 150);
    close(fig_dsi);
    fprintf('  Saved: %s\n', fname_dsi);

    % --- Aspect ratio comparison ---
    fig_ar = plot_aspect_ratio_comparison(ds.results);
    % Update title
    ax = findall(fig_ar, 'Type', 'axes');
    if ~isempty(ax)
        title(ax(1), sprintf('Tuning Aspect Ratio — %s', ds_label));
    end

    fname_ar = sprintf('comb_aspect_ratio_%s.png', ds.tag);
    exportgraphics(fig_ar, fullfile(preview_dir, fname_ar), 'Resolution', 150);
    close(fig_ar);
    fprintf('  Saved: %s\n', fname_ar);

    fprintf('\n');
end

fprintf('=== Done. All outputs in: %s ===\n', preview_dir);


%% ========================= Helper Functions ============================

function r_out = harmonize_results(r_in, batch_label)
% HARMONIZE_RESULTS  Extract common fields into a uniform struct array.
%   Adds a .batch field ('early' or 'late') for traceability.

    n = numel(r_in);
    r_out = struct([]);

    for k = 1:n
        s.date_str        = r_in(k).date_str;
        s.is_on           = r_in(k).is_on;
        s.is_ttl          = r_in(k).is_ttl;
        s.group           = r_in(k).group;
        s.max_v_aligned   = r_in(k).max_v_aligned;
        s.dsi_vector      = r_in(k).dsi_vector;
        s.dsi_pdnd        = r_in(k).dsi_pdnd;
        s.dir_tuning_fwhm = r_in(k).dir_tuning_fwhm;
        s.batch           = batch_label;

        if isempty(r_out)
            r_out = s;
        else
            r_out(end + 1) = s; %#ok<AGROW>
        end
    end

end


function combined = harmonize_and_merge(early, late)
% HARMONIZE_AND_MERGE  Merge early and late results into one struct array.

    early_h = harmonize_results(early, 'early');
    late_h  = harmonize_results(late, 'late');

    combined = [early_h, late_h];

end

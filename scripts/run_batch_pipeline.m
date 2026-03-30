%% run_batch_pipeline.m — Run the population analysis on all 25 experiments
% Calls batch_analyze_1DRF, then writes a diagnostic summary to disk.

close all; clear all; %#ok<CLALL>

addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF/';
save_dir  = fullfile(data_root, 'population_results');

opts = struct();
opts.save_figs   = true;
opts.save_dir    = save_dir;
opts.stat_method = 'median_mad';

fprintf('=== Starting batch pipeline ===\n');
tic;

try
    results = batch_analyze_1DRF(data_root, opts);
catch ME
    fprintf('\n!!! batch_analyze_1DRF FAILED: %s\n', ME.message);
    fprintf('    at %s line %d\n', ME.stack(1).name, ME.stack(1).line);
    for k = 2:min(5, numel(ME.stack))
        fprintf('    -> %s line %d\n', ME.stack(k).name, ME.stack(k).line);
    end
    return;
end

elapsed = toc;
fprintf('\nTotal elapsed time: %.1f seconds\n', elapsed);

%% Write diagnostic summary
summary_path = fullfile(save_dir, 'batch_summary.txt');
fid = fopen(summary_path, 'w');

fprintf(fid, 'Batch analysis summary\n');
fprintf(fid, '======================\n');
fprintf(fid, 'Date: %s\n', datetime('now'));
fprintf(fid, 'Data root: %s\n', data_root);
fprintf(fid, 'Stat method: %s\n', opts.stat_method);
fprintf(fid, 'Elapsed: %.1f s\n\n', elapsed);
fprintf(fid, 'Total cells: %d\n\n', numel(results));

% Group counts
groups = {results.group};
for g = ["on_control", "on_ttl", "off_control", "off_ttl"]
    fprintf(fid, '  %s: %d\n', g, sum(strcmp(groups, g)));
end

% Per-cell details
fprintf(fid, '\n\n--- Per-cell details ---\n');
for i = 1:numel(results)
    r = results(i);
    fprintf(fid, '[%2d] %-28s  group=%-12s  PD=%6.1f deg  strain=%s  frame=%d\n', ...
        i, r.folder, r.group, r.pd_direction, r.strain, r.frame);
end

% PD comparison: LUT-based vs circular stats
fprintf(fid, '\n\n--- PD direction comparison (LUT vs circular stats) ---\n');
fprintf(fid, '%-28s  PD_lut  PD_circ  diff\n', 'folder');
for i = 1:numel(results)
    r = results(i);
    angles = r.max_v_aligned(:,1);
    responses = r.max_v_aligned(:,2);
    [~, max_idx] = max(responses);
    pd_circ = rad2deg(angles(max_idx));
    fprintf(fid, '%-28s  %6.1f  %6.1f  %+7.1f\n', ...
        r.folder, r.pd_direction, pd_circ, pd_circ - r.pd_direction);
end

fclose(fid);
fprintf('\nDiagnostic summary -> %s\n', summary_path);

% Also save the full results + summary as .mat
save(fullfile(save_dir, 'batch_results_full.mat'), 'results', 'opts', 'elapsed');

fprintf('=== Done ===\n');

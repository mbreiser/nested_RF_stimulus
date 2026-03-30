%% run_batch_chunks.m — Process all 25 experiments and save results
% Designed to be run as a single script with local functions.

close all; clear all; %#ok<CLALL>

addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF/';
save_dir  = fullfile(data_root, 'population_results');
if ~isfolder(save_dir), mkdir(save_dir); end

%% Options
opts = struct();
opts.lut_path       = fullfile(fileparts(which('batch_analyze_1DRF')), 'bar_lut.mat');
opts.plot_order     = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
opts.baseline_range = [1000 9000];
opts.stim_trim_end  = 7000;
opts.percentile     = 98;
opts.flash_baseline = 1:5000;
opts.flash_ylim     = [-15 35];
opts.pattern_offset = 2;
opts.on_threshold   = 129;
opts.stat_method    = 'median_mad';
opts.prop_int       = 0.5;

%% Load LUT
S_lut = load(opts.lut_path, 'Tbl');
Tbl = S_lut.Tbl;

%% Discover experiments
d = dir(data_root);
d = d([d.isdir]);
d = d(~startsWith({d.name}, '.'));
valid = false(numel(d), 1);
for i = 1:numel(d)
    valid(i) = isfile(fullfile(data_root, d(i).name, 'currentExp.mat'));
end
d = d(valid);
n_exp = numel(d);
fprintf('Found %d experiment folders\n', n_exp);

%% Process all experiments
results = struct([]);
tic;

for exp_idx = 1:n_exp
    folder = d(exp_idx).name;
    exp_folder = fullfile(data_root, folder);
    fprintf('[%2d/%d] %-25s ', exp_idx, n_exp, folder);

    try
        r = process_one(exp_folder, Tbl, opts);
        r.folder = folder;

        if isempty(results)
            results = r;
        else
            results(end + 1) = r; %#ok<AGROW>
        end

        fprintf('-> %-12s  PD=%6.1f°  %s (frame=%d)\n', ...
            r.group, r.pd_direction, r.strain, r.frame);

    catch ME
        fprintf('ERROR: %s\n', ME.message);
    end
end

elapsed = toc;
fprintf('\n=== Processed %d cells in %.1f s ===\n', numel(results), elapsed);

%% Group counts
groups = {results.group};
for g = ["on_control", "on_ttl", "off_control", "off_ttl"]
    fprintf('  %s: %d\n', g, sum(strcmp(groups, g)));
end

%% Save results
save(fullfile(save_dir, 'batch_results.mat'), 'results', 'opts');
fprintf('\nResults saved.\n');

%% PD comparison table
fprintf('\n--- PD direction: LUT vs circular-stats ---\n');
fprintf('%-25s  PD_lut  PD_circ  diff\n', 'folder');
for i = 1:numel(results)
    r = results(i);
    angles_deg = rad2deg(r.max_v_aligned(:,1));
    responses  = r.max_v_aligned(:,2);
    [~, mi]    = max(responses);
    pd_circ    = angles_deg(mi);
    fprintf('%-25s  %6.1f  %6.1f  %+7.1f\n', ...
        r.folder, r.pd_direction, pd_circ, pd_circ - r.pd_direction);
end

%% Write summary file
fid = fopen(fullfile(save_dir, 'batch_summary.txt'), 'w');
fprintf(fid, 'Batch pipeline summary  %s\n', datetime('now'));
fprintf(fid, '=========================================\n');
fprintf(fid, 'n_cells = %d   elapsed = %.1f s\n', numel(results), elapsed);
fprintf(fid, 'prop_int = %.2f   stat_method = %s\n\n', opts.prop_int, opts.stat_method);
groups2 = {results.group};
for g = ["on_control", "on_ttl", "off_control", "off_ttl"]
    fprintf(fid, '  %s: %d\n', g, sum(strcmp(groups2, g)));
end
fprintf(fid, '\n--- Per-cell ---\n');
for i = 1:numel(results)
    r = results(i);
    fprintf(fid, '[%2d] %-25s  %-12s  PD=%6.1f  %s  frame=%d\n', ...
        i, r.folder, r.group, r.pd_direction, r.strain, r.frame);
end
fprintf(fid, '\n--- PD comparison ---\n');
for i = 1:numel(results)
    r = results(i);
    angles_deg = rad2deg(r.max_v_aligned(:,1));
    responses  = r.max_v_aligned(:,2);
    [~, mi]    = max(responses);
    pd_circ    = angles_deg(mi);
    fprintf(fid, '%-25s  LUT=%6.1f  Circ=%6.1f  diff=%+.1f\n', ...
        r.folder, r.pd_direction, pd_circ, pd_circ - r.pd_direction);
end
fclose(fid);
fprintf('\nSummary -> %s\n', fullfile(save_dir, 'batch_summary.txt'));


%% ========================= Local Functions ============================

function r = process_one(exp_folder, Tbl, opts)
    S = load(fullfile(exp_folder, 'currentExp.mat'), 'currentExp');
    ex = S.currentExp;
    Log = ex.Log;
    v_data = Log.ADC.Volts(2,:) * 10;
    f_data = Log.Frames;

    r.date_str = ex.date;
    r.strain   = ex.fly.Strain;
    r.frame    = ex.fly.Frame;
    r.is_on    = r.frame > opts.on_threshold;
    r.is_ttl   = contains(lower(r.strain), 'ttl');

    if r.is_on && ~r.is_ttl,      r.group = 'on_control';
    elseif r.is_on && r.is_ttl,   r.group = 'on_ttl';
    elseif ~r.is_on && ~r.is_ttl, r.group = 'off_control';
    else,                          r.group = 'off_ttl';
    end

    % Bar sweep responses
    [max_v, ~] = compute_bar_sweep_responses(v_data, f_data, ...
        opts.baseline_range, opts.stim_trim_end, opts.percentile);
    [~, ~, r.max_v_aligned] = find_PD_and_order_idx(max_v);

    % PD from LUT
    pd_info = find_pd_from_lut(max_v, Tbl.Direction, Tbl.Orientation, ...
        Tbl.PatNum, Tbl.FnNum, opts.plot_order, Tbl, opts.pattern_offset);
    r.pd_direction   = pd_info.pd_direction;
    r.pd_orientation = pd_info.pd_orientation;

    % Bar flash
    [~, ~, mean_slow_bf, ~] = parse_bar_flash_data(f_data, v_data, opts.prop_int);

    bl = opts.flash_baseline;
    r.pd_flash_bl    = extract_traces(mean_slow_bf, pd_info.bar_flash_col, pd_info.pos_order, bl);
    r.ortho_flash_bl = extract_traces(mean_slow_bf, pd_info.ortho_flash_col, pd_info.pos_order, bl);
end


function traces = extract_traces(mean_slow_bf, flash_col, pos_order, bl_samples)
    n_pos = numel(pos_order);
    traces = [];
    for p = 1:n_pos
        row  = pos_order(p);
        data = mean_slow_bf{row, flash_col};
        data = data - mean(data(bl_samples));
        if isempty(traces)
            traces = zeros(n_pos, numel(data));
        end
        L = min(numel(data), size(traces, 2));
        traces(p, 1:L) = data(1:L);
    end
end

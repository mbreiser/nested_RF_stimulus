% PLOT_THREE_EPOCH_VOLTAGE  Voltage distributions across 3 non-overlapping epochs.
%
%   Generates two figures:
%     1. Bar Sweeps (all 48 cells: 25 late + 23 early)
%     2. Bar Flashes (25 late cells only)
%
%   Each figure: 2x3 panels — Row 1 pooled histograms, Row 2 box-whisker+dots.
%   Columns: pre-stimulus, during stimulus, post-stimulus.
%
%   Usage:
%     run('scripts/plot_three_epoch_voltage.m')

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root   = '/Users/reiserm/Documents/ttl_1DRF';
preview_dir = fullfile(data_root, 'figure_previews');
if ~isfolder(preview_dir), mkdir(preview_dir); end

%% Histogram bins
bin_edges = -80:0.5:-10;
n_bins = numel(bin_edges) - 1;

%% Load batch metadata
fprintf('Loading batch results for metadata...\n');

S_late = load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
late_results = S_late.results;
n_late = numel(late_results);
fprintf('  Late batch: %d cells\n', n_late);

S_early = load(fullfile(data_root, 'pre-bar-flash', 'population_results', ...
    'batch_results_pre_bf.mat'), 'results');
early_results = S_early.results;
n_early = numel(early_results);
fprintf('  Early batch: %d cells\n', n_early);

n_total = n_late + n_early;
fprintf('  Total: %d cells\n\n', n_total);

%% Build combined cell list
groups_list = {'on_control', 'on_ttl', 'off_control', 'off_ttl'};
group_to_idx = containers.Map(groups_list, {1, 2, 3, 4});

cell_list = struct('folder', {}, 'group', {}, 'group_idx', {}, ...
    'is_late', {}, 'batch_idx', {});

for k = 1:n_late
    r = late_results(k);
    c.folder    = fullfile(data_root, r.folder);
    c.group     = r.group;
    c.group_idx = group_to_idx(r.group);
    c.is_late   = true;
    c.batch_idx = k;
    if isempty(cell_list)
        cell_list = c;
    else
        cell_list(end+1) = c; %#ok<SAGROW>
    end
end

for k = 1:n_early
    r = early_results(k);
    c.folder    = r.folder;  % already full path
    c.group     = r.group;
    c.group_idx = group_to_idx(r.group);
    c.is_late   = false;
    c.batch_idx = k;
    if isempty(cell_list)
        cell_list = c;
    else
        cell_list(end+1) = c; %#ok<SAGROW>
    end
end

%% Initialize accumulators

% Sweep: 4 groups x 3 epochs
sweep_hist = cell(4, 3);
for g = 1:4
    for e = 1:3
        sweep_hist{g,e} = zeros(1, n_bins);
    end
end

% Flash: 4 groups x 3 epochs (late cells only)
flash_hist = cell(4, 3);
for g = 1:4
    for e = 1:3
        flash_hist{g,e} = zeros(1, n_bins);
    end
end

% Per-cell stats
sweep_stats = struct('group', {}, 'epoch_median', {});
flash_stats = struct('group', {}, 'epoch_median', {});

%% ========================= MAIN LOOP ====================================
fprintf('Processing %d cells...\n', n_total);

for ci = 1:n_total
    cl = cell_list(ci);
    fprintf('[%2d/%d] %s  (%s, %s) ...', ...
        ci, n_total, cl.folder, cl.group, ternary(cl.is_late, 'late', 'early'));

    orig_dir = pwd;
    try
        % Load raw data
        [~, ~, Log, ~, ~] = load_protocol2_data(cl.folder);
        cd(orig_dir);  % restore immediately

        f_data = Log.ADC.Volts(1, :);
        v_data = Log.ADC.Volts(2, :) * 10;

        %% ---------- BAR SWEEP EPOCHS ----------
        if cl.is_late
            bar_data = parse_bar_data(f_data, v_data);
        else
            bar_data = parse_bar_data_pre_bf(f_data, v_data);
        end

        % Pool voltage from all 16 slow-bar mean traces (column 4)
        n_dirs = size(bar_data, 1);
        pre_v_all  = [];
        stim_v_all = [];
        post_v_all = [];

        for d = 1:n_dirs
            trace = bar_data{d, 4};
            if isempty(trace), continue; end
            trace = trace(:)';  % ensure row vector
            tlen = numel(trace);

            % Pre-stimulus: 1000:9000 (skip first 100ms for edge artifacts)
            pre_end = min(9000, tlen);
            pre_v_all = [pre_v_all, trace(1000:pre_end)]; %#ok<AGROW>

            % During stimulus: 9001:(end-7000)
            stim_end = tlen - 7000;
            if stim_end > 9001
                stim_v_all = [stim_v_all, trace(9001:stim_end)]; %#ok<AGROW>
            end

            % Post-stimulus: (end-7000+1):(end-1000)
            post_start = stim_end + 1;
            post_end   = tlen - 1000;
            if post_end > post_start
                post_v_all = [post_v_all, trace(post_start:post_end)]; %#ok<AGROW>
            end
        end

        g = cl.group_idx;
        sw.group = cl.group;
        sw.epoch_median = [median(pre_v_all), median(stim_v_all), median(post_v_all)];

        if isempty(sweep_stats)
            sweep_stats = sw;
        else
            sweep_stats(end+1) = sw; %#ok<SAGROW>
        end

        % Accumulate histograms
        sweep_hist{g,1} = sweep_hist{g,1} + histcounts(pre_v_all,  bin_edges);
        sweep_hist{g,2} = sweep_hist{g,2} + histcounts(stim_v_all, bin_edges);
        sweep_hist{g,3} = sweep_hist{g,3} + histcounts(post_v_all, bin_edges);

        fprintf(' sweep OK (%dk/%dk/%dk)', ...
            round(numel(pre_v_all)/1000), ...
            round(numel(stim_v_all)/1000), ...
            round(numel(post_v_all)/1000));

        %% ---------- BAR FLASH EPOCHS (late batch only) ----------
        if cl.is_late
            [~, ~, mean_slow, ~] = parse_bar_flash_data(f_data, v_data, 0.5);

            flash_pre_all  = [];
            flash_stim_all = [];
            flash_post_all = [];

            [n_pos, n_ori] = size(mean_slow);
            for ori = 1:n_ori
                for pos = 1:n_pos
                    trace = mean_slow{pos, ori};
                    if isempty(trace), continue; end
                    trace = trace(:)';
                    tlen = numel(trace);

                    % Pre-stimulus: 1:5000
                    pre_end_f = min(5000, tlen);
                    flash_pre_all = [flash_pre_all, trace(1:pre_end_f)]; %#ok<AGROW>

                    % During stimulus: 5001:5801
                    if tlen >= 5801
                        flash_stim_all = [flash_stim_all, trace(5001:5801)]; %#ok<AGROW>
                    end

                    % Post-stimulus: 5802:end
                    if tlen > 5802
                        flash_post_all = [flash_post_all, trace(5802:end)]; %#ok<AGROW>
                    end
                end
            end

            fl.group = cl.group;
            fl.epoch_median = [median(flash_pre_all), ...
                               median(flash_stim_all), ...
                               median(flash_post_all)];

            if isempty(flash_stats)
                flash_stats = fl;
            else
                flash_stats(end+1) = fl; %#ok<SAGROW>
            end

            flash_hist{g,1} = flash_hist{g,1} + histcounts(flash_pre_all,  bin_edges);
            flash_hist{g,2} = flash_hist{g,2} + histcounts(flash_stim_all, bin_edges);
            flash_hist{g,3} = flash_hist{g,3} + histcounts(flash_post_all, bin_edges);

            fprintf(' | flash OK (%dk/%dk/%dk)', ...
                round(numel(flash_pre_all)/1000), ...
                round(numel(flash_stim_all)/1000), ...
                round(numel(flash_post_all)/1000));
        end

        fprintf('\n');

    catch ME
        cd(orig_dir);
        fprintf(' ERROR: %s\n', ME.message);
        continue;
    end
end

fprintf('\n=== Data extraction complete ===\n\n');

%% Print summary
fprintf('SWEEP STATS (n=%d):\n', numel(sweep_stats));
for g = 1:4
    mask = strcmp({sweep_stats.group}, groups_list{g});
    if ~any(mask), continue; end
    meds = vertcat(sweep_stats(mask).epoch_median);
    fprintf('  %s (n=%d):  Pre=%.1f  Stim=%.1f  Post=%.1f mV\n', ...
        groups_list{g}, sum(mask), ...
        mean(meds(:,1)), mean(meds(:,2)), mean(meds(:,3)));
end

fprintf('\nFLASH STATS (n=%d):\n', numel(flash_stats));
for g = 1:4
    mask = strcmp({flash_stats.group}, groups_list{g});
    if ~any(mask), continue; end
    meds = vertcat(flash_stats(mask).epoch_median);
    fprintf('  %s (n=%d):  Pre=%.1f  Stim=%.1f  Post=%.1f mV\n', ...
        groups_list{g}, sum(mask), ...
        mean(meds(:,1)), mean(meds(:,2)), mean(meds(:,3)));
end

%% ========================= PLOT ==========================================

fprintf('\nGenerating figures...\n');

% --- Sweep figure ---
fig1 = plot_three_epoch_distributions(sweep_hist, sweep_stats, bin_edges, ...
    sprintf('Bar Sweep Voltage — 3 Epochs (n=%d cells)', numel(sweep_stats)));
png1 = fullfile(preview_dir, 'pop_voltage_epochs_sweeps.png');
exportgraphics(fig1, png1, 'Resolution', 150);
fprintf('  Saved: %s\n', png1);

% --- Flash figure ---
fig2 = plot_three_epoch_distributions(flash_hist, flash_stats, bin_edges, ...
    sprintf('Bar Flash Voltage — 3 Epochs (n=%d cells)', numel(flash_stats)));
png2 = fullfile(preview_dir, 'pop_voltage_epochs_flashes.png');
exportgraphics(fig2, png2, 'Resolution', 150);
fprintf('  Saved: %s\n', png2);

fprintf('\n=== Done ===\n');


%% ========================= Helpers ======================================

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

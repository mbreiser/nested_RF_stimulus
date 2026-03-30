% COMPARE_M2_VS_M5_STATS  Compare significance counts between M2 and M5 alignment.
%
%   Runs both classic (per-position raw mV Wilcoxon) and Gruntman
%   (normalized sliding-mean Wilcoxon) statistical tests using M2 vs M5
%   alignment and prints a summary comparison table.
%
%   No figures generated — stats only.
%
%   Usage:
%     run('scripts/compare_m2_vs_m5_stats.m')

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
out_dir   = fullfile(data_root, 'manuscript_figures');

%% Load batch results
fprintf('Loading batch results...\n');
S = load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
results = S.results;
fprintf('  %d cells loaded\n', numel(results));

% Group masks
on_mask   = [results.is_on];
off_mask  = ~on_mask;
ctrl_mask = ~[results.is_ttl];
ttl_mask  = [results.is_ttl];

fprintf('  ON ctrl=%d, ON TTL=%d, OFF ctrl=%d, OFF TTL=%d\n', ...
    sum(on_mask & ctrl_mask), sum(on_mask & ttl_mask), ...
    sum(off_mask & ctrl_mask), sum(off_mask & ttl_mask));

%% Parameters
STIM_ON   = 5001;
RESP_END  = 6551;
DEP_PCTILE  = 99.5;
HYP_PCTILE  = 0.5;
REJECT_THRESH = 0.5;  % mV, for classic only
DEP_SD_MULT  = 3;
HYP_SD_MULT  = 2;

alignments = {'peak', 'm5'};  % peak = M2, m5 = M5
align_labels = {'M2', 'M5'};
cell_types = {'ON', 'OFF'};

positions = -5:5;

%% Open output file
stats_file = fullfile(out_dir, 'compare_m2_vs_m5_stats.txt');
fid = fopen(stats_file, 'w');
fprintf(fid, 'M2 vs M5 Alignment — Significance Comparison\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

%% Collect results
% Structure: sig_counts(alignment, cell_type, method, polarity)
summary = {};  % will build a table

for a = 1:2
    align = alignments{a};
    align_lab = align_labels{a};

    if strcmp(align, 'peak')
        pd_field = 'pd_flash_peak_aligned';
    else
        pd_field = 'pd_flash_m5_aligned';
    end

    fprintf(fid, '================================================================\n');
    fprintf(fid, '  ALIGNMENT: %s (%s)\n', align_lab, pd_field);
    fprintf(fid, '================================================================\n\n');

    for t = 1:2
        type_str = cell_types{t};
        if t == 1
            mask = on_mask;
            type_long = 'T4 (ON)';
        else
            mask = off_mask;
            type_long = 'T5 (OFF)';
        end

        ctrl_sel = mask & ctrl_mask;
        ttl_sel  = mask & ttl_mask;
        nc = sum(ctrl_sel);
        nt = sum(ttl_sel);

        fprintf(fid, '--- %s: ctrl n=%d, tutl- n=%d ---\n\n', type_long, nc, nt);

        %% ===== Classic: D/E trace-level stats (PD + ortho axis) =====
        % PD axis
        pd_ctrl  = {results(ctrl_sel).(pd_field)};
        pd_ttl   = {results(ttl_sel).(pd_field)};

        % Ortho axis
        if strcmp(align, 'peak')
            ortho_field = 'ortho_flash_peak_aligned';
        else
            ortho_field = 'ortho_flash_m5_aligned';
        end
        ortho_ctrl = {results(ctrl_sel).(ortho_field)};
        ortho_ttl  = {results(ttl_sel).(ortho_field)};

        % Trace-level per-position Wilcoxon (D/E style)
        pd_trace_stats    = compute_position_trace_stats(pd_ctrl, pd_ttl, STIM_ON, RESP_END);
        ortho_trace_stats = compute_position_trace_stats(ortho_ctrl, ortho_ttl, STIM_ON, RESP_END);

        n_pd_sig    = sum([pd_trace_stats.sig]);
        n_ortho_sig = sum([ortho_trace_stats.sig]);

        fprintf(fid, 'Classic D/E (trace-level 99.5%% pctile, per-position Wilcoxon):\n');
        write_trace_stats_table(fid, 'PD axis', pd_trace_stats, positions);
        write_trace_stats_table(fid, 'Ortho axis', ortho_trace_stats, positions);

        %% ===== Classic: F/G amplitude stats (depol + hyperpol) =====
        [ctrl_dep, ctrl_hyp] = extract_robust_amps(results(ctrl_sel), pd_field, ...
            [STIM_ON, RESP_END], [STIM_ON, Inf], DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
        [ttl_dep, ttl_hyp] = extract_robust_amps(results(ttl_sel), pd_field, ...
            [STIM_ON, RESP_END], [STIM_ON, Inf], DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

        dep_amp_stats = compute_amp_stats(ctrl_dep, ttl_dep);
        hyp_amp_stats = compute_amp_stats(ctrl_hyp, ttl_hyp);

        n_dep_amp_sig = sum([dep_amp_stats.sig]);
        n_hyp_amp_sig = sum([hyp_amp_stats.sig]);

        fprintf(fid, 'Classic F/G (raw mV amplitudes, per-position Wilcoxon):\n');
        write_amp_stats_table(fid, 'Depolarization', dep_amp_stats, positions);
        write_amp_stats_table(fid, 'Hyperpolarization', hyp_amp_stats, positions);

        %% ===== Gruntman: normalized + sliding-mean =====
        [ctrl_dep_raw, ctrl_hyp_raw, cd_zeroed, ch_zeroed] = ...
            extract_position_profiles(results(ctrl_sel), pd_field, ...
            STIM_ON, RESP_END, DEP_PCTILE, HYP_PCTILE, DEP_SD_MULT, HYP_SD_MULT);
        [ttl_dep_raw, ttl_hyp_raw, td_zeroed, th_zeroed] = ...
            extract_position_profiles(results(ttl_sel), pd_field, ...
            STIM_ON, RESP_END, DEP_PCTILE, HYP_PCTILE, DEP_SD_MULT, HYP_SD_MULT);

        % Normalize
        ctrl_dep_norm = normalize_profiles(ctrl_dep_raw, 'depol');
        ctrl_hyp_norm = normalize_profiles(ctrl_hyp_raw, 'hyperpol');
        ttl_dep_norm  = normalize_profiles(ttl_dep_raw, 'depol');
        ttl_hyp_norm  = normalize_profiles(ttl_hyp_raw, 'hyperpol');

        % Sliding-mean Wilcoxon (with FDR)
        dep_sm_stats = sliding_mean_wilcoxon(ctrl_dep_norm, ttl_dep_norm);
        hyp_sm_stats = sliding_mean_wilcoxon(ctrl_hyp_norm, ttl_hyp_norm);

        % Count raw p<0.05 (ignoring FDR)
        n_dep_sm_raw = sum([dep_sm_stats.p] < 0.05);
        n_hyp_sm_raw = sum([hyp_sm_stats.p] < 0.05);

        % Count FDR-corrected
        n_dep_sm_fdr = sum([dep_sm_stats.sig]);
        n_hyp_sm_fdr = sum([hyp_sm_stats.sig]);

        fprintf(fid, 'Gruntman (normalized, sliding-mean Wilcoxon, 3-pos window):\n');
        fprintf(fid, '  SD zeroed: ctrl dep=%d, ctrl hyp=%d, TTL dep=%d, TTL hyp=%d\n', ...
            cd_zeroed, ch_zeroed, td_zeroed, th_zeroed);
        write_sm_stats_table(fid, 'Depol sliding-mean', dep_sm_stats);
        write_sm_stats_table(fid, 'Hyperpol sliding-mean', hyp_sm_stats);

        %% ===== Also run per-position Wilcoxon on normalized data (no sliding mean) =====
        dep_norm_pp_stats = compute_amp_stats(ctrl_dep_norm, ttl_dep_norm);
        hyp_norm_pp_stats = compute_amp_stats(ctrl_hyp_norm, ttl_hyp_norm);
        n_dep_norm_pp = sum([dep_norm_pp_stats.sig]);
        n_hyp_norm_pp = sum([hyp_norm_pp_stats.sig]);

        fprintf(fid, 'Normalized per-position Wilcoxon (no sliding mean):\n');
        write_amp_stats_table(fid, 'Depol (norm)', dep_norm_pp_stats, positions);
        write_amp_stats_table(fid, 'Hyperpol (norm)', hyp_norm_pp_stats, positions);

        %% Collect summary row
        summary{end+1} = struct( ...
            'align', align_lab, 'type', type_str, ...
            'pd_trace', n_pd_sig, 'ortho_trace', n_ortho_sig, ...
            'dep_amp', n_dep_amp_sig, 'hyp_amp', n_hyp_amp_sig, ...
            'dep_sm_raw', n_dep_sm_raw, 'hyp_sm_raw', n_hyp_sm_raw, ...
            'dep_sm_fdr', n_dep_sm_fdr, 'hyp_sm_fdr', n_hyp_sm_fdr, ...
            'dep_norm_pp', n_dep_norm_pp, 'hyp_norm_pp', n_hyp_norm_pp); %#ok<SAGROW>

        fprintf(fid, '\n');
    end
end

%% Print summary table
fprintf(fid, '\n================================================================\n');
fprintf(fid, '  SUMMARY: Number of positions with raw p < 0.05\n');
fprintf(fid, '================================================================\n\n');

fprintf(fid, '%-6s %-5s | D/E PD  D/E ort | F/G dep F/G hyp | SM dep  SM hyp  | Norm dep Norm hyp\n', ...
    'Align', 'Type');
fprintf(fid, '%s\n', repmat('-', 1, 95));

for i = 1:numel(summary)
    s = summary{i};
    fprintf(fid, '%-6s %-5s |   %d       %d     |   %d       %d     |   %d       %d       |   %d        %d\n', ...
        s.align, s.type, ...
        s.pd_trace, s.ortho_trace, ...
        s.dep_amp, s.hyp_amp, ...
        s.dep_sm_raw, s.hyp_sm_raw, ...
        s.dep_norm_pp, s.hyp_norm_pp);
end

fprintf(fid, '\nColumn key:\n');
fprintf(fid, '  D/E PD    = Classic trace-level, PD axis, per-position Wilcoxon (uncorrected)\n');
fprintf(fid, '  D/E ort   = Classic trace-level, ortho axis, per-position Wilcoxon (uncorrected)\n');
fprintf(fid, '  F/G dep   = Classic raw mV depol amplitude, per-position Wilcoxon (uncorrected)\n');
fprintf(fid, '  F/G hyp   = Classic raw mV hyperpol amplitude, per-position Wilcoxon (uncorrected)\n');
fprintf(fid, '  SM dep    = Gruntman normalized, 3-pos sliding-mean Wilcoxon (raw p<0.05)\n');
fprintf(fid, '  SM hyp    = Gruntman normalized, 3-pos sliding-mean Wilcoxon (raw p<0.05)\n');
fprintf(fid, '  Norm dep  = Gruntman normalized, per-position Wilcoxon (no sliding mean, uncorrected)\n');
fprintf(fid, '  Norm hyp  = Gruntman normalized, per-position Wilcoxon (no sliding mean, uncorrected)\n');

% Also print FDR-corrected counts
fprintf(fid, '\n%-6s %-5s | SM dep(FDR) SM hyp(FDR)\n', 'Align', 'Type');
fprintf(fid, '%s\n', repmat('-', 1, 40));
for i = 1:numel(summary)
    s = summary{i};
    fprintf(fid, '%-6s %-5s |   %d           %d\n', s.align, s.type, s.dep_sm_fdr, s.hyp_sm_fdr);
end

fclose(fid);

% Also print to console
fprintf('\n');
type(stats_file);

fprintf('\nSaved: %s\n', stats_file);


%% ========================= Local Functions ===============================

function stats = compute_position_trace_stats(traces_ctrl, traces_ttl, stim_on, resp_end)
% Per-position Wilcoxon rank-sum on 99.5th pctile of trace data.
    n_pos = 11;
    stats = struct('pos', num2cell(1:n_pos), 'n_ctrl', 0, 'n_ttl', 0, ...
        'mean_ctrl', NaN, 'mean_ttl', NaN, 'p', NaN, 'sig', false);

    for p = 1:n_pos
        ctrl_vals = extract_peak_at_pos(traces_ctrl, p, stim_on, resp_end);
        ttl_vals  = extract_peak_at_pos(traces_ttl, p, stim_on, resp_end);
        stats(p).n_ctrl = numel(ctrl_vals);
        stats(p).n_ttl  = numel(ttl_vals);
        if ~isempty(ctrl_vals), stats(p).mean_ctrl = mean(ctrl_vals); end
        if ~isempty(ttl_vals),  stats(p).mean_ttl  = mean(ttl_vals); end

        if numel(ctrl_vals) >= 2 && numel(ttl_vals) >= 2
            stats(p).p = ranksum(ctrl_vals, ttl_vals);
            stats(p).sig = stats(p).p < 0.05;
        end
    end
end


function vals = extract_peak_at_pos(traces, pos, stim_on, resp_end)
% Extract 99.5th percentile peak for each cell at a given position.
    vals = [];
    for c = 1:numel(traces)
        mat = traces{c};
        if isempty(mat) || pos > size(mat, 1), continue; end
        row = mat(pos, :);
        if all(isnan(row)), continue; end
        win_end = min(resp_end, numel(row));
        if stim_on > numel(row), continue; end
        win = row(stim_on:win_end);
        win = win(~isnan(win));
        if isempty(win), continue; end
        vals(end+1) = prctile(win, 99.5); %#ok<AGROW>
    end
end


function [dep_mat, hyp_mat] = extract_robust_amps(results_sub, trace_field, ...
    dep_win, hyp_win, dep_pct, hyp_pct, reject_thresh)
% Extract per-cell robust amplitude metrics at each of 11 positions.
    n = numel(results_sub);
    dep_mat = NaN(n, 11);
    hyp_mat = NaN(n, 11);

    for k = 1:n
        if ~isfield(results_sub(k), trace_field), continue; end
        traces = results_sub(k).(trace_field);
        if isempty(traces), continue; end
        n_samp = size(traces, 2);

        for pos = 1:11
            row = traces(pos, :);
            if all(isnan(row)), continue; end

            stim_start = dep_win(1);
            if stim_start > n_samp, continue; end
            mean_resp = mean(row(stim_start:end), 'omitnan');
            if abs(mean_resp) < reject_thresh, continue; end

            % Depolarization
            de = min(dep_win(2), n_samp);
            if ~isinf(de)
                dep_mat(k, pos) = prctile(row(dep_win(1):de), dep_pct);
            else
                dep_mat(k, pos) = prctile(row(dep_win(1):end), dep_pct);
            end

            % Hyperpolarization
            he = hyp_win(2);
            if isinf(he), he = n_samp; end
            he = min(he, n_samp);
            hyp_mat(k, pos) = prctile(row(hyp_win(1):he), hyp_pct);
        end
    end
end


function stats = compute_amp_stats(ctrl_data, ttl_data)
% Per-position Wilcoxon rank-sum on amplitude matrices.
    n_pos = size(ctrl_data, 2);
    stats = struct('pos', num2cell(1:n_pos), 'n_ctrl', 0, 'n_ttl', 0, ...
        'mean_ctrl', NaN, 'mean_ttl', NaN, 'p', NaN, 'sig', false);

    for i = 1:n_pos
        cv = ctrl_data(:, i); cv = cv(~isnan(cv));
        tv = ttl_data(:, i);  tv = tv(~isnan(tv));
        stats(i).n_ctrl = numel(cv);
        stats(i).n_ttl  = numel(tv);
        if ~isempty(cv), stats(i).mean_ctrl = mean(cv); end
        if ~isempty(tv), stats(i).mean_ttl  = mean(tv); end

        if numel(cv) >= 2 && numel(tv) >= 2
            stats(i).p = ranksum(cv, tv);
            stats(i).sig = stats(i).p < 0.05;
        end
    end
end


function [dep_mat, hyp_mat, n_dep_zeroed, n_hyp_zeroed] = extract_position_profiles( ...
    results_sub, trace_field, stim_on, resp_end, dep_pct, hyp_pct, dep_sd_mult, hyp_sd_mult)
% Extract per-cell amplitude profiles with Gruntman SD thresholding.
    n = numel(results_sub);
    dep_mat = NaN(n, 11);
    hyp_mat = NaN(n, 11);
    n_dep_zeroed = 0;
    n_hyp_zeroed = 0;

    for k = 1:n
        if ~isfield(results_sub(k), trace_field), continue; end
        traces = results_sub(k).(trace_field);
        if isempty(traces), continue; end
        n_samp = size(traces, 2);
        if stim_on > n_samp, continue; end

        % Pooled baseline SD across all 11 positions
        baseline_samples = [];
        for pos = 1:11
            row = traces(pos, :);
            if all(isnan(row)), continue; end
            bl = row(1:stim_on-1);
            bl = bl(~isnan(bl));
            baseline_samples = [baseline_samples, bl]; %#ok<AGROW>
        end
        if numel(baseline_samples) < 100, continue; end
        baseline_sd = std(baseline_samples);

        for pos = 1:11
            row = traces(pos, :);
            if all(isnan(row)), continue; end

            de = min(resp_end, n_samp);
            win_dep = row(stim_on:de);
            win_dep = win_dep(~isnan(win_dep));
            if ~isempty(win_dep)
                val = prctile(win_dep, dep_pct);
                if val < dep_sd_mult * baseline_sd
                    dep_mat(k, pos) = 0;
                    n_dep_zeroed = n_dep_zeroed + 1;
                else
                    dep_mat(k, pos) = val;
                end
            end

            win_hyp = row(stim_on:end);
            win_hyp = win_hyp(~isnan(win_hyp));
            if ~isempty(win_hyp)
                val = prctile(win_hyp, hyp_pct);
                if abs(val) < hyp_sd_mult * baseline_sd
                    hyp_mat(k, pos) = 0;
                    n_hyp_zeroed = n_hyp_zeroed + 1;
                else
                    hyp_mat(k, pos) = val;
                end
            end
        end
    end
end


function norm_mat = normalize_profiles(raw_mat, polarity)
% Normalize per-cell profiles.
    [n, p] = size(raw_mat);
    norm_mat = NaN(n, p);

    for k = 1:n
        profile = raw_mat(k, :);
        valid = ~isnan(profile);
        if sum(valid) < 3, continue; end

        if strcmp(polarity, 'depol')
            profile(profile < 0) = 0;
            pk = max(profile);
            if pk > 0
                norm_mat(k, :) = profile / pk;
            end
        else
            profile(profile > 0) = 0;
            mn = min(profile);
            if mn < 0
                norm_mat(k, :) = profile / abs(mn);
            end
        end
    end
end


function stats = sliding_mean_wilcoxon(ctrl_norm, ttl_norm)
% 3-position sliding-mean Wilcoxon rank-sum with BH FDR q=0.05.
    centers = -4:4;
    n_centers = numel(centers);
    stats = struct('center', num2cell(centers), ...
        'n_ctrl', 0, 'n_ttl', 0, ...
        'mean_ctrl', NaN, 'mean_ttl', NaN, ...
        'p', NaN, 'sig', false);

    raw_p = NaN(1, n_centers);
    for i = 1:n_centers
        ci = centers(i) + 6;
        cols = (ci-1):(ci+1);

        ctrl_vals = nanmean(ctrl_norm(:, cols), 2); %#ok<NANMEAN>
        ctrl_vals = ctrl_vals(~isnan(ctrl_vals));
        ttl_vals = nanmean(ttl_norm(:, cols), 2); %#ok<NANMEAN>
        ttl_vals = ttl_vals(~isnan(ttl_vals));

        stats(i).n_ctrl = numel(ctrl_vals);
        stats(i).n_ttl  = numel(ttl_vals);
        if ~isempty(ctrl_vals), stats(i).mean_ctrl = mean(ctrl_vals); end
        if ~isempty(ttl_vals),  stats(i).mean_ttl  = mean(ttl_vals); end

        if numel(ctrl_vals) >= 2 && numel(ttl_vals) >= 2
            raw_p(i) = ranksum(ctrl_vals, ttl_vals);
            stats(i).p = raw_p(i);
        end
    end

    % BH FDR correction
    valid_idx = find(~isnan(raw_p));
    if ~isempty(valid_idx)
        [sorted_p, sort_order] = sort(raw_p(valid_idx));
        m = numel(sorted_p);
        q = 0.05;
        threshold = (1:m)' / m * q;
        sig_mask = sorted_p(:) <= threshold;
        last_sig = find(sig_mask, 1, 'last');
        if ~isempty(last_sig)
            sig_indices = sort_order(1:last_sig);
            for j = sig_indices(:)'
                stats(valid_idx(j)).sig = true;
            end
        end
    end
end


function write_trace_stats_table(fid, label, stats, positions)
    fprintf(fid, '  %s:\n', label);
    fprintf(fid, '  %-5s  %5s  %5s  %8s  %8s  %8s  %s\n', ...
        'Pos', 'nCtrl', 'nTTL', 'MnCtrl', 'MnTTL', 'p', 'Sig');
    fprintf(fid, '  %s\n', repmat('-', 1, 55));
    pos_labels = make_pos_labels(positions);
    for i = 1:numel(stats)
        sig_str = '';
        if stats(i).sig
            if stats(i).p < 0.001, sig_str = '***';
            elseif stats(i).p < 0.01, sig_str = '**';
            else, sig_str = '*';
            end
        end
        fprintf(fid, '  %-5s  %5d  %5d  %8.2f  %8.2f  %8.4f  %s\n', ...
            pos_labels{i}, stats(i).n_ctrl, stats(i).n_ttl, ...
            stats(i).mean_ctrl, stats(i).mean_ttl, stats(i).p, sig_str);
    end
    fprintf(fid, '\n');
end


function write_amp_stats_table(fid, label, stats, positions)
    fprintf(fid, '  %s:\n', label);
    fprintf(fid, '  %-5s  %5s  %5s  %8s  %8s  %8s  %s\n', ...
        'Pos', 'nCtrl', 'nTTL', 'MnCtrl', 'MnTTL', 'p', 'Sig');
    fprintf(fid, '  %s\n', repmat('-', 1, 55));
    for i = 1:numel(stats)
        sig_str = '';
        if stats(i).sig
            if stats(i).p < 0.001, sig_str = '***';
            elseif stats(i).p < 0.01, sig_str = '**';
            else, sig_str = '*';
            end
        end
        fprintf(fid, '  %+3d    %5d  %5d  %8.3f  %8.3f  %8.4f  %s\n', ...
            positions(i), stats(i).n_ctrl, stats(i).n_ttl, ...
            stats(i).mean_ctrl, stats(i).mean_ttl, stats(i).p, sig_str);
    end
    fprintf(fid, '\n');
end


function write_sm_stats_table(fid, label, sm_stats)
    fprintf(fid, '  %s:\n', label);
    fprintf(fid, '  %-7s  %5s  %5s  %8s  %8s  %8s  %s\n', ...
        'Center', 'nCtrl', 'nTTL', 'MnCtrl', 'MnTTL', 'p', 'Sig(FDR)');
    fprintf(fid, '  %s\n', repmat('-', 1, 60));
    for i = 1:numel(sm_stats)
        sig_str = '';
        raw_str = '';
        if sm_stats(i).sig
            sig_str = ' FDR*';
        end
        if ~isnan(sm_stats(i).p) && sm_stats(i).p < 0.05
            raw_str = '*';
        end
        fprintf(fid, '  %+3d      %5d  %5d  %8.3f  %8.3f  %8.4f  %s%s\n', ...
            sm_stats(i).center, sm_stats(i).n_ctrl, sm_stats(i).n_ttl, ...
            sm_stats(i).mean_ctrl, sm_stats(i).mean_ttl, ...
            sm_stats(i).p, raw_str, sig_str);
    end
    fprintf(fid, '\n');
end


function labels = make_pos_labels(positions)
    labels = cell(1, numel(positions));
    for i = 1:numel(positions)
        p = positions(i);
        if p == -5,      labels{i} = 'ND';
        elseif p == 5,   labels{i} = 'PD';
        elseif p == 0,   labels{i} = '0';
        else,            labels{i} = sprintf('%+d', p);
        end
    end
end

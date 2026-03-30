% EXPORT_ALL_STATS_TO_XLSX  Extract all statistical results from both
% manuscript figure scripts and write to a multi-sheet Excel file.
%
%   Output: manuscript_figures/figure_statistics.xlsx
%
%   Sheets:
%     1. T4_Direction_Wilcoxon — 16 directions, per-direction Wilcoxon (Panel A)
%     2. T5_Direction_Wilcoxon — 16 directions, per-direction Wilcoxon (Panel B)
%     3. PanelC_DSI_AR         — DSI and Aspect Ratio pairwise Wilcoxon (Panel C)
%     4. EFGH_Pooled_Ranksum   — Pooled 3-position Wilcoxon (Panels F/H)
%     5. EFGH_FWHM             — Per-cell FWHM Wilcoxon (Panels F/H)

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
out_file  = fullfile(data_root, 'manuscript_figures', 'figure_statistics.xlsx');

% Delete existing file to avoid stale sheets
if isfile(out_file), delete(out_file); end

%% ========================================================================
%  SHEETS 1 & 2: Per-direction Wilcoxon (from ABC cache)
%  ========================================================================

cache_file = fullfile(data_root, 'population_results', 'ring_of_traces_cache_gauss_999.mat');
fprintf('Loading ABC cache: %s\n', cache_file);
C = load(cache_file);
all_cells         = C.all_cells;
groups_final      = C.groups_final;
pd_aligned_angles = C.pd_aligned_angles;

% T4 (ON): groups_final(1) = ON ctrl, groups_final(2) = ON tutl-
[on_stats, ~] = compute_direction_stats_local( ...
    all_cells, groups_final(1).indices, groups_final(2).indices, pd_aligned_angles);

% T5 (OFF): groups_final(3) = OFF ctrl, groups_final(4) = OFF tutl-
[off_stats, ~] = compute_direction_stats_local( ...
    all_cells, groups_final(3).indices, groups_final(4).indices, pd_aligned_angles);

% Build tables
T4_tbl = build_direction_table(on_stats);
T5_tbl = build_direction_table(off_stats);

writetable(T4_tbl, out_file, 'Sheet', 'T4_Direction_Wilcoxon');
fprintf('  Wrote sheet: T4_Direction_Wilcoxon\n');
writetable(T5_tbl, out_file, 'Sheet', 'T5_Direction_Wilcoxon');
fprintf('  Wrote sheet: T5_Direction_Wilcoxon\n');

%% ========================================================================
%  SHEET 3: Panel C — DSI and Aspect Ratio pairwise Wilcoxon
%  ========================================================================

POLAR_PERCENTILE = 99.9;

% Build combined struct from all_cells
combined = struct([]);
for ci = 1:numel(all_cells)
    c = all_cells(ci);
    s.is_on  = c.is_on;
    s.is_ttl = c.is_ttl;
    d = c.peak_amps;  % 16x2 [angles, peaks]
    s.dsi_pdnd = (d(5,2) - d(13,2)) / (d(5,2) + d(13,2));
    s.max_v_aligned = d;
    combined = [combined, s]; %#ok<AGROW>
end

% DSI values per group
group_masks = {
    [combined.is_on] & ~[combined.is_ttl], ...   % T4 ctrl
    [combined.is_on] &  [combined.is_ttl], ...    % T4 tutl-
    ~[combined.is_on] & ~[combined.is_ttl], ...   % T5 ctrl
    ~[combined.is_on] &  [combined.is_ttl]        % T5 tutl-
};
group_names = {'T4_ctrl', 'T4_tutl-', 'T5_ctrl', 'T5_tutl-'};

dsi_groups = cell(1, 4);
ar_groups  = cell(1, 4);
for g = 1:4
    mask = group_masks{g};
    dsi_groups{g} = [combined(mask).dsi_pdnd];

    ar_vals = NaN(1, sum(mask));
    idx = find(mask);
    for k = 1:numel(idx)
        d = combined(idx(k)).max_v_aligned;
        if isnumeric(d) && size(d,1) == 16 && size(d,2) == 2
            ar_vals(k) = compute_ar_local(d);
        end
    end
    ar_groups{g} = ar_vals(~isnan(ar_vals));
end

% 3 planned comparisons per metric
comparisons = {
    'T4 ctrl vs T4 tutl-', 1, 2;
    'T5 ctrl vs T5 tutl-', 3, 4;
    'T4 ctrl vs T5 ctrl',  1, 3;
};

comp_labels = cell(6, 1);
metric_col  = cell(6, 1);
n1_col = NaN(6, 1); n2_col = NaN(6, 1);
median1 = NaN(6, 1); median2 = NaN(6, 1);
mean1 = NaN(6, 1); mean2 = NaN(6, 1);
p_col = NaN(6, 1);
sig_col = cell(6, 1);

row = 0;
for ci = 1:3
    row = row + 1;
    g1 = comparisons{ci, 2}; g2 = comparisons{ci, 3};
    comp_labels{row} = comparisons{ci, 1};
    metric_col{row}  = 'DSI (PD-ND)';
    v1 = dsi_groups{g1}; v2 = dsi_groups{g2};
    n1_col(row) = numel(v1); n2_col(row) = numel(v2);
    median1(row) = median(v1); median2(row) = median(v2);
    mean1(row) = mean(v1); mean2(row) = mean(v2);
    if numel(v1) >= 2 && numel(v2) >= 2
        p_col(row) = ranksum(v1, v2);
    end
    if ~isnan(p_col(row))
        if p_col(row) < 0.001, sig_col{row} = '***';
        elseif p_col(row) < 0.01, sig_col{row} = '**';
        elseif p_col(row) < 0.05, sig_col{row} = '*';
        else, sig_col{row} = 'n.s.';
        end
    else
        sig_col{row} = '';
    end
end
for ci = 1:3
    row = row + 1;
    g1 = comparisons{ci, 2}; g2 = comparisons{ci, 3};
    comp_labels{row} = comparisons{ci, 1};
    metric_col{row}  = 'Aspect Ratio';
    v1 = ar_groups{g1}; v2 = ar_groups{g2};
    n1_col(row) = numel(v1); n2_col(row) = numel(v2);
    median1(row) = median(v1); median2(row) = median(v2);
    mean1(row) = mean(v1); mean2(row) = mean(v2);
    if numel(v1) >= 2 && numel(v2) >= 2
        p_col(row) = ranksum(v1, v2);
    end
    if ~isnan(p_col(row))
        if p_col(row) < 0.001, sig_col{row} = '***';
        elseif p_col(row) < 0.01, sig_col{row} = '**';
        elseif p_col(row) < 0.05, sig_col{row} = '*';
        else, sig_col{row} = 'n.s.';
        end
    else
        sig_col{row} = '';
    end
end

panelC_tbl = table(metric_col, comp_labels, n1_col, n2_col, ...
    median1, median2, mean1, mean2, p_col, sig_col, ...
    'VariableNames', {'Metric', 'Comparison', 'n_Group1', 'n_Group2', ...
    'Median_Group1', 'Median_Group2', 'Mean_Group1', 'Mean_Group2', ...
    'p_Wilcoxon', 'Significance'});

writetable(panelC_tbl, out_file, 'Sheet', 'PanelC_DSI_AR');
fprintf('  Wrote sheet: PanelC_DSI_AR\n');

%% ========================================================================
%  SHEET 4: EFGH Pooled 3-position Wilcoxon rank-sum
%  ========================================================================

res_file = fullfile(data_root, 'population_results', 'batch_results.mat');
fprintf('Loading EFGH batch results: %s\n', res_file);
S = load(res_file, 'results');
results = S.results;

% Group masks
on_ctrl  = [results.is_on] & ~[results.is_ttl];
on_ttl   = [results.is_on] &  [results.is_ttl];
off_ctrl = ~[results.is_on] & ~[results.is_ttl];
off_ttl  = ~[results.is_on] &  [results.is_ttl];

% Constants (matching EFGH script)
STIM_ONSET  = 5001;
DEP_WINDOW  = [STIM_ONSET, STIM_ONSET + 2000 - 1];
HYP_WINDOW  = [STIM_ONSET + 1000, Inf];
DEP_PCTILE  = 99.9;
HYP_PCTILE  = 0.1;
REJECT_THRESH = 0.5;

% Extract amplitude matrices
[pd_on_c_dep, pd_on_c_hyp]   = extract_robust_amps_local(results(on_ctrl),  'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[pd_on_t_dep, pd_on_t_hyp]   = extract_robust_amps_local(results(on_ttl),   'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[pd_off_c_dep, pd_off_c_hyp] = extract_robust_amps_local(results(off_ctrl), 'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[pd_off_t_dep, pd_off_t_hyp] = extract_robust_amps_local(results(off_ttl),  'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

[ort_on_c_dep, ort_on_c_hyp]   = extract_robust_amps_local(results(on_ctrl),  'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[ort_on_t_dep, ort_on_t_hyp]   = extract_robust_amps_local(results(on_ttl),   'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[ort_off_c_dep, ort_off_c_hyp] = extract_robust_amps_local(results(off_ctrl), 'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[ort_off_t_dep, ort_off_t_hyp] = extract_robust_amps_local(results(off_ttl),  'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

% Compute pooled rank-sum for all 8 panels
pool_centers = -4:4;
[pd_on_dep_pp, ~]  = compute_pooled_ranksum_local(pd_on_c_dep,  pd_on_t_dep);
[pd_on_hyp_pp, ~]  = compute_pooled_ranksum_local(pd_on_c_hyp,  pd_on_t_hyp);
[pd_off_dep_pp, ~] = compute_pooled_ranksum_local(pd_off_c_dep, pd_off_t_dep);
[pd_off_hyp_pp, ~] = compute_pooled_ranksum_local(pd_off_c_hyp, pd_off_t_hyp);
[ort_on_dep_pp, ~]  = compute_pooled_ranksum_local(ort_on_c_dep,  ort_on_t_dep);
[ort_on_hyp_pp, ~]  = compute_pooled_ranksum_local(ort_on_c_hyp,  ort_on_t_hyp);
[ort_off_dep_pp, ~] = compute_pooled_ranksum_local(ort_off_c_dep, ort_off_t_dep);
[ort_off_hyp_pp, ~] = compute_pooled_ranksum_local(ort_off_c_hyp, ort_off_t_hyp);

% Build table: 8 panels x 9 pool centers
n_panels = 8;
panel_labels = {'PD T4 Depol', 'PD T5 Depol', 'PD T4 Hyperpol', 'PD T5 Hyperpol', ...
                'Ort T4 Depol', 'Ort T5 Depol', 'Ort T4 Hyperpol', 'Ort T5 Hyperpol'};
all_pp = [pd_on_dep_pp; pd_off_dep_pp; pd_on_hyp_pp; pd_off_hyp_pp; ...
          ort_on_dep_pp; ort_off_dep_pp; ort_on_hyp_pp; ort_off_hyp_pp];

rows = n_panels * 9;
panel_col = cell(rows, 1);
center_col = NaN(rows, 1);
pool_range_col = cell(rows, 1);
p_pool_col = NaN(rows, 1);
sig_pool_col = cell(rows, 1);
n_sig_col = cell(rows, 1);

idx = 0;
for pi = 1:n_panels
    for ci = 1:9
        idx = idx + 1;
        panel_col{idx} = panel_labels{pi};
        center_col(idx) = pool_centers(ci);
        pool_range_col{idx} = sprintf('[%d, %d, %d]', pool_centers(ci), pool_centers(ci)+1, pool_centers(ci)+2);
        p_pool_col(idx) = all_pp(pi, ci);
        if ~isnan(p_pool_col(idx)) && p_pool_col(idx) < 0.05
            if p_pool_col(idx) < 0.001, sig_pool_col{idx} = '***';
            elseif p_pool_col(idx) < 0.01, sig_pool_col{idx} = '**';
            else, sig_pool_col{idx} = '*';
            end
        else
            sig_pool_col{idx} = '';
        end
    end
    n_sig_col{(pi-1)*9 + 1} = sprintf('%d / 9', sum(all_pp(pi,:) < 0.05, 'omitnan'));
end

pool_tbl = table(panel_col, center_col, pool_range_col, p_pool_col, sig_pool_col, ...
    'VariableNames', {'Panel', 'Pool_Center', 'Positions_Pooled', ...
    'p_Wilcoxon_Pooled', 'Significance'});

writetable(pool_tbl, out_file, 'Sheet', 'EFGH_Pooled_Ranksum');
fprintf('  Wrote sheet: EFGH_Pooled_Ranksum\n');

%% ========================================================================
%  SHEET 5: EFGH FWHM rank-sum
%  ========================================================================

positions = -5:5;
fwhm_panels = {'PD T4 Depol', 'PD T5 Depol', 'Ort T4 Depol', 'Ort T5 Depol'};
ctrl_dep_mats = {pd_on_c_dep, pd_off_c_dep, ort_on_c_dep, ort_off_c_dep};
ttl_dep_mats  = {pd_on_t_dep, pd_off_t_dep, ort_on_t_dep, ort_off_t_dep};

fwhm_panel_col = cell(4, 1);
fwhm_ctrl_mean_col = NaN(4, 1);
fwhm_ttl_mean_col  = NaN(4, 1);
fwhm_ctrl_group_col = NaN(4, 1);  % group mean FWHM
fwhm_ttl_group_col  = NaN(4, 1);
n_ctrl_fwhm = NaN(4, 1);
n_ttl_fwhm  = NaN(4, 1);
p_fwhm_col  = NaN(4, 1);
sig_fwhm_col = cell(4, 1);

for fi = 1:4
    fwhm_panel_col{fi} = fwhm_panels{fi};

    % Group mean FWHM
    ctrl_mn = mean(ctrl_dep_mats{fi}, 1, 'omitnan');
    ttl_mn  = mean(ttl_dep_mats{fi},  1, 'omitnan');
    [fw_c, ~, ~] = compute_fwhm_local(positions, ctrl_mn);
    [fw_t, ~, ~] = compute_fwhm_local(positions, ttl_mn);
    fwhm_ctrl_group_col(fi) = fw_c;
    fwhm_ttl_group_col(fi)  = fw_t;

    % Per-cell FWHM
    cv_fw = compute_percell_fwhm_local(positions, ctrl_dep_mats{fi});
    tv_fw = compute_percell_fwhm_local(positions, ttl_dep_mats{fi});
    cv_fw = cv_fw(~isnan(cv_fw)); tv_fw = tv_fw(~isnan(tv_fw));

    n_ctrl_fwhm(fi) = numel(cv_fw);
    n_ttl_fwhm(fi)  = numel(tv_fw);
    fwhm_ctrl_mean_col(fi) = mean(cv_fw);
    fwhm_ttl_mean_col(fi)  = mean(tv_fw);

    if numel(cv_fw) >= 2 && numel(tv_fw) >= 2
        p_fwhm_col(fi) = ranksum(cv_fw, tv_fw);
    end
    if ~isnan(p_fwhm_col(fi))
        if p_fwhm_col(fi) < 0.001, sig_fwhm_col{fi} = '***';
        elseif p_fwhm_col(fi) < 0.01, sig_fwhm_col{fi} = '**';
        elseif p_fwhm_col(fi) < 0.05, sig_fwhm_col{fi} = '*';
        else, sig_fwhm_col{fi} = 'n.s.';
        end
    else
        sig_fwhm_col{fi} = '';
    end
end

fwhm_tbl = table(fwhm_panel_col, fwhm_ctrl_group_col, fwhm_ttl_group_col, ...
    n_ctrl_fwhm, n_ttl_fwhm, fwhm_ctrl_mean_col, fwhm_ttl_mean_col, ...
    p_fwhm_col, sig_fwhm_col, ...
    'VariableNames', {'Panel', 'FWHM_GroupMean_Ctrl', 'FWHM_GroupMean_Tutl', ...
    'n_Ctrl', 'n_Tutl', 'FWHM_PerCell_Mean_Ctrl', 'FWHM_PerCell_Mean_Tutl', ...
    'p_Wilcoxon', 'Significance'});

writetable(fwhm_tbl, out_file, 'Sheet', 'EFGH_FWHM');
fprintf('  Wrote sheet: EFGH_FWHM\n');

%% Done
fprintf('\n=== All statistics exported to: %s ===\n', out_file);


%% =========================================================================
%%                          LOCAL FUNCTIONS
%% =========================================================================

function [stats, tbl] = compute_direction_stats_local(all_cells, ctrl_idx, ttl_idx, pd_aligned_angles)
    n_ctrl = numel(ctrl_idx);
    n_ttl  = numel(ttl_idx);
    n_dirs = 16;

    stats = struct('angle_deg', cell(n_dirs,1), 'n_ctrl', cell(n_dirs,1), ...
        'n_ttl', cell(n_dirs,1), 'median_ctrl', cell(n_dirs,1), ...
        'median_ttl', cell(n_dirs,1), 'mean_ctrl', cell(n_dirs,1), ...
        'mean_ttl', cell(n_dirs,1), 'p_raw', cell(n_dirs,1), ...
        'p_fdr', cell(n_dirs,1), 'sig_fdr', cell(n_dirs,1));

    p_raw_all = NaN(n_dirs, 1);

    for di = 1:n_dirs
        ctrl_vals = NaN(n_ctrl, 1);
        for k = 1:n_ctrl
            pa = all_cells(ctrl_idx(k)).peak_amps;
            if ~isempty(pa) && size(pa, 1) >= di
                ctrl_vals(k) = pa(di, 2);
            end
        end
        ttl_vals = NaN(n_ttl, 1);
        for k = 1:n_ttl
            pa = all_cells(ttl_idx(k)).peak_amps;
            if ~isempty(pa) && size(pa, 1) >= di
                ttl_vals(k) = pa(di, 2);
            end
        end
        ctrl_vals = ctrl_vals(~isnan(ctrl_vals));
        ttl_vals  = ttl_vals(~isnan(ttl_vals));

        stats(di).angle_deg   = pd_aligned_angles(di);
        stats(di).n_ctrl      = numel(ctrl_vals);
        stats(di).n_ttl       = numel(ttl_vals);
        stats(di).median_ctrl = median(ctrl_vals);
        stats(di).median_ttl  = median(ttl_vals);
        stats(di).mean_ctrl   = mean(ctrl_vals);
        stats(di).mean_ttl    = mean(ttl_vals);

        if numel(ctrl_vals) >= 2 && numel(ttl_vals) >= 2
            p_raw_all(di) = ranksum(ctrl_vals, ttl_vals);
        end
        stats(di).p_raw = p_raw_all(di);
    end

    % BH-FDR correction
    valid = ~isnan(p_raw_all);
    p_valid = p_raw_all(valid);
    n_valid = numel(p_valid);
    valid_idx = find(valid);

    [p_sorted, sort_order] = sort(p_valid);
    p_fdr_sorted = p_sorted;
    for i = 1:n_valid
        p_fdr_sorted(i) = p_sorted(i) * n_valid / i;
    end
    for i = n_valid-1:-1:1
        p_fdr_sorted(i) = min(p_fdr_sorted(i), p_fdr_sorted(i+1));
    end
    p_fdr_sorted = min(p_fdr_sorted, 1);

    p_fdr_all = NaN(n_dirs, 1);
    p_fdr_unsorted = NaN(n_valid, 1);
    p_fdr_unsorted(sort_order) = p_fdr_sorted;
    for i = 1:n_valid
        p_fdr_all(valid_idx(i)) = p_fdr_unsorted(i);
    end

    for di = 1:n_dirs
        stats(di).p_fdr   = p_fdr_all(di);
        stats(di).sig_fdr = ~isnan(p_fdr_all(di)) && p_fdr_all(di) < 0.05;
    end

    tbl = stats;  % return same struct
end


function T = build_direction_table(stats)
    n = numel(stats);
    angle_deg = NaN(n, 1);
    n_ctrl = NaN(n, 1); n_ttl = NaN(n, 1);
    median_ctrl = NaN(n, 1); median_ttl = NaN(n, 1);
    mean_ctrl = NaN(n, 1); mean_ttl = NaN(n, 1);
    p_raw = NaN(n, 1); p_fdr = NaN(n, 1);
    sig = cell(n, 1);

    for i = 1:n
        angle_deg(i)   = stats(i).angle_deg;
        n_ctrl(i)      = stats(i).n_ctrl;
        n_ttl(i)       = stats(i).n_ttl;
        median_ctrl(i) = stats(i).median_ctrl;
        median_ttl(i)  = stats(i).median_ttl;
        mean_ctrl(i)   = stats(i).mean_ctrl;
        mean_ttl(i)    = stats(i).mean_ttl;
        p_raw(i)       = stats(i).p_raw;
        p_fdr(i)       = stats(i).p_fdr;
        if stats(i).sig_fdr
            if p_fdr(i) < 0.001, sig{i} = '***';
            elseif p_fdr(i) < 0.01, sig{i} = '**';
            else, sig{i} = '*';
            end
        else
            sig{i} = '';
        end
    end

    T = table(angle_deg, n_ctrl, n_ttl, median_ctrl, median_ttl, ...
        mean_ctrl, mean_ttl, p_raw, p_fdr, sig, ...
        'VariableNames', {'Direction_deg', 'n_Ctrl', 'n_Tutl', ...
        'Median_Ctrl_mV', 'Median_Tutl_mV', 'Mean_Ctrl_mV', 'Mean_Tutl_mV', ...
        'p_Raw', 'p_FDR', 'Significant_FDR005'});
end


function [dep_mat, hyp_mat] = extract_robust_amps_local(results_sub, trace_field, ...
    dep_win, hyp_win, dep_pct, hyp_pct, ~)
    n = numel(results_sub);
    dep_mat = NaN(n, 11); hyp_mat = NaN(n, 11);
    for k = 1:n
        if ~isfield(results_sub(k), trace_field), continue; end
        traces = results_sub(k).(trace_field);
        if isempty(traces), continue; end
        n_samp = size(traces, 2);
        for pos = 1:11
            row = traces(pos, :);
            if all(isnan(row)), continue; end
            de = min(dep_win(2), n_samp);
            if ~isinf(de)
                dep_mat(k, pos) = prctile(row(dep_win(1):de), dep_pct);
            else
                dep_mat(k, pos) = prctile(row(dep_win(1):end), dep_pct);
            end
            he = hyp_win(2); if isinf(he), he = n_samp; end
            he = min(he, n_samp);
            hval = prctile(row(hyp_win(1):he), hyp_pct);
            hyp_mat(k, pos) = min(hval, 0);
        end
    end
end


function [pool_pvals, pool_centers] = compute_pooled_ranksum_local(ctrl_mat, ttl_mat)
    pool_centers = -4:4;
    pool_pvals = NaN(1, 9);
    for p = 1:9
        cols = p:p+2;
        c = ctrl_mat(:, cols); c = c(:); c = c(~isnan(c));
        t = ttl_mat(:, cols);  t = t(:); t = t(~isnan(t));
        if numel(c) >= 2 && numel(t) >= 2
            pool_pvals(p) = ranksum(c, t);
        end
    end
end


function ar = compute_ar_local(d_aligned)
    angles = d_aligned(:, 1);
    resps  = d_aligned(:, 2);
    [~, pd_idx]     = min(abs(angles - pi/2));
    [~, ortho1_idx] = min(abs(angles - 0));
    [~, ortho2_idx] = min(abs(angles - pi));
    pd_resp    = resps(pd_idx);
    ortho_mean = mean([resps(ortho1_idx), resps(ortho2_idx)]);
    if ortho_mean > 0
        ar = pd_resp / ortho_mean;
    else
        ar = NaN;
    end
end


function [fw, left_x, right_x] = compute_fwhm_local(positions, mean_vals)
    fw = NaN; left_x = NaN; right_x = NaN;
    valid = ~isnan(mean_vals);
    if sum(valid) < 3, return; end
    pos_v = positions(valid); val_v = mean_vals(valid);
    [pk, pk_idx] = max(val_v);
    if pk <= 0, return; end
    half_max = pk / 2;
    for j = pk_idx:-1:2
        if val_v(j-1) <= half_max
            frac = (half_max - val_v(j-1)) / (val_v(j) - val_v(j-1));
            left_x = pos_v(j-1) + frac * (pos_v(j) - pos_v(j-1)); break;
        end
    end
    if isnan(left_x), left_x = pos_v(1); end
    for j = pk_idx:numel(val_v)-1
        if val_v(j+1) <= half_max
            frac = (half_max - val_v(j+1)) / (val_v(j) - val_v(j+1));
            right_x = pos_v(j+1) - frac * (pos_v(j+1) - pos_v(j)); break;
        end
    end
    if isnan(right_x), right_x = pos_v(end); end
    fw = right_x - left_x;
end


function fwhm_vec = compute_percell_fwhm_local(positions, amp_mat)
    n = size(amp_mat, 1); fwhm_vec = NaN(n, 1);
    for k = 1:n
        [fw, ~, ~] = compute_fwhm_local(positions, amp_mat(k, :));
        fwhm_vec(k) = fw;
    end
end

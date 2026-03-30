%% Quick diagnostic: check pooled rank-sum p-values for all 8 amplitude panels
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

S = load('/Users/reiserm/Documents/ttl_1DRF/population_results/batch_results.mat', 'results');
results = S.results;

on_ctrl  = [results.is_on] & ~[results.is_ttl];
on_ttl   = [results.is_on] &  [results.is_ttl];
off_ctrl = ~[results.is_on] & ~[results.is_ttl];
off_ttl  = ~[results.is_on] &  [results.is_ttl];

DEP_WINDOW = [5001, 7000]; HYP_WINDOW = [6001, Inf];
DEP_PCTILE = 99.9; HYP_PCTILE = 0.1; REJECT_THRESH = 0.5;

% PD axis
[pd_on_c_dep, pd_on_c_hyp]   = extract_robust_amps_diag(results(on_ctrl),  'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[pd_on_t_dep, pd_on_t_hyp]   = extract_robust_amps_diag(results(on_ttl),   'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[pd_off_c_dep, pd_off_c_hyp] = extract_robust_amps_diag(results(off_ctrl), 'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[pd_off_t_dep, pd_off_t_hyp] = extract_robust_amps_diag(results(off_ttl),  'pd_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

% Ortho axis
[ort_on_c_dep, ort_on_c_hyp]   = extract_robust_amps_diag(results(on_ctrl),  'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[ort_on_t_dep, ort_on_t_hyp]   = extract_robust_amps_diag(results(on_ttl),   'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[ort_off_c_dep, ort_off_c_hyp] = extract_robust_amps_diag(results(off_ctrl), 'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);
[ort_off_t_dep, ort_off_t_hyp] = extract_robust_amps_diag(results(off_ttl),  'ortho_flash_m6_aligned', DEP_WINDOW, HYP_WINDOW, DEP_PCTILE, HYP_PCTILE, REJECT_THRESH);

fprintf('\n=== Pooled rank-sum p-values (asterisks at p<0.05) ===\n');
fprintf('Pool centers: -4    -3    -2    -1     0    +1    +2    +3    +4\n');

panels = {'PD T4 dep', 'PD T5 dep', 'PD T4 hyp', 'PD T5 hyp', ...
           'Ort T4 dep', 'Ort T5 dep', 'Ort T4 hyp', 'Ort T5 hyp'};
ctrl_mats = {pd_on_c_dep, pd_off_c_dep, pd_on_c_hyp, pd_off_c_hyp, ...
             ort_on_c_dep, ort_off_c_dep, ort_on_c_hyp, ort_off_c_hyp};
ttl_mats  = {pd_on_t_dep, pd_off_t_dep, pd_on_t_hyp, pd_off_t_hyp, ...
             ort_on_t_dep, ort_off_t_dep, ort_on_t_hyp, ort_off_t_hyp};

for i = 1:8
    [pp, pc] = compute_pooled_ranksum_diag(ctrl_mats{i}, ttl_mats{i});
    sig = sum(pp < 0.05);
    fprintf('%-12s', panels{i});
    for j = 1:9
        if pp(j) < 0.05
            fprintf(' %.3f*', pp(j));
        else
            fprintf(' %.3f ', pp(j));
        end
    end
    fprintf('  sig=%d\n', sig);
end

% Also print data counts per position for debugging
fprintf('\n=== Non-NaN counts per position (ctrl / ttl) ===\n');
for i = 1:8
    cm = ctrl_mats{i}; tm = ttl_mats{i};
    fprintf('%-12s ctrl:', panels{i});
    for j = 1:11, fprintf(' %d', sum(~isnan(cm(:,j)))); end
    fprintf('  ttl:');
    for j = 1:11, fprintf(' %d', sum(~isnan(tm(:,j)))); end
    fprintf('\n');
end


function [pool_pvals, pool_centers] = compute_pooled_ranksum_diag(ctrl_mat, ttl_mat)
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

function [dep_mat, hyp_mat] = extract_robust_amps_diag(results_sub, trace_field, dep_win, hyp_win, dep_pct, hyp_pct, reject_thresh)
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
            mean_resp = mean(row(min(dep_win(1),n_samp):end), 'omitnan');
            if abs(mean_resp) < reject_thresh, continue; end
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

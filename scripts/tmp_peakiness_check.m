% Quick peakiness comparison: is fast flash RF peakier than slow?
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';
S = load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
results = S.results;
LUT = load('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src/analysis/protocol2/bar_lut.mat');
plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

n_cells = numel(results);
peak_ratio_slow = NaN(n_cells,1); peak_ratio_fast = NaN(n_cells,1);
kurtosis_slow = NaN(n_cells,1); kurtosis_fast = NaN(n_cells,1);
m6_slow_c = NaN(n_cells,1); m6_fast_c = NaN(n_cells,1);

for ci = 1:n_cells
    r = results(ci);
    exp_folder = fullfile(data_root, r.folder);
    orig_dir = pwd;
    [~, ~, Log, ~, ~] = load_protocol2_data(exp_folder); cd(orig_dir);
    f_data = Log.ADC.Volts(1,:); v_data = Log.ADC.Volts(2,:)*10;
    [~, ~, mean_slow, mean_fast] = parse_bar_flash_data(f_data, v_data, 0.5);

    ce = load(fullfile(exp_folder, 'currentExp.mat'), 'pattern_order', 'func_order');
    [ld, lo, lp, lf] = verify_lut_directions(LUT.Tbl, ce.pattern_order, ce.func_order, plot_order);
    bar_data = parse_bar_data(f_data, v_data);
    so.baseline_range=[1000 9000]; so.stim_trim_end=7000; so.percentile=98;
    max_v = compute_bar_sweep_responses(bar_data, plot_order, so);
    pd_info = find_pd_from_lut(max_v, ld, lo, lp, lf, plot_order, LUT.Tbl, 2);

    A_s = zeros(1,11); A_f = zeros(1,11);
    for pi = 1:11
        fp = pd_info.pos_order(pi);
        ts = mean_slow{fp, pd_info.bar_flash_col};
        if ~isempty(ts)
            bl = mean(ts(1:min(5000,numel(ts))));
            A_s(pi) = prctile(ts(5001:min(6551,numel(ts))) - bl, 99.5);
        end
        ts = mean_fast{fp, pd_info.bar_flash_col};
        if ~isempty(ts)
            bl = mean(ts(1:min(2500,numel(ts))));
            A_f(pi) = prctile(ts(2501:min(4051,numel(ts))) - bl, 99.5);
        end
    end

    % Normalize to peak=1
    As = max(A_s,0); Af = max(A_f,0);
    if max(As)>0, As = As/max(As); end
    if max(Af)>0, Af = Af/max(Af); end

    pos_s = As(As>0); pos_f = Af(Af>0);
    if ~isempty(pos_s), peak_ratio_slow(ci) = 1/mean(pos_s); end
    if ~isempty(pos_f), peak_ratio_fast(ci) = 1/mean(pos_f); end
    kurtosis_slow(ci) = kurtosis(As);
    kurtosis_fast(ci) = kurtosis(Af);

    m6_slow_c(ci) = compute_m6_local(max(A_s,0));
    m6_fast_c(ci) = compute_m6_local(max(A_f,0));
end

% Group analysis
is_on = [results.is_on]; is_ttl = [results.is_ttl];
groups = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
masks = {is_on&~is_ttl, is_on&is_ttl, ~is_on&~is_ttl, ~is_on&is_ttl};

fprintf('\n=== Peakiness: Fast vs Slow (normalized profiles) ===\n');
fprintf('Higher peak_ratio = more peaked/concentrated\n\n');
fprintf('%-10s  n  PkRatio_S PkRatio_F  Kurt_S  Kurt_F  Fast_peakier%%\n', 'Group');
for g = 1:4
    m = masks{g}; if sum(m)==0, continue; end
    fp = peak_ratio_fast(m) > peak_ratio_slow(m);
    fprintf('%-10s %2d  %5.2f     %5.2f     %5.2f   %5.2f    %3.0f%%\n', ...
        groups{g}, sum(m), mean(peak_ratio_slow(m)), mean(peak_ratio_fast(m)), ...
        mean(kurtosis_slow(m)), mean(kurtosis_fast(m)), 100*mean(fp));
end
fp_all = peak_ratio_fast > peak_ratio_slow;
fprintf('%-10s %2d  %5.2f     %5.2f     %5.2f   %5.2f    %3.0f%%\n', ...
    'ALL', n_cells, mean(peak_ratio_slow), mean(peak_ratio_fast), ...
    mean(kurtosis_slow), mean(kurtosis_fast), 100*mean(fp_all));

[p_pk,~] = signrank(peak_ratio_fast, peak_ratio_slow);
[p_ku,~] = signrank(kurtosis_fast, kurtosis_slow);
fprintf('\nSigned-rank p: peak_ratio p=%.4f, kurtosis p=%.4f\n', p_pk, p_ku);

% Directional bias
delta = m6_fast_c - m6_slow_c;
fprintf('\n=== Directional bias: fast_M6 - slow_M6 ===\n');
fprintf('  Mean shift: %+.2f (positive = fast shifts toward higher pos = PD)\n', mean(delta));
fprintf('  Median shift: %+.2f\n', median(delta));
fprintf('  Positive (fast > slow): %d/%d\n', sum(delta>0), n_cells);
fprintf('  Negative (fast < slow): %d/%d\n', sum(delta<0), n_cells);
fprintf('  Equal: %d/%d\n', sum(abs(delta)<0.01), n_cells);
[p_dir,~] = signrank(delta);
fprintf('  Signed-rank p (shift != 0): %.4f\n', p_dir);

fprintf('\nBy group:\n');
for g = 1:4
    m = masks{g}; if sum(m)==0, continue; end
    d = delta(m);
    fprintf('  %-10s: mean=%+.2f, %d/%d positive, %d/%d negative\n', ...
        groups{g}, mean(d), sum(d>0), sum(m), sum(d<0), sum(m));
end

fprintf('\n=== TTL-specific: is fast M6 systematically different for TTL cells? ===\n');
ttl_mask = is_ttl;
ctrl_mask = ~is_ttl;
fprintf('  TTL cells (n=%d): mean shift = %+.2f\n', sum(ttl_mask), mean(delta(ttl_mask)));
fprintf('  Ctrl cells (n=%d): mean shift = %+.2f\n', sum(ctrl_mask), mean(delta(ctrl_mask)));
[p_grp, ~] = ranksum(delta(ttl_mask), delta(ctrl_mask));
fprintf('  Rank-sum p (TTL vs ctrl shift): %.4f\n', p_grp);

function c = compute_m6_local(A)
    [pv, pp] = max(A); ta = sum(A);
    if pv<=0||ta<=0, c=6; return; end
    target=0.68*ta; l=pp; r=pp; cum=A(pp);
    while cum<target && (l>1||r<11)
        cl=(l>1); cr=(r<11);
        if cl&&cr
            if A(l-1)>=A(r+1), l=l-1; cum=cum+A(l); else, r=r+1; cum=cum+A(r); end
        elseif cl, l=l-1; cum=cum+A(l); else, r=r+1; cum=cum+A(r); end
    end
    bump=l:r; Ab=A(bump); c=sum(Ab.*bump)/sum(Ab);
end

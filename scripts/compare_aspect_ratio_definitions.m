% COMPARE_ASPECT_RATIO_DEFINITIONS  Phase 0C: Compare AR_pd vs AR_pdnd.
%
%   Compares two aspect ratio definitions on the combined dataset (n=48):
%     AR_pd   = PD / mean(ortho_CW, ortho_CCW)           — current definition
%     AR_pdnd = (PD + ND) / (ortho_CW + ortho_CCW)       — full-axis definition
%
%   Both computed from max_v_aligned (16×2), where PD = 90° (row 5),
%   ND = 270° (row 13), ortho_CW = 0° (row 1), ortho_CCW = 180° (row 9).
%
%   Outputs:
%     diag_aspect_ratio_comparison.png — 1×2 figure: AR_pd and AR_pdnd box+dot
%     diag_aspect_ratio_stats.txt     — group stats + rank-sum p-values
%
%   Usage:
%     run('scripts/compare_aspect_ratio_definitions.m')

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
fprintf('  Late batch: %d cells\n', numel(late));

S_early = load(fullfile(data_root, 'pre-bar-flash', 'population_results', ...
    'batch_results_pre_bf.mat'), 'results');
early = S_early.results;
fprintf('  Early batch: %d cells\n', numel(early));

%% Compute both AR definitions for all cells
n_late  = numel(late);
n_early = numel(early);
n_total = n_late + n_early;

ar_pd    = NaN(n_total, 1);
ar_pdnd  = NaN(n_total, 1);
is_on    = false(n_total, 1);
is_ttl   = false(n_total, 1);
groups   = cell(n_total, 1);

for k = 1:n_total
    if k <= n_late
        d = late(k).max_v_aligned;
        is_on(k)  = late(k).is_on;
        is_ttl(k) = late(k).is_ttl;
        groups{k}  = late(k).group;
    else
        ei = k - n_late;
        d = early(ei).max_v_aligned;
        is_on(k)  = early(ei).is_on;
        is_ttl(k) = early(ei).is_ttl;
        groups{k}  = early(ei).group;
    end

    angles = d(:, 1);
    resps  = d(:, 2);

    % Identify key directions (PD at π/2, ND at 3π/2, ortho at 0 and π)
    [~, pd_idx]     = min(abs(angles - pi/2));
    [~, nd_idx]     = min(abs(angles - 3*pi/2));
    [~, ortho1_idx] = min(abs(angles - 0));
    [~, ortho2_idx] = min(abs(angles - pi));

    pd_resp    = resps(pd_idx);
    nd_resp    = resps(nd_idx);
    ortho1     = resps(ortho1_idx);
    ortho2     = resps(ortho2_idx);
    ortho_mean = mean([ortho1, ortho2]);
    ortho_sum  = ortho1 + ortho2;

    % AR_pd = PD / mean(ortho)
    if ortho_mean > 0
        ar_pd(k) = pd_resp / ortho_mean;
    end

    % AR_pdnd = (PD + ND) / (ortho_CW + ortho_CCW)
    if ortho_sum > 0
        ar_pdnd(k) = (pd_resp + nd_resp) / ortho_sum;
    end
end

%% Group indices
on_ctrl  = is_on & ~is_ttl;
on_ttl   = is_on &  is_ttl;
off_ctrl = ~is_on & ~is_ttl;
off_ttl  = ~is_on &  is_ttl;

group_names = {'ON ctrl', 'ON TTL', 'OFF ctrl', 'OFF TTL'};
group_masks = {on_ctrl, on_ttl, off_ctrl, off_ttl};
group_colors = {[0.3 0.3 0.3], [0.85 0.2 0.2], [0.3 0.3 0.3], [0.85 0.2 0.2]};

%% Print summary stats
fprintf('\n=== Aspect Ratio Summary ===\n');
fprintf('%-12s  %5s  %8s  %8s  %8s  %8s\n', ...
    'Group', 'n', 'AR_pd', '±SEM', 'AR_pdnd', '±SEM');
fprintf('%s\n', repmat('-', 1, 60));

for g = 1:4
    mask = group_masks{g};
    n_g = sum(mask);
    fprintf('%-12s  %5d  %8.3f  %8.3f  %8.3f  %8.3f\n', ...
        group_names{g}, n_g, ...
        mean(ar_pd(mask), 'omitnan'), std(ar_pd(mask), 'omitnan')/sqrt(n_g), ...
        mean(ar_pdnd(mask), 'omitnan'), std(ar_pdnd(mask), 'omitnan')/sqrt(n_g));
end

%% Statistical tests
fprintf('\n=== Wilcoxon Rank-Sum Tests ===\n');

comparisons = {
    'ON ctrl vs TTL',   on_ctrl, on_ttl;
    'OFF ctrl vs TTL',  off_ctrl, off_ttl;
    'ON vs OFF (ctrl)', on_ctrl, off_ctrl;
    'ON vs OFF (TTL)',  on_ttl, off_ttl
};

p_values = struct('comparison', {}, 'ar_pd_p', {}, 'ar_pdnd_p', {}, ...
    'ar_pd_rbc', {}, 'ar_pdnd_rbc', {});

fprintf('%-22s  %10s  %10s  %10s  %10s\n', ...
    'Comparison', 'AR_pd p', 'AR_pd rbc', 'AR_pdnd p', 'AR_pdnd rbc');
fprintf('%s\n', repmat('-', 1, 70));

for c = 1:size(comparisons, 1)
    m1 = comparisons{c, 2};
    m2 = comparisons{c, 3};

    v1_pd = ar_pd(m1);   v2_pd = ar_pd(m2);
    v1_pn = ar_pdnd(m1);  v2_pn = ar_pdnd(m2);

    % Remove NaN
    v1_pd = v1_pd(~isnan(v1_pd)); v2_pd = v2_pd(~isnan(v2_pd));
    v1_pn = v1_pn(~isnan(v1_pn)); v2_pn = v2_pn(~isnan(v2_pn));

    if numel(v1_pd) >= 2 && numel(v2_pd) >= 2
        [p_pd, ~, stats_pd] = ranksum(v1_pd, v2_pd);
        n1 = numel(v1_pd); n2 = numel(v2_pd);
        U_pd = stats_pd.ranksum - n1*(n1+1)/2;
        rbc_pd = 1 - 2*U_pd / (n1*n2);
    else
        p_pd = NaN; rbc_pd = NaN;
    end

    if numel(v1_pn) >= 2 && numel(v2_pn) >= 2
        [p_pn, ~, stats_pn] = ranksum(v1_pn, v2_pn);
        n1 = numel(v1_pn); n2 = numel(v2_pn);
        U_pn = stats_pn.ranksum - n1*(n1+1)/2;
        rbc_pn = 1 - 2*U_pn / (n1*n2);
    else
        p_pn = NaN; rbc_pn = NaN;
    end

    fprintf('%-22s  %10.4f  %10.3f  %10.4f  %10.3f\n', ...
        comparisons{c, 1}, p_pd, rbc_pd, p_pn, rbc_pn);

    p_values(c).comparison = comparisons{c, 1};
    p_values(c).ar_pd_p    = p_pd;
    p_values(c).ar_pdnd_p  = p_pn;
    p_values(c).ar_pd_rbc  = rbc_pd;
    p_values(c).ar_pdnd_rbc = rbc_pn;
end

%% Plot: 1×2 figure — AR_pd (left), AR_pdnd (right)
fig = figure('Position', [50 50 900 400], 'Color', 'w');
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ar_labels = {'AR_{pd}', 'AR_{pdnd}'};
ar_data   = {ar_pd, ar_pdnd};

for panel = 1:2
    ax = nexttile(t);
    hold(ax, 'on');

    data = ar_data{panel};
    x_positions = [1, 2, 3.5, 4.5];  % spacing between ON/OFF groups
    jitter_width = 0.15;

    for g = 1:4
        mask = group_masks{g};
        vals = data(mask);
        vals = vals(~isnan(vals));

        if isempty(vals), continue; end

        x = x_positions(g);

        % Box plot manually
        q1 = prctile(vals, 25);
        q3 = prctile(vals, 75);
        med = median(vals);
        iqr_val = q3 - q1;
        whisker_lo = max(min(vals), q1 - 1.5*iqr_val);
        whisker_hi = min(max(vals), q3 + 1.5*iqr_val);

        % Box
        fill(ax, [x-0.2, x+0.2, x+0.2, x-0.2], [q1 q1 q3 q3], ...
            group_colors{g}, 'FaceAlpha', 0.2, 'EdgeColor', group_colors{g});
        plot(ax, [x-0.2, x+0.2], [med med], 'Color', group_colors{g}, 'LineWidth', 2);

        % Whiskers
        plot(ax, [x x], [whisker_lo q1], 'Color', group_colors{g}, 'LineWidth', 0.5);
        plot(ax, [x x], [q3 whisker_hi], 'Color', group_colors{g}, 'LineWidth', 0.5);

        % Jittered dots
        x_jitter = x + jitter_width * (rand(size(vals)) - 0.5);
        scatter(ax, x_jitter, vals, 20, group_colors{g}, 'filled', ...
            'MarkerFaceAlpha', 0.6);
    end

    % Add p-value brackets for ctrl vs TTL
    for pair = 1:2
        if pair == 1
            x1 = x_positions(1); x2 = x_positions(2);
            ar_names = {'pd', 'pdnd'};
            p_val = p_values(1).(sprintf('ar_%s_p', ar_names{panel}));
        else
            x1 = x_positions(3); x2 = x_positions(4);
            ar_names = {'pd', 'pdnd'};
            p_val = p_values(2).(sprintf('ar_%s_p', ar_names{panel}));
        end

        if ~isnan(p_val)
            y_top = max(data(~isnan(data))) * 1.05;
            plot(ax, [x1 x2], [y_top y_top], 'k-', 'LineWidth', 0.5);
            if p_val < 0.001
                p_str = 'p<0.001';
            else
                p_str = sprintf('p=%.3f', p_val);
            end
            text(ax, (x1+x2)/2, y_top * 1.03, p_str, ...
                'HorizontalAlignment', 'center', 'FontSize', 7);
        end
    end

    set(ax, 'XTick', x_positions, 'XTickLabel', group_names, ...
        'FontSize', 8, 'XTickLabelRotation', 30);
    ylabel(ax, ar_labels{panel}, 'FontSize', 10);
    title(ax, sprintf('%s (n=%d)', ar_labels{panel}, sum(~isnan(data))), 'FontSize', 10);

    % Reference line at 1 (circular = no DS)
    yline(ax, 1, ':', 'Color', [0.5 0.5 0.5]);
end

title(t, 'Aspect Ratio Definition Comparison — Combined Dataset', 'FontSize', 12);

fname = 'diag_aspect_ratio_comparison.png';
exportgraphics(fig, fullfile(preview_dir, fname), 'Resolution', 150);
fprintf('\nSaved: %s\n', fname);
close(fig);

%% Save stats to text file
fid = fopen(fullfile(preview_dir, 'diag_aspect_ratio_stats.txt'), 'w');
fprintf(fid, 'Aspect Ratio Definition Comparison\n');
fprintf(fid, '===================================\n\n');

fprintf(fid, 'Definitions:\n');
fprintf(fid, '  AR_pd   = PD / mean(ortho_CW, ortho_CCW)\n');
fprintf(fid, '  AR_pdnd = (PD + ND) / (ortho_CW + ortho_CCW)\n\n');

fprintf(fid, '%-12s  %5s  %12s  %12s\n', 'Group', 'n', 'AR_pd', 'AR_pdnd');
fprintf(fid, '%s\n', repmat('-', 1, 50));
for g = 1:4
    mask = group_masks{g};
    n_g = sum(mask);
    fprintf(fid, '%-12s  %5d  %5.2f ± %4.2f  %5.2f ± %4.2f\n', ...
        group_names{g}, n_g, ...
        mean(ar_pd(mask), 'omitnan'), std(ar_pd(mask), 'omitnan')/sqrt(n_g), ...
        mean(ar_pdnd(mask), 'omitnan'), std(ar_pdnd(mask), 'omitnan')/sqrt(n_g));
end

fprintf(fid, '\nWilcoxon Rank-Sum Tests:\n');
fprintf(fid, '%-22s  %12s  %12s\n', 'Comparison', 'AR_pd p', 'AR_pdnd p');
fprintf(fid, '%s\n', repmat('-', 1, 50));
for c = 1:numel(p_values)
    fprintf(fid, '%-22s  %12.4f  %12.4f\n', ...
        p_values(c).comparison, p_values(c).ar_pd_p, p_values(c).ar_pdnd_p);
end

fprintf(fid, '\nEffect sizes (rank-biserial correlation):\n');
fprintf(fid, '%-22s  %12s  %12s\n', 'Comparison', 'AR_pd rbc', 'AR_pdnd rbc');
fprintf(fid, '%s\n', repmat('-', 1, 50));
for c = 1:numel(p_values)
    fprintf(fid, '%-22s  %12.3f  %12.3f\n', ...
        p_values(c).comparison, p_values(c).ar_pd_rbc, p_values(c).ar_pdnd_rbc);
end

fclose(fid);
fprintf('Saved: diag_aspect_ratio_stats.txt\n');

fprintf('\n=== Phase 0C Complete ===\n');

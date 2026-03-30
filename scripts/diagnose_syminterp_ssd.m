% DIAGNOSE_SYMINTERP_SSD  Inspect the SSD landscape for a flagged cell.
%
%   Loads a specific flagged cell's tuning data, computes the SSD at every
%   candidate axis angle, and plots the landscape to understand why SymInterp
%   picks an unexpected axis.

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF';

%% Load the flagged cell: 2025_11_05 OFF TTL (cell #39 or #41 from table)
S_late = load(fullfile(data_root, 'population_results', 'batch_results.mat'), 'results');
late = S_late.results;

% Find the 2025_11_05 OFF TTL cells
for k = 1:numel(late)
    if strcmp(late(k).date_str, '2025_11_05') && ~late(k).is_on && late(k).is_ttl
        fprintf('Found: cell %d — %s, ON=%d, TTL=%d, group=%s\n', ...
            k, late(k).date_str, late(k).is_on, late(k).is_ttl, late(k).group);

        d = late(k).max_v_aligned;
        angles_16    = d(:, 1);
        responses_16 = d(:, 2);

        % Print the raw data
        fprintf('\n  Raw 16-point tuning curve:\n');
        fprintf('  %8s  %10s\n', 'Deg', 'Response');
        for j = 1:16
            fprintf('  %8.1f  %10.2f\n', rad2deg(angles_16(j)), responses_16(j));
        end

        %% Interpolate to fine grid
        angles_wrap = [angles_16; angles_16(1) + 2*pi];
        resp_wrap   = [responses_16; responses_16(1)];
        fine_angles = linspace(0, 2*pi, 145)';
        fine_angles = fine_angles(1:end-1);
        fine_resp   = interp1(angles_wrap, resp_wrap, fine_angles, 'linear');

        %% Compute vector sum PD
        x = sum(responses_16 .* cos(angles_16));
        y = sum(responses_16 .* sin(angles_16));
        vec_pd = atan2(y, x);
        if vec_pd < 0, vec_pd = vec_pd + 2*pi; end
        fprintf('\n  Vector sum PD: %.1f deg\n', rad2deg(vec_pd));

        %% Compute SSD at every fine-grid candidate axis (unconstrained)
        n_fine = numel(fine_angles);
        ssd_all     = NaN(n_fine, 1);
        ssd_weighted = NaN(n_fine, 1);

        fine_wrap = [fine_angles; fine_angles(1) + 2*pi];
        fresp_wrap = [fine_resp; fine_resp(1)];

        for c = 1:n_fine
            ax = fine_angles(c);
            mirror = mod(2*ax - fine_angles, 2*pi);
            mirror_resp = interp1(fine_wrap, fresp_wrap, mirror, 'linear');

            % Unweighted SSD (current method)
            ssd_all(c) = sum((fine_resp - mirror_resp).^2);

            % Response-weighted SSD: weight by mean of pair
            w = (fine_resp + mirror_resp) / 2;
            w = max(w, 0);  % no negative weights
            ssd_weighted(c) = sum(w .* (fine_resp - mirror_resp).^2);
        end

        %% Find minima
        % Constrained range: ±45° around vec_pd
        cand_dist = abs(angle(exp(1i*fine_angles) ./ exp(1i*vec_pd)));
        cand_dist_anti = abs(angle(exp(1i*fine_angles) ./ exp(1i*mod(vec_pd+pi, 2*pi))));
        in_window = (cand_dist <= pi/4) | (cand_dist_anti <= pi/4);

        % Current method: unweighted, constrained
        ssd_constrained = ssd_all;
        ssd_constrained(~in_window) = Inf;
        [~, idx_uw] = min(ssd_constrained);

        % Proposed: weighted, constrained
        ssd_w_constrained = ssd_weighted;
        ssd_w_constrained(~in_window) = Inf;
        [~, idx_w] = min(ssd_w_constrained);

        % Also unconstrained minima
        [~, idx_uw_free] = min(ssd_all);
        [~, idx_w_free] = min(ssd_weighted);

        fprintf('\n  SSD minima:\n');
        fprintf('    Unweighted unconstrained: %.1f deg (SSD=%.1f)\n', ...
            rad2deg(fine_angles(idx_uw_free)), ssd_all(idx_uw_free));
        fprintf('    Unweighted constrained:   %.1f deg (SSD=%.1f)\n', ...
            rad2deg(fine_angles(idx_uw)), ssd_all(idx_uw));
        fprintf('    Weighted unconstrained:   %.1f deg (SSD=%.1f)\n', ...
            rad2deg(fine_angles(idx_w_free)), ssd_weighted(idx_w_free));
        fprintf('    Weighted constrained:     %.1f deg (SSD=%.1f)\n', ...
            rad2deg(fine_angles(idx_w)), ssd_weighted(idx_w));

        % SSD at 90° for comparison
        [~, idx_90] = min(abs(fine_angles - pi/2));
        fprintf('    Unweighted at 90°:        SSD=%.1f\n', ssd_all(idx_90));
        fprintf('    Weighted at 90°:          SSD=%.1f\n', ssd_weighted(idx_90));

        %% Plot SSD landscape
        fig = figure('Position', [50 100 1200 800]);
        tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        % Panel 1: Tuning curve
        nexttile; hold on;
        polarplot([angles_16; angles_16(1)], [responses_16; responses_16(1)], ...
            '-o', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'MarkerFaceColor', [0.5 0.5 0.5]);
        title(sprintf('%s — OFF TTL tuning curve', late(k).date_str));

        % Panel 2: Unweighted SSD landscape
        nexttile; hold on;
        deg_axis = rad2deg(fine_angles);
        plot(deg_axis, ssd_all, 'b-', 'LineWidth', 1.5);
        xline(rad2deg(vec_pd), 'g--', 'LineWidth', 1.5);
        xline(rad2deg(fine_angles(idx_uw)), 'r-', 'LineWidth', 2);
        xline(90, 'k:', 'LineWidth', 1);
        % Shade constraint window
        yl = ylim;
        fill_x = [rad2deg(vec_pd)-45, rad2deg(vec_pd)+45, rad2deg(vec_pd)+45, rad2deg(vec_pd)-45];
        fill_y = [yl(1) yl(1) yl(2) yl(2)];
        patch(fill_x, fill_y, [0.9 0.95 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        % Also shade antipodal window
        fill_x2 = fill_x + 180;
        patch(fill_x2, fill_y, [0.9 0.95 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        plot(deg_axis, ssd_all, 'b-', 'LineWidth', 1.5);  % replot on top
        xline(rad2deg(fine_angles(idx_uw)), 'r-', 'LineWidth', 2);
        xlabel('Candidate axis (deg)');
        ylabel('SSD (unweighted)');
        title(sprintf('Unweighted SSD — min at %.1f°', rad2deg(fine_angles(idx_uw))));
        legend('SSD', 'Vec PD', sprintf('Min (%.1f°)', rad2deg(fine_angles(idx_uw))), ...
            '90°', 'Location', 'best');
        xlim([0 360]);

        % Panel 3: Weighted SSD landscape
        nexttile; hold on;
        patch(fill_x, [0 0 max(ssd_weighted)*1.1 max(ssd_weighted)*1.1], ...
            [0.9 0.95 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        patch(fill_x2, [0 0 max(ssd_weighted)*1.1 max(ssd_weighted)*1.1], ...
            [0.9 0.95 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        plot(deg_axis, ssd_weighted, 'b-', 'LineWidth', 1.5);
        xline(rad2deg(vec_pd), 'g--', 'LineWidth', 1.5);
        xline(rad2deg(fine_angles(idx_w)), 'r-', 'LineWidth', 2);
        xline(90, 'k:', 'LineWidth', 1);
        xlabel('Candidate axis (deg)');
        ylabel('SSD (response-weighted)');
        title(sprintf('Weighted SSD — min at %.1f°', rad2deg(fine_angles(idx_w))));
        xlim([0 360]);

        % Panel 4: Overlay of both (normalized)
        nexttile; hold on;
        ssd_norm = ssd_all / max(ssd_all);
        ssd_w_norm = ssd_weighted / max(ssd_weighted);
        plot(deg_axis, ssd_norm, 'b-', 'LineWidth', 1.5);
        plot(deg_axis, ssd_w_norm, 'r-', 'LineWidth', 1.5);
        xline(90, 'k:', 'LineWidth', 1);
        xline(rad2deg(fine_angles(idx_uw)), 'b--', 'LineWidth', 1.5);
        xline(rad2deg(fine_angles(idx_w)), 'r--', 'LineWidth', 1.5);
        xlabel('Candidate axis (deg)');
        ylabel('Normalized SSD');
        title('Normalized: unweighted (blue) vs weighted (red)');
        legend('Unweighted', 'Weighted', '90°', ...
            sprintf('UW min (%.1f°)', rad2deg(fine_angles(idx_uw))), ...
            sprintf('W min (%.1f°)', rad2deg(fine_angles(idx_w))), ...
            'Location', 'best');
        xlim([0 360]);

        preview_dir = fullfile(data_root, 'figure_previews');
        exportgraphics(fig, fullfile(preview_dir, ...
            sprintf('diag_syminterp_ssd_%s.png', late(k).date_str)), ...
            'Resolution', 150);
        close(fig);

        fprintf('\n  Diagnostic plot saved.\n\n');
        break;  % just do the first matching cell
    end
end

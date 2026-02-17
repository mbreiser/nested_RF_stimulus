function fig = plot_voltage_histograms(v_data, f_data, median_v, title_str)
% PLOT_VOLTAGE_HISTOGRAMS  Voltage distribution: full recording vs stimulus periods.
%
%   FIG = PLOT_VOLTAGE_HISTOGRAMS(V_DATA, F_DATA, MEDIAN_V, TITLE_STR)
%   creates a two-panel figure showing the distribution of membrane voltage
%   for the entire recording (top) and during stimulus presentation only
%   (bottom). This helps visualise what the baseline calculations are
%   based on.
%
%   INPUTS:
%     v_data    - 1xN voltage trace in mV (already scaled)
%     f_data    - 1xN frame signal; f_data ~= 0 marks stimulus periods
%     median_v  - Median voltage of the full recording (mV)
%     title_str - Figure title string
%
%   OUTPUT:
%     fig - Figure handle
%
%   FIGURE LAYOUT:
%     Top panel:  Histogram of the full recording voltage with median
%                 marked by a red dashed line.
%     Bottom panel: Histogram of voltage during stimulus periods only
%                 (samples where f_data ~= 0), with both the global
%                 median (red) and the stimulus-period mean (green)
%                 marked.
%
%   PURPOSE:
%     The analysis pipeline uses the global median voltage as a baseline
%     reference for bar sweep responses, and the local pre-flash mean for
%     bar flash responses. This figure shows the full distribution so you
%     can judge whether the baseline is representative and whether the
%     recording is stable.
%
%   See also ANALYZE_SINGLE_EXPERIMENT_MR

    %% Compute stimulus-period statistics
    stim_mask = f_data ~= 0;
    v_stim = v_data(stim_mask);
    mean_v_stim = mean(v_stim);

    %% Create figure
    fig = figure('Name', 'Voltage Distributions', ...
        'Position', [100 200 800 600]);
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    %% Top panel: Full recording
    ax1 = nexttile;
    histogram(v_data, 200, 'FaceColor', [0.5 0.5 0.7], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold on;
    xline(median_v, 'r--', sprintf('median = %.1f mV', median_v), ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'top', ...
        'LabelHorizontalAlignment', 'left', 'FontSize', 10);
    xlabel('Vm (mV)');
    ylabel('Count');
    title('Full Recording');
    box off;
    ax1.TickDir = 'out';
    ax1.FontSize = 11;

    %% Bottom panel: Stimulus periods only
    ax2 = nexttile;
    histogram(v_stim, 200, 'FaceColor', [0.5 0.7 0.5], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold on;
    xline(median_v, 'r--', sprintf('median = %.1f mV', median_v), ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'top', ...
        'LabelHorizontalAlignment', 'left', 'FontSize', 10);
    xline(mean_v_stim, '--', sprintf('stim mean = %.1f mV', mean_v_stim), ...
        'Color', [0.2 0.7 0.2], 'LineWidth', 1.5, ...
        'LabelVerticalAlignment', 'top', ...
        'LabelHorizontalAlignment', 'right', 'FontSize', 10);
    xlabel('Vm (mV)');
    ylabel('Count');
    title(sprintf('During Stimulus Only (%.1f%% of recording)', ...
        100 * sum(stim_mask) / numel(stim_mask)));
    box off;
    ax2.TickDir = 'out';
    ax2.FontSize = 11;

    %% Link x-axes for comparison
    linkaxes([ax1, ax2], 'x');

    sgtitle(title_str, 'FontSize', 14);

end

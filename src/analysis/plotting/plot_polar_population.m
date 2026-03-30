function fig = plot_polar_population(aligned_data_ctrl, aligned_data_ttl, title_str, opts)
% PLOT_POLAR_POPULATION  PD-aligned polar tuning curve with shaded spread.
%
%   FIG = PLOT_POLAR_POPULATION(ALIGNED_DATA_CTRL, ALIGNED_DATA_TTL, TITLE_STR)
%   creates a polar plot showing the population directional tuning curve
%   for two groups (control and ttl). Each cell's tuning data has been
%   rotated so that the preferred direction aligns to pi/2 (90 degrees,
%   pointing up). The center line and shaded spread band are shown.
%
%   FIG = PLOT_POLAR_POPULATION(..., OPTS) uses the options structure to
%   override default parameters.
%
%   INPUTS:
%     aligned_data_ctrl - Cell array of 16x2 matrices, one per control cell.
%                         Each matrix has columns [aligned_angle_rad, response].
%                         PD is at pi/2 after alignment by find_PD_and_order_idx.
%     aligned_data_ttl  - Cell array of 16x2 matrices, one per ttl cell.
%                         Same format as aligned_data_ctrl.
%     title_str         - Figure title string
%     opts              - (Optional) structure with fields:
%                           .stat_method - 'median_mad' (default) or 'mean_sem'
%                           .group1_label - label for group 1 (default: 'control')
%                           .group2_label - label for group 2 (default: 'ttl')
%                           .group1_line_color - [r g b] line color for group 1
%                           .group1_fill_color - [r g b] fill color for group 1
%                           .group2_line_color - [r g b] line color for group 2
%                           .group2_fill_color - [r g b] fill color for group 2
%
%   OUTPUT:
%     fig - Figure handle
%
%   FIGURE LAYOUT:
%     Polar axes showing directional tuning. Control group in black with
%     gray shading. TTL group in red with light red shading. The shaded
%     band represents MAD (default) or SEM depending on stat_method.
%     Legend shows group names and sample sizes.
%
%   APPROACH:
%     Uses a dual-axis technique: polar axes provide the grid and ticks,
%     with a transparent Cartesian axes overlay for both the data lines
%     and spread patches (MATLAB's polaraxes do not natively support patch
%     objects). Drawing both on the same Cartesian axes guarantees perfect
%     alignment between the center line and shaded band.
%
%   See also FIND_PD_AND_ORDER_IDX, POLAR_MEAN_BY_STRAIN,
%            BATCH_ANALYZE_1DRF
% ________________________________________________________________________

    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'stat_method'), opts.stat_method = 'median_mad'; end
    if ~isfield(opts, 'group1_label'), opts.group1_label = 'control'; end
    if ~isfield(opts, 'group2_label'), opts.group2_label = 'ttl'; end

    % Compute group statistics
    stats_ctrl = compute_group_stats(aligned_data_ctrl, opts.stat_method);
    stats_ttl  = compute_group_stats(aligned_data_ttl, opts.stat_method);

    % Get canonical theta from whichever group has data
    if ~isempty(stats_ctrl.center)
        theta = stats_ctrl.theta;
    elseif ~isempty(stats_ttl.center)
        theta = stats_ttl.theta;
    else
        error('No valid aligned data in either group.');
    end

    % Colours — use custom if provided, otherwise defaults
    if isfield(opts, 'group1_line_color')
        col_ctrl_line = opts.group1_line_color;
    else
        col_ctrl_line = [0 0 0];            % black
    end
    if isfield(opts, 'group1_fill_color')
        col_ctrl_fill = opts.group1_fill_color;
    else
        col_ctrl_fill = [0.80 0.80 0.80];   % light gray
    end
    if isfield(opts, 'group2_line_color')
        col_ttl_line = opts.group2_line_color;
    else
        col_ttl_line  = [1 0 0];            % red
    end
    if isfield(opts, 'group2_fill_color')
        col_ttl_fill = opts.group2_fill_color;
    else
        col_ttl_fill  = [1 0.70 0.70];      % light red
    end

    % Create figure with polar axes
    fig = figure('Name', title_str);
    axPolar = polaraxes;
    hold(axPolar, 'on');

    % Set r-limits based on data (center + spread)
    rmax = estimate_rmax(stats_ctrl, stats_ttl);
    if ~isnan(rmax) && rmax > 0
        rlim(axPolar, [0 rmax * 1.1]);
    end

    % Create overlay Cartesian axes for shaded patches
    axFill = axes('Position', axPolar.Position, 'Color', 'none', ...
                  'XColor', 'none', 'YColor', 'none', 'HitTest', 'off');
    axis(axFill, 'equal');
    hold(axFill, 'on');

    % Initial sync of Cartesian overlay limits to polar r-limits
    rL = rlim(axPolar);
    set(axFill, 'XLim', [-rL(2) rL(2)], 'YLim', [-rL(2) rL(2)]);

    % Plot control group
    hCtrl = gobjects(1, 1);
    if ~isempty(stats_ctrl.center)
        hCtrl = plot_with_shade(axFill, theta, ...
            stats_ctrl.center, stats_ctrl.spread, ...
            col_ctrl_line, col_ctrl_fill, 0.40, 2);
    end

    % Plot ttl group
    hTtl = gobjects(1, 1);
    if ~isempty(stats_ttl.center)
        hTtl = plot_with_shade(axFill, theta, ...
            stats_ttl.center, stats_ttl.spread, ...
            col_ttl_line, col_ttl_fill, 0.40, 2);
    end

    % Re-sync Cartesian overlay limits after plotting (polar rlim may
    % have auto-adjusted), ensuring patch and line positions match exactly
    rL = rlim(axPolar);
    set(axFill, 'XLim', [-rL(2) rL(2)], 'YLim', [-rL(2) rL(2)]);
    set(axFill, 'Position', axPolar.Position);

    % Keep polar axes underneath for grid/ticks, overlay on top for data
    axPolar.Color = 'none';
    uistack(axFill, 'top');

    % Legend on the Cartesian overlay (where the visible lines are)
    hs = []; labs = {};
    if ~isempty(stats_ctrl.center)
        hs(end+1)   = hCtrl;
        labs{end+1}  = sprintf('%s (n=%d)', opts.group1_label, stats_ctrl.n);
    end
    if ~isempty(stats_ttl.center)
        hs(end+1)   = hTtl;
        labs{end+1}  = sprintf('%s (n=%d)', opts.group2_label, stats_ttl.n);
    end
    if ~isempty(hs)
        legend(axFill, hs, labs, 'Location', 'best');
    end

    axPolar.FontSize = 12;
    title(axPolar, title_str);

end


%% ========================= Local Functions ============================

function stats = compute_group_stats(aligned_data, stat_method)
% COMPUTE_GROUP_STATS  Stack aligned responses and compute center/spread.
%   stat_method: 'median_mad' (default) or 'mean_sem'.

    stats.theta  = [];
    stats.center = [];
    stats.spread = [];
    stats.n      = 0;

    if isempty(aligned_data)
        return;
    end

    vals = [];
    for k = 1:numel(aligned_data)
        d = aligned_data{k};
        if isnumeric(d) && size(d, 2) == 2 && size(d, 1) >= 2
            if isempty(stats.theta)
                stats.theta = d(:, 1);
            end
            vals = [vals, d(:, 2)]; %#ok<AGROW>
        end
    end

    if isempty(vals)
        return;
    end

    stats.n = size(vals, 2);

    if strcmpi(stat_method, 'mean_sem')
        stats.center = mean(vals, 2, 'omitnan');
        sd           = std(vals, 0, 2, 'omitnan');
        stats.spread = sd ./ sqrt(stats.n);
    else  % 'median_mad'
        stats.center = median(vals, 2, 'omitnan');
        stats.spread = mad(vals, 1, 2);  % median absolute deviation (flag=1)
    end

end


function rmax = estimate_rmax(sc, st)
% ESTIMATE_RMAX  Find maximum r-value across both groups for axis scaling.

    allVals = [];
    if ~isempty(sc.center)
        allVals = [allVals; sc.center(:) + sc.spread(:)];
    end
    if ~isempty(st.center)
        allVals = [allVals; st.center(:) + st.spread(:)];
    end
    if isempty(allVals)
        rmax = NaN;
    else
        rmax = max(allVals, [], 'omitnan');
        if ~isfinite(rmax) || rmax <= 0
            rmax = NaN;
        end
    end

end


function hLine = plot_with_shade(axFill, theta, centerVals, bandVals, ...
    lineColor, fillColor, alphaVal, lw)
% PLOT_WITH_SHADE  Draw center line and spread patch on the Cartesian overlay.
%   Both the shaded polygon and the center line are drawn on the same
%   Cartesian axes (axFill) using pol2cart, guaranteeing perfect alignment.
%   Returns the line handle for use in the legend.

    th = theta(:);
    m  = centerVals(:);
    b  = bandVals(:);

    % Ensure non-negative radii
    upper = max(m + b, 0);
    lower = max(m - b, 0);

    % Close the loop by appending the first point
    th_closed    = [th; th(1) + 2*pi];
    upper_closed = [upper; upper(1)];
    lower_closed = [lower; lower(1)];
    m_closed     = [m; m(1)];

    % Build closed polygon: outer ring forward, inner ring backward
    % No interpolation — straight segments match the center line exactly
    th_poly = [th_closed; flipud(th_closed)];
    r_poly  = [upper_closed; flipud(lower_closed)];
    [x_poly, y_poly] = pol2cart(th_poly, r_poly);

    % Shaded spread patch
    patch('XData', x_poly, 'YData', y_poly, ...
          'FaceColor', fillColor, 'FaceAlpha', alphaVal, ...
          'EdgeColor', 'none', 'Parent', axFill, 'HitTest', 'off');

    % Center line (same axes = guaranteed overlap with patch)
    [x_line, y_line] = pol2cart(th_closed, m_closed);
    hLine = plot(axFill, x_line, y_line, '-', 'Color', lineColor, 'LineWidth', lw);

end

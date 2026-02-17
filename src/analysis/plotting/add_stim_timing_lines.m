function add_stim_timing_lines(fig, stim_onset, stim_offset, line_color, line_style)
% ADD_STIM_TIMING_LINES  Overlay vertical lines marking stimulus ON/OFF on all axes.
%
%   ADD_STIM_TIMING_LINES(FIG, STIM_ONSET, STIM_OFFSET) draws vertical
%   solid green lines at the stimulus onset and offset sample indices on
%   every axes object in the figure whose x-limits encompass those values.
%   The lines are drawn behind existing data traces.
%
%   ADD_STIM_TIMING_LINES(..., LINE_COLOR) specifies the line colour as
%   an RGB triplet (default: [0.2 0.7 0.2], green).
%
%   ADD_STIM_TIMING_LINES(..., LINE_COLOR, LINE_STYLE) specifies the line
%   style (default: '-').
%
%   INPUTS:
%     fig          - Figure handle
%     stim_onset   - Sample index where stimulus turns ON (e.g. 5001)
%     stim_offset  - Sample index where stimulus turns OFF
%     line_color   - (Optional) RGB triplet, default [0.2 0.7 0.2]
%     line_style   - (Optional) Line style string, default '-'
%
%   NOTES:
%     Only axes whose current x-limits include stim_onset are annotated.
%     This means polar axes or axes with different x-ranges (e.g. the
%     central polar plot in the bar sweep figure) are automatically skipped.
%
%     Lines are drawn as plain plot() objects and moved to the bottom of
%     the axes children stack so that data traces render on top.
%
%   See also ANALYZE_SINGLE_EXPERIMENT_MR, PLOT_BAR_FLASH_1X11,
%            PLOT_BAR_FLASH_HEATMAP

    if nargin < 4 || isempty(line_color)
        line_color = [0.2 0.7 0.2];
    end
    if nargin < 5 || isempty(line_style)
        line_style = '-';
    end

    all_axes = findobj(fig, 'Type', 'axes');

    for k = 1:numel(all_axes)
        ax = all_axes(k);

        % Skip polar axes (they don't have linear x-limits)
        if isa(ax, 'matlab.graphics.axis.PolarAxes')
            continue;
        end

        xl = get(ax, 'XLim');
        yl = get(ax, 'YLim');

        % Only annotate axes whose x-range includes the stimulus timing
        if stim_onset >= xl(1) && stim_onset <= xl(2)
            hold(ax, 'on');

            % Draw as plain plot objects (full-height vertical lines)
            h1 = plot(ax, [stim_onset stim_onset], yl, line_style, ...
                'Color', [line_color 0.8], 'LineWidth', 1.2);
            h2 = plot(ax, [stim_offset stim_offset], yl, line_style, ...
                'Color', [line_color 0.8], 'LineWidth', 1.2);

            % Move lines to the bottom of the axes children stack
            % so data traces render on top
            ch = get(ax, 'Children');
            % h1 and h2 are at the top of the stack (most recent)
            % Move them to the bottom
            is_timing = ismember(ch, [h1; h2]);
            set(ax, 'Children', [ch(~is_timing); ch(is_timing)]);
        end
    end

end

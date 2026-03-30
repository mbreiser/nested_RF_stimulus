function fig = plot_flash_1x11_individual_rows(traces_cell, labels, title_str, opts)
% PLOT_FLASH_1X11_INDIVIDUAL_ROWS  Multi-row bar flash figure, one row per cell.
%
%   FIG = PLOT_FLASH_1X11_INDIVIDUAL_ROWS(TRACES_CELL, LABELS, TITLE_STR)
%   creates an N_cells × N_pos tiled layout where each row shows one cell's
%   bar flash traces across spatial positions. Each tile shows the actual
%   voltage trace at that position.
%
%   FIG = PLOT_FLASH_1X11_INDIVIDUAL_ROWS(..., OPTS) uses the options
%   structure to override defaults.
%
%   INPUTS:
%     traces_cell - Cell array (one per cell). Each element is an N_pos×N
%                   matrix of baseline-subtracted mean flash traces.
%                   N_pos can be 11 (original) or 15 (shifted canvas).
%                   NaN rows are rendered as blank tiles.
%     labels      - Cell array of row labels (one per cell), e.g.
%                   '2025-11-10-10-17 (c=5.7, s=+0)'
%     title_str   - Figure title string
%     opts        - (Optional) structure with fields:
%                     .y_limits   - [ymin ymax] in mV (default: [-15 35])
%                     .col_labels - Cell array of column header strings.
%                                   If empty, auto-generates degree labels
%                                   centered on the middle column.
%                     .center_col - Column index for the aligned center
%                                   (default: ceil(n_pos/2)). Used for
%                                   auto-generating degree labels and for
%                                   highlighting the center column.
%                     .line_color - N_cells×3 colour matrix, one row per
%                                   cell (default: MATLAB lines colormap)
%
%   OUTPUT:
%     fig - Figure handle
%
%   See also SHIFT_TRACES_TO_CANVAS, PLOT_FLASH_1X11_POPULATION,
%            ANALYZE_SINGLE_EXPERIMENT_MR
% ________________________________________________________________________

    %% Defaults
    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'y_limits'),   opts.y_limits   = [-15 35]; end
    if ~isfield(opts, 'col_labels'), opts.col_labels  = {}; end
    if ~isfield(opts, 'center_col'), opts.center_col  = []; end
    if ~isfield(opts, 'line_color'), opts.line_color   = []; end

    n_cells = numel(traces_cell);

    % Derive n_pos from first non-empty trace
    n_pos = 0;
    for kk = 1:n_cells
        if ~isempty(traces_cell{kk})
            n_pos = size(traces_cell{kk}, 1);
            break;
        end
    end
    if n_pos == 0
        fig = figure('Name', title_str);
        return;
    end

    % Center column
    if isempty(opts.center_col)
        opts.center_col = ceil(n_pos / 2);
    end
    center_col = opts.center_col;

    % Auto-generate degree labels if not provided
    if isempty(opts.col_labels)
        opts.col_labels = cell(1, n_pos);
        for p = 1:n_pos
            deg_val = (p - center_col) * 2.5;
            if p == 1
                opts.col_labels{p} = sprintf('%.1f° (Lead)', deg_val);
            elseif p == center_col
                opts.col_labels{p} = '0° (Center)';
            elseif p == n_pos
                opts.col_labels{p} = sprintf('+%.1f° (Trail)', deg_val);
            elseif deg_val > 0
                opts.col_labels{p} = sprintf('+%.1f°', deg_val);
            else
                opts.col_labels{p} = sprintf('%.1f°', deg_val);
            end
        end
    end

    % Colours
    if isempty(opts.line_color)
        opts.line_color = lines(n_cells);
    end

    %% Create figure
    fig = figure('Name', title_str);
    tiledlayout(n_cells, n_pos, 'TileSpacing', 'tight', 'Padding', 'compact');

    for cell_idx = 1:n_cells
        mat = traces_cell{cell_idx};  % n_pos × N timepoints
        col = opts.line_color(cell_idx, :);

        for pos_idx = 1:n_pos
            nexttile;
            hold on;

            has_data = false;
            if ~isempty(mat) && pos_idx <= size(mat, 1)
                trace = mat(pos_idx, :);
                if ~all(isnan(trace))
                    x = 1:numel(trace);
                    plot(x, trace, '-', 'Color', col, 'LineWidth', 1.0);
                    has_data = true;
                end
            end

            % Gray background for tiles with no data
            if ~has_data
                set(gca, 'Color', [0.95 0.95 0.95]);
            end

            ylim(opts.y_limits);
            set(gca, 'XTick', []);

            % Highlight center column with subtle left-edge line
            if pos_idx == center_col
                xline(1, '-', 'Color', [0.6 0 0], 'LineWidth', 1.5, 'Alpha', 0.4);
            end

            % Y-axis label on first column only (show row label)
            if pos_idx == 1
                ylabel(labels{cell_idx}, 'FontSize', 6, 'Interpreter', 'none');
            else
                set(gca, 'YTickLabel', []);
            end

            % Column headers on first row only
            if cell_idx == 1
                title(opts.col_labels{pos_idx}, 'FontSize', 7);
            end

            box off;
            ax = gca;
            ax.LineWidth = 0.6;
            ax.TickDir = 'out';
            ax.TickLength = [0.02 0.02];
            ax.FontSize = 6;
        end
    end

    sgtitle(title_str, 'FontSize', 13, 'FontWeight', 'bold');

    % Figure size: adaptive
    fig_width  = max(1200, 120 * n_pos);
    fig_height = max(400, 180 * n_cells);
    set(fig, 'Position', [50 50 fig_width fig_height]);

end

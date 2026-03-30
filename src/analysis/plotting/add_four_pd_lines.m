function add_four_pd_lines(fig, pd_stim_deg_4, max_rho, method_labels)
% ADD_FOUR_PD_LINES  Draw 4 PD direction lines on the polar plot.
%
%   ADD_FOUR_PD_LINES(FIG, PD_STIM_DEG_4, MAX_RHO) draws 4 diametric lines
%   on the PolarAxes in FIG, one per PD alignment method.
%
%   ADD_FOUR_PD_LINES(FIG, PD_STIM_DEG_4, MAX_RHO, METHOD_LABELS) uses
%   custom legend labels (default: {'VecNN','SymNN','VecInterp','SymInterp'}).
%
%   INPUTS:
%     fig            - Figure handle containing a PolarAxes
%     pd_stim_deg_4  - 4x1 vector of PD directions in stimulus degrees
%                      [VecNN, SymNN, VecInterp, SymInterp]
%     max_rho        - Maximum rho value (for line extent)
%     method_labels  - (Optional) 4x1 cell of legend labels
%
%   Line styles:
%     VecNN:     solid black,   LineWidth 2.0
%     SymNN:     dashed green,  LineWidth 1.5
%     VecInterp: dotted blue,   LineWidth 1.5
%     SymInterp: solid magenta, LineWidth 2.0
%
%   See also ADD_POLAR_DIRECTION_LINES, COMPUTE_PD_FOUR_METHODS

    if nargin < 4
        method_labels = {'VecNN', 'SymNN', 'VecInterp', 'SymInterp'};
    end

    pax = findobj(fig, 'Type', 'PolarAxes');
    if isempty(pax), return; end

    r_max = max_rho * 1.05;
    hold(pax, 'on');

    % Style definitions: {Color, LineStyle, LineWidth}
    styles = {
        [0    0    0   ], '-',  2.0;   % VecNN: black solid
        [0.15 0.65 0.15], '--', 1.5;   % SymNN: green dashed
        [0.15 0.35 0.70], ':',  1.5;   % VecInterp: blue dotted
        [0.75 0.00 0.65], '-',  2.0;   % SymInterp: magenta solid
    };

    h = gobjects(4, 1);
    for m = 1:4
        pd_rad = deg2rad(pd_stim_deg_4(m));
        h(m) = polarplot(pax, [pd_rad, pd_rad + pi], [r_max, r_max], ...
            styles{m, 2}, 'Color', styles{m, 1}, 'LineWidth', styles{m, 3});
    end

    % Push all 4 lines behind the tuning curve / arrow
    ch = get(pax, 'Children');
    is_guide = ismember(ch, h);
    set(pax, 'Children', [ch(~is_guide); ch(is_guide)]);

    % Legend
    legend(h, method_labels, 'Location', 'southeast', 'FontSize', 7, ...
        'Box', 'off');

end

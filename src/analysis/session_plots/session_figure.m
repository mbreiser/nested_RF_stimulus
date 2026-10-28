function fig = session_figure(name, visible, pos)
%SESSION_FIGURE  The figure window of a session plot: reused by name, white, compact default fonts.
%   FIG = SESSION_FIGURE(NAME, VISIBLE, POS). 28 October 2026.
fig = findobj(groot, 'Type', 'figure', 'Name', name);
if isempty(fig)
    fig = figure('Name', name, 'NumberTitle', 'off', 'Color', 'w', 'Position', pos, 'Visible', tern(visible, 'on', 'off'));
else
    fig = fig(1); clf(fig); fig.Visible = tern(visible, 'on', 'off');
end
set(fig, 'DefaultAxesFontSize', 8, 'DefaultAxesTitleFontWeight', 'normal', 'DefaultAxesTickDir', 'out', ...
    'DefaultAxesBox', 'off', 'DefaultLineLineWidth', 1);
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end

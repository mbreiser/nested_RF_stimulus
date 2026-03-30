% OPEN_MANUSCRIPT_FIG  Open a saved .fig file for interactive layout tweaking.
%
%   Opens the manuscript figure in MATLAB for manual adjustment of panel
%   positions, font sizes, and other properties. Use the Property Inspector
%   (double-click elements) or command line to tweak, then re-export:
%
%     exportgraphics(gcf, 'fig_ds_panels_ABC.pdf', 'ContentType', 'vector');
%     exportgraphics(gcf, 'fig_ds_panels_ABC.png', 'Resolution', 300);
%
%   Usage:
%     run('scripts/open_manuscript_fig.m')

fig_file = '/Users/reiserm/Documents/ttl_1DRF/manuscript_figures/fig_ds_panels_ABC.fig';

if ~isfile(fig_file)
    error('Fig file not found: %s\nRun generate_manuscript_fig_ds.m first.', fig_file);
end

fig = openfig(fig_file, 'visible');
fprintf('Opened: %s\n', fig_file);
fprintf('\nTips:\n');
fprintf('  - Drag/resize axes interactively\n');
fprintf('  - Use findall(fig, ''Type'', ''axes'') to list all axes\n');
fprintf('  - set(ax, ''Position'', [x y w h]) to reposition\n');
fprintf('  - exportgraphics(fig, ''output.pdf'', ''ContentType'', ''vector'') to re-export\n');
fprintf('  - savefig(fig, ''%s'') to save changes back\n', fig_file);

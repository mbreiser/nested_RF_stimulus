% REGENERATE_ALL_CENTROID_FIGURES  Re-run centering on all 25 cells and
% regenerate population alignment figures using the updated bump-based
% centroid algorithm.
%
% Steps:
%   1. Run ON cells (10) with visualize_everything → new centering summary
%   2. Run OFF cells (15) with visualize_everything → new centering summary
%   3. Run population alignment figure generation

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

fprintf('========================================\n');
fprintf('  Step 1: ON cells (10)\n');
fprintf('========================================\n');
run('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/scripts/run_on_cells_mr.m');

fprintf('\n\n========================================\n');
fprintf('  Step 2: OFF cells (15)\n');
fprintf('========================================\n');
run('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/scripts/run_off_cells_mr.m');

fprintf('\n\n========================================\n');
fprintf('  Step 3: Population alignment figures\n');
fprintf('========================================\n');
run('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/scripts/plot_centroid_alignment_figures.m');

fprintf('\n\nAll done.\n');

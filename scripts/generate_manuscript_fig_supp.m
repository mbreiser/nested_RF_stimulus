% GENERATE_MANUSCRIPT_FIG_SUPP  Build the supplementary manuscript figure.
%   Thin wrapper around generate_manuscript_fig('supp').

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));
generate_manuscript_fig('supp');

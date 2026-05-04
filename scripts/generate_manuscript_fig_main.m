% GENERATE_MANUSCRIPT_FIG_MAIN  Build the main manuscript figure.
%   Thin wrapper around generate_manuscript_fig('main').

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));
generate_manuscript_fig('main');

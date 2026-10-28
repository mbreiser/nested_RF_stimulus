function [pattern_order, func_order, trial_dur] = create_protocol2(exp_folder, metadata)
% CREATE_PROTOCOL2  Assemble currentExp.mat from patterns and functions.
%
%   [PATTERN_ORDER, FUNC_ORDER, TRIAL_DUR] = CREATE_PROTOCOL2(EXP_FOLDER, METADATA)
%   reads all pattern and position function files from the experiment folder,
%   organizes them into the correct presentation order for the 1DRF protocol,
%   and saves the experiment configuration to currentExp.mat.
%
%   INPUTS:
%     exp_folder - char | string
%                  Full path to the experiment directory containing
%                  'Patterns/' and 'Functions/' subdirectories.
%     metadata   - struct
%                  Structure from GET_INPUT_PARAMETERS containing:
%                    .Frame  (numeric) - peak frame from Protocol 1
%                    .Side   (char)    - 'L' or 'R'
%                    .Age    (char)    - fly age in days
%                    .Strain (char)    - fly genotype string
%
%   Kept in step with the rig's nested_RF_stimulus (branch setup-pharma,
%   commits e99cc82, 5f269a7, 59a57d5), which makes Jin Yong's Protocol 2
%   packages: five sweep speeds, a first "mock" sweep, no square-flash
%   condition (28 Oct 2026).
%
%   OUTPUTS:
%     pattern_order - 1xN double
%                     Pattern presentation order indices. Structure:
%                       [grey, mock sweep (bar_patt(1), 28 dps FLIP; excluded
%                        from the analysis),
%                        bar_patt x2 (fwd+flip) x n_bar_patt at 28 dps, grey,
%                        ... at 56, 168, 250 and 500 dps, each followed by grey,
%                        bar_flash_slow, grey, bar_flash_fast]
%     func_order    - 1xN double
%                     Position function indices matched to pattern_order.
%     trial_dur     - 1xN double
%                     Duration of each trial in seconds.
%
%   PATTERN ASSIGNMENT:
%     Pattern 1       : 4px flash pattern (shown on its frame 1 as the grey;
%                       the flashes are not played)
%     Patterns 2..N-1 : bar sweep patterns (8 orientations)
%     Pattern N       : bar flash pattern (FLASHES)
%
%   FUNCTION ASSIGNMENT:
%     Function 1      : 4px flash position function (not played)
%     Functions 2-3   : 28 dps bar fwd/flip
%     Functions 4-5   : 56 dps bar fwd/flip
%     Functions 6-7   : 168 dps bar fwd/flip
%     Functions 8-9   : 250 dps bar fwd/flip
%     Functions 10-11 : 500 dps bar fwd/flip
%     Functions 12-14 : bar flash slow (80 ms, Rep1-3)
%     Functions 15-17 : bar flash fast (14 ms, Rep1-3)
%     Function N-1    : 10s static (initial grey baseline)
%     Function N      : 3s static (inter-block gap)
%
%   SAVED FILE:
%     currentExp.mat - contains currentExp struct with pattern/function
%                      metadata, plus pattern_order, func_order, trial_dur,
%                      and metadata arrays.
%
%   See also GENERATE_PROTOCOL2, RUN_PROTOCOL2, GENERATE_STATIC_FUNCTION

%% check for correct folders
pattern_folder = fullfile(exp_folder, 'Patterns');
function_folder = fullfile(exp_folder, 'Functions');

assert(exist(pattern_folder,'dir')==7,'did not detect pattern folder')
assert(exist(function_folder,'dir')==7,'did not detect position function folder')

%% read patterns
matinfo = dir(fullfile(pattern_folder, "*.mat"));
patinfo = dir(fullfile(pattern_folder, '*.pat'));
% num_files is the number of patterns.
num_files = length({patinfo.name});

flash_patt = 1;
bar_patt = flash_patt(end)+1:num_files-1;
bar_flash_patt = num_files; % Last pattern is the bar flashes.
n_bar_patt = numel(bar_patt);
n_dir = 2; % 2 reps = forward and reverse direction
n_speeds = 5; % 5 reps = n_speeds
pattern_order = [1, bar_patt(1), repmat([repelem(bar_patt, n_dir), 1], [1, n_speeds]), bar_flash_patt, 1, bar_flash_patt];

% two "flash_patt" at the beginning - one for grey presentation.
n_patts = numel(pattern_order);

for p = 1:n_patts
    f = pattern_order(p); % Index of pattern to load
    currentExp.pattern.pattNames{p} = matinfo(f).name;
    currentExp.pattern.patternList{p} = patinfo(f).name;

    patternIN = load(fullfile(pattern_folder, matinfo(f).name));
    currentExp.pattern.x_num(p) = patternIN.pattern.x_num; % Matlab now has its own built in class "pattern" which conflicts here.
    currentExp.pattern.y_num(p) = patternIN.pattern.y_num;
    currentExp.pattern.gs_val(p) = patternIN.pattern.gs_val;
    currentExp.pattern.arena_pitch(p) = round(rad2deg(patternIN.pattern.param.arena_pitch));
end

currentExp.pattern.num_patterns = num_files;

%% read position functions
generate_static_function(function_folder, 10)
generate_static_function(function_folder, 3)

matinfo = dir(fullfile(function_folder, "*.mat"));
pfninfo = dir(fullfile(function_folder, "*.pfn"));
num_files = length({pfninfo.name});
total_dur = 0;
trial_dur = nan([1, n_patts]);

fns_28dps = [2,3];
fns_56dps = [4,5];
% fns_112dps = [6,7];
fns_168dps = [6,7];
fns_250dps = [8,9];
fns_500dps = [10,11];

fn_bar_flash_slow = 12; % run_protocol2 plays the -Rep2 / -Rep3 functions of the same name on repeats 2 and 3
fn_bar_flash_fast = 15;
fn_static_10 = num_files-1;
fn_static_3 = num_files;

% Add a first bar sweep in 180 direction - exclude this from the analysis.
func_order = [fn_static_10,...
    fns_28dps(2), ...
    repmat(fns_28dps, [1,n_bar_patt]), fn_static_3,...
    repmat(fns_56dps, [1,n_bar_patt]), fn_static_3,...
    repmat(fns_168dps, [1,n_bar_patt]), fn_static_3,...
    repmat(fns_250dps, [1,n_bar_patt]), fn_static_3,...
    repmat(fns_500dps, [1,n_bar_patt]), fn_static_3,...
    fn_bar_flash_slow, fn_static_3,...
    fn_bar_flash_fast];

    % repmat(fns_112dps, [1,n_bar_patt]), fn_static_3,...

for p = 1:n_patts
    f = func_order(p);
    currentExp.function.functionName{p} = matinfo(f).name;
    currentExp.function.functionList{p} = pfninfo(f).name;
    load(fullfile(function_folder, matinfo(f).name))
    currentExp.function.functionSize(p) = pfnparam.size;
    currentExp.function.trial_dur(p) = sum(pfnparam.dur);
    total_dur = total_dur + sum(pfnparam.dur);
    trial_dur(p) = sum(pfnparam.dur);
end
currentExp.totalDuration = total_dur;
currentExp.function.numFunc = num_files;

%% save currentExp structure
save(fullfile(exp_folder, 'currentExp.mat'),'currentExp', 'pattern_order', 'func_order', 'trial_dur', 'metadata');

end
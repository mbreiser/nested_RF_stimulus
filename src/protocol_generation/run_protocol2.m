function run_protocol2(exp_folder, pattern_order, func_order, trial_dur, n_reps)
% RUN_PROTOCOL2  Execute Protocol 2 on the G4 LED arena.
%
%   RUN_PROTOCOL2(EXP_FOLDER, PATTERN_ORDER, FUNC_ORDER, TRIAL_DUR)
%   connects to the G4 Panels controller and runs the complete Protocol 2
%   experiment, presenting all patterns with their position functions and
%   logging voltage and frame timing data.
%
%   INPUTS:
%     exp_folder    - Full path to experiment directory containing
%                     currentExp.mat, Patterns/, and Functions/
%     pattern_order - 1xN array of pattern IDs in presentation order
%     func_order    - 1xN array of function IDs matching pattern_order
%     trial_dur     - 1xN array of trial durations in seconds
%
%   EXPERIMENT EXECUTION:
%     - Runs 3 repetitions of all conditions (n_reps = 3)
%     - Control mode 1 (position function control)
%     - Logs data to TDMS files in 'Log Files/' subdirectory
%     - Displays progress to command window
%     - Presents grey screen at end of experiment
%     - Converts TDMS files to .mat format after completion
%
%   PREREQUISITES:
%     - G4 Display Tools installed (PanelsController class available)
%     - userSettings.m configured with correct paths
%     - Arena connected and powered on
%     - No existing log files in experiment folder
%
%   OUTPUT:
%     Creates in exp_folder/Log Files/:
%       - G4_TDMS_Log_*.mat files with voltage and frame timing data
%
%   DURATION:
%     Estimated time is printed to console before experiment starts.
%     Typical Protocol 2 duration: ~45-60 minutes for 3 repetitions.
%
%   See also CREATE_PROTOCOL2, GENERATE_PROTOCOL2, G4_TDMS_FOLDER2STRUCT

    %% Experiment metadata from user input:
    % n_reps = 3;
    tic;
    %% set up for experiment
    %Load configuration and start G4 Host
    % Check the use of this
    userSettings;

    load(fullfile(exp_folder,'currentExp.mat'));
    num_conditions = numel(pattern_order);
    log_folder = fullfile(exp_folder,'Log Files');
    if ~exist(log_folder,'dir')
        mkdir(log_folder);
    end

    %% Position function per rep and condition
    % The bar flash pattern is the one whose file name holds FLASHES; it is
    % played every rep. Its position functions come in -Rep1/-Rep2/-Rep3
    % versions (a different flash order each rep): func_order names the Rep1
    % function, rep r uses the Rep<r> function of the same name.
    % IDs are the 4-digit prefixes of the files in Patterns\ and Functions\.
    pat_files = dir(fullfile(exp_folder, 'Patterns', '*.mat'));
    pat_ids = cellfun(@(n) str2double(n(1:4)), {pat_files.name});
    bar_flash_patterns = pat_ids(contains({pat_files.name}, 'FLASHES'));
    func_files = dir(fullfile(exp_folder, 'Functions', '*.mat'));
    func_ids = cellfun(@(n) str2double(n(1:4)), {func_files.name});
    func_names = cellfun(@(n) n(5:end), {func_files.name}, 'UniformOutput', false);

    rep_func = repmat(func_order(1,:), n_reps, 1);
    for c = 1:num_conditions
        if ~ismember(pattern_order(1,c), bar_flash_patterns)
            continue
        end
        f1 = func_names{func_ids == func_order(1,c)};
        if isempty(regexp(f1, 'Rep1(?!\d)', 'once'))
            continue
        end
        for r = 2:n_reps
            fr = regexprep(f1, 'Rep1(?!\d)', sprintf('Rep%d', r));
            k = find(strcmp(func_names, fr));
            if isscalar(k)
                rep_func(r, c) = func_ids(k);
            else
                fprintf('WARNING: no function %s for rep %d of cond %d; using %s\n', fr, r, c, f1);
            end
        end
    end
    fprintf('Bar flash pattern(s): %s\n', mat2str(bar_flash_patterns));
    for c = find(any(rep_func ~= func_order(1,:), 1))
        fprintf('Cond %d (pattern %d): functions %s for reps 1-%d\n', c, pattern_order(1,c), mat2str(rep_func(:,c).'), n_reps);
    end

    %% Open new Panels controller instance
    ctlr = PanelsController();
    ctlr.open(true);
    
    %% Change to root directory
    ctlr.setRootDirectory(exp_folder);
    
    %check if log files already present for this experiment
    assert(~exist(fullfile(exp_folder, 'Log Files\*'),'file'),'unsorted log files present in save folder, remove before restarting experiment\n');
    % assert(~exist(fullfile(exp_folder, 'Results\', fly_name),'dir'),'Results folder already exists with that fly name\n');
    
    %finish setting up for experiment
    exp_seconds = currentExp.totalDuration;
    fprintf(['Estimated experiment duration: ' num2str((exp_seconds*n_reps)/60) ' minutes\n']);
    
    %% start experiment
    ctlr.startLog(); %starts logging data in .tdms files

    for r = 1:n_reps
        for c = 1:num_conditions
            %trial portion
            ctlr.setControlMode(1);
            ctlr.setPatternID(pattern_order(1,c));
            % Use different position functions for bar flashes each rep.
            ctlr.setPatternFunctionID(rep_func(r,c));
            trial_t = trial_dur(1, c);
            fprintf(['Rep ', num2str(r), ' of ', num2str(n_reps), ', cond ' num2str(c) ' of ' num2str(num_conditions) ': ' strjoin(currentExp.pattern.pattNames(pattern_order(1,c))) '\n']);
            ctlr.startDisplay(ceil(trial_t*10)-1); %duration expected in 100ms units
        end
    end 

    % ADD Greyscale at the end - only at the very end - not at the end of
    % every rep. 
    ctlr.setControlMode(1);
    ctlr.setPatternID(pattern_order(1,1));
    ctlr.setPatternFunctionID(func_order(1,1));
    trial_t = trial_dur(1, 1);
    fprintf('Greyscale interval at the end.');
    ctlr.startDisplay(ceil(trial_t*10)-1); %duration expected in 100ms units

    %rename/move results folder
    ctlr.stopLog('showTimeoutMessage', true);
    ctlr.stopDisplay()
    ctlr.close()
    disp('finished');

    toc;
    % % Convert TDMS files to mat file - current issues.
    G4_TDMS_folder2struct(log_folder)
    

end 
function generate_bar_flash_pos_fns(save_dir, n_flank, n_reps)

% 80ms flashes + 920ms interval - - - match 28 dps moving bars. 
% For each bar pattern, generate a random order 

n_orientations = 8;
n_frames_orient = (n_flank*2 + 1);
n_frames_pattern = n_frames_orient * n_orientations;

for sp = [1,2]
    
    if sp == 1
        flash_dur = 0.08;
        flash_str = "80";
        int_dur = 0.92; 
        int_str = "920";
    elseif sp == 2
        flash_dur = 0.014;
        flash_str = "14";
        int_dur = 0.486; 
        int_str = "486";
    end 
    
    % Make 3 different position functions, one for each rep.
    for r = 1:n_reps
    
        % Random order in which to present flashes (n_flank 5: randperm(88), as before 30 Oct 2026).
        rand_order = randperm(n_frames_pattern);
        % Swap the value of 1 for the last frame (89 for n_flank 5) - first frame is blank
        rand_order(rand_order==1)=n_frames_pattern+1;
    
        value_array = zeros(1, 2 * numel(rand_order)); 
        value_array(1:2:end) = 1; 
        value_array(2:2:end) = rand_order;
        n_vals = numel(value_array);
    
        % user-defined function parameters
        pfnparam.type = 'pfn'; 
        pfnparam.frames = n_frames_pattern+1; %number of frames in pattern
        pfnparam.gs_val = 4; %brightness bits in pattern
        pfnparam.section = repmat({'static'}, 1, n_vals); 
        pfnparam.dur = repmat([int_dur flash_dur], 1, n_frames_pattern); %section duration (in s)
        pfnparam.val = value_array; %function value for static sections
        pfnparam.high = ones(1, n_vals)*(n_frames_pattern+1);
        pfnparam.low = ones(1, n_vals);
        pfnparam.freq = ones(1, n_vals);
        pfnparam.flip = zeros(1, n_vals);
        
        %% generate function
        func = G4_Function_Generator(pfnparam);
        
        %% save function
        pfnparam.ID = get_function_ID('pfn',save_dir);
        filename = strcat("bar_flashes_", flash_str, 'ms_', int_str, 'ms_', string(n_frames_pattern), 'frames-Rep', string(r));
        save_function_G4(func, pfnparam, save_dir, filename);
    
    end 

end

end 


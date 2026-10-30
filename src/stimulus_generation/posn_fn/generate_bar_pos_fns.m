function generate_bar_pos_fns(save_dir, px_crop)
% GENERATE_BAR_POS_FNS  Create position functions for moving bar stimuli.
%
%   GENERATE_BAR_POS_FNS(SAVE_DIR, PX_CROP) generates sawtooth position
%   functions for moving bar stimuli at 5 speeds (28, 56, 168, 250, 500 dps)
%   in both forward and flip (reverse) directions.
%
%   Kept in step with the rig's nested_RF_stimulus (branch setup-pharma,
%   commit e99cc82 "Add extra bar sweep speeds"), which makes Jin Yong's
%   Protocol 2 packages (28 Oct 2026).
%
%   INPUTS:
%     save_dir - char | string
%                Directory to save the .pfn and .mat function files.
%     px_crop  - double (30 or 45)
%                Size of the cropped display region in pixels. Determines
%                the frame range and timing parameters for the sawtooth.
%
%   SPEEDS GENERATED (px_crop = 30; static first frame, then the sweep):
%     28 dps  - 1 s static,    2.328 s sweep
%     56 dps  - 1 s static,    1.175 s sweep
%     168 dps - 0.75 s static, 0.40 s sweep
%     250 dps - 0.5 s static,  0.266 s sweep
%     500 dps - 0.5 s static,  0.135 s sweep
%     (px_crop = 45 defines 28 and 56 dps only, each after 1 s static)
%
%   DIRECTIONS:
%     For each speed, two functions are generated:
%       Forward : frames ramp from low to high
%       Flip    : frames ramp from high to low (opposite direction)
%
%   FRAME RANGES:
%     px_crop = 45 : frames 1-68  (full visible bar movement)
%     px_crop = 30 : frames 11-62 (cropped to visible region)
%
%   TIMING CALCULATION:
%     Duration = (frame_range / degrees_per_frame) / dps
%     Frequency = 1 / duration (for sawtooth waveform generation)
%
%   OUTPUT FILES (10 total = 5 speeds x 2 directions):
%     For each speed/direction combination:
%       <dps>dps_[FLIP_]<dur>s_<low>-<high>_288fr_1-25step_3s_static.pfn
%       Corresponding .mat file with pfnparam structure
%
%   See also GENERATE_BAR_STIMULUS_XY, G4_FUNCTION_GENERATOR, CREATE_PROTOCOL2

% July 2025 - added 1s of static first frame before each
% sawtooth.

    for dps = [28, 56, 168, 250, 500]

        if px_crop == 45

            dur_static = 1;
            if dps == 14
                dur_sweep = 6.08;
                dur_str = '6-08';
                freq = 0.16471;
            elseif dps == 28
                dur_sweep = 3.04;
                dur_str = '3-04';
                freq = 0.32941;
            elseif dps == 56
                dur_sweep = 1.52;
                dur_str = '1-52';
                freq = 0.65882;
            end

            fr_low = 1;
            fr_high = 68;

        elseif px_crop == 30

            if dps == 28
                dur_static = 1;
                dur_sweep = 2.328;
                dur_str = '2-328';
                freq = 0.43077;
            elseif dps == 56
                dur_static = 1;
                dur_sweep = 1.175;
                dur_str = '1-175';
                freq = 0.86154;
            elseif dps == 112
                dur_static = 0.75;
                dur_sweep = 0.6;
                dur_str = '0-6';
                freq = 1.7231;
            elseif dps == 168
                dur_static = 0.75;
                dur_sweep = 0.40;
                dur_str = '0-40';
                freq = 2.5846;
            elseif dps == 250
                dur_static = 0.5;
                dur_sweep = 0.266;
                dur_str = '0-266';
                freq = 3.8462;
            elseif dps == 500
                dur_static = 0.5;
                dur_sweep = 0.135;
                dur_str = '0-135';
                freq = 7.6923;
            end

            fr_low = 11;
            fr_high = 62;

        end

        for flip = [0, 1]
            if flip == 1
                flip_str = 'FLIP_';
            else
                flip_str = '';
            end

            % user-defined function parameters
            pfnparam.type = 'pfn'; %number of frames in pattern
            pfnparam.frames = 288; %number of frames in pattern
            pfnparam.gs_val = 4; %brightness bits in pattern
            pfnparam.section = {'static' 'sawtooth'}; %static, sawtooth, traingle, sine, cosine, or square
            pfnparam.dur = [dur_static dur_sweep]; %section duration (in s)
            pfnparam.val = [1 1]; %function value for static sections
            pfnparam.high = [62 fr_high]; %high end of function range {for non-static sections}
            pfnparam.low = [1 fr_low ]; %low end of function range {for non-static sections}
            pfnparam.freq = [1 freq ]; %frequency of section {for non-static sections}
            pfnparam.size_speed_ratio = [40 40]; %size/speed ratio {for looms}
            pfnparam.flip = [0 flip]; %flip the range of values of function {for non-static sections}

            %% generate function
            func = G4_Function_Generator(pfnparam);

            %% save function
            pfnparam.ID = get_function_ID('pfn',save_dir);
            filename = strcat(string(dps), 'dps_', flip_str, dur_str, 's_', string(fr_low), '-', string(fr_high), '_288fr_1-25step_3s_static');
            save_function_G4(func, pfnparam, save_dir, filename);
        end
    end

end

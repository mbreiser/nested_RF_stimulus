function generate_bar_flash_stimulus_xy(x, y, px_intensity, px_crop, on_off, exp_folder, n_flank)
% Generate inset bar patterns centred on [X, Y]
    
% From an array of the "central" frame for each pattern. Choose "n_flank"
% frames either side of the central frame. 
% If "n_flank" = 5, then there will be 11 frames per bar orientation. 
% There are 8 orientations in total so if "n_flank" = 5, then 11 x 8 = 88
% frames (n_flank 7: 15 x 8 = 120 frames, bars out to +/-14 px; the
% full-field patterns have those frames, 30 Oct 2026).

% These will then be randomly selected. 

    centre_frames = [30, 33, 33, 34, 34, 34, 34, 35, 34, 34, 32, 32, 32, 34, 34, 34]; 

    save_dir = strcat(exp_folder, '\patterns');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    
    % Path to the folder with the full field bar stimuli.
    full_bar_patterns_dir = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\results\patterns\protocol2\full_field_bars4'; % July 2025 - changed to 4bkg
    cd(full_bar_patterns_dir)
    
    bar_patts = dir('*.mat');
    
    if on_off == "on"
        patt_ids = 1:1:8;
    elseif on_off == "off"
        patt_ids = 9:1:16;
    elseif on_off == "both"
        patt_ids = 1:1:16;
    end 

    % Make empty pattern - this will include the frames for each
    % orientation.
    n_orientations = 8;
    n_frames_orient = (n_flank*2 + 1);
    n_frames_pattern = n_frames_orient * n_orientations;
    bkg_color = px_intensity(1);
    arena_px_h = 48;
    arena_px_w = 192;
    
    % Make empty pattern array purely with grey background.
    Pats = ones(arena_px_h, arena_px_w, n_frames_pattern+1)*bkg_color;

    frames_filled = 0;
    
    % Generate 'cropped' bar stimuli centred on the [x,y] coordinates.
    for p = patt_ids

        % Load the appropriate bar pattern:
        fname = bar_patts(p).name;
        patt_file = fullfile(bar_patts(p).folder, fname);
        load(patt_file, 'pattern');

        % FInd where on the screen to position the cropped stimulus
        [disp_y1, disp_y2, disp_x1, disp_x2] = centeredSquare(x, y, px_crop);

        % crop around the centre of the arena screen [24, 96];
        arena_px_centre_h = arena_px_h/2;
        arena_px_centre_w = arena_px_w/2;
        px_side = (px_crop-1)/2;
        crop_h_st = int8(arena_px_centre_h-px_side); 
        crop_h_end = int8(arena_px_centre_h+px_side);
        crop_w_st = int8(arena_px_centre_w-px_side);
        crop_w_end = int8(arena_px_centre_w+px_side);
        px_rng = [crop_h_st, crop_h_end, crop_w_st, crop_w_end];
        
        overall_frame_range = (1:n_frames_orient)+frames_filled;

        cf_of_pattern = centre_frames(p); % Which frame number has the bar in the centre of the arena.
        frames_from_patt = cf_of_pattern-(n_flank*2): 2 : cf_of_pattern+(n_flank*2);
        assert(frames_from_patt(1) >= 1 && frames_from_patt(end) <= size(pattern.Pats, 3), ...
            'generate_bar_flash_stimulus_xy: n_flank %d needs frames %d..%d of %s (%d frames)', n_flank, frames_from_patt(1), frames_from_patt(end), fname, size(pattern.Pats, 3));

        % Fill in the location with the central part of pattern. 
        for i = 1:n_frames_orient

            % Index in which to find the frame
            frame_id = frames_from_patt(i);

            % index in which to store the frame
            overall_id = overall_frame_range(i);

            Pats(int16(disp_y1):int16(disp_y2), int16(disp_x1):int16(disp_x2), overall_id+1) = pattern.Pats(crop_h_st:crop_h_end, crop_w_st:crop_w_end, frame_id);
        end

        frames_filled = frames_filled + n_frames_orient;

    end 

    patName = strcat(fname(6:end-14), 'CROP', string(px_crop), '_X', string(x), '_Y', string(y), '_FLASHES');
    param = pattern.param;
    param.stretch = zeros(size(Pats, 3), 1);
    param.gs_val = 4;
    param.arena_pitch = 0;
    param.px_rng = px_rng;
    
    param.ID = get_pattern_ID(save_dir);
        
    save_pattern_G4(Pats, param, save_dir, patName);

    %test
    % figure
    % for i = 1:size(Pats, 3)
    % aa = Pats(:, :, i);
    % imagesc(aa)
    % pause(0.2)
    % end 

end 
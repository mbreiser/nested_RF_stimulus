% PLOT_PRE_BF_SWEEPS_INDIVIDUAL  Polar ring-of-traces figures for each
% cell in the pre-bar-flash dataset.
%
%   Uses the same plot_slow_bar_sweep_polar function as the main dataset,
%   producing a central polar tuning curve surrounded by 16 radially
%   arranged time-series traces (gray = per-cycle reps, blue = mean).
%
%   Also adds PD/ortho direction lines and DSI annotation, matching the
%   enhanced style of analyze_single_experiment_mr.m Figure 1.
%
%   Output: One PNG per cell in <data_root>/figure_previews/sweeps/
%           Named: sweep_<group>_<date_str>.png

%% Setup
addpath(genpath('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/src'));
addpath('/Users/reiserm/HHMI Dropbox/Michael Reiser/Matlab_work/CircStat2012a');

data_root = '/Users/reiserm/Documents/ttl_1DRF/pre-bar-flash';
save_dir  = fullfile(data_root, 'figure_previews', 'sweeps');
if ~isfolder(save_dir), mkdir(save_dir); end

% Standard plot order and sweep response options (same as batch pipeline)
plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
sweep_opts.baseline_range = [1000 9000];
sweep_opts.stim_trim_end  = 7000;
sweep_opts.percentile     = 98;

% Load LUT
lut_path = fullfile(fileparts(mfilename('fullpath')), ...
    '..', 'src', 'analysis', 'protocol2', 'bar_lut.mat');
lut_path = char(java.io.File(lut_path).getCanonicalPath());
S_lut = load(lut_path, 'Tbl');
Tbl = S_lut.Tbl;

%% Discover experiments (same logic as batch pipeline)
exp_list = struct('folder', {}, 'date_str', {}, ...
    'treatment', {}, 'cell_type', {});

for treatment = ["control", "ttl"]
    for cell_type = ["ON", "OFF"]
        subdir = fullfile(data_root, treatment, cell_type);
        if ~isfolder(subdir), continue; end

        dd = dir(subdir);
        dd = dd([dd.isdir]);
        dd = dd(~startsWith({dd.name}, '.'));

        for k = 1:numel(dd)
            exp_folder = fullfile(subdir, dd(k).name);
            if ~isfile(fullfile(exp_folder, 'currentExp.mat'))
                continue;
            end
            ei.folder    = exp_folder;
            ei.date_str  = dd(k).name;
            ei.treatment = char(treatment);
            ei.cell_type = char(cell_type);
            exp_list(end + 1) = ei; %#ok<SAGROW>
        end
    end
end

fprintf('Found %d experiments\n\n', numel(exp_list));

%% Process each experiment
orig_dir = pwd;

for exp_idx = 1:numel(exp_list)
    ei = exp_list(exp_idx);

    % Build group label
    if strcmpi(ei.cell_type, 'ON') && strcmpi(ei.treatment, 'control')
        group = 'on_control';
    elseif strcmpi(ei.cell_type, 'ON') && strcmpi(ei.treatment, 'ttl')
        group = 'on_ttl';
    elseif strcmpi(ei.cell_type, 'OFF') && strcmpi(ei.treatment, 'control')
        group = 'off_control';
    else
        group = 'off_ttl';
    end

    fprintf('[%d/%d] %s  (%s) ...', exp_idx, numel(exp_list), ei.date_str, group);

    try
        % Load data
        cd(ei.folder);
        [~, ~, Log, ~, ~] = load_protocol2_data(ei.folder);
        cd(orig_dir);

        f_data = Log.ADC.Volts(1, :);
        v_data = Log.ADC.Volts(2, :) * 10;
        median_v = median(v_data);

        % Load currentExp for LUT verification
        ce = load(fullfile(ei.folder, 'currentExp.mat'), ...
            'pattern_order', 'func_order', 'metadata');

        % Get LUT directions
        [lut_directions, ~, ~, ~] = ...
            verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, plot_order);

        % Parse bar data
        bar_data = parse_bar_data_pre_bf(f_data, v_data);

        % Correct dark-bar polarity swap for pre-Oct-15 OFF cells
        is_off = strcmpi(ei.cell_type, 'OFF');
        [bar_data, polarity_corrected] = correct_off_polarity_swap( ...
            bar_data, ce.pattern_order, ce.func_order, ei.date_str, is_off);
        if polarity_corrected
            fprintf(' [polarity-corrected]');
        end

        % Compute bar sweep responses (max_v for polar plot)
        max_v = compute_bar_sweep_responses(bar_data, plot_order, sweep_opts);

        % --- Compute DSI for annotation ---
        lut_dirs_ordered = lut_directions(plot_order);
        [~, sort_idx] = sort(lut_dirs_ordered);
        max_v_sorted = max_v(sort_idx);
        max_v_polar = [max_v_sorted; max_v_sorted(1)];
        [d_aligned, ~, ~, ~, dir_fwhm, ~, ~, ~] = ...
            find_PD_and_order_idx(max_v_polar, 0);
        [~, dsi_vector, dsi_pdnd] = compute_bar_response_metrics(d_aligned);

        % PD direction from vector sum
        theta_rad = deg2rad(lut_dirs_ordered);
        responses = max_v(:);
        angle_rad = atan2(sum(responses .* sin(theta_rad(:))), ...
                          sum(responses .* cos(theta_rad(:))));
        pd_direction = mod(rad2deg(angle_rad), 360);

        % --- Build title ---
        title_str = sprintf('28 dps — %s — %s / %s', ...
            strrep(ei.date_str, '_', '-'), ei.cell_type, ei.treatment);

        % --- Create polar ring-of-traces figure ---
        fig = plot_slow_bar_sweep_polar(bar_data, max_v, lut_directions, ...
            plot_order, median_v, title_str);
        set(fig, 'Visible', 'off');

        % --- Annotate: PD line (red) on central polar plot ---
        pax = findobj(fig, 'Type', 'polaraxes');
        if ~isempty(pax)
            hold(pax, 'on');
            pd_rad = deg2rad(pd_direction);
            r_max = max(max_v) * 0.95;
            polarplot(pax, [pd_rad, pd_rad + pi], [r_max, r_max], ...
                '-', 'Color', [0.8 0 0], 'LineWidth', 1.5);
        end

        % --- Annotate: DSI + FWHM text box ---
        lines = {};
        lines{end+1} = sprintf('DSI_{vec} = %.2f', dsi_vector);
        lines{end+1} = sprintf('DSI_{pdnd} = %.2f', dsi_pdnd);
        lines{end+1} = sprintf('FWHM = %.0f°', dir_fwhm);
        annotation(fig, 'textbox', [0.01 0.01 0.25 0.12], ...
            'String', strjoin(lines, '\n'), ...
            'FitBoxToText', 'on', 'FontSize', 9, 'FontWeight', 'bold', ...
            'EdgeColor', [0.7 0.7 0.9], 'BackgroundColor', [0.95 0.95 1]);

        % --- Save ---
        fname = sprintf('sweep_%s_%s.png', group, ei.date_str);
        exportgraphics(fig, fullfile(save_dir, fname), 'Resolution', 150);
        close(fig);
        fprintf(' DSI=%.2f FWHM=%.0f°  saved\n', dsi_vector, dir_fwhm);

    catch ME
        cd(orig_dir);
        fprintf(' ERROR: %s\n', ME.message);
    end
end

cd(orig_dir);
fprintf('\nDone. PNGs saved to: %s\n', save_dir);

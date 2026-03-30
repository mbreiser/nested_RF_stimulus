function methods = compute_pd_four_methods(angles_16, responses_16)
% COMPUTE_PD_FOUR_METHODS  Compare 4 PD extraction/alignment methods.
%
%   METHODS = COMPUTE_PD_FOUR_METHODS(ANGLES_16, RESPONSES_16) computes the
%   preferred direction and aligned tuning curve using four methods:
%
%     vec_nn     — Vector sum, snap to nearest 22.5 deg grid (current method)
%     sym_nn     — Symmetry axis from 16 candidates, snap to 22.5 deg grid
%     vec_interp — Vector sum (continuous PD), interpolate to 2.5 deg grid
%     sym_interp — Interpolate to 2.5 deg, brute-force symmetry axis
%
%   INPUTS:
%     angles_16    - 16x1 vector of stimulus directions in radians (0 to 2*pi)
%     responses_16 - 16x1 vector of peak depolarization responses (mV)
%
%   OUTPUT:
%     methods - struct with fields: vec_nn, sym_nn, vec_interp, sym_interp
%               Each field is a struct with:
%                 .pd_angle   - PD direction in radians (0 to 2*pi)
%                 .pd_deg     - PD direction in degrees (0 to 360)
%                 .aligned    - Nx2 [aligned_angles, responses], PD at pi/2
%                 .fwhm_deg   - FWHM of aligned tuning curve (degrees)
%                 .sym_score  - symmetry score (0 = perfect, normalized SSD)
%
%   ALGORITHM DETAILS:
%
%   VecNN (Method 1):
%     Standard vector sum: PD = atan2(sum(R*sin(th)), sum(R*cos(th))).
%     Snap to nearest of 16 evenly spaced grid points (22.5 deg spacing).
%     Rotate all angles so snapped PD maps to pi/2.
%     Output: 16x2 on the 22.5 deg grid.
%
%   SymNN (Method 2):
%     Test each of 16 directions as a symmetry axis. Fold the tuning curve:
%     pair each direction with its reflection across the axis. Compute SSD
%     between paired responses. Select axis with minimum SSD.
%     Resolve 180 deg ambiguity by choosing the end closer to the vector
%     sum PD from Method 1.
%     Snap and rotate identically to Method 1. Output: 16x2.
%
%   VecInterp (Method 3):
%     Interpolate 16-point curve to 144-point (2.5 deg) using circular
%     linear interpolation. Compute vector sum on original 16 points for
%     continuous PD angle. Rotate the 144-point curve so PD maps to pi/2.
%     Output: 144x2 on the 2.5 deg grid.
%
%   SymInterp (Method 4):
%     Interpolate to 144-point. Search all 144 candidate axes for the one
%     that minimizes fold SSD. Resolve 180 deg ambiguity with vector sum.
%     Rotate 144-point curve. Output: 144x2 on the 2.5 deg grid.
%
%   See also FIND_PD_AND_ORDER_IDX, COMPUTE_FWHM, COMPUTE_BAR_RESPONSE_METRICS

    angles_16 = angles_16(:);
    responses_16 = responses_16(:);

    %% --- Common: vector sum PD (continuous) ---
    vec_pd = compute_vector_sum_pd(angles_16, responses_16);

    %% --- Common: interpolated curve ---
    [fine_angles, fine_responses] = interpolate_tuning_curve(angles_16, responses_16);

    %% --- Method 1: VecNN ---
    snapped_pd = snap_to_grid(vec_pd, angles_16);
    aligned_16 = rotate_and_align(angles_16, responses_16, snapped_pd, pi/2);
    methods.vec_nn.pd_angle  = snapped_pd;
    methods.vec_nn.pd_deg    = rad2deg(snapped_pd);
    methods.vec_nn.aligned   = aligned_16;
    methods.vec_nn.fwhm_deg      = compute_fwhm_from_aligned(aligned_16);
    methods.vec_nn.sym_score     = compute_symmetry_score(aligned_16);
    methods.vec_nn.dsi_pdnd      = compute_dsi_pdnd(aligned_16);
    methods.vec_nn.aspect_ratio  = compute_aspect_ratio(aligned_16);

    %% --- Method 2: SymNN ---
    sym_axis = find_symmetry_axis(angles_16, responses_16, angles_16);
    sym_pd   = resolve_180_ambiguity(sym_axis, vec_pd);
    snapped_sym_pd = snap_to_grid(sym_pd, angles_16);
    aligned_sym = rotate_and_align(angles_16, responses_16, snapped_sym_pd, pi/2);
    methods.sym_nn.pd_angle  = snapped_sym_pd;
    methods.sym_nn.pd_deg    = rad2deg(snapped_sym_pd);
    methods.sym_nn.aligned   = aligned_sym;
    methods.sym_nn.fwhm_deg      = compute_fwhm_from_aligned(aligned_sym);
    methods.sym_nn.sym_score     = compute_symmetry_score(aligned_sym);
    methods.sym_nn.dsi_pdnd      = compute_dsi_pdnd(aligned_sym);
    methods.sym_nn.aspect_ratio  = compute_aspect_ratio(aligned_sym);

    %% --- Method 3: VecInterp ---
    aligned_fine_vec = rotate_and_align(fine_angles, fine_responses, vec_pd, pi/2);
    methods.vec_interp.pd_angle  = vec_pd;
    methods.vec_interp.pd_deg    = rad2deg(vec_pd);
    methods.vec_interp.aligned   = aligned_fine_vec;
    methods.vec_interp.fwhm_deg      = compute_fwhm_from_aligned(aligned_fine_vec);
    methods.vec_interp.sym_score     = compute_symmetry_score(aligned_fine_vec);
    methods.vec_interp.dsi_pdnd      = compute_dsi_pdnd(aligned_fine_vec);
    methods.vec_interp.aspect_ratio  = compute_aspect_ratio(aligned_fine_vec);

    %% --- Method 4: SymInterp ---
    % Search all fine-grid angles for best symmetry axis (unconstrained).
    % Previous ±45° constraint was a workaround for an interpolation NaN bug
    % (now fixed by both-end wrapping in interpolate_tuning_curve).
    sym_axis_fine = find_symmetry_axis(fine_angles, fine_responses, fine_angles);
    sym_pd_fine   = resolve_180_ambiguity(sym_axis_fine, vec_pd);
    aligned_fine_sym = rotate_and_align(fine_angles, fine_responses, sym_pd_fine, pi/2);
    methods.sym_interp.pd_angle  = sym_pd_fine;
    methods.sym_interp.pd_deg    = rad2deg(sym_pd_fine);
    methods.sym_interp.aligned   = aligned_fine_sym;
    methods.sym_interp.fwhm_deg      = compute_fwhm_from_aligned(aligned_fine_sym);
    methods.sym_interp.sym_score     = compute_symmetry_score(aligned_fine_sym);
    methods.sym_interp.dsi_pdnd      = compute_dsi_pdnd(aligned_fine_sym);
    methods.sym_interp.aspect_ratio  = compute_aspect_ratio(aligned_fine_sym);

end


%% ========================= Helper Functions ============================

function pd = compute_vector_sum_pd(angles, responses)
% COMPUTE_VECTOR_SUM_PD  Vector sum preferred direction (continuous).

    x = sum(responses .* cos(angles));
    y = sum(responses .* sin(angles));
    pd = atan2(y, x);
    if pd < 0
        pd = pd + 2*pi;
    end

end


function snapped = snap_to_grid(angle_rad, grid_angles)
% SNAP_TO_GRID  Find nearest grid angle to the target angle.

    d = abs(circ_dist_local(grid_angles, angle_rad));
    [~, idx] = min(d);
    snapped = grid_angles(idx);

end


function d = circ_dist_local(a, b)
% CIRC_DIST_LOCAL  Signed circular distance between angles.

    d = angle(exp(1i*a) ./ exp(1i*b));

end


function [fine_angles, fine_responses] = interpolate_tuning_curve(angles, responses)
% INTERPOLATE_TUNING_CURVE  Circular linear interpolation to 2.5 deg grid.

    % Wrap both ends to ensure full [0, 2*pi) coverage for interpolation.
    % Without the prepended point, cells whose sorted angles start at 22.5°
    % (not 0°) leave query points 0°-20° outside the data range → NaN.
    angles_wrap    = [angles(end) - 2*pi; angles; angles(1) + 2*pi];
    responses_wrap = [responses(end);     responses; responses(1)];

    % Fine grid: 0 to 2*pi - step, with step = 2.5 deg = pi/72
    n_fine = 144;
    fine_angles = linspace(0, 2*pi, n_fine + 1)';
    fine_angles = fine_angles(1:end-1);  % 0, 2.5, 5, ... 357.5 deg

    fine_responses = interp1(angles_wrap, responses_wrap, fine_angles, 'linear');

end


function aligned = rotate_and_align(angles, responses, pd_angle, target)
% ROTATE_AND_ALIGN  Rotate angles so PD maps to target, sort by angle.

    offset = target - pd_angle;
    aligned_angles = mod(angles + offset, 2*pi);
    aligned = sortrows([aligned_angles, responses], 1);

end


function axis_angle = find_symmetry_axis(angles, responses, candidate_axes)
% FIND_SYMMETRY_AXIS  Find the axis of symmetry that minimizes fold SSD.
%
%   Tests each candidate axis. For each, reflects every angle across the
%   axis and looks up the mirror response via interpolation. Computes SSD.

    n_cand = numel(candidate_axes);
    ssd = NaN(n_cand, 1);

    % Build interpolation function (circular, both-end wrap)
    angles_wrap    = [angles(end) - 2*pi; angles; angles(1) + 2*pi];
    responses_wrap = [responses(end);     responses; responses(1)];

    for c = 1:n_cand
        ax = candidate_axes(c);
        % Reflect each angle across the axis: mirror = 2*ax - angle
        mirror_angles = mod(2*ax - angles, 2*pi);
        % Look up mirror responses via interpolation
        mirror_responses = interp1(angles_wrap, responses_wrap, mirror_angles, 'linear');
        % SSD between original and mirrored
        ssd(c) = sum((responses - mirror_responses).^2);
    end

    [~, best_idx] = min(ssd);
    axis_angle = candidate_axes(best_idx);

end


function pd = resolve_180_ambiguity(axis_angle, vec_pd)
% RESOLVE_180_AMBIGUITY  Pick the end of the symmetry axis closer to vec PD.

    candidate_1 = axis_angle;
    candidate_2 = mod(axis_angle + pi, 2*pi);

    d1 = abs(circ_dist_local(candidate_1, vec_pd));
    d2 = abs(circ_dist_local(candidate_2, vec_pd));

    if d1 <= d2
        pd = candidate_1;
    else
        pd = candidate_2;
    end

end


function fwhm = compute_fwhm_from_aligned(aligned)
% COMPUTE_FWHM_FROM_ALIGNED  FWHM from aligned tuning curve (any grid).

    angles    = aligned(:, 1);
    responses = aligned(:, 2);

    R_max = max(responses);
    if R_max <= 0
        fwhm = NaN;
        return;
    end

    half_max = R_max / 2;
    above = responses >= half_max;

    if ~any(above)
        fwhm = NaN;
        return;
    end

    idx = find(above);
    ang1 = angles(idx(1));
    ang2 = angles(idx(end));

    fwhm = rad2deg(ang2 - ang1);
    if fwhm < 0
        fwhm = fwhm + 360;
    end

end


function score = compute_symmetry_score(aligned)
% COMPUTE_SYMMETRY_SCORE  Normalized SSD of fold around PD (at pi/2).
%   Returns 0 for perfectly symmetric, values approaching 1 for asymmetric.

    angles    = aligned(:, 1);
    responses = aligned(:, 2);

    % Build circular interpolation for mirror lookup
    angles_wrap    = [angles(end) - 2*pi; angles; angles(1) + 2*pi];
    responses_wrap = [responses(end);     responses; responses(1)];

    % Fold around pi/2: mirror angle = 2*(pi/2) - angle = pi - angle
    mirror_angles = mod(pi - angles, 2*pi);
    mirror_resp   = interp1(angles_wrap, responses_wrap, mirror_angles, 'linear');

    ssd = sum((responses - mirror_resp).^2);
    ss_total = sum(responses.^2);

    if ss_total > 0
        score = ssd / ss_total;
    else
        score = NaN;
    end

end


function dsi = compute_dsi_pdnd(aligned)
% COMPUTE_DSI_PDND  DSI from aligned curve: (R_PD - R_ND)/(R_PD + R_ND).
%   PD is at pi/2, ND is at 3*pi/2. Uses interpolation for fine grids.

    angles    = aligned(:, 1);
    responses = aligned(:, 2);

    % Build circular interpolation (both-end wrap)
    angles_wrap    = [angles(end) - 2*pi; angles; angles(1) + 2*pi];
    responses_wrap = [responses(end);     responses; responses(1)];

    r_pd = interp1(angles_wrap, responses_wrap, pi/2, 'linear');
    r_nd = interp1(angles_wrap, responses_wrap, 3*pi/2, 'linear');

    denom = r_pd + r_nd;
    if denom > 0
        dsi = (r_pd - r_nd) / denom;
    else
        dsi = NaN;
    end

end


function ar = compute_aspect_ratio(aligned)
% COMPUTE_ASPECT_RATIO  Aspect ratio from aligned curve: PD / mean(ortho).
%   PD at pi/2; orthogonal at 0 and pi. Uses interpolation for fine grids.

    angles    = aligned(:, 1);
    responses = aligned(:, 2);

    % Build circular interpolation (both-end wrap)
    angles_wrap    = [angles(end) - 2*pi; angles; angles(1) + 2*pi];
    responses_wrap = [responses(end);     responses; responses(1)];

    r_pd     = interp1(angles_wrap, responses_wrap, pi/2, 'linear');
    r_ortho1 = interp1(angles_wrap, responses_wrap, 0, 'linear');
    r_ortho2 = interp1(angles_wrap, responses_wrap, pi, 'linear');

    ortho_mean = mean([r_ortho1, r_ortho2]);
    if ortho_mean > 0
        ar = r_pd / ortho_mean;
    else
        ar = NaN;
    end

end

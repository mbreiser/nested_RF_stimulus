function m6 = compute_m6_centroid(A, area_fraction)
% COMPUTE_M6_CENTROID  68%-area bump centroid (M6 method).
%
%   M6 = COMPUTE_M6_CENTROID(A, AREA_FRACTION) finds the contiguous
%   region of length-N positive amplitude vector A that contains at
%   least AREA_FRACTION of the total positive area, by greedy expansion
%   from the peak (at each step add the adjacent position with the
%   larger amplitude).  Returns a struct M6 with fields:
%       centroid     - fractional centroid of the bump (1..N)
%       centroid_int - rounded integer centroid (1..N)
%       bump_range   - [left, right] inclusive integer bounds
%       bump_width   - right - left + 1
%
%   Used by BATCH_ANALYZE_1DRF to align flash trace columns by the
%   spatial centroid of each cell's depolarization profile.

    n_pos = numel(A);
    [peak_val, peak_pos] = max(A);
    total_area = sum(A);

    if peak_val <= 0 || total_area <= 0
        m6.centroid     = 6;
        m6.centroid_int = 6;
        m6.bump_range   = [1, n_pos];
        m6.bump_width   = n_pos;
        return;
    end

    target = area_fraction * total_area;
    left  = peak_pos;
    right = peak_pos;
    cum   = A(peak_pos);

    while cum < target && (left > 1 || right < n_pos)
        can_left  = (left > 1);
        can_right = (right < n_pos);
        if can_left && can_right
            if A(left - 1) >= A(right + 1)
                left = left - 1;  cum = cum + A(left);
            else
                right = right + 1; cum = cum + A(right);
            end
        elseif can_left
            left = left - 1;  cum = cum + A(left);
        else
            right = right + 1; cum = cum + A(right);
        end
    end

    bump   = left:right;
    A_bump = A(bump);
    m6.centroid     = sum(A_bump(:) .* bump(:)) / sum(A_bump);
    m6.centroid_int = round(m6.centroid);
    m6.centroid_int = max(1, min(n_pos, m6.centroid_int));
    m6.bump_range   = [left, right];
    m6.bump_width   = right - left + 1;

end

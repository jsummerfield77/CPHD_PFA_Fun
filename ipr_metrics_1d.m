function metrics = ipr_metrics_1d(mag_prof, axis_m, n_sidelobes)
% mag_prof : 1D magnitude (NOT dB), column or row
% axis_m   : same length, position in meters
% n_sidelobes : how many sidelobes to integrate on EACH side

if nargin < 3
    n_sidelobes = 5;
end

mag_prof = mag_prof(:).';     % make row
axis_m   = axis_m(:).';       % make row

% Peak
[peak_val, peak_idx] = max(mag_prof);
peak_pow = peak_val.^2;

% ----- 3 dB (half-power) width -----
hp_level = peak_val / sqrt(2);    % -3 dB in amplitude
% left crossing
left_idx = peak_idx;
while left_idx > 1 && mag_prof(left_idx) > hp_level
    left_idx = left_idx - 1;
end
% right crossing
right_idx = peak_idx;
while right_idx < numel(mag_prof) && mag_prof(right_idx) > hp_level
    right_idx = right_idx + 1;
end

% linear interpolate boundaries
if left_idx == 1
    x_left = axis_m(1);
else
    x1 = axis_m(left_idx); x2 = axis_m(left_idx+1);
    y1 = mag_prof(left_idx); y2 = mag_prof(left_idx+1);
    t = (hp_level - y1) / (y2 - y1);
    x_left = x1 + t*(x2 - x1);
end

if right_idx == numel(mag_prof)
    x_right = axis_m(end);
else
    x1 = axis_m(right_idx-1); x2 = axis_m(right_idx);
    y1 = mag_prof(right_idx-1); y2 = mag_prof(right_idx);
    t = (hp_level - y1) / (y2 - y1);
    x_right = x1 + t*(x2 - x1);
end

hpbw_m = x_right - x_left;

% ----- find sidelobes (peak-based) -----
% define mainlobe samples: everything between x_left and x_right
main_mask = axis_m >= x_left & axis_m <= x_right;
main_power = sum( mag_prof(main_mask).^2 );

% we’ll scan both sides for local maxima
pwr = mag_prof.^2;
idxs = 1:numel(pwr);

% left side (toward smaller indices)
left_sig = pwr(1:peak_idx-1);
left_lobes = [];
for k = peak_idx-1:-1:2
    if pwr(k) > pwr(k-1) && pwr(k) > pwr(k+1)
        left_lobes = [left_lobes, k];
    end
end

% right side
right_sig = pwr(peak_idx+1:end);
right_lobes = [];
for k = peak_idx+1:numel(pwr)-1
    if pwr(k) > pwr(k-1) && pwr(k) > pwr(k+1)
        right_lobes = [right_lobes, k];
    end
end

% PSLR: strongest sidelobe (either side)
all_lobes = [left_lobes, right_lobes];
if isempty(all_lobes)
    pslr_db = -Inf;
else
    [~, i_max_sl] = max(pwr(all_lobes));
    sl_pow = pwr(all_lobes(i_max_sl));
    pslr_db = 10*log10(sl_pow / peak_pow);
end

% ----- ISLR over first N sidelobes each side -----
% For ISLR we just sum power in first N peaks each side (if there aren’t N, we take what we have)
left_lobes  = fliplr(left_lobes);   % closest first
right_lobes = right_lobes;          % already from center out

left_lobes  = left_lobes(1:min(n_sidelobes, numel(left_lobes)));
right_lobes = right_lobes(1:min(n_sidelobes, numel(right_lobes)));

% For each lobe, integrate halfway to neighbors
isl_power = 0;

% helper to get integration bounds between two indices
    function [i1,i2] = bounds_around(k)
        % k is index of a local max
        % integrate from midpoint to left neighbor to midpoint to right neighbor
        i1 = k; i2 = k;
        % left bound
        if k > 1
            i1 = floor((k + (k-1))/2);
        end
        % right bound
        if k < numel(pwr)
            i2 = ceil((k + (k+1))/2);
        end
    end

for k = left_lobes
    [i1,i2] = bounds_around(k);
    isl_power = isl_power + sum(pwr(i1:i2));
end
for k = right_lobes
    [i1,i2] = bounds_around(k);
    isl_power = isl_power + sum(pwr(i1:i2));
end

islr_db = 10*log10( isl_power / main_power );

% pack
metrics = struct();
metrics.HPBW_m  = hpbw_m;
metrics.PSLR_dB = pslr_db;
metrics.ISLR_dB = islr_db;
metrics.peakVal = peak_val;
metrics.peakIdx = peak_idx;
end

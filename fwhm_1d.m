function bw = fwhm_1d(x, y)
%FWHM_1D  Full width at half max (3-dB) for a 1-D pattern in dB.
%
% x : vector of angle parameter (e.g., direction cosine)
% y : vector of pattern values in dB
%
% bw : full width (same units as x); NaN if cannot find 3-dB points.

    x = x(:).';
    y = y(:).';

    [ymax, imax] = max(y);
    target = ymax - 3;  % 3-dB down

    % Search on positive side only (assumes symmetry about x = 0)
    % Use index of xmax as "center"
    % Right side
    ir = find(y(imax:end) <= target, 1, 'first');
    if isempty(ir) || imax+ir-1 <= 1
        bw = NaN;
        return;
    end
    ir = imax + ir - 1;

    % Left side
    il = find(y(1:imax) <= target, 1, 'last');
    if isempty(il) || il >= numel(x)
        bw = NaN;
        return;
    end

    % Interpolate for more accurate crossing locations
    xr = interp1(y([ir-1 ir]), x([ir-1 ir]), target, 'linear');
    xl = interp1(y([il il+1]), x([il il+1]), target, 'linear');

    bw = xr - xl;
end

function y_interp = sinc_interp_Ntap(y_original, x_original, x_interp, N, use_window)
% N-tap, 1-D sinc interpolation for UNIFORMLY spaced x_original.
% Safely handles edges and out-of-band queries.

    if nargin < 5
        use_window = true;
    end

    y_original = y_original(:).';   % row
    x_original = x_original(:).';   % row

    M  = numel(y_original);
    dx = x_original(2) - x_original(1);
    N  = floor(N/2)*2;              % force even
    H  = N/2;

    % quick uniformity check (tolerant)
    if max(abs(diff(x_original) - dx)) > 1e-6*abs(dx)
        % fall back to interp1 if grid not uniform enough
        y_interp = interp1(x_original, y_original, x_interp, 'linear', 0);
        return;
    end

    y_interp = zeros(size(x_interp));
    if use_window
        base_win = hann(N+1).';
    end

    for k = 1:numel(x_interp)
        % fractional index (1-based)
        t = (x_interp(k) - x_original(1))/dx + 1;

        % if clearly out of band, return 0
        if t < 1 - H || t > M + H
            y_interp(k) = 0;
            continue;
        end

        c = floor(t);

        i_start = c - H;
        i_end   = c + H;

        i_start_clipped = max(1, i_start);
        i_end_clipped   = min(M, i_end);

        idx = i_start_clipped:i_end_clipped;
        x_samp = x_original(idx);

        u = (x_interp(k) - x_samp) / dx;
        w = sinc(u);

        if use_window
            L = numel(idx);
            if L == (N+1)
                win = base_win;
            else
                win = hann(L).';
            end
            w = w .* win;
        end

        s = sum(w);
        if s ~= 0
            w = w / s;
        end

        y_interp(k) = sum(y_original(idx) .* w);
    end
end

function [im_corr, phi_hist] = pga_autofocus(im_in, Nit, opts)
% PGA_AUTOFOCUS  Phase Gradient Autofocus for complex SAR images
%
%   [im_corr, phi_hist] = pga_autofocus(im_in, Nit, opts)
%
% INPUTS
%   im_in : complex SAR image, size [Naz x Nr]  (rows = azimuth, cols = range)
%   Nit   : number of PGA iterations (e.g., 4–10)
%   opts  : struct of options (all optional)
%       .Ktargets  : number of bright, isolated targets to use per iter (default 64)
%       .patch     : odd chip height (az) in pixels for each target (default 33)
%       .width     : chip width (range) in pixels for each target (default 33)
%       .minSep    : minimum separation (L_inf) between targets (default 16)
%       .apodize   : logical, apply Hann along azimuth within patch (default true)
%       .taperFull : logical, apply Hann along full image azimuth before FFT (default false)
%       .robustW   : exponent for magnitude weighting of gradients (default 1.0)
%       .recalcPeaksEachIter : logical, re-pick peaks each iter (default true)
%       .show      : logical, plot per-iter diagnostics (default false)
%
% OUTPUTS
%   im_corr : autofocus-corrected complex image (same size as input)
%   phi_hist: [Naz x Nit] phase error estimates (per-iteration, over azimuth freq bins)
%
% NOTES
%  - Implements a classical image-domain PGA: estimate a global 1-D azimuth
%    phase error by averaging phase-gradient estimates from K bright,
%    isolated point-like chips; then apply a frequency-domain phase
%    correction shared by all columns (range bins). Iterate to convergence.
%  - Removes zero-mean and linear (group delay) components from φ(f) each
%    iteration to avoid trivial shift/tilt terms.
%
% Reference (conceptual): Jakowatz et al., Spotlight-Mode SAR (PGA concept).

% ---------------- Defaults ----------------
if nargin < 2 || isempty(Nit),  Nit = 6; end
if nargin < 3,  opts = struct; end
def.Ktargets  = 64;
def.patch     = 33;      % azimuth chip height (odd)
def.width     = 33;      % range chip width (odd)
def.minSep    = 16;      % nonmax suppression radius (L_inf)
def.apodize   = true;
def.taperFull = false;
def.robustW   = 1.0;     % magnitude^robustW
def.recalcPeaksEachIter = true;
def.show      = false;
opts = setdefaults(opts, def);

% ---------------- Shapes ----------------
[Naz, Nr] = size(im_in);
im_corr   = im_in;
phi_hist  = zeros(Naz, Nit, 'double');  % store per-iter phase errors in az-freq domain

% Precompute azimuth FFT frequency index mapping (Naz bins)
% We keep phi on FFT-shifted frequency ordering to match fftshift usage.
azFreqIdx = (-floor(Naz/2)):(ceil(Naz/2)-1);
azFreqIdx = azFreqIdx(:);

% Optional full-image azimuth taper to reduce spectral leakage
if opts.taperFull
    han_full = hann(Naz, 'periodic');
    im_corr  = (han_full .* ones(Naz,1)) .* im_corr;
end

% ---------------- Iterations ----------------
for it = 1:Nit
    % Choose targets (on current image magnitude)
    if it == 1 || opts.recalcPeaksEachIter
        K = opts.Ktargets;
        pCenters = pick_isolated_peaks(abs(im_corr), K, opts.minSep);
    end

    % Estimate azimuth-phase gradient across frequency bins
    grad_sum   = zeros(Naz-1, 1);  % gradient length = Naz-1 (diff along freq)
    weight_sum = zeros(Naz-1, 1);

    for k = 1:size(pCenters,1)
        r0 = pCenters(k,1);
        c0 = pCenters(k,2);

        % Extract chip (odd sizes)
        halfA = floor(opts.patch/2);
        halfR = floor(opts.width/2);

        r1 = max(1,     r0 - halfA);
        r2 = min(Naz,   r0 + halfA);
        c1 = max(1,     c0 - halfR);
        c2 = min(Nr,    c0 + halfR);

        chip = im_corr(r1:r2, c1:c2);      % [Ha x Wa], Ha ~ patch height
        Ha   = size(chip,1);

        % Apodize along azimuth (rows) inside chip
        if opts.apodize
            win = hann(Ha, 'periodic');
            chip = (win .* ones(Ha,1)) .* chip;
        end

        % FFT along azimuth, shift to centered frequency
        ChipSpec = fftshift(fft(chip, [], 1), 1);           % [Ha x Wa]
        magSpec  = abs(ChipSpec);
        angSpec  = unwrap(angle(ChipSpec), [], 1);          % unwrap along azimuth frequency

        % Differentiate phase along azimuth frequency (row-wise diff)
        % Use magnitude^robustW as weights; then average across range columns
        dphi    = diff(angSpec, 1, 1);                      % [(Ha-1) x Wa]
        W       = max(magSpec(1:end-1,:), eps).^opts.robustW;

        g_est   = sum(W .* dphi, 2) ./ (sum(W,2) + eps);    % [(Ha-1) x 1], weighted average across columns

        % Accumulate into central Naz-1 support (embed Ha-1 into Naz-1 centered)
        % Map chip frequency bins (Ha) centered into global (Naz)
        idxGlob = center_embed_indices(Ha, Naz-1);          % indices in [1, Naz-1]
        grad_sum(idxGlob)   = grad_sum(idxGlob)   + g_est;
        weight_sum(idxGlob) = weight_sum(idxGlob) + ones(numel(g_est),1);
    end

    % Average gradient across chips
    g_avg = zeros(Naz-1,1);
    msk   = weight_sum > 0;
    g_avg(msk) = grad_sum(msk) ./ weight_sum(msk);

    % Integrate gradient to get φ(f) up to an additive constant
    phi_f = [0; cumsum(g_avg)];              % length Naz, on FFT-shifted frequency axis (conceptually)

    % Remove DC and linear terms (best-fit a + b*f and subtract)
    phi_f = detrend_linear_phase(phi_f, azFreqIdx);

    % Save and apply
    phi_hist(:,it) = phi_f(:);

    % Apply correction to whole image: FFT along azimuth, multiply by exp(-j*phi)
    IMspec  = fftshift(fft(im_corr, [], 1), 1);          % [Naz x Nr]
    corrVec = exp(-1j * phi_f);                          % [Naz x 1]
    IMspecC = corrVec .* IMspec;
    im_corr = ifft(ifftshift(IMspecC, 1), [], 1);

    if opts.show
        fprintf('PGA iter %d: peak phase std = %.3f rad\n', it, std(phi_f));
        drawnow;
    end
end

end % pga_autofocus

% =====================================================================
%                               Helpers
% =====================================================================

function out = setdefaults(opts, def)
    out = def;
    f = fieldnames(opts);
    for i = 1:numel(f)
        out.(f{i}) = opts.(f{i});
    end
end

function centers = pick_isolated_peaks(magImg, K, minSep)
% Greedy non-maximum suppression without toolboxes
% magImg: real, >=0; K: number to pick; minSep: L_inf radius to suppress
    [Naz, Nr] = size(magImg);
    centers   = zeros(0,2);
    taken     = false(Naz, Nr);

    % Flatten and sort by magnitude
    [~, idx] = sort(magImg(:), 'descend');

    for ii = 1:numel(idx)
        if size(centers,1) >= K, break; end
        [r, c] = ind2sub([Naz, Nr], idx(ii));
        if taken(r,c), continue; end

        % Accept this peak
        centers(end+1,:) = [r, c]; %#ok<AGROW>

        % Suppress L_inf neighborhood
        r1 = max(1, r - minSep); r2 = min(Naz, r + minSep);
        c1 = max(1, c - minSep); c2 = min(Nr,  c + minSep);
        taken(r1:r2, c1:c2) = true;
    end
end

function idxGlob = center_embed_indices(Ha, NzMinus1)
% Map local (Ha-1) gradient bins into global (Naz-1) bins, centered.
% Both are FFT-shift style sizes. We align centers and trim/pad symmetrically.
    if Ha < 2
        idxGlob = []; return;
    end
    Hm1 = Ha - 1;
    % centers
    cG = floor(NzMinus1/2) + 1;             % 1-based center index of global grad array
    cL = floor(Hm1/2) + 1;                   % center of local grad array
    % We embed local indices 1:HM1 around global center
    offset = (1: Hm1) - cL;
    idxGlob = cG + offset;
    % clip to valid range
    idxGlob(idxGlob < 1)        = [];
    idxGlob(idxGlob > NzMinus1) = [];
end

function phi_out = detrend_linear_phase(phi_in, freqIdx)
% Remove a + b*f trend from phi_in
    A = [ones(numel(freqIdx),1), freqIdx(:)];
    x = A \ phi_in(:);
    phi_out = phi_in(:) - A*x;
end

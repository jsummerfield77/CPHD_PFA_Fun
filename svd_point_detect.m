function [az_pk, r_pk, score_pk] = svd_point_detect(fil_sar, r, az, varargin)
%SVD_POINT_DETECT  Bright point detection in SAR images using SVD metric.
%
% [az_pk, r_pk, score_pk] = svd_point_detect(fil_sar, r, az)
% [ ... ] = svd_point_detect(..., 'Patch',64, 'Stride',16, 'MinSep',32, ...
%                              'Kmax',300, 'PctThr',99.7, 'NormPct',99.9)
%
% INPUTS
%   fil_sar : speckle-filtered complex (or magnitude) SAR image, size Naz x Nr
%   r       : 1 x Nr vector of column axis positions (range)
%   az      : 1 x Naz vector of row axis positions (azimuth)
%
% OUTPUTS
%   az_pk    : azimuth-axis positions (same units as az)
%   r_pk     : range-axis positions (same units as r)
%   score_pk : SVD scores (s(1)/sum(s)) for returned detections

% -------------------- defaults --------------------
p = inputParser;
p.addParameter('Patch',   64);
p.addParameter('Stride',  16);
p.addParameter('MinSep',  32);
p.addParameter('Kmax',    300);
p.addParameter('PctThr',  99.7);
p.addParameter('NormPct', 99.9);
p.parse(varargin{:});
patch   = p.Results.Patch;
stride  = p.Results.Stride;
minSep  = p.Results.MinSep;
Kmax    = p.Results.Kmax;
pctThr  = p.Results.PctThr;
normPct = p.Results.NormPct;

% -------------------- prep ------------------------
[Naz, Nr] = size(fil_sar);
lin = abs(fil_sar);
lin_scale = prctile(lin(:), normPct);
lin = lin ./ max(lin_scale, eps('single'));

r_starts = 1:stride:(Naz - patch + 1);
c_starts = 1:stride:(Nr  - patch + 1);
R = numel(r_starts);
C = numel(c_starts);

svd_metric = zeros(R, C, 'single');
bright_row = zeros(R, C, 'uint32');
bright_col = zeros(R, C, 'uint32');

% -------------------- sliding SVD + brightest px -------------------------
for ii = 1:R
    r0 = r_starts(ii);
    slab = lin(r0:r0+patch-1, :);
    for jj = 1:C
        c0 = c_starts(jj);
        P = slab(:, c0:c0+patch-1);

        s = svd(P, 'econ');
        svd_metric(ii, jj) = s(1) / max(sum(s), eps('single'));

        [~, idxMax] = max(P(:));
        [rLocal, cLocal] = ind2sub([patch, patch], idxMax);
        bright_row(ii, jj) = uint32(r0 + rLocal - 1);
        bright_col(ii, jj) = uint32(c0 + cLocal - 1);
    end
end

% -------------------- CFAR-like thresholding -----------------------------
tau = prctile(svd_metric(:), pctThr);
cand = svd_metric >= tau;
[ii_win, jj_win] = find(cand);
if isempty(ii_win)
    az_pk = []; r_pk = []; score_pk = [];
    return;
end

scores = svd_metric(cand);
[score_sorted, order] = sort(scores, 'descend');
ii_win = ii_win(order);
jj_win = jj_win(order);

row_c = double(bright_row(sub2ind([R C], ii_win, jj_win)));
col_c = double(bright_col(sub2ind([R C], ii_win, jj_win)));

% -------------------- non-maximum suppression ----------------------------
keep = false(numel(row_c),1);
for k = 1:numel(row_c)
    if ~keep(k)
        keep(k) = true;
        if nnz(keep) >= Kmax, break; end
        dr = row_c - row_c(k);
        dc = col_c - col_c(k);
        near = (dr.^2 + dc.^2) <= (minSep^2);
        near(k) = false;
        keep(near) = false;
    end
end

row_pk   = row_c(keep);
col_pk   = col_c(keep);
score_pk = score_sorted(keep);

% -------------------- convert to axis coords -----------------------------
row_i = max(1, min(Naz, round(row_pk)));
col_i = max(1, min(Nr, round(col_pk)));
az_pk = az(row_i);
r_pk  = r(col_i);
end

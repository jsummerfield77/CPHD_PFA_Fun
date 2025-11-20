function G_dB = array_gain_dB(dcx, dcy, G0, A)
%ARRAY_GAIN_DB  Evaluate CPHD polynomial array gain pattern (in dB)
%
%   G_dB = array_gain_dB(dcx, dcy, G0, A)
%
% Inputs:
%   dcx : direction cosine in X (scalar or vector)
%   dcy : direction cosine in Y (scalar or vector)
%   G0  : boresight gain in dB (meta.Antenna.AntPattern.GainZero)
%   A   : GainPoly matrix (meta.Antenna.AntPattern.Array.GainPoly)
%
% The polynomial is:
%   Gain(dcx,dcy) = G0 + sum_{m,n} A(m+1,n+1)*dcx^m*dcy^n
%
% Output:
%   G_dB : gain in dB for each (dcx,dcy) pair.

dcx = dcx(:);      % column vector
dcy = dcy(:).';    % row vector

% We'll compute a full 2-D grid of gain values
[M, N] = size(A);

Gpoly = zeros(length(dcx), length(dcy));

for m = 0:M-1
    for n = 0:N-1
        coeff = A(m+1, n+1);
        if coeff ~= 0
            Gpoly = Gpoly + coeff .* (dcx.^m) .* (dcy.^n);
        end
    end
end

G_dB = G0 + Gpoly;
end

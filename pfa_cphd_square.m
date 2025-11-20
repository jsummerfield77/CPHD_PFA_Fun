function [sar, dx_r, dx_az, res_r, res_az, oversample_ratio] = pfa_cphd_square(cphd_filename,oversample_ratio)
% PFA image formation from CPHD (spotlight) with an inscribed SQUARE passband.

% -------- Tunables -------------------------------------------------------
% oversample_ratio = 1.2;      % BW_pix / BW_res
Ntaps            = 12;         % taps for sinc_interp_Ntap
use_single       = true;       % use single precision for big arrays


% -------- Requirements ---------------------------------------------------
req_funcs = {'read_cphd_data','open_cphd_reader','ecef2lla','ecef2ned','sinc_interp_Ntap'};
missing = req_funcs(cellfun(@(f) exist(f,'file')==0, req_funcs));
if ~isempty(missing), error('Missing required functions on path: %s', strjoin(missing, ', ')); end

c_sol = physconst('lightspeed');

% ===== STEP 0: metadata (Antenna info) ==================================
%I need to get some Antenna information to determin the extent of the
%image scene. This will be a Nyquist check. 

% Use a CPHD reader object
rdr = open_cphd_reader(cphd_filename);

%This will open 1 of 2 types of meta data. This meta data describes the
%collection as a whole.  Later, pulse-by-pulse meta data will descxribe
%waveform frequency and platform position information.
meta = rdr.get_meta();

ap = meta.Antenna.AntPattern;
G0 = ap.GainZero;           % dB at boresight
A  = ap.Array.GainPoly;     % 2-D polynomial coefficients

% Nested function to evaluate array gain in dB at given dcx, dcy
eval_gain_dB = @(dcx, dcy) array_gain_dB(dcx, dcy, G0, A);

% ---- Azimuth cut (vary DCX, DCY = 0) ----
dcx_max = 0.05;  % direction cosine span (~ up to ~3 deg); adjust if needed
dcx = linspace(-dcx_max, dcx_max, 2001);
dcy = 0;
G_az = eval_gain_dB(dcx, dcy);   % 1 x N

% ---- Elevation cut (vary DCY, DCX = 0) ----
dcy_max = 0.05;
dcy_vec = linspace(-dcy_max, dcy_max, 2001);
dcx0 = 0;
G_el = eval_gain_dB(dcx0, dcy_vec);  % 1 x N

% ---- Compute FWHM in direction-cosine space ----
dcx_bw = fwhm_1d(dcx, G_az);      % full width in DCX
dcy_bw = fwhm_1d(dcy_vec, G_el);  % full width in DCY

% ---- Convert to angles in degrees (using asin) ----
% Half-width in DC, then convert that to angle
dcx_hw = dcx_bw / 2;
dcy_hw = dcy_bw / 2;

beam_az_FWHP = 2*asin(dcx_hw);
beam_el_FWHP = 2*asin(dcy_hw);


% ===== STEP 1: metadata (all pulses) ====================================

% [~, nbdata] = read_cphd_data(cphd_filename);
[~, nbdata] = rdr.read_cphd(1:meta.Data.Channel.NumVectors);   % [Nf x Nt]

if isempty(strfind(cphd_filename,'SP'))
    warning('Filename does not contain "SP" — ensure spotlight-like data.');
end

% ===== STEP 2: LOS in NED ===============================================
TxPos  = nbdata.TxPos;  SRPPos = nbdata.SRPPos;  Npulses = size(TxPos,1);
Pos_diff = SRPPos - TxPos;  TxRange  = sqrt(dot(Pos_diff, Pos_diff, 2));
uLOS_ecef = [Pos_diff(:,1)./TxRange, Pos_diff(:,2)./TxRange, Pos_diff(:,3)./TxRange];
srp_ecef = SRPPos(1,:);  srp_lla = ecef2lla(srp_ecef);  ref_lla = srp_lla;
uLOS_ned = ecef2ned(uLOS_ecef, ref_lla);
uLOS_ned = [uLOS_ned(:,1)./sqrt(dot(uLOS_ned,uLOS_ned,2)),uLOS_ned(:,2)./sqrt(dot(uLOS_ned,uLOS_ned,2)),uLOS_ned(:,3)./sqrt(dot(uLOS_ned,uLOS_ned,2))];
 
% ===== STEP 3: pulses needed for square passband =========================
BW_avg   = mean(nbdata.FX2 - nbdata.FX1);
fmin_avg = mean(nbdata.FX1);
grazing_ang_rad = atan2(uLOS_ned(:,3), hypot(uLOS_ned(:,1),uLOS_ned(:,2)));
Az_ang_rad = unwrap(atan2(uLOS_ned(:,2), uLOS_ned(:,1)));
delta_Az_min = BW_avg / fmin_avg;
delta_Az     = abs(Az_ang_rad - Az_ang_rad(1));
idx_use      = find(delta_Az >= delta_Az_min, 1, 'first');  if isempty(idx_use), idx_use = Npulses; end
N_used = idx_use;

%Don't need it, save some RAM
clear nbdata;

% ===== STEP 4: load only used pulses ====================================
[data_keep, nb_k] = rdr.read_cphd(1:N_used);   % [Nf x Nt]
TxPos  = nb_k.TxPos;  SRPPos = nb_k.SRPPos;
Pos_diff = SRPPos - TxPos;  TxRange  = sqrt(dot(Pos_diff, Pos_diff, 2));
uLOS_ecef = [Pos_diff(:,1)./ TxRange,Pos_diff(:,2)./ TxRange,Pos_diff(:,3)./ TxRange];
srp_ecef = SRPPos(1,:);  srp_lla  = ecef2lla(srp_ecef);  ref_lla = srp_lla;
uLOS_ned = ecef2ned(uLOS_ecef, ref_lla);  
uLOS_ned = [uLOS_ned(:,1)./ sqrt(dot(uLOS_ned,uLOS_ned,2)),uLOS_ned(:,2)./ sqrt(dot(uLOS_ned,uLOS_ned,2)),uLOS_ned(:,3)./ sqrt(dot(uLOS_ned,uLOS_ned,2))];
TxVel_ned = ecef2ned(nb_k.TxVel, ref_lla);
temp = [TxVel_ned(:,1)./ TxRange,TxVel_ned(:,2)./ TxRange,TxVel_ned(:,2)./ TxRange];
d_LOS_ned = -cross(uLOS_ned, cross(temp, uLOS_ned, 2), 2);

% ===== STEP 5: ground-plane K-space =====================================
Down_vec = [0 0 1];
R_vec  = sum(uLOS_ned,1);  R_vec  = R_vec  / norm(R_vec);
R_vec_GP = cross(Down_vec, cross(R_vec, Down_vec)); R_vec_GP = R_vec_GP / norm(R_vec_GP);
D_vec  = sum(d_LOS_ned,1); D_vec  = D_vec  / norm(D_vec);
D_vec_GP = cross(Down_vec, cross(D_vec, Down_vec)); D_vec_GP = D_vec_GP / norm(D_vec_GP);

Nf = size(data_keep,1);  Nt = size(data_keep,2);
freq = linspace(nb_k.FX1(1), nb_k.FX2(1), Nf);
grazing_ang_rad = atan2(uLOS_ned(:,3), hypot(uLOS_ned(:,1),uLOS_ned(:,2)));
Az_ang_rad = unwrap(atan2(uLOS_ned(:,2), uLOS_ned(:,1)));  Az_ang_rad = Az_ang_rad - mean(Az_ang_rad);

K_r  = (4*pi/c_sol) * (cos(Az_ang_rad) .* cos(grazing_ang_rad)) * freq;  % [Nt x Nf]
K_az = (4*pi/c_sol) * (sin(Az_ang_rad) .* cos(grazing_ang_rad)) * freq;  % [Nt x Nf]

% ----- Inscribed square passband -----
kr_min = max(K_r(:,1));  kr_max = min(K_r(:,end));  kr_bw = kr_max - kr_min;  kr_c = 0.5*(kr_max + kr_min);
idx1 = find(K_r(1,:)   <= kr_min, 1, 'first');
idx2 = find(K_r(end,:) <= kr_min, 1, 'first');
if K_az(1,idx1) < K_az(end,idx2), kaz_min = K_az(1,idx1); kaz_max = K_az(end,idx2);
else,                           kaz_min = K_az(end,idx2); kaz_max = K_az(1,idx1);
end
kaz_bw = kaz_max - kaz_min; kaz_c = 0.5*(kaz_max + kaz_min);
K_bw = 0.9 * min(kaz_bw, kr_bw);   % margin against edge artifacts

kr_rect_max  = kr_c + 0.5*K_bw;  kr_rect_min  = kr_c - 0.5*K_bw;
kaz_rect_max = kaz_c + 0.5*K_bw;  kaz_rect_min = kaz_c - 0.5*K_bw;
kr_os_max  = kr_c + 0.5*oversample_ratio*K_bw;  kr_os_min  = kr_c - 0.5*oversample_ratio*K_bw;
kaz_os_max = kaz_c + 0.5*oversample_ratio*K_bw; kaz_os_min = kaz_c - 0.5*oversample_ratio*K_bw;

% ===== STEP 6: grid size from spans =====================================
t_start = min(nb_k.TOAE1);  t_end = max(nb_k.TOAE2);
slant_rng_extent  = 0.5 * c_sol * (t_end - t_start);
%Range scene length based on time of arrival (this may be based on range gating)
ground_rng_extent = slant_rng_extent * cos(mean(grazing_ang_rad));

%now I want to check of the ground illumination ellipse reduces the ground_rng_extent
R0 = mean(TxRange);
% ground_rng = R0*cos(mean(grazing_ang_rad));
%equivalent height assuming a flat earth model
%I know it is not flat, this is an approximation
h = R0*sin(mean(grazing_ang_rad));
ellipse_near = h/tan(mean(grazing_ang_rad)+.5*beam_el_FWHP);
ellipse_far = h/tan(mean(grazing_ang_rad)-.5*beam_el_FWHP);
ellipse_extent = ellipse_far-ellipse_near; %Lengrh in m

% The range extent could be limited by range gating or by ground
% illumination ellipse.  The smallest determins the extent.
ground_rng_extent = min([ground_rng_extent ellipse_extent]);

% The ground illumination ellipse determins the az width of the scene.
%The 3dB width is R0*beam_az_FWHP but I find that 2x works better as you
%can still image why below the 3dB point.  
az_extent = R0 * 2*beam_az_FWHP; % added 100% margin


%This is the point when we calculate the Nyquist requirments for step sizes
%in K-space.  Max step size is 2*pi/extent
dk_r_max  = 2*pi / ground_rng_extent;  
dk_az_max = 2*pi / az_extent;
%We are goinf to do square pixels so we will set the step size to the min 
% of the two. 
dk_max    = min(dk_r_max, dk_az_max);


%Determine the next power of 2 that will 
Nx = 2^ceil(log(oversample_ratio * K_bw / dk_max) / log(2));  

%Uniform K-space sample spacing (These will be the IFFT positions)
Kr  = linspace(kr_os_min,  kr_os_max,  Nx);
Kaz = linspace(kaz_os_min, kaz_os_max, Nx);

if use_single, Kr = single(Kr); Kaz = single(Kaz); end
dKr = Kr(2)-Kr(1);
x_extent = 2*pi/dKr;  dx = x_extent/(Nx-1);
dx_r = dx; dx_az = dx;
res_r = 2*pi / K_bw;  res_az = res_r;

% ===== STEP 7: Interp #1 (Kr)  — PARFOR-SAFE =============================
if use_single
    data1 = complex(zeros(Nx, Nt, 'single'));
    Kaz1  = zeros(Nx, Nt, 'single');
    template_row = zeros(1, Nx, 'single');   % <— avoid reading data1 in parfor
else
    data1 = complex(zeros(Nx, Nt, 'double'));
    Kaz1  = zeros(Nx, Nt, 'double');
    template_row = zeros(1, Nx, 'double');
end

parfor ii = 1:Nt
    kr_src = K_r(ii,:);           % 1 x Nf, increasing
    s_src  = data_keep(:,ii).';   % 1 x Nf

    in = (Kr >= kr_src(1)) & (Kr <= kr_src(end));

    % PARFOR-SAFE: do NOT use 'like',data1 here
    re = template_row;  im = template_row;
    re(in) = sinc_interp_Ntap(real(s_src), kr_src, Kr(in), Ntaps, false);
    im(in) = sinc_interp_Ntap(imag(s_src), kr_src, Kr(in), Ntaps, false);

    data1(:,ii) = complex(re, im).';      % sliced write ok
    Kaz1(:,ii)  = tan(Az_ang_rad(ii)) * Kr;
end

% NaN cleanup (scalar expansion; no 'like' reads)
nan_idx = isnan(data1);
if any(nan_idx(:)), data1(nan_idx) = complex(0,0); end

% ===== STEP 8: Interp #2 (Kaz) — PARFOR-SAFE =============================
if use_single
    data2 = complex(zeros(Nx, Nx, 'single'));
else
    data2 = complex(zeros(Nx, Nx, 'double'));
end

parfor ii = 1:Nx
    min_Kaz = min(Kaz1(ii,:));
    max_Kaz = max(Kaz1(ii,:));
    mask_ii = (Kaz >= min_Kaz) & (Kaz <= max_Kaz) & (Kaz >= kaz_rect_min) & (Kaz <= kaz_rect_max);

    % No 'like' allocations, no data2 reads here
    re_ii = interp1(Kaz1(ii,:), real(data1(ii,:)), Kaz, 'pchip', 0);
    im_ii = interp1(Kaz1(ii,:), imag(data1(ii,:)), Kaz, 'pchip', 0);

    data2(:,ii) = (mask_ii .* (re_ii + 1i*im_ii)).';
end

nan_idx = isnan(data2);
if any(nan_idx(:)), data2(nan_idx) = complex(0,0); end

% ===== STEP 9: mask & IFFT ===============================================
mask_az = double((Kaz >= kaz_rect_min) & (Kaz <= kaz_rect_max));
mask_r  = double((Kr  >= kr_rect_min)  & (Kr  <= kr_rect_max));
S_k = single(mask_az.' * mask_r) .* data2;
sar = ifftshift(ifft2(S_k));
sar = sar ./ max(abs(sar(:)) + eps);

end

% rng_ReMoComp.m
function [rng_shift_sig,Ka_Range_Map,Ka_Range_Map2,chip2] = rng_ReMoComp(chip,r_c_idx,az_c_idx,Nr,Naz,dx_r,dx_az,n_chip_r,n_chip_az,res_r,res_az,zp_factor)
%chip is a point-like target within a SAR image
%r_c_idx is the rng index of the center of the chip - in SAR rng pixels
%az_c_idx is the az index of the center of the chip - in SAR az pixels
%Nr is the # of SAR rng pixels 
%Naz is the # of SAR az pixels 
%dx_r is the rng spacing between SAR pixels   
%dx_az is the az spacing between SAR pixels   
%n_chip_r is the # of Chip rng pixels
%n_chip_az is the # of Chip az pixels
%res_r is the rng resolution of SAR image based on K-space width
%res_az is the az resolution of SAR image based on K-space width

%Naz & Nr must be even
az = [-Naz/2:Naz/2-1]*dx_az;
r = [-Nr/2:Nr/2-1]*dx_r;

az_idx = [-n_chip_az/2:n_chip_az/2-1]+az_c_idx;

d_kr_chip = 2*pi/(n_chip_r*zp_factor*dx_r);
kr_chip = [-n_chip_r*zp_factor/2:n_chip_r*zp_factor/2-1]*d_kr_chip;

%K-space spacing - note not the same a resolution, and we are oversampled
d_Kaz = 2*pi/(Naz*dx_az);
BW_az = 2*pi/res_az;
%K-space index
Kaz = [-Naz/2:Naz/2-1]*d_Kaz;

d_Kr_chip = 2*pi/(n_chip_r*dx_r);


K_keep_mask = (Kaz >= -BW_az/2) & (Kaz <= BW_az/2);
K_keep_idx = find(K_keep_mask==1);
%%%%%%%%%%%%%%%%%
%Chip with zero padded in az direction
%We are going to use the entire SAR az width so that we can look at
%re-mocomping the range to a bright point.
sar_isolate = zeros(Naz,n_chip_r);
sar_isolate(az_idx,:) = chip;

%Now we are going to do a 2D FFT, but I want to zero-pad in the range
%direction.  Note: we are only using the chip's range ext (much smaller than full SAR image) 
K_2D_isolate = zeros(Naz,n_chip_r*zp_factor);
K_2D_isolate(:,[-n_chip_r/2:n_chip_r/2-1]+n_chip_r*zp_factor/2) = fft2(sar_isolate);

%This is just to display a before an after 
K_2D_isolate2 = zeros(Naz,n_chip_r*zp_factor);

%Now that we have zero padded, we want a 1D iFFT, this will result in
%Ka/Range image

Ka_Range_Map = ifft(K_2D_isolate,[],2);
Ka_Range_Map = Ka_Range_Map/max(abs(Ka_Range_Map(:)));
Ka_Range_Map2 = zeros(size(Ka_Range_Map)); 

%we want to measure the shit it will take to move the max value to the
%zero range bin of the chip
rng_shift_sig = zeros(1,Naz);
% length(K_keep_idx)
zero_idx = n_chip_r*zp_factor/2+1;
%Do this for each Ka bin
for jj=min(K_keep_idx):max(K_keep_idx),

    %Note - hand wave- we are assuming that there is only a single
    %point-target in the chip.  This can be wrong.  I hope when we average
    %the estimate over a large number of chips that most chips truly did
    %have a single point target or at least a point-target that dominated
    %all other scatterers.
    idx = find( abs(Ka_Range_Map(jj,:)) == max(abs(Ka_Range_Map(jj,:))));
    %uppsampled offset from zero range bin
    rng_shift_sig(jj) = (zero_idx-idx(1))*dx_r/zp_factor;

    %apply a linear phase shift in Kr direction
    K_2D_isolate2(jj,:) = K_2D_isolate(jj,:).*exp(-j*kr_chip*rng_shift_sig(jj));


end
Ka_Range_Map2 = ifft(K_2D_isolate2,[],2);
Ka_Range_Map2 = Ka_Range_Map2/max(abs(Ka_Range_Map2(:)));

sar_isolate2 = ifft2(K_2D_isolate2(:,[-n_chip_r/2:n_chip_r/2-1]+n_chip_r*zp_factor/2));
chip2 =  sar_isolate2(az_idx,:);
chip2 = chip2/max(abs(chip2(:)));

%wrapper.m
%This wrapper is a test case to try an ISAR-like image focus technique. SAR
%images formed from airbone platforms typically require auto-focus.  Phase
% gradient auto-focus (PGAF), is the standard.  It is robust.  The problem
% is that it is only a slow-time (Kaz to be precise) phase correction as it
% is correcting for Mo-meas error in platform position, a slow-time
% changing parameter.  Here we are dealing wih some posible jitter in the
% range gate positions due to some bad timing.  This means there can be
% shifts in fast-time (Kr to be precise).  This is normal for an ISAR image.
% ISAR processing requires a re-MoComp.  We will attempt this process here.
%
% Step: 1) Form a SAR image using CPHD data.  We are going to use a
% Vanilla polar format image formation algorithm (PFA).  It is Vanilla
% because it will not have any of the phase corrections used to extend the
% scene size of PFA. PFA is based on a plane wave model, this model is
% accurate for small scene sizes.  PFA benifits from some phase corrections
% that correct for wave curvature when dealing with large scenes.  Note:
% these phase corrections are going to obfuscate the range shifts I'm
% trying to correct.
%
% Step: 2) Point target detection.  I want to identify some point-like
% targets in the scene.  I'm going to use a method that is based on
% singular value decomposition (SVD).  This will grab point-targets, it
% will also grab some non-point-targets and thus there is some false-alarms
% here (deal with it).
%
% step: 3) For each of these targets, I'm going to try the ISAR range
% correction.  In ISAR you only focus a localized area (think ship in a
% large SAR image).  ISAR does not care if you are treating the imaging
% process as spatially variant or invariant.  PFA is treating the entire
% scene as spatially invariant (that is the plane wave assumption).  This
% might be where this attempt fails.  I'm going to average the range
% coreections (fast-time phase corecctions) over all of the spatialy
% diverse point targets.  I'm highly confident I can focus a localized
% target, the question is: can I focus an entire image be averaging the
% corrections of spatially diverse targets.
%
%
%Step: 4) Last touch will be a PGAF.  I going to try a Vanilla PGAF, but I
% may want to try using the phase errors from only these detected points.

%This makes sure my PATH includes NGA's SAR MATLAB Toolbox
addpath(genpath('~/Matlab_files/Global_functions/'));

%CPHD file: This is going to be very, very big.  Bring a computer with a
%lot of RAM.
cphd_filename = '/home/jsummerfield/Data/cphd/CAPELLA_C15_SP_CPHD_HH_20241002003230_20241002003257.cphd';


%I'm going to save some images to this directory
output_dir = './Images/C15_20241002';
if exist(output_dir)==0,
    mkdir(output_dir);
end

csvFile = [output_dir '/ipr_metrics.csv'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 1 - PFA SAR image

%I'm going to oversample the pixel positions relative to Nyquist.  This
%also means that in K-space, I'm going to have some zero-padded region outside
% of the square passband region.
oversample = 1.05;

%PFA is CPU intensive, the rest of the code is graphics card intensive.
% I'm getting fustrated with how slow generating images is when using a
% computer over a VNC connection. (softwareGL = Bad, WebGL=Bad).
%
% I have added a flag method that lets me optionally load a pre-processed
% SAR image.  This lets me generate a SAR image on one of the big servers
% with lots of RAM but lacks use of a GPU (Note this is all stupid software
% openGL stuff,  we have big GPUs that let me run CUDA but not openGL :( ).

%Flag - generate a new image or used saved data
gen_new_data_flag = 1;
data_filename = 'PFA_Square_Pixel_SAR.mat';

if(gen_new_data_flag), %New image-> bring a computer with lots of RAM
    [sar, dx_r, dx_az, res_r, res_az, osr] = pfa_cphd_square(cphd_filename,oversample);

    %save the data to a struct.
    data.sar = sar;
    data.dx_r = dx_r;
    data.dx_az = dx_az;
    data.res_r = res_r;
    data.res_az = res_az;
    data.osr = osr;
    %save the struct to a .mat file.  Use v7.3.  This is more than 1 GB
    save(data_filename, 'data', '-v7.3', '-nocompression');
    %clear the struct to save mem
    clear data;
else %-> load saved data.
    %load the data from a .mat file
    load(data_filename);
    %transfer the data from the struct to the original variable names
    sar     = data.sar;
    dx_r    = data.dx_r;
    dx_az   = data.dx_az;
    res_r   = data.res_r;
    res_az  = data.res_az;
    osr     = data.osr;
    %clear the struct to save mem
    clear data;

    fid = fopen(csvFile, 'w');
    % ---- Write header row ----
    header = [ ...
        'TargetName,' ...
        'rng_HPBW_m,' ...
        'rng_PSLR_dB,' ...
        'rng_ISLR_dB,' ...
        'rng_peakVal,' ...
        'rng_peakIdx,' ...
        'rng_RMS_PhaseErr_rad,' ...
        'az_HPBW_m,' ...
        'az_PSLR_dB,' ...
        'az_ISLR_dB,' ...
        'az_peakVal,' ...
        'az_peakIdx,' ...
        'az_RMS_PhaseErr_rad' ...
        ];

    fprintf(fid, '%s\n', header);

    %from here on continue as usual.  Best to have hardware openGL (ie laptop
    %better than big server to make a pretty picture).

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Stop complaining about computers and get to work
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %Normalize the SAR image.  When plotting max value will be 0 dB
    sar = sar/max(abs(sar(:)));

    %Size of the SAR image
    [Nr,Naz] = size(sar);

    %Range pixels
    r = [-Nr/2:Nr/2-1]*dx_r;
    %Az (X-Range) pixels
    az = [-Naz/2:Naz/2-1]*dx_az;

    %Display the SAR image with 80 dB of dynamic range.  This image will have a
    %lot of speckle (noise-like but not truly noise).
    fig1 = figure(1);
    set(fig1, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    imagesc(r,az,20*log10(abs(sar)),[-80 0]);
    colormap gray;
    axis xy equal tight;
    colorbar;
    xlabel('Range m');
    ylabel('X-Range m');
    title('SAR image dB');

    %Save the figure as a .png file
    filename = fullfile(output_dir, 'SAR_image.png');
    print(fig1, filename, '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 2 - Point target detection

    %Trim-
    %Note the SAR image for formed with square K-space and sqare pixels
    %The image is going to be limited by the range swath and zero outside it.
    %The point target detection can get fooled by the range swath edge as it is
    %edge-like (big gradient). I'm only going to let the SVD point-target
    %detection run on the interior of the image.  Probably fine in the az
    % direction but I'm going to keep it square shaped for now
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %square area +/- 5 km in Range and az
    mask_r = r>=-5e3 & r<5e3;
    mask_r_idx = find(mask_r ==1);

    mask_az = az>=-5e3 & az<5e3;
    mask_az_idx = find(mask_az ==1);


    %speckle filter -
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %May or may not use it.  Can impact SVD performace: less noise (not really noise)
    %but can distory IPR shapes.  The unfiltered images have "pure" IPR shapes
    %if the target is truly point-like
    fil_sar = speckle_filter(sar, 'radius', 2, 'ENL', 1, 'filter', 'Lee');
    fil_sar = fil_sar/max(abs(fil_sar(:)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Display the de-speckle SAR image.  Save it as a .png file.
    %These always look better.
    fig2 = figure(2);
    set(fig2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    imagesc(r,az,20*log10(abs(fil_sar)),[-80 0]);
    colormap gray;
    axis xy equal tight;
    colorbar;
    xlabel('Range m');
    ylabel('X-Range m');
    title('Speckle Filtered SAR image dB');

    filename = fullfile(output_dir, 'Speckled_Filtered_SAR_image.png');
    print(fig2, filename, '-dpng', '-r300');

    % SVD point-Target detection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %SVD on Speckled Filtered SAR image
    % [az_pk, r_pk, score_pk] = svd_point_detect(fil_sar(mask_az_idx,mask_r_idx),r(mask_r_idx),az(mask_az_idx));

    %seems to work better on the non-filtered image.  (filtered looks better though)
    [az_pk, r_pk, score_pk] = svd_point_detect(sar(mask_az_idx,mask_r_idx),r(mask_r_idx),az(mask_az_idx));


    %Overlay the detection positions ontop of the de-speckle SAR image. Save as
    %a .png file.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig3 = figure(3);
    set(fig3, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    imagesc(r,az,20*log10(abs(fil_sar)),[-80 0]);
    colormap gray;
    hold on;
    scatter(r_pk, az_pk, 20, score_pk, 'filled');
    axis xy equal tight;
    colorbar;
    xlabel('Range m');
    ylabel('X-Range m');
    title('SVD Point Detections');

    filename = fullfile(output_dir, 'SVD_Detections_SAR_image.png');
    print(fig3, filename, '-dpng', '-r300');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 3 - Range Correction.  This is the start of the fun.
    %We are going to chip out square regions near each detection.  I'm going to
    %do a non-symetric 2D FFT to get the data in the Ka/Kr domain.  I'm using
    %full image number of pixels in the Kaz direction but only a small number
    %of pixels in the Kr directions.  I zero pad in the range direction (zero
    %pad relative to the chip size not the full image size), then I do a 1D IFFT
    % in the range/Kr direction.  For each Kaz bin, I'm going to find shift the
    % brightest point to the zero-range bin.  This is the re-MoComp operation.

    %How many of the detected point targets are we going to averag over
    N_chips = 10; %For now as we debug, use like 50 when we do this for real

    %This is the size of the chip.  I'm going to make the length relative to
    %the number of sidelobes.  Remember this is square K-space so should be
    %squrea resolution
    N_SLs = 8;


    %length in meters
    chip_r_length = 2*N_SLs*res_r;
    chip_az_length = 2*N_SLs*res_az;

    %length in num of pixels
    n_chip_r  = 2*round(.5*chip_r_length/dx_r);
    n_chip_az = 2*round(.5*chip_az_length/dx_az);

    %chip pixels
    chip_r_m = [-n_chip_r/2:n_chip_r/2-1]*dx_r;
    chip_az_m = [-n_chip_az/2:n_chip_az/2-1]*dx_az;

    %zero Pad factor
    zp_factor = 2^3;

    %zero pad chip pixels
    chip_r_zp_m = [-n_chip_r*zp_factor/2:n_chip_r*zp_factor/2-1]*dx_r/zp_factor;
    chip_az_zp_m = [-n_chip_az*zp_factor/2:n_chip_az*zp_factor/2-1]*dx_az/zp_factor;

    %pixel pos diff in range
    d_kr_chip = 2*pi/(max(chip_r_zp_m)-min(chip_r_zp_m));
    %Kr K-space
    kr_chip = [-n_chip_r*zp_factor/2:n_chip_r*zp_factor/2-1]*d_kr_chip;
    %pixel pos diff in az
    d_Kaz = 2*pi/(max(az)-min(az));
    %Kaz K-space
    Kaz = [-Naz/2:Naz/2-1]*d_Kaz;
    %spatial Bandwidth in az, I'll use this for the zero area outside Ka
    %passband
    BW_az = 2*pi/res_az;

    %Use this to only process Kaz within the passband
    K_keep_mask = (Kaz >= -BW_az/2) & (Kaz <= BW_az/2);
    K_keep_idx = find(K_keep_mask==1);

    %Logic
    % max(Kaz)-min(Kaz) = oversample*BW_az
    avg_rng_shift_sig = zeros(1,Naz);



    fig4 = figure(4);
    set(fig4, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    %Loop through select number of detected point-targets
    for ii=1:N_chips,

        %Center location of the point-targets
        r_c_idx  = find(abs(r-r_pk(ii)) == min(abs(r-r_pk(ii))));
        az_c_idx = find(abs(az-az_pk(ii)) == min(abs(az-az_pk(ii))));

        %Center location of the point-targets
        r_idx = [-n_chip_r/2:n_chip_r/2-1]+r_c_idx;
        az_idx = [-n_chip_az/2:n_chip_az/2-1]+az_c_idx;

        chip_ii = sar(az_idx,r_idx);
        chip_ii = chip_ii/max(abs(chip_ii(:)));


        CHIP2 = zeros(n_chip_az*zp_factor,n_chip_r*zp_factor);
        CHIP2([-n_chip_az/2:n_chip_az/2-1]+n_chip_az*zp_factor/2,[-n_chip_r/2:n_chip_r/2-1]+n_chip_r*zp_factor/2)=...
            fftshift(fft2(chip_ii));

        chip_zp_ii = ifft2(ifftshift(CHIP2));
        chip_zp_ii = chip_zp_ii/max(abs(chip_zp_ii(:)));

        %%%%%%%%%%%%%%%%%%
        figure(4);
        clf;

        %%%%%%%%%%%%%%%%%%
        subplot(2,2,1);
        imagesc(chip_r_zp_m,chip_az_zp_m,20*log10(abs(chip_zp_ii)),[-40 0]);
        axis equal xy tight;
        xlabel('Range m');
        ylabel('X-Range m');
        title(sprintf('2D IPR %i',ii));

        zero_r_idx  = find(chip_r_zp_m==0);
        zero_az_idx = find(chip_az_zp_m==0);

        rng_ipr = chip_zp_ii(zero_az_idx,:);
        az_ipr = chip_zp_ii(:,zero_r_idx);

        rng_metrics = ipr_metrics_1d(abs(rng_ipr), chip_r_zp_m, N_SLs );
        az_metrics = ipr_metrics_1d(abs(az_ipr), chip_az_zp_m, N_SLs );

        subplot(2,2,2);
        plot(chip_r_zp_m,abs(rng_ipr));
        axis square tight;
        grid on;
        xlabel('Range m');
        ylabel('Rng IPR abs');
        title(sprintf('Range IPR %i',ii));

        subplot(2,2,3);
        plot(chip_az_zp_m,abs(az_ipr));
        axis square tight;
        grid on;
        xlabel('X-Range m');
        ylabel('X-Rng IPR abs');
        title(sprintf('X-Range IPR %i',ii));


        RNG_K = fftshift(fft(rng_ipr));
        rng_ppp = RNG_K(2:end).*RNG_K(1:end-1);
        rng_ang_ppp = unwrap(angle(rng_ppp));
        rng_ang_ppp = rng_ang_ppp - mean(rng_ang_ppp);
        rng_pe = cumsum(rng_ang_ppp);
        rms_rng_pe = sqrt(mean(rng_pe.^2));

        AZ_K = fftshift(fft(az_ipr));
        az_ppp = AZ_K(2:end).*AZ_K(1:end-1);
        az_ang_ppp = unwrap(angle(az_ppp));
        az_ang_ppp = az_ang_ppp - mean(az_ang_ppp);
        az_pe = cumsum(az_ang_ppp);
        az_rms_pe = sqrt(mean(az_pe.^2));

        subplot(2,2,4);
        axis off;   % Hide axes for clean text display

        % Convert struct fields to cell arrays of strings
        % Convert struct fields to cell arrays of strings
        rng_lines = {
            sprintf('HPBW_m:   %.4f', rng_metrics.HPBW_m)
            sprintf('PSLR_d_B:  %.4f', rng_metrics.PSLR_dB)
            sprintf('ISLR_d_B:  %.4f', rng_metrics.ISLR_dB)
            sprintf('peakVal:  %.4f', rng_metrics.peakVal)
            sprintf('peakIdx:  %d',    rng_metrics.peakIdx)
            sprintf('RMS Phase Err_r_a_d:  %.4f', rms_rng_pe)
            };

        az_lines = {
            sprintf('HPBW_m:   %.4f', az_metrics.HPBW_m)
            sprintf('PSLR_d_B:  %.4f', az_metrics.PSLR_dB)
            sprintf('ISLR_d_B:  %.4f', az_metrics.ISLR_dB)
            sprintf('peakVal:  %.4f', az_metrics.peakVal)
            sprintf('peakIdx:  %d',    az_metrics.peakIdx)
            sprintf('RMS Phase Err_r_a_d:  %.4f', az_rms_pe)
            };

        % --- R2014a-compatible join for column cell arrays ---
        rng_text = sprintf('%s\n', rng_lines{:});   % concatenates each line with newline
        az_text  = sprintf('%s\n', az_lines{:});

        % Display left column (range metrics)
        text(0.05, 0.95, 'Range Metrics', ...
            'FontWeight','bold', 'FontSize', 12);

        text(0.05, 0.85, rng_text, ...
            'VerticalAlignment','top', 'FontSize', 11);


        % Display left column (range metrics)
        text(0.55, 0.95, 'Az Metrics', ...
            'FontWeight','bold', 'FontSize', 12);

        text(0.55, 0.85, az_text, ...
            'VerticalAlignment','top', 'FontSize', 11);

        filename = fullfile(output_dir, sprintf('Chip_%i.png',ii));
        print(fig4, filename, '-dpng', '-r300');

        % ---- Construct line for this target ----
        fprintf(fid, ...
            ['PntTgt_%d,' ...                % Target name
            '%.6f,%.6f,%.6f,' ...           % Range HPBW, PSLR, ISLR
            '%.6f,%d,%.6f,' ...             % Range peakVal, peakIdx, RMS PE
            '%.6f,%.6f,%.6f,' ...           % Az HPBW, PSLR, ISLR
            '%.6f,%d,%.6f\n'], ...          % Az peakVal, peakIdx, RMS PE
            ...
            ii, ...
            rng_metrics.HPBW_m, rng_metrics.PSLR_dB, rng_metrics.ISLR_dB, ...
            rng_metrics.peakVal, rng_metrics.peakIdx, rms_rng_pe, ...
            az_metrics.HPBW_m,  az_metrics.PSLR_dB,  az_metrics.ISLR_dB, ...
            az_metrics.peakVal, az_metrics.peakIdx, az_rms_pe ...
            );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Full slow-time width but limited in range to chip index

        [rng_shift_sig,Ka_Range_Map,Ka_Range_Map2,chip2] = rng_ReMoComp(chip_ii,r_c_idx,az_c_idx,Nr,Naz,dx_r,dx_az,n_chip_r,n_chip_az,res_r,res_az,zp_factor);

        %accumulate the range shift signatures for each detect
        avg_rng_shift_sig  = avg_rng_shift_sig  + rng_shift_sig;

        figure(4);
        clf;
        subplot(221);
        imagesc(abs(Ka_Range_Map));
        colormap gray;
        subplot(222);
        imagesc(abs(Ka_Range_Map2));
        colormap gray;

        subplot(223)
        imagesc(chip_r_zp_m,chip_az_zp_m,20*log10(abs(chip_ii)),[-40 0]);
        axis equal xy tight;
        xlabel('Range m');
        ylabel('X-Range m');

        subplot(224)
        imagesc(chip_r_zp_m,chip_az_zp_m,20*log10(abs(chip2)),[-40 0]);
        axis equal xy tight;
        xlabel('Range m');
        ylabel('X-Range m');

        filename = fullfile(output_dir, sprintf('Rng_Corr_Chip_%i.png',ii));
        print(fig4, filename, '-dpng', '-r300');

    end

    %average range shift signatures is accumulation divided by num detects used
    avg_rng_shift_sig = avg_rng_shift_sig/N_chips;

    %apply the range shift correction
    %1D Range to Kr FFT
    Kr_az_map = fft(sar,[],2);
    d_kr = 2*pi/(max(r)-min(r));
    kr = [-Nr/2:Nr/2-1]*d_kr;
    %This double casting is needed for older versions of MATLAB
    %I'm trying matrix operations
    %negative sign for conjugate linear phase in the Kr direction
    phase_corr = -single(double(kr).'*double(avg_rng_shift_sig));

    %range shift corrected SAR image is 1D FFT of the Kr correction times Kr/az
    %map
    sar_rng_corr = ifft(Kr_az_map.*exp(j*phase_corr),[],2);

    fig5 = figure(5);
    set(fig5, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    imagesc(r,az,20*log10(abs(sar_rng_corr)),[-80 0]);
    colormap gray;
    axis xy equal tight;
    colorbar;
    xlabel('Range m');
    ylabel('X-Range m');
    title('Rng Corrected SAR image dB');

    %Save the figure as a .png file
    filename = fullfile(output_dir, 'Rng_corr_SAR_image.png');
    print(fig5, filename, '-dpng', '-r300');

    %Run SVD detection a second time
    [az_pk2, r_pk2, score_pk2] = svd_point_detect(sar_rng_corr(mask_az_idx,mask_r_idx),r(mask_r_idx),az(mask_az_idx));

    fig6 = figure(6);
    set(fig6, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    %Loop through select number of detected point-targets
    for ii=1:N_chips,

        %Center location of the point-targets
        r_c_idx  = find(abs(r-r_pk2(ii)) == min(abs(r-r_pk2(ii))));
        az_c_idx = find(abs(az-az_pk2(ii)) == min(abs(az-az_pk2(ii))));

        %Center location of the point-targets
        r_idx = [-n_chip_r/2:n_chip_r/2-1]+r_c_idx;
        az_idx = [-n_chip_az/2:n_chip_az/2-1]+az_c_idx;

        % The range correction can cause shifts, find the brightest point and re-center
        chip_ii = sar_rng_corr(az_idx,r_idx);
        chip_ii = chip_ii/max(abs(chip_ii(:)));

        CHIP2 = zeros(n_chip_az*zp_factor,n_chip_r*zp_factor);
        CHIP2([-n_chip_az/2:n_chip_az/2-1]+n_chip_az*zp_factor/2,[-n_chip_r/2:n_chip_r/2-1]+n_chip_r*zp_factor/2)=...
            fftshift(fft2(chip_ii));

        chip_zp_ii = ifft2(ifftshift(CHIP2));
        chip_zp_ii = chip_zp_ii/max(abs(chip_zp_ii(:)));

        %%%%%%%%%%%%%%%%%%
        figure(6);
        clf;

        %%%%%%%%%%%%%%%%%%
        subplot(2,2,1);
        imagesc(chip_r_zp_m,chip_az_zp_m,20*log10(abs(chip_zp_ii)),[-40 0]);
        axis equal xy tight;
        xlabel('Range m');
        ylabel('X-Range m');
        title(sprintf('Rng Corr 2D IPR %i',ii));

        zero_r_idx  = find(chip_r_zp_m==0);
        zero_az_idx = find(chip_az_zp_m==0);

        rng_ipr = chip_zp_ii(zero_az_idx,:);
        az_ipr = chip_zp_ii(:,zero_r_idx);

        rng_metrics = ipr_metrics_1d(abs(rng_ipr), chip_r_zp_m, N_SLs );
        az_metrics = ipr_metrics_1d(abs(az_ipr), chip_az_zp_m, N_SLs );

        subplot(2,2,2);
        plot(chip_r_zp_m,abs(rng_ipr));
        axis square tight;
        grid on;
        xlabel('Range m');
        ylabel('Rng IPR abs');
        title(sprintf('Rng Corr Range IPR %i',ii));

        subplot(2,2,3);
        plot(chip_az_zp_m,abs(az_ipr));
        axis square tight;
        grid on;
        xlabel('X-Range m');
        ylabel('X-Rng IPR abs');
        title(sprintf('Rng Corr X-Range IPR %i',ii));


        RNG_K = fftshift(fft(rng_ipr));
        rng_ppp = RNG_K(2:end).*RNG_K(1:end-1);
        rng_ang_ppp = unwrap(angle(rng_ppp));
        rng_ang_ppp = rng_ang_ppp - mean(rng_ang_ppp);
        rng_pe = cumsum(rng_ang_ppp);
        rms_rng_pe = sqrt(mean(rng_pe.^2));

        AZ_K = fftshift(fft(az_ipr));
        az_ppp = AZ_K(2:end).*AZ_K(1:end-1);
        az_ang_ppp = unwrap(angle(az_ppp));
        az_ang_ppp = az_ang_ppp - mean(az_ang_ppp);
        az_pe = cumsum(az_ang_ppp);
        az_rms_pe = sqrt(mean(az_pe.^2));

        subplot(2,2,4);
        axis off;   % Hide axes for clean text display

        % Convert struct fields to cell arrays of strings
        % Convert struct fields to cell arrays of strings
        rng_lines = {
            sprintf('HPBW_m:   %.4f', rng_metrics.HPBW_m)
            sprintf('PSLR_d_B:  %.4f', rng_metrics.PSLR_dB)
            sprintf('ISLR_d_B:  %.4f', rng_metrics.ISLR_dB)
            sprintf('peakVal:  %.4f', rng_metrics.peakVal)
            sprintf('peakIdx:  %d',    rng_metrics.peakIdx)
            sprintf('RMS Phase Err_r_a_d:  %.4f', rms_rng_pe)
            };

        az_lines = {
            sprintf('HPBW_m:   %.4f', az_metrics.HPBW_m)
            sprintf('PSLR_d_B:  %.4f', az_metrics.PSLR_dB)
            sprintf('ISLR_d_B:  %.4f', az_metrics.ISLR_dB)
            sprintf('peakVal:  %.4f', az_metrics.peakVal)
            sprintf('peakIdx:  %d',    az_metrics.peakIdx)
            sprintf('RMS Phase Err_r_a_d:  %.4f', az_rms_pe)
            };

        % --- R2014a-compatible join for column cell arrays ---
        rng_text = sprintf('%s\n', rng_lines{:});   % concatenates each line with newline
        az_text  = sprintf('%s\n', az_lines{:});

        % Display left column (range metrics)
        text(0.05, 0.95, 'Range Metrics', ...
            'FontWeight','bold', 'FontSize', 12);

        text(0.05, 0.85, rng_text, ...
            'VerticalAlignment','top', 'FontSize', 11);


        % Display left column (range metrics)
        text(0.55, 0.95, 'Az Metrics', ...
            'FontWeight','bold', 'FontSize', 12);

        text(0.55, 0.85, az_text, ...
            'VerticalAlignment','top', 'FontSize', 11);

        filename = fullfile(output_dir, sprintf('Rng_Corr_Chip_%i.png',ii));
        print(fig6, filename, '-dpng', '-r300');
        % ---- Construct line for this target ----
        fprintf(fid, ...
            ['Rng_Corr_PntTgt_%d,' ...                % Target name
            '%.6f,%.6f,%.6f,' ...           % Range HPBW, PSLR, ISLR
            '%.6f,%d,%.6f,' ...             % Range peakVal, peakIdx, RMS PE
            '%.6f,%.6f,%.6f,' ...           % Az HPBW, PSLR, ISLR
            '%.6f,%d,%.6f\n'], ...          % Az peakVal, peakIdx, RMS PE
            ...
            ii, ...
            rng_metrics.HPBW_m, rng_metrics.PSLR_dB, rng_metrics.ISLR_dB, ...
            rng_metrics.peakVal, rng_metrics.peakIdx, rms_rng_pe, ...
            az_metrics.HPBW_m,  az_metrics.PSLR_dB,  az_metrics.ISLR_dB, ...
            az_metrics.peakVal, az_metrics.peakIdx, az_rms_pe ...
            );
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 4 - Hit this with a stock/Vanilla PGAF - we can do better than this
    % if we add a pulse pair product to the loop over each chip.  But baby steps

    [sar_PGAF, ~] = pga_autofocus(sar_rng_corr, 10);


    %Run SVD detection a 3rd time
    [az_pk3, r_pk3, score_pk3] = svd_point_detect(sar_PGAF(mask_az_idx,mask_r_idx),r(mask_r_idx),az(mask_az_idx));

    fig7 = figure(7);
    set(fig7, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    imagesc(r,az,20*log10(abs(sar_PGAF)),[-80 0]);

    colormap gray;
    axis xy equal tight;
    colorbar;
    xlabel('Range m');
    ylabel('X-Range m');
    title('Rng Corrected PGAF SAR image dB');

    %Save the figure as a .png file
    filename = fullfile(output_dir, 'Rng_corr_PGAF_SAR_image.png');
    print(fig7, filename, '-dpng', '-r300');


    fig8 = figure(8);
    set(fig8, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    %Loop through select number of detected point-targets
    for ii=1:N_chips,

        %Center location of the point-targets
        r_c_idx  = find(abs(r-r_pk3(ii)) == min(abs(r-r_pk3(ii))));
        az_c_idx = find(abs(az-az_pk3(ii)) == min(abs(az-az_pk3(ii))));

        %Center location of the point-targets
        r_idx = [-n_chip_r/2:n_chip_r/2-1]+r_c_idx;
        az_idx = [-n_chip_az/2:n_chip_az/2-1]+az_c_idx;

        chip_ii = sar_PGAF(az_idx,r_idx);
        chip_ii = chip_ii/max(abs(chip_ii(:)));

        CHIP2 = zeros(n_chip_az*zp_factor,n_chip_r*zp_factor);
        CHIP2([-n_chip_az/2:n_chip_az/2-1]+n_chip_az*zp_factor/2,[-n_chip_r/2:n_chip_r/2-1]+n_chip_r*zp_factor/2)=...
            fftshift(fft2(chip_ii));

        chip_zp_ii = ifft2(ifftshift(CHIP2));
        chip_zp_ii = chip_zp_ii/max(abs(chip_zp_ii(:)));

        %%%%%%%%%%%%%%%%%%
        figure(8);
        clf;

        %%%%%%%%%%%%%%%%%%
        subplot(2,2,1);
        imagesc(chip_r_zp_m,chip_az_zp_m,20*log10(abs(chip_zp_ii)),[-40 0]);
        axis equal xy tight;
        xlabel('Range m');
        ylabel('X-Range m');
        title(sprintf('PGAF 2D IPR %i',ii));

        zero_r_idx  = find(chip_r_zp_m==0);
        zero_az_idx = find(chip_az_zp_m==0);

        rng_ipr = chip_zp_ii(zero_az_idx,:);
        az_ipr = chip_zp_ii(:,zero_r_idx);

        rng_metrics = ipr_metrics_1d(abs(rng_ipr), chip_r_zp_m, N_SLs );
        az_metrics = ipr_metrics_1d(abs(az_ipr), chip_az_zp_m, N_SLs );

        subplot(2,2,2);
        plot(chip_r_zp_m,abs(rng_ipr));
        axis square tight;
        grid on;
        xlabel('Range m');
        ylabel('Rng IPR abs');
        title(sprintf('PGAF Range IPR %i',ii));

        subplot(2,2,3);
        plot(chip_az_zp_m,abs(az_ipr));
        axis square tight;
        grid on;
        xlabel('X-Range m');
        ylabel('X-Rng IPR abs');
        title(sprintf('PGAF X-Range IPR %i',ii));


        RNG_K = fftshift(fft(rng_ipr));
        rng_ppp = RNG_K(2:end).*RNG_K(1:end-1);
        rng_ang_ppp = unwrap(angle(rng_ppp));
        rng_ang_ppp = rng_ang_ppp - mean(rng_ang_ppp);
        rng_pe = cumsum(rng_ang_ppp);
        rms_rng_pe = sqrt(mean(rng_pe.^2));

        AZ_K = fftshift(fft(az_ipr));
        az_ppp = AZ_K(2:end).*AZ_K(1:end-1);
        az_ang_ppp = unwrap(angle(az_ppp));
        az_ang_ppp = az_ang_ppp - mean(az_ang_ppp);
        az_pe = cumsum(az_ang_ppp);
        az_rms_pe = sqrt(mean(az_pe.^2));

        subplot(2,2,4);
        axis off;   % Hide axes for clean text display

        % Convert struct fields to cell arrays of strings
        % Convert struct fields to cell arrays of strings
        rng_lines = {
            sprintf('HPBW_m:   %.4f', rng_metrics.HPBW_m)
            sprintf('PSLR_d_B:  %.4f', rng_metrics.PSLR_dB)
            sprintf('ISLR_d_B:  %.4f', rng_metrics.ISLR_dB)
            sprintf('peakVal:  %.4f', rng_metrics.peakVal)
            sprintf('peakIdx:  %d',    rng_metrics.peakIdx)
            sprintf('RMS Phase Err_r_a_d:  %.4f', rms_rng_pe)
            };

        az_lines = {
            sprintf('HPBW_m:   %.4f', az_metrics.HPBW_m)
            sprintf('PSLR_d_B:  %.4f', az_metrics.PSLR_dB)
            sprintf('ISLR_d_B:  %.4f', az_metrics.ISLR_dB)
            sprintf('peakVal:  %.4f', az_metrics.peakVal)
            sprintf('peakIdx:  %d',    az_metrics.peakIdx)
            sprintf('RMS Phase Err_r_a_d:  %.4f', az_rms_pe)
            };

        % --- R2014a-compatible join for column cell arrays ---
        rng_text = sprintf('%s\n', rng_lines{:});   % concatenates each line with newline
        az_text  = sprintf('%s\n', az_lines{:});

        % Display left column (range metrics)
        text(0.05, 0.95, 'Range Metrics', ...
            'FontWeight','bold', 'FontSize', 12);

        text(0.05, 0.85, rng_text, ...
            'VerticalAlignment','top', 'FontSize', 11);


        % Display left column (range metrics)
        text(0.55, 0.95, 'Az Metrics', ...
            'FontWeight','bold', 'FontSize', 12);

        text(0.55, 0.85, az_text, ...
            'VerticalAlignment','top', 'FontSize', 11);

        filename = fullfile(output_dir, sprintf('PGAF_Chip_%i.png',ii));
        print(fig8, filename, '-dpng', '-r300');
        % ---- Construct line for this target ----
        fprintf(fid, ...
            ['PGAF_PntTgt_%d,' ...                % Target name
            '%.6f,%.6f,%.6f,' ...           % Range HPBW, PSLR, ISLR
            '%.6f,%d,%.6f,' ...             % Range peakVal, peakIdx, RMS PE
            '%.6f,%.6f,%.6f,' ...           % Az HPBW, PSLR, ISLR
            '%.6f,%d,%.6f\n'], ...          % Az peakVal, peakIdx, RMS PE
            ...
            ii, ...
            rng_metrics.HPBW_m, rng_metrics.PSLR_dB, rng_metrics.ISLR_dB, ...
            rng_metrics.peakVal, rng_metrics.peakIdx, rms_rng_pe, ...
            az_metrics.HPBW_m,  az_metrics.PSLR_dB,  az_metrics.ISLR_dB, ...
            az_metrics.peakVal, az_metrics.peakIdx, az_rms_pe ...
            );
    end

    % ---- Close file ----
    fclose(fid);

    fprintf('CSV written: %s\n', csvFile);

end


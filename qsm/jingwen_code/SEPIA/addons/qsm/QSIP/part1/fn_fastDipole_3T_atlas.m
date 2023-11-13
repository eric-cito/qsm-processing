%%
function [nfm_final, nfm_dipole, atlas, msk, x, scales, num_iter] = fn_fastDipole_3T_atlas(subjnum_sph)
%
% subjnum_sph : is a string. ie. '0011_12192007'

%% Inputs for all subjects except '0011_12192007' (all inputs are Radiological)
% atlas = read_avw('FM_atlas');
% 
% threshold atlas to remove inferior region where there was CT
% observations from only 1 subject ( zubal CT ) 
% ind = find(atlas < 0.96);
% atlas(ind) = 0;
% 
% atlas = imrotate(atlas,90);
% 
% filename = strcat(subjnum_sph,'_brain_mask_ero');
% msk = read_avw(filename);  
% msk = imrotate(double(msk),90);
% 
% filename = strcat(subjnum_sph,'_sphase_mm.img');
% magn = load_nii(filename);
% magn = imrotate(double(magn.img),90);
% magn = magn / max(magn(:));
% 
% filename = strcat(subjnum_sph,'_unwphs.img');
% uphase = load_nii(filename);
% uphase = imrotate(double(uphase.img),90);

%% Inputs for subject 0011_12192007
%atlas = read_avw('FMPRC_atlas');
atlas = read_avw('FM_atlas');
atlas = atlas(:,:,:,1);

%threshold atlas to remove inferior region where there was CT
%observations from only 1 subject ( zubal CT ) 
ind = find(atlas < 0.96);
atlas(ind) = 0;

%atlas = imrotate(atlas,90);

%msk = read_avw('FMPRC_mag_cmb_brain_mask2_ero');   
%msk = read_avw('FMPRC_mag_cmb_brain_mask'); 
%msk = read_avw('mag_cmb_brain_mask');
%msk = read_avw('FM_mag_cmb_brain_mask2_ero');
msk = read_avw('FM_mag_cmb_brain_mask2_ero_sphere20x3');
%msk = imrotate(msk,90);


%magn = read_avw('FMPRC_mag_cmb');
[magn, grot, scales] = read_avw('mag_cmb_brain');
%magn = imrotate(magn, 90);
magn = magn / max(magn(:));

%uphase = read_avw('FMPRC_phase_cmb_rescale3');
uphase = read_avw('phase_cmb');
%uphase = imrotate(uphase,90);
%%

uphase = uphase .* msk;
magn = magn .* msk;

for t = 1:size(msk,3)
    figure(1), subplot(1,4,1), imagesc( magn(:,:,t), [0,.4] ), colormap(gray), axis image off, title(num2str(t)), 
    figure(1), subplot(1,4,2), imagesc( uphase(:,:,t) ), colormap(gray), axis image off, title(num2str(t)), 
    figure(1), subplot(1,4,3), imagesc( msk(:,:,t) ), colormap(gray), axis image off, title(num2str(t)),
    figure(1), subplot(1,4,4), imagesc( atlas(:,:,t) ), colormap(gray), axis image off, title(num2str(t)), pause(.1)
end


%%
%scales = [0.8594, 0.8594, 2.0];
%scales = [0.42969, 0.42969, 2.0];  

N = size(uphase);
%FOV = N .* scales;  % (in milimeters)
FOV = N .* scales(1:3)';  % (in milimeters)

B0 = 3.0;   % Tesla
gyro = 2*pi*42.58;
TE = 28e-3;  % second

chi0 = 0.4;     % unitless, in ppm
delta = -9.5;   % unitless, in ppm

nfm = -uphase / (B0 * gyro * TE);   % normalized field map, unitless, in ppm


% pad the image to twice the size and use kernel corresponding to this FOV
nfm_pad = padarray(nfm,N/2);
msk_pad = padarray(msk,N/2);


D_spatial = ifft3c(get_aniso_kernel( FOV, N, N/2+1 ));

D = fft3c(padarray(D_spatial,N/2));

%% fast CG dipole computation

atlas_chi = chi0 * ones(N(1),N(2),N(3)) + delta * atlas;
ind_brain = find(msk == 1);
atlas_chi(ind_brain) = chi0 + delta;
%x0 = padarray(atlas_chi,N/2);
x0 = 0 * nfm_pad;

num_iter = 200;

x = CG_dipole(nfm_pad, x0, D, msk_pad, num_iter);


%% 

x_dipole = x .* (1-msk_pad);
nfm_dipole = ifft3c(D .* fft3c(x_dipole));

nfm_final = msk_pad.*(nfm_pad - nfm_dipole);
%nfm_final = nfm_final(257-128:256+128, 257-128:256+128, 65-32:64+32);
nfm_final = nfm_final((N(1)+1)-(N(1)/2):N(2)+(N(1)/2),(N(1)+1)-(N(1)/2):N(2)+(N(1)/2), (N(3)+1)-(N(3)/2):N(3)+(N(3)/2));

for t = 1:N(3)
    figure(1), imagesc( nfm_final(:,:,t), [-.1,.1] ), axis image off, colormap(gray), pause(.1)
end


%% Outputs
% nfm_final = imrotate(nfm_final,-90);
% nfm_dipole = imrotate(nfm_dipole,-90);
% x = imrotate(x,-90);
% msk = imrotate(msk,-90);


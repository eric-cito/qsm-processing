%%
function [nfm_final, nfm_dipole, atlas, mask, x, scales, num_iter, residual] = ...
    fn_fastDipole_atlas_new(atlas, mask, magni, phase_uw, params)
%
% subjnum_sph : is a string. 

%% Inputs for subjects

atlas = atlas(:,:,:,1);

% threshold atlas to remove inferior region where there was CT
% observations from only 1 subject ( zubal CT ) 

atlas(atlas < params.atlas_thr) = 0;
scales = params.voxelsize;
magni = magni / max(magni(:));

phase_uw = phase_uw .* mask;
magni = magni .* mask;

%% 

N = size(phase_uw);
FOV = N .* scales(1:3);  % (in milimeters)

B0 = params.bo;   % Tesla
gyro = 2*pi*42.58;
TE = params.te;  % second

chi0 = params.chi0;     % unitless, in ppm
delta = params.delta;   % unitless, in ppm

nfm = -phase_uw / (B0 * gyro * TE);   % normalized field map, unitless, in ppm


% pad the image to twice the size and use kernel corresponding to this FOV
nfm_pad = padarray(nfm,N/2);
msk_pad = padarray(mask,N/2);


D_spatial = ifft3c(get_aniso_kernel( FOV, N, N/2+1 ));

D = fft3c(padarray(D_spatial,N/2));

%% fast CG dipole computation

atlas_chi = chi0 * ones(N(1),N(2),N(3)) + delta * atlas;
atlas_chi(mask == 1) = chi0 + delta;

if params.atlas_flag == 1
    x0 = padarray(atlas_chi,N/2);
    disp('x0=atlas')
else 
    x0 = 0 * nfm_pad;
    disp('x0=0')
end

num_iter = params.num_iter;

%[x,residual] = CG_dipole2(nfm_pad, x0, D, msk_pad, num_iter);
[x,residual] = CG_dipole3(nfm_pad, x0, D, msk_pad, num_iter);

%% 

x_dipole = x .* (1-msk_pad);
nfm_dipole = ifft3c(D .* fft3c(x_dipole));

nfm_final = msk_pad.*(nfm_pad - nfm_dipole);

nfm_final = nfm_final((N(1)+1)-(N(1)/2):N(1)+(N(1)/2), (N(2)+1)-(N(2)/2):N(2)+(N(2)/2), (N(3)+1)-(N(3)/2):N(3)+(N(3)/2));




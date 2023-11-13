%%
function [nfm_final, nfm_dipole, atlas, msk, x, scales, num_iter, residual] = fn_fastDipole_7T_HD_atlas(fname_atlas, fname_mask, fname_mag, fname_uphase, params)
%
% subjnum_sph : is a string. 

%% Inputs for subjects

atlas = read_avw(fname_atlas);
atlas = atlas(:,:,:,1);

%threshold atlas to remove inferior region where there was CT
%observations from only 1 subject ( zubal CT ) 
ind = find(atlas < params.atlas_thr);
atlas(ind) = 0;

msk = read_avw(fname_mask);

[magn, grot, scales] = read_avw(fname_mag);
magn = magn / max(magn(:));

uphase = read_avw(fname_uphase);
%%

uphase = uphase .* msk;
magn = magn .* msk;

% for t = 1:size(msk,3)
%     figure(1), subplot(1,4,1), imagesc( magn(:,:,t), [0,.4] ), colormap(gray), axis image off, title(num2str(t)), 
%     figure(1), subplot(1,4,2), imagesc( uphase(:,:,t) ), colormap(gray), axis image off, title(num2str(t)), 
%     figure(1), subplot(1,4,3), imagesc( msk(:,:,t) ), colormap(gray), axis image off, title(num2str(t)),
%     figure(1), subplot(1,4,4), imagesc( atlas(:,:,t) ), colormap(gray), axis image off, title(num2str(t)), pause(.1)
% end


%% 

N = size(uphase);
FOV = N .* scales(1:3)';  % (in milimeters)

B0 = params.bo;   % Tesla
gyro = 2*pi*42.58;
TE = params.te;  % second

chi0 = params.chi0;     % unitless, in ppm
delta = params.delta;   % unitless, in ppm

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

nfm_final = nfm_final((N(1)+1)-(N(1)/2):N(2)+(N(1)/2),(N(1)+1)-(N(1)/2):N(2)+(N(1)/2), (N(3)+1)-(N(3)/2):N(3)+(N(3)/2));

% for t = 1:N(3)
%     figure(1), imagesc( nfm_final(:,:,t), [-.1,.1] ), axis image off, colormap(gray), pause(.1)
% end




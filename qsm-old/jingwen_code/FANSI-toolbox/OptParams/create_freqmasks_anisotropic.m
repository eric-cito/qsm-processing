function [ m1, m2, m3 ] = create_freqmasks_anisotropic( spatial_res,kernel )
% Calculates the masks that define Regions of Interest in the Fourier domain of
% reconstructed QSM images. These masks are defined by boundares that depend on:
% 1) Dipole kernel coefficients.
% 2) Absolute frequency radial range.
% Default values are those suggested by:
%
% Milovic C, et al. Comparison of Parameter Optimization Methods for Quantitative 
% Susceptibility Mapping. Magn Reson Med. 2020. DOI: 10.1002/MRM.28435 
%
% Parameters:
% spatial_res: voxel size in mm.
% kernel:  dipole kernel matrix, as calculated by the dipole_kernel function.
%
% Output:
% m1, m2, m3: binary masks defining three regions in the Fourier domain.
%             by default m1 covers the magic cone, whereas m2 and m3 do not.
%             For regularization weight optimization purposes, m2 and m3 are recommended.
%
% Use this function if voxels are highly asymmetric, with large thickness in the z axis.
% Please fine-tune the ranges in other scenarios.
%
% Last Modified by Carlos Milovic, 08.07.2020

N = size(kernel);

center = 1+floor(N/2);

rin = max(spatial_res)*0.65/2; % Radial boundaries of the masks are defined as absolute 
rout = max(spatial_res)*0.95/2;% frequency values [1/mm]. Please modify if needed.

kx = 1:N(1);
ky = 1:N(2);
kz = 1:N(3);

kx = single(kx) - center(1);
ky = single(ky) - center(2);
kz = single(kz) - center(3);

delta_kx = spatial_res(1)/N(1);
delta_ky = spatial_res(2)/N(2);
delta_kz = spatial_res(3)/N(3);

kx = kx * delta_kx;
ky = ky * delta_ky;
kz = kz * delta_kz;

kx = reshape(kx,[length(kx),1,1]);
ky = reshape(ky,[1,length(ky),1]);
kz = reshape(kz,[1,1,length(kz)]);

kx = repmat(kx,[1,N(2),N(3)]);
ky = repmat(ky,[N(1),1,N(3)]);
kz = repmat(kz,[N(1),N(2),1]);

k2 = kx.^2 + ky.^2 + kz.^2;

m0 = ones(N);
m0(k2>rout^2) = 0.0;
m0(k2<rin^2) = 0.0;
m0 = fftshift(m0); % Internal mask defining the radial range.

m1 = (m0).*(1-(abs(kernel)>0.080)); % Fine tune limits if necesary.
m2 = (m0).*(((kernel)>0.1) - ((kernel)>0.2) );
m3 = (m0).*(((kernel)>0.22) - ((kernel)>0.32) );

end


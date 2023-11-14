function [ kernel ] = dipole_kernel_angulated( N, spatial_res, B0_dir )
% This function calculates the dipole kernel used in QSM (single orientation)
% that models the susceptibility-to-field convolution.
% Use this function for angulated acquisitions, i.e. the main field direction is not along the z axis.
%
% This function uses an extension to the continuous kernel proposed by Salomir, et al. 2003.
%
% Parameters:
% N: array size
% spatial_res: voxel size in mm.
% B0_dir: main field direction, e.g. [0 0 1]
%
% Output:
% kernel: dipole kernel in the frequency space
%
% Created by Carlos Milovic, 30.03.2017
% Modified by Julio Acosta-Cabronero, 26.05.2017
% Last Modified by Carlos Milovic, 06.07.2020

[ky,kx,kz] = (meshgrid(-floor(N(2)/2):ceil(N(2)/2)-1, -floor(N(1)/2):ceil(N(1)/2)-1, -floor(N(3)/2):ceil(N(3)/2)-1));

kx = (single(kx) / max(abs(single(kx(:))))) / spatial_res(1);
ky = (single(ky) / max(abs(single(ky(:))))) / spatial_res(2);
kz = (single(kz) / max(abs(single(kz(:))))) / spatial_res(3);


k2 = kx.^2 + ky.^2 + kz.^2;
k2(k2==0) = eps;

%R_tot = eye(3); % Original formulation with a rotation matrix
%kernel = fftshift( 1/3 - (kx * R_tot(3,1) + ky * R_tot(3,2) + kz * R_tot(3,3)).^2 ./ k2 ); 

% JAC
kernel = fftshift( 1/3 - (kx*B0_dir(1) + ky*B0_dir(2) + kz*B0_dir(3)).^2 ./ k2 );    
kernel(1,1,1) = 0.0;

end


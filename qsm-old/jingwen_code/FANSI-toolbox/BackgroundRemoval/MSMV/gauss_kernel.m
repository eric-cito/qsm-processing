function [ kernel ] = gauss_kernel( N, spatial_res, sigma )
% Creates a Gaussian convolution kernel, in the frequency domain
%
% Parameters:
% N: voxel dimensions.
% spatial_res: voxel size, in mm.
% sigma: standard deviation, in mm.
%
% Output:
% kernel: Gaussian blur operator, in the frequency domain.
%
%  Last Modified by Carlos Milovic, 12.07.2020


[ky,kx,kz] = meshgrid(-floor(N(2)/2):ceil(N(2)/2)-1, -floor(N(1)/2):ceil(N(1)/2)-1, -floor(N(3)/2):ceil(N(3)/2)-1);
kx = single(kx)*spatial_res(1);
ky = single(ky)*spatial_res(2);
kz = single(kz)*spatial_res(3);
k2 = kx.^2 + ky.^2 + kz.^2;


kernel = exp( -k2/(2*sigma*sigma) ); 

DC = sum( kernel(:) );

kernel = fftn(fftshift(kernel/DC));
    
end


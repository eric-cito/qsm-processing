function chi_tik = tik( phase_use, mask_use, kernel, beta, N )
% Direct QSM solver with a Tikhonov regularization of the solution.
% This promotes smooth solutions with low energy.
%
% Parameters:
% phase_use: local field map
% mask_use: binary 3D image that defines the ROI.
% kernel: dipole kernel in the frequency space
% beta: regularization parameter
% N: array size
%
% output:
% chi_L2 - susceptibility map.
%
% Last modified by Carlos Milovic in 2020.07.07
%

%tic
K2 = abs(kernel).^2;

    chi_L2 = real( ifftn(conj(kernel) .* fftn(phase_use) ./ (K2 + beta)) ) .* mask_use;
%toc


end

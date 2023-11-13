function chi_tkd = tkd( phase_use, mask_use, kernel, kthre, N )
% QSM through the Truncation of the Dipole Kernel method (Shmueli et al 2009)
%
% Parameters:
% phase_use: local field map
% mask_use: binary 3D image that defines the ROI.
% kernel: dipole kernel in the frequency space
% kthre: threshold in the frequency space to truncate the kernel
% N: array size
%
% output:
% chi_tkd: susceptibility map.
%
% Please cite:
%Shmueli K., de Zwart J. a, van Gelderen P., Li T.-Q., Dodd S. J., Duyn J. H. 
%Magnetic susceptibility mapping of brain tissue in vivo using MRI phase data. 
%Magn Reson Med. 2009;62, 1510–1522. DOI:10.1002/mrm.22135
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Created by Carlos Milovic in 2017.03.30
% Last modified by Carlos Milovic in 2020.07.07

tic
%kernel_inv = zeros(N); % This formulation is similar to Wharton et al 2010
%kernel_inv( abs(kernel) > kthre ) = 1 ./ kernel(abs(kernel) > kthre);

kernel_inv = ones(N)/kthre; 
kernel_inv( abs(kernel) > kthre ) = 1 ./ kernel(abs(kernel) > kthre);

chi_tkd = real( ifftn( kernel_inv.* fftn(phase_use) ) ) .* mask_use; 


toc


end

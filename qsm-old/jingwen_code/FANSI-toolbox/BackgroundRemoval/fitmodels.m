function [ background, kappa ] = fitmodels( img, mask, weight, phi_e )
%   This function fits a background field model to acquired data in a least squares sense.
%   Use this to fit arbitrary background models (i.e. analytic models) to the actual background.
%   The fitting model also includes linear gradients in each direction.
%
%  Parameters:
%  img: acquired or total field data.
%  mask: binary mask defining the region of interest.
%  weight: reliability map. This might be an estimation of the SNR or noise whitening matrix.
%  phi_e: background model to be fitted to the data.
%
%  Output:
%  background: Linearly rescaled background model.
%  kappa: vector containing the scaling factors
%         kappa(1): offset/pedestal.
%         kappa(2): linear scaling X direction.
%         kappa(3): linear scaling Y direction.
%         kappa(4): linear scaling Z direction.
%         kappa(5): linear scaling of phi_e.
%
%  Last Modified by Carlos Milovic, 07.07.2020


N = size(img);
%se = strel('sphere',3);
%mask3=imerode(mask,se);
mask3 = mask;

% Deprecated code. Simple linear fitting.
%    kappa = 0;
%     phi0 = 0;
% 
%     a1 = sum( weight(:).*mask3(:).*phi_e(:) );
%     a2 = sum( weight(:).*mask3(:).*phi_e(:).*phi_e(:) );
%     a3 = sum( weight(:).*mask3(:) );
%     f1 = sum( weight(:).*mask3(:).*phi_e(:).*z2(:) );
%     f2 = sum( weight(:).*mask3(:).*z2(:) );
%     det = a1*a1-a2*a3;
%     
%     phi0 = (a1*f1 - a2*f2)/det;
%     kappa = (-a3*f1 + a1*f2)/det;
 
% Create linear gradients
[ky,kx,kz] = meshgrid(-floor(N(2)/2):ceil(N(2)/2)-1, -floor(N(1)/2):ceil(N(1)/2)-1, -floor(N(3)/2):ceil(N(3)/2)-1);
kx = (single(kx) / max(abs(single(kx(:))))) / spatial_res(1);
ky = (single(ky) / max(abs(single(ky(:))))) / spatial_res(2);
kz = (single(kz) / max(abs(single(kz(:))))) / spatial_res(3);
    
% Calculate the coefficients of the linear system    
A(1,1) = sum( weight(:).*mask3(:) );
A(2,2) = sum( weight(:).*mask3(:).*kx(:).*kx(:) );
A(3,3) = sum( weight(:).*mask3(:).*ky(:).*ky(:) );
A(4,4) = sum( weight(:).*mask3(:).*kz(:).*kz(:) );
A(5,5) = sum( weight(:).*mask3(:).*phi_e(:).*phi_e(:) );

A(1,2) = sum( weight(:).*mask3(:).*kx(:) );
A(2,1) = A(1,2);
A(1,3) = sum( weight(:).*mask3(:).*ky(:) );
A(3,1) = A(1,3);
A(1,4) = sum( weight(:).*mask3(:).*kz(:) );
A(4,1) = A(1,4);
A(1,5) = sum( weight(:).*mask3(:).*phi_e(:) );
A(5,1) = A(1,5);

A(2,3) = sum( weight(:).*mask3(:).*ky(:).*kx(:) );
A(3,2) = A(2,3);
A(2,4) = sum( weight(:).*mask3(:).*kz(:).*kx(:) );
A(4,2) = A(2,4);
A(2,5) = sum( weight(:).*mask3(:).*phi_e(:).*kx(:) );
A(5,2) = A(2,5);

A(3,4) = sum( weight(:).*mask3(:).*kz(:).*ky(:) );
A(4,3) = A(3,4);
A(3,5) = sum( weight(:).*mask3(:).*phi_e(:).*ky(:) );
A(5,3) = A(3,5);

A(4,5) = sum( weight(:).*mask3(:).*phi_e(:).*kz(:) );
A(5,4) = A(4,5);

F(1) = sum( weight(:).*mask3(:).*img(:) );
F(2) = sum( weight(:).*mask3(:).*img(:).*kx(:) );
F(3) = sum( weight(:).*mask3(:).*img(:).*ky(:) );
F(4) = sum( weight(:).*mask3(:).*img(:).*kz(:) );
F(5) = sum( weight(:).*mask3(:).*img(:).*phi_e(:) );

% Invert linear system to find scaling factors
kappa = A\F';

% Create fitted model
background = kappa(5)*phi_e+kappa(4)*kz+kappa(3)*ky+kappa(2)*ky+kappa(1);

end


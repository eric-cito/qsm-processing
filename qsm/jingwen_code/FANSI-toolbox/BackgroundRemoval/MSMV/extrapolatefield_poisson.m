function [ phi_out ] = extrapolatefield_poisson( phi, mask )
% This function extrapolates the phase values outside of a user-defined
% region of interest. Typically, this is done to "inpaint" regions without 
% MRI signal, such as air cavities or external fields.
%
% This is done by directly inverting a Poisson equation where Lap(phi_out) = modified-Lap(phi).
% The field outside the mask is imposed to be harmonic by modifiying the laplacian
% of the input phase and setting it to zero. Von Neuman boundary conditions are
% imposed at the boundary of the ROI to better fit external fields.
%
% Parameters:
% phi: total field or phase map.
% mask: ROI mask that corresponds to the tissue.
%
% Output:
% phi_out: extrapolated field map
%
% Last Modified by Carlos Milovic, 12.07.2020
% TO DO: Replace with operators in the frequency domain


% Calculate the mean phase value inside the ROI mask
meanIn = sum(phi(:).*mask(:))/sum(mask(:));

% Backward differentiation
dpx = phi-phi([1,1:(end-1)],:,:);
dpy = phi-phi(:,[1,1:(end-1)],:);
dpz = phi-phi(:,:,[1,1:(end-1)]);

dmx = mask-mask([1,1:(end-1)],:,:);
dmy = mask-mask(:,[1,1:(end-1)],:);
dmz = mask-mask(:,:,[1,1:(end-1)]);

% Forward differentiation
tx = phi([2:end,end],:,:)-phi;
ty = phi(:,[2:end,end],:)-phi;
tz = phi(:,:,[2:end,end])-phi;

% Replace the backward differentiation with the forward differentiation for all 
% voxels in the boundary where external voxels were included.
% This imposes Von Neuman boundary conditions to the ROI.
dpx( dmx > 0 ) = tx( dmx > 0 );
dpy( dmy > 0 ) = ty( dmy > 0 );
dpz( dmz > 0 ) = tz( dmz > 0 );

% Use forward differentiation to obtain the second derivative
d2px = dpx([2:end,end],:,:)-dpx;
d2py = dpy([2:end,end],:,:)-dpy;
d2pz = dpz([2:end,end],:,:)-dpz;

% Recalculate the gradient of the mask with forward operator.
dmx = mask([2:end,end],:,:)-mask;
dmy = mask(:,[2:end,end],:)-mask;
dmz = mask(:,:,[2:end,end])-mask;

% Impose Von Neuman boundary conditions, as before.
d2px( dmx < 0 ) = 0.0;
d2py( dmy < 0 ) = 0.0;
d2pz( dmz < 0 ) = 0.0;

% Calculate the Laplacian of the phase map. Set to zero values outside of the mask.
% This means external fields are harmonic.
Lp = mask.*(d2px+d2py+d2pz);
FLp = fftn(Lp);
FLp(1,1,1) = 0.0; % Ensure that it is zero meaned.

% Create Laplacian operator in the space domain and then transform it.
N = size(phi);
Lap = zeros(N);
Lap(1,1,1) = -6.0;
Lap(2,1,1) = 1.0;
Lap(end,1,1) = 1.0;
Lap(1,2,1) = 1.0;
Lap(1,end,1) = 1.0;
Lap(1,1,2) = 1.0;
Lap(1,1,end) = 1.0;
L = fftn(Lap);

% Invert the Poisson equation
phi_out = ifftn(FLp./(fftn(Lap)+eps));

% Demean inside the ROI mask and fit input data
meanOut = sum(phi_out(:).*mask(:))/sum(mask(:));
phi_out = phi_out-meanOut+meanIn;



end

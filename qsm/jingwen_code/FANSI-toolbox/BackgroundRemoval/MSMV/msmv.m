function [ local ] = msmv( p0, mask_use, spatial_res, layers, showresults )
% Perform the background field removal step using a Multiscale Spherical
% Mean Value algorithm.
% This algorithms uses Gaussian kernels to build a Laplacian pyramid. Since each
% layer contains onli high-passed information by a spherical kernel, no background
% fields are expected. Discarding the last (residual) layer yields only local
% fields. 
% A deconvolution step may be needed to recover large scale features.
%
% External fields are extrapolated as harmonic functions to avoid corrupting the
% data at each layer.
%
% Parameters:
% p0: total field or phase map.
% mask_use: ROI mask that corresponds to the local tissue.
% spatial_res: voxel size, in mm.
% layers: number of layers in a dyadic sequence of scales (recommended = 6, max = 8)
%
% Output:
% local: local field map
%
% Please cite:
% Milovic C, Langkammer C, Uribe S, Irarrazaval P, Acosta-Cabronero J, and Tejos C. 
% Multiscale Spherical Mean Value based background field removal method for 
% Quantitative Susceptibility Mapping. 27th International Conference of the ISMRM, 
% Montreal, Canada, 2019.
%
% Created by Carlos Milovic, 31.10.2018
% Last Modified by Carlos Milovic, 14.07.2020


if nargin < 4
    layers = 6;
end


N = size(p0);

if nargin < 5
showresults = false; % Set to true if you want to see and evaluate each layer.
end

% Create the Gaussian 3D convolutional kernels
g0 = gauss_kernel( N, spatial_res, 1/2 );
g1 = gauss_kernel( N, spatial_res, 1 );
g2 = gauss_kernel( N, spatial_res, 2 );

% Create the spherical morphological element used to erode the mask at each layer.
se =  strel('sphere',1);
    
    
% First layer    
    mt = mask_use;
    pt = p0;
    gt = g0; % Use filter with smaller standard deviation.
    mte = imerode(mt,se);
% Extrapolate external fields to avoid corrupting the internal data with noise.    
phi_out  = extrapolatefield_iter( pt, mt,1e2,100 ); % Iterate to refine the inner solution
phi_out(mt > 0) = pt(mt>0);
phi_out  = extrapolatefield_iter( phi_out, mt,1e1,100 );
phi_out(mt > 0) = pt(mt>0);
phi_out  = extrapolatefield_iter( phi_out, mt,1e0,100 );
phi_out(mt > 0) = pt(mt>0);
phi_out  = extrapolatefield_iter( phi_out, mt,1e-1,100 );
phi_out(mte > 0) = pt(mte>0);

b0 = ifftn( gt.*fftn(phi_out) ); % blurred layer
l0 = (phi_out-b0).*(mt+mte)/2; % Laplacian layer
b0 = phi_out-l0; % update, leaving external fields unchanged

if showresults == true
    imagesc3d2(l0, N/2, 10, [90,90,-90], [-1,1], [], 'l0')
end

if layers > 1
% Second layer, initialize variables with previous layer
    mt = mte; % Use previously eroded mask
    pt = b0; % Input data is the blurred first layer.
    gt = g1; % Duplicate the standard deviation of the filter.
    mte = imerode(mt,se);

% Refine again the extrapolation of the fields, with the new mask    
phi_out  = extrapolatefield_iter( pt, mt,2e-1,100 );
phi_out(mte > 0) = pt(mte>0);
b1 = ifftn( gt.*fftn(phi_out) );
    mte2 = imerode(mte,se);
l1 = (phi_out-b1).* mte.*ifftn( gt.*fftn(mte2)); % Second Laplacian layer
b1 = phi_out-l1;
    
if showresults == true
    imagesc3d2(l1, N/2, 11, [90,90,-90], [-1,1], [], 'l1')
end
    
if layers > 2
% Third layer, initialize variables with previous layer
    mt = mte; % Use previously eroded mask
    pt = b1; % Input data is the blurred first layer.
    gt = g2; % Duplicate the standard deviation of the filter.
    mte = imerode(mt,se);

% Refine again the extrapolation of the fields, with the new mask    
phi_out  = extrapolatefield_iter( pt, mt,2e-1,100 );
phi_out(mte > 0) = pt(mte>0);
b2 = ifftn( gt.*fftn(phi_out) );
    mte2 = imerode(mte,se);
l2 = (phi_out-b1).* mte.*ifftn( gt.*fftn(mte2)); % Second Laplacian layer
b2 = phi_out-l1;
    
if showresults == true
    imagesc3d2(l2, N/2, 12, [90,90,-90], [-1,1], [], 'l2')
end

if layers > 3
% Fourth layer    
% Here, instead of changing the filter size, we downsample the data to make the convolution more efficient.
% This allows padding this layer to avoid large scale artifacts.
    mt = zeros(N);
    mt( (round(N(1)/4)+1):(round(N(1)/2)+round(N(1)/4)),(round(N(2)/4)+1):round((N(2)/2)+round(N(2)/4)),(round(N(3)/4)+1):(round(N(3)/2)+round(N(3)/4)) ) = imresizen( single(mte), 0.5 ); %
    pt = zeros(N);
    pt( (round(N(1)/4)+1):(round(N(1)/2)+round(N(1)/4)),(round(N(2)/4)+1):round((N(2)/2)+round(N(2)/4)),(round(N(3)/4)+1):(round(N(3)/2)+round(N(3)/4)) ) = imresizen( b2, 0.5 ); 
    mte = imerode(mt,se);
%gt = gauss_kernel( size(mt), spatial_res, 2 );
    
phi_out  = extrapolatefield_iter( pt, mt,2e-1,250 );
phi_out(mte > 0) = pt(mte>0);
b3 = ifftn( gt.*fftn(phi_out) );
    mte2 = imerode(mte,se);
l3 = (phi_out-b3).* mte.*ifftn( gt.*fftn(mte2));
b3 = phi_out-l3;
    
if showresults == true
    imagesc3d2(l3, N/2, 13, [90,90,-90], [-1,1], [], 'l3')
end
    
if layers > 4    
% Fifth layer
% Hereforth, data is subsampled without padding to increase efficiency.
    mt = imresizen( single(mte), 0.5 );
    pt = imresizen( b3, 0.5 );
    mte = imerode(mt,se);
    gt = gauss_kernel( N/2, spatial_res, 2 );
    
phi_out  = extrapolatefield_iter( pt, mt,1e2,100 );
phi_out(mt > 0) = pt(mt>0);
phi_out  = extrapolatefield_iter( phi_out, mt,1e1,100 );
phi_out(mt > 0) = pt(mt>0);
phi_out  = extrapolatefield_iter( phi_out, mt,1e0,100 );
phi_out(mt > 0) = pt(mt>0);
phi_out  = extrapolatefield_iter( phi_out, mt,1e-1,100 );
phi_out(mte > 0) = pt(mte>0);

b4 = ifftn( gt.*fftn(phi_out) );
    mte2 = imerode(mte,se);
l4 = (phi_out-b4).* mte.*ifftn( gt.*fftn(mte2));
b4 = phi_out-l4;
    
if showresults == true
    imagesc3d2(l4, N/4, 14, [90,90,-90], [-1,1], [], 'l4')
end
    
    
if layers > 5    
% Sixth layer
    mt = imresizen( single(mte), 0.5 );
    pt = imresizen( b4, 0.5 );
    mte = imerode(mt,se);
    gt = gauss_kernel( N/4, spatial_res, 2 );
    
phi_out  = extrapolatefield_iter( pt, mt,2e-1,100 );
phi_out(mte > 0) = pt(mte>0);
b5 = ifftn( gt.*fftn(phi_out) );
    mte2 = imerode(mte,se);
l5 = (phi_out-b5).* mte.*ifftn( gt.*fftn(mte2)).^1;
b5 = phi_out-l5;
    
if showresults == true
    imagesc3d2(l5, N/8, 15, [90,90,-90], [-1,1], [], 'l5')
end
    
if layers > 6    
% Seventh layer
    mt = imresizen( single(mte), 0.5 );
    pt = imresizen( b5, 0.5 );
    mte = imerode(mt,se);
    gt = gauss_kernel( N/8, spatial_res, 2 );
    
phi_out  = extrapolatefield_iter( pt, mt,2e-1,100 );
phi_out(mte > 0) = pt(mte>0);
b6 = ifftn( gt.*fftn(phi_out) );
    mte2 = imerode(mte,se);
l6 = (phi_out-b6).* mte.*ifftn( gt.*fftn(mte2)).^1;
b6 = phi_out-l6;
    
if showresults == true
    imagesc3d2(l6, N/16, 16, [90,90,-90], [-1,1], [], 'l6')
end
    % This layer represents structures around 32 voxels in size. Typically stop here
    % to avoid corruption due to regions without signal, and limit QSM underestimation of 
    % large structures. If needed, consider using a convolution step with sigma = 16 or similar.
    
    
if layers > 7    
% Eigth layer
    mt = imresizen( single(mte), 0.5 );
    pt = imresizen( b6, 0.5 );
    mte = imerode(mt,se);
    gt = gauss_kernel( N/16, spatial_res, 2 );
    
phi_out  = extrapolatefield_iter( pt, mt,2e-1,100 );
phi_out(mte > 0) = pt(mte>0);
b7 = ifftn( gt.*fftn(phi_out) );
    mte2 = imerode(mte,se);
l7 = (phi_out-b7).* mte.*ifftn( gt.*fftn(mte2)).^1;
b7 = phi_out-l7;
    
if showresults == true
    imagesc3d2(l7, N/32, 17, [90,90,-90], [-1,1], [], 'l7')
    imagesc3d2(b7, N/32, 18, [90,90,-90], [-1,1], [], 'b7')
end

else
l7 = zeros( round(size(l6)/2) );  
if showresults == true
    imagesc3d2(b6, N/32, 17, [90,90,-90], [-1,1], [], 'b6')
end  
end % 8th layer

else
l6 = zeros( round(size(l5)/2) ); 
l7 = zeros( round(size(l6)/2) );  
if showresults == true
    imagesc3d2(b5, N/16, 16, [90,90,-90], [-1,1], [], 'b5')
end  
end % 7th layer

else
l5 = zeros( round(size(l4)/2 )); 
l6 = zeros( round(size(l5)/2 )); 
l7 = zeros( round(size(l6)/2 ));  
if showresults == true
    imagesc3d2(b4, N/8, 15, [90,90,-90], [-1,1], [], 'b4')
end  
end % 6th layer

else
l4 = zeros( round(size(l3)/2 )); 
l5 = zeros( round(size(l4)/2 )); 
l6 = zeros( round(size(l5)/2 )); 
l7 = zeros( round(size(l6)/2 ));  
if showresults == true
    imagesc3d2(b3, N/4, 14, [90,90,-90], [-1,1], [], 'b3')
end  
end % 5th layer

else
l3 = zeros( round(size(l2)/2 )); 
l4 = zeros( round(size(l3)/2 )); 
l5 = zeros( round(size(l4)/2 )); 
l6 = zeros( round(size(l5)/2 )); 
l7 = zeros( round(size(l6)/2 ));  
if showresults == true
    imagesc3d2(b2, N/2, 13, [90,90,-90], [-1,1], [], 'b2')
end  
end % 4th layer

else
l2 = zeros( round(size(l1)) ); 
l3 = zeros( round(size(l2)/2 )); 
l4 = zeros( round(size(l3)/2) ); 
l5 = zeros( round(size(l4)/2 )); 
l6 = zeros( round(size(l5)/2) ); 
l7 = zeros( round(size(l6)/2) );  
if showresults == true
    imagesc3d2(b1, N/2, 12, [90,90,-90], [-1,1], [], 'b1')
end  
end % 3rd layer

else
l1 = zeros( round(size(l0) )); 
l2 = zeros( round(size(l1) )); 
l3 = zeros( round(size(l2)/2 )); 
l4 = zeros( round(size(l3)/2 )); 
l5 = zeros( round(size(l4)/2 )); 
l6 = zeros( round(size(l5)/2 )); 
l7 = zeros( round(size(l6)/2) );  
if showresults == true
    imagesc3d2(b0, N/2, 11, [90,90,-90], [-1,1], [], 'b0')
end  
end % 2nd layer

    % Upscale and aggregate detail layers
    
    cl6 = imresizen( l7, size(l6)./size(l7),'spline' ) + l6;
    cl5 = imresizen( cl6, size(l5)./size(l6),'spline' ) + l5;
    cl4 = imresizen( cl5, size(l4)./size(l5),'spline' ) + l4;
    cl3 = imresizen( cl4, size(l3)./size(l4),'spline' ) + l3;
    cropl3 = cl3( (round(N(1)/4)+1):(N(1)/2+round(N(1)/4)),(round(N(2)/4)+1):(N(2)/2+round(N(2)/4)),(round(N(3)/4)+1):(N(3)/2+round(N(3)/4)) ); % Crop
    cl2 = imresizen( cropl3,  size(l2)./size(cropl3) ) + l2;
    cl1 = cl2 + l1;
    local = cl1 + l0;
    
if showresults == true
    imagesc3d2(local, N/2, 19, [90,90,-90], [-1,1], [], 'local')
end

end


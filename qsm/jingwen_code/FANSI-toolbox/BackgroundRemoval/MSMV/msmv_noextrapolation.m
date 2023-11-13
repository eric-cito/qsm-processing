function [ local ] = msmv_noextrapolation( p0, mask_use, spatial_res, layers, showresults )
% Perform the background field removal step using a Multiscale Spherical
% Mean Value algorithm.
% This algorithms uses Gaussian kernels to build a Laplacian pyramid. Since each
% layer contains onli high-passed information by a spherical kernel, no background
% fields are expected. Discarding the last (residual) layer yields only local
% fields. 
% A deconvolution step may be needed to recover large scale features.
%
% External fields are not extrapolated to speed up the calculation.
%
% Parameters:
% p0: total field or phase map.
% mask_use: ROI mask that corresponds to the local tissue.
% spatial_res: voxel size, in mm.
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

if nargin < 5
showresults = false; % Set to true if you want to see and evaluate each layer.
end

N = size(p0);

% Create the Gaussian 3D convolutional kernels
    g0 = gauss_kernel( N, spatial_res, 1/2 );
    g1 = gauss_kernel( N, spatial_res, 1 );
    g2 = gauss_kernel( N, spatial_res, 2 );

% Create the spherical morphological element used to erode the mask at each layer.
    se =  strel('sphere',1);
    
    
    
% First layer    
    mt = mask_use;
    pt = p0;
    gt = g0;% Use filter with smaller standard deviation.
    mte = imerode(mt,se);
    
b0 = ifftn( gt.*fftn(pt) );
l0 = (pt-b0).*(mt);
%b0 = pt-l0;
if showresults == true
    imagesc3d2(l0, N/2, 10, [90,90,-90], [-1,1], [], 'l0')
end

if layers > 1
% Second layer, initialize variables with previous layer
    mt = mte;
    pt = b0;
    gt = g1;
    mte = imerode(mt,se);
b1 = ifftn( gt.*fftn(pt) );
    mte2 = imerode(mte,se);
l1 = (pt-b1).* mte.*ifftn( g0.*g1.*fftn(mte2));
%b1 = pt-l1;
    
if showresults == true
    imagesc3d2(l1, N/2, 11, [90,90,-90], [-1,1], [], 'l1')
end
    
    
if layers > 2
% Third layer    
    mt = mte;
    pt = b1;
    gt = g2;
    mte = imerode(mt,se);
b2 = ifftn( gt.*fftn(pt) );
    mte2 = imerode(mte,se);
    mte2 = imerode(mte2,se);
l2 = (pt-b2).* mte.*ifftn( g0.*g1.*g2.*fftn(mte2));
%b2 = pt-l2;
    
if showresults == true
    imagesc3d2(l2, N/2, 12, [90,90,-90], [-1,1], [], 'l2')
end
    
    
if layers > 3
% Fourth layer
    %sz = size(mte);
    %dsz2 = mod(sz,2);    
    pt = b2;%padarray( b2,dsz2,'replicate','post');  % Padding to even size
    
    mt = imresizen( single(mte), 0.5 );
    pt = imresizen( pt, 0.5,'spline' );
    mte = imerode(mt,se);
    mte = imerode(mte,se);
    gt = gauss_kernel( size(pt), spatial_res, 2 );

b3 = ifftn( gt.*fftn(pt) );
    mte2 = imerode(mte,se);
    mte2 = imerode(mte2,se);
l3 = (pt-b3).* mte.*ifftn( gt.*fftn(mte2));
%b3 = pt-l3;
    
if showresults == true
    imagesc3d2(l3, N/4, 13, [90,90,-90], [-1,1], [], 'l3')
end
    
    
    
    
if layers > 4
% Fifth layer
    %sz = size(mte);
    %dsz3 = mod(sz,2);    
    pt = b3;%padarray( b3,dsz3,'replicate','post'); % Padding to even size
    
    mt = imresizen( single(mte), 0.5 );
    pt = imresizen( pt, 0.5,'spline' );
    mte = imerode(mt,se);
    mte = imerode(mte,se);
    gt = gauss_kernel( size(pt), spatial_res, 2 );
b4 = ifftn( gt.*fftn(pt) );
    mte2 = imerode(mte,se);
    mte2 = imerode(mte2,se);
l4 = (pt-b4).* mte.*ifftn( gt.*fftn(mte2));
%b4 = pt-l4;
    
if showresults == true
    imagesc3d2(l4, N/8, 14, [90,90,-90], [-1,1], [], 'l4')
end
    
    
    
    
if layers > 5
% Sixth layer
    %sz = size(mte);
    %dsz4 = mod(sz,2);    
    pt = b4;%padarray( b4,dsz4,'replicate','post'); % Padding to even size
    
    mt = imresizen( single(mte), 0.5 );
    pt = imresizen( pt, 0.5,'spline' );
    mte = imerode(mt,se);
    mte = imerode(mte,se);
    gt = gauss_kernel( size(pt), spatial_res, 2 );
b5 = ifftn( gt.*fftn(pt) );
    mte2 = imerode(mte,se);
    mte2 = imerode(mte2,se);
l5 = (pt-b5).* mte.*ifftn( gt.*fftn(mte2)).^1;
%b5 = pt-l5;

if showresults == true
    imagesc3d2(l5, N/16, 15, [90,90,-90], [-1,1], [], 'l5')
end
    
    
    
if layers > 6   
% Seventh layer
    %sz = size(mte);
    %dsz5 = mod(sz,2);    
    pt = b5;%padarray( b5,dsz5,'replicate','post'); % Padding to even size
    
    mt = imresizen( single(mte), 0.5 );
    pt = imresizen( pt, 0.5,'spline' );
    mte = imerode(mt,se);
    mte = imerode(mte,se);
    gt = gauss_kernel( size(pt), spatial_res, 2 );
b6 = ifftn( gt.*fftn(pt) );
    mte2 = imerode(mte,se);
    mte2 = imerode(mte2,se);
l6 = (pt-b6).* mte.*ifftn( gt.*fftn(mte2)).^1;
%b6 = pt-l6;
    
if showresults == true
    imagesc3d2(l6, N/32, 16, [90,90,-90], [-1,1], [], 'l6')
end
    % This layer represents structures around 32 voxels in size. Typically stop here
    % to avoid corruption due to regions without signal, and limit QSM underestimation of 
    % large structures. If needed, consider using a convolution step with sigma = 16 or similar.
    
    
    
if layers > 7    
% Eight layer
    %sz = size(mte);
    %dsz5 = mod(sz,2);    
    pt = b6;%padarray( b5,dsz5,'replicate','post'); % Padding to even size
    
    mt = imresizen( single(mte), 0.5 );
    pt = imresizen( pt, 0.5,'spline' );
    mte = imerode(mt,se);
    mte = imerode(mte,se);
    gt = gauss_kernel( size(pt), spatial_res, 2 );
b7 = ifftn( gt.*fftn(pt) );
    mte2 = imerode(mte,se);
    mte2 = imerode(mte2,se);
l7 = (pt-b7).* mte.*ifftn( gt.*fftn(mte2)).^1;
%b6 = pt-l6;
    
if showresults == true
    imagesc3d2(l7, N/64, 17, [90,90,-90], [-1,1], [], 'l7')
end
if showresults == true
    imagesc3d2(b7, N/64, 18, [90,90,-90], [-1,1], [], 'b7')
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
    cl2 = imresizen( cl3, size(l2)./size(l3),'spline' ) + l2;
    cl1 = cl2 + l1;
    local = (cl1 + l0).*mask_use;
    
if showresults == true
    imagesc3d2(local, N/2, 19, [90,90,-90], [-1,1], [], 'local')
end

end


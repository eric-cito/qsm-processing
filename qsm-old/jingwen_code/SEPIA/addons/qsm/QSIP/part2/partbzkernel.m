function [F] = partbzkernel(n,dx2,dy2,dz2,scales)
%
%
% AUTHOR: Clare Poynton
%
% PURPOSE:  Compute F(x) in eq. 16 of MJ's paper (see citation). 
%           The zeroth order term, B0(0) = (0,0,1), for a normalized (constant) field
%           along the z-direction. 
%
% INPUTS: 
%
% n = vector of image dimensions. This should be the same dimensions as chi map to be estimated later.
%
% scales = vector containing voxel sizes. ie. [0.9375, 0.9375, 2.5]
%
% dx2, dy2, dz2 = voxel sizes / 2.  
%
% CITATION:
%
% Jenkinson M, Wilson JL, Jezzard P. Perturbation method for magnetic field
% calculations of nonconductive objects. Magn Reson Med. 2004 Sep; 52(3):
% 471-7. 
% 
%%


nx=n(1);
ny=n(2);
nz=n(3);

dx = scales(1);  
dy = scales(2);  
dz = scales(3);

xoff = dx2;
yoff = dy2;
zoff = dz2;

z1 = [-(nz-1):1:(nz-1)];
y1 = [-(ny-1):1:(ny-1)];
x1 = [-(nx-1):1:(nx-1)];

kernel = zeros(2*nx-1,2*ny-1,2*nz-1);

for z1_ind = 1:size(z1,2) 
    for y1_ind = 1:size(y1,2) 
        for x1_ind = 1:size(x1,2) 
          
           x = x1(x1_ind)*dx+xoff; 
           y = y1(y1_ind)*dy+yoff; 
           z = z1(z1_ind)*dz+zoff;
           
           r = sqrt(x*x + y*y + z*z);
           
           xy_zr = (x*y)/(z*r);
           
           kernel(x1(x1_ind)+nx,y1(y1_ind)+ny,z1(z1_ind)+nz) = atan(xy_zr);
           
        end
    end
end

F = kernel * (1/(4*pi));



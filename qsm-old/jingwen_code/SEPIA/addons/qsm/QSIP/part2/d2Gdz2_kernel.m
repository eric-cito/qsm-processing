function [d2G_dz2_kernel] = d2Gdz2_kernel(n,scales)
%
% AUTHOR: Clare Poynton
%
% PURPOSE:  Compute d2G_dz2_kernel = H(x) in eq. 15 of MJ's paper (see citation). 
%           The zeroth order term, B0(0) = (0,0,1), for a normalized (constant) field
%           along the z-direction. 
%
% INPUTS: 
%
% n = vector of image dimensions. This should be the same dimensions as chi map to be estimated later.
%
% scales = vector containing voxel sizes. ie. [0.9375, 0.9375, 2.5]
%
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

dx2 = dx/2.0;  
dy2 = dy/2.0;  
dz2 = dz/2.0;

kernel = partbzkernel(n,dx2,dy2,dz2,scales);
kernel = kernel - partbzkernel(n,-dx2,dy2,dz2,scales);
kernel = kernel - partbzkernel(n,dx2,-dy2,dz2,scales);
kernel = kernel - partbzkernel(n,dx2,dy2,-dz2,scales);
kernel = kernel + partbzkernel(n,dx2,-dy2,-dz2,scales);
kernel = kernel + partbzkernel(n,-dx2,dy2,-dz2,scales);
kernel = kernel + partbzkernel(n,-dx2,-dy2,dz2,scales);
d2G_dz2_kernel = kernel - partbzkernel(n,-dx2,-dy2,-dz2,scales);





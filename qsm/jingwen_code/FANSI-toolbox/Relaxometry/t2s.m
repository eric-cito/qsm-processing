function [ M0, R2s, T2s ] = t2s( MV, TE, weight )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tic
W = weight.*weight;
lMV = W.*log(MV+eps);

N = size(MV(:,:,:,1));
f1 = zeros(N,'double');
f2 = zeros(N,'double');

a11 = zeros(N,'double');
a22 = zeros(N,'double');
a21 = zeros(N,'double');
a12 = zeros(N,'double');

for c = 1:length(TE)
    a11 = a11 + TE(c)*TE(c)*W(:,:,:,c) ;
    a21 = a21 - TE(c)*W(:,:,:,c);    
    a22 = a22 - W(:,:,:,c);

    f1 = f1-TE(c)*lMV(:,:,:,c);
    f2 = f2-lMV(:,:,:,c);
end
a12 = -a21;
dA = eps+a11.*a22 - a21.*a12;

R2s = ((a22.*f1-a21.*f2)./(dA));
M0 = exp( (-a12.*f1+a11.*f2)./dA );
T2s = 1.0./(R2s+eps);

toc

end
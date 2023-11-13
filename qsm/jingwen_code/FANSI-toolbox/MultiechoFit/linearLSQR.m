function [ phi0, deltaB ] = linearLSQR( Phase, TE, weight )
%This function performs a linear least squared error fit to multiecho GRE phase data. 
%Please remember to use unwrapped phase data!
%
% Required fields:
% Phase: all echoes phase container (4th dimension echo index), in radians, and unwrapped
% TE: echo times, in ms.
% weight: weights data container (all echoes), rescaled to the [0,1] range. Magnitude or Magnitude*TE is suggested.
%
% Last modified by Carlos Milovic in 2020.10.20

tic
W = weight.*weight;
wp = W.*Phase;

N = size(Phase(:,:,:,1));
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

    f1 = f1-TE(c)*wp(:,:,:,c);
    f2 = f2-wp(:,:,:,c);
end
a12 = -a21;
dA = eps+a11.*a22 - a21.*a12;

deltaB = -((a22.*f1-a21.*f2)./(dA));
phi0 =  (-a12.*f1+a11.*f2)./dA ;

toc

end

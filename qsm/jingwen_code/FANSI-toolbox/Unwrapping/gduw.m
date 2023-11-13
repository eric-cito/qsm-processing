function [ kout uout ] = gduw( img, mask, vsize )
% Direct gradient domain unwrapping. 
% A is found by displacing the wrapping artifacts by offsetting the input phase.
% The gradient of several shifted phases are computed and the median value is calculated.
% A solution is found for the image that generates the modified gradient by solving this
% problem in a least squared sense. This is equivalent to solve a Poisson equation
% Here a rapid solver based on the FFT is used.
% Von Neuman boundary conditions are set at the mask boundary.
%
% Parameters:
% img: wrapped image
% mask: ROI mask
% vsize: voxel size.
%
% Output:
% kout: unwrapped image by 2pi integer additive factors
% uout: approximate solution to the Poisson equation
%
% Modified by Carlos Milovic in 2017.11.23
% Last modified by Carlos Milovic in 2020.07.14

N = size(img);

% Padd the images to avoid artifacts due to boundary conditions
padd = round(N(1)/8);
img2 = zeros(N+2*padd);
img2((1+padd):(padd+N(1)),(1+padd):(padd+N(2)),(1+padd):(padd+N(3))) = img;
mask2 = zeros(N+2*padd);
mask2((1+padd):(padd+N(1)),(1+padd):(padd+N(2)),(1+padd):(padd+N(3))) = mask;

N2 = size(img2);
GM = zeros( [N2(1) N2(2) N2(3) 5] );
GM(:,:,:,1) = (mask2([2:end,end],:,:) - mask2)/vsize(1);
GM(:,:,:,2) = (mask2(:,[2:end,end],:) - mask2)/vsize(2);
GM(:,:,:,3) = (mask2(:,:,[2:end,end]) - mask2)/vsize(3);
GM = 1-single(abs(GM) > 0.0);


% Calculate the gradient of the input data
GX(:,:,:,1) = GM(:,:,:,1).*mask2.*(img2([2:end,end],:,:) - img2)/vsize(1);
GY(:,:,:,1) = GM(:,:,:,2).*mask2.*(img2(:,[2:end,end],:) - img2)/vsize(2);
GZ(:,:,:,1) = GM(:,:,:,3).*mask2.*(img2(:,:,[2:end,end]) - img2)/vsize(3);


% Shift the phase data, and recalculate the gradients
temp = img2 + pi/3;
temp = angle(exp(1i*temp));
GX(:,:,:,2) = GM(:,:,:,1).*mask2.*(temp([2:end,end],:,:) - temp)/vsize(1);
GY(:,:,:,2) = GM(:,:,:,2).*mask2.*(temp(:,[2:end,end],:) - temp)/vsize(2);
GZ(:,:,:,2) = GM(:,:,:,3).*mask2.*(temp(:,:,[2:end,end]) - temp)/vsize(3);

% Shift the phase data again, and recalculate the gradients
temp = img2 + 2*pi/3;
temp = angle(exp(1i*temp));
GX(:,:,:,3) = GM(:,:,:,1).*mask2.*(temp([2:end,end],:,:) - temp)/vsize(1);
GY(:,:,:,3) = GM(:,:,:,2).*mask2.*(temp(:,[2:end,end],:) - temp)/vsize(2);
GZ(:,:,:,3) = GM(:,:,:,3).*mask2.*(temp(:,:,[2:end,end]) - temp)/vsize(3);


temp = img2 + 1*pi/4;
temp = angle(exp(1i*temp));
GX(:,:,:,4) = GM(:,:,:,1).*mask2.*(temp([2:end,end],:,:) - temp)/vsize(1);
GY(:,:,:,4) = GM(:,:,:,2).*mask2.*(temp(:,[2:end,end],:) - temp)/vsize(2);
GZ(:,:,:,4) = GM(:,:,:,3).*mask2.*(temp(:,:,[2:end,end]) - temp)/vsize(3);

temp = img2 + 3*pi/4;
temp = angle(exp(1i*temp));
GX(:,:,:,5) = GM(:,:,:,1).*mask2.*(temp([2:end,end],:,:) - temp)/vsize(1);
GY(:,:,:,5) = GM(:,:,:,2).*mask2.*(temp(:,[2:end,end],:) - temp)/vsize(2);
GZ(:,:,:,5) = GM(:,:,:,3).*mask2.*(temp(:,:,[2:end,end]) - temp)/vsize(3);


% Calculate the median value
v(:,:,:,1) = median(GX,4);
v(:,:,:,2) = median(GY,4);
v(:,:,:,3) = median(GZ,4);

% Don't allow large gradients.
for i=1:3
v(v>pi) = v(v>pi)-2*pi;
v(v<-pi) = v(v<-pi)+2*pi;
end


% Calculate the divergence of such gradient data.
u = zeros( [N2(1) N2(2) N2(3)] );
u = (v(:,:,:,3) - v(:,:,[1,1:(end-1)],3))/vsize(3) + (v(:,:,:,2) - v(:,[1,1:(end-1)],:,2))/vsize(2) + (v(:,:,:,1) - v([1,1:(end-1)],:,:,1))/vsize(1);
u = mask2.*u;


% Define the operators
FOV = N2.*vsize;
center = 1+N2/2;
kx = 1:N2(1);
ky = 1:N2(2);
kz = 1:N2(3);

kx = kx - center(1);
ky = ky - center(2);
kz = kz - center(3);

delta_kx = 1/FOV(1);
delta_ky = 1/FOV(2);
delta_kz = 1/FOV(3);


kx = kx * delta_kx;
ky = ky * delta_ky;
kz = kz * delta_kz;

kx = reshape(kx,[length(kx),1,1]);
ky = reshape(ky,[1,length(ky),1]);
kz = reshape(kz,[1,1,length(kz)]);

kx = repmat(kx,[1,N2(2),N2(3)]);
ky = repmat(ky,[N2(1),1,N2(3)]);
kz = repmat(kz,[N2(1),N2(2),1]);

k2 = -3+cos(2*pi*kx)+cos(2*pi*ky)+cos(2*pi*kz); % discretized Laplacian
k2 = 2*k2;
k2(k2==0) = eps;
kernel = 1.0 ./ k2;
DC = (kx==0) & (ky==0) & (kz==0);
kernel(DC==1) = 0;
kernel = fftshift(kernel);

% Solve Poisson equation
u = real(ifftn(kernel.*fftn(u)));

% Crop the solution
uout = u((1+padd):(padd+N(1)),(1+padd):(padd+N(2)),(1+padd):(padd+N(3)));
uout = uout+img(round(N(1)/2),round(N(2)/2),round(N(3)/2))-uout(round(N(1)/2),round(N(2)/2),round(N(3)/2));

% Find the best 2pi integers that unwraps the image to match the previous solution
kout = img;
for i = 1:150
    out_old = kout;
kout = kout + 2*pi*round( (uout - kout)/(2*pi) );

if sum(abs(out_old(:)-kout(:))) < 1
    break;
end

end
% Output the approximate and integer solutions.

end


function [ kernel ] = dipole_kernel_fansi( N, spatial_res, mode )
% This function calculates the dipole kernel used in QSM (single orientation)
% that models the susceptibility-to-field convolution.
% Use this function if the main field direction is along the z axis.
% In addition to the standard dipole kernel function, additional approximations are provided.
% See description below for details.
%
% Parameters:
% N: array size
% spatial_res: voxel size in mm.
% mode:  0 for the continuous kernel proposed by Salomir, et al. 2003.
%        1 for the discrete kernel proposed by Milovic, et al. 2017.
%        2 for the Integrated Green function proposed by Jenkinson, et al. 2004
%
% Output:
% kernel: dipole kernel in the frequency space
%
% Created by Carlos Milovic, 30.03.2017
% Last Modified by Carlos Milovic, 06.07.2020



if mode == 0 % Continuous kernel

[ky,kx,kz] = (meshgrid(-floor(N(2)/2):ceil(N(2)/2)-1, -floor(N(1)/2):ceil(N(1)/2)-1, -floor(N(3)/2):ceil(N(3)/2)-1));

kx = (single(kx) / max(abs(single(kx(:))))) / spatial_res(1);
ky = (single(ky) / max(abs(single(ky(:))))) / spatial_res(2);
kz = (single(kz) / max(abs(single(kz(:))))) / spatial_res(3);

k2 = kx.^2 + ky.^2 + kz.^2;
k2(k2==0) = eps;

kernel = fftshift( 1/3 - (kz.^2) ./ k2 ); % Keep the center of the frequency domain at [1,1,1]
%kernel13 = fftshift( kz.* kx  ./ k2 );      % STI anisotropic component
%kernel23 = fftshift( kz.* ky  ./ k2 );      % STI anisotropic component
kernel(1,1,1) = 0.0;

elseif mode == 1 % Discrete kernel

    
FOV = N.*spatial_res;
center = 1+floor(N/2);
kx = 1:N(1);
ky = 1:N(2);
kz = 1:N(3);

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

kx = repmat(kx,[1,N(2),N(3)]);
ky = repmat(ky,[N(1),1,N(3)]);
kz = repmat(kz,[N(1),N(2),1]);

k2 = -3+cos(2*pi*kx)+cos(2*pi*ky)+cos(2*pi*kz);
k2(k2==0) = eps;

kernel = 1/3 - (-1+cos(2*pi*kz)) ./ k2;

kernel = fftshift(kernel);  % Keep the center of the frequency domain at [1,1,1]
kernel(1,1,1) = 0.0;
    



elseif mode == 2 % Integrated Green function
    
FOV = N.*spatial_res;
center = 1+floor(N/2);
    
kx = 1:N(1);
ky = 1:N(2);
kz = 1:N(3);

kx = kx - center(1);
ky = ky - center(2);
kz = kz - center(3);


% determine the step sizes delta_kx, delta_ky, delta_kz in k-space
delta_kx = N(1)/FOV(1);
delta_ky = N(2)/FOV(2);
delta_kz = N(3)/FOV(3);


kx = kx * delta_kx;
ky = ky * delta_ky;
kz = kz * delta_kz;

kx = reshape(kx,[length(kx),1,1]);
ky = reshape(ky,[1,length(ky),1]);
kz = reshape(kz,[1,1,length(kz)]);

kx = repmat(kx,[1,N(2),N(3)]);
ky = repmat(ky,[N(1),1,N(3)]);
kz = repmat(kz,[N(1),N(2),1]);

spatial_kernel = zeros(N);
for i=1:N(1)
    for j=1:N(2)
        for k=1:N(3)
            x(1) = kx(i,j,k)-0.5*delta_kx;
            x(2) = kx(i,j,k)+0.5*delta_kx;
            y(1) = ky(i,j,k)-0.5*delta_ky;
            y(2) = ky(i,j,k)+0.5*delta_ky;
            z(1) = kz(i,j,k)-0.5*delta_kz;
            z(2) = kz(i,j,k)+0.5*delta_kz;
            
            m = zeros([2 2 2]);
            for di = 1:2
                for dj = 1:2
                    for dk = 1:2
                        m(di,dj,dk) = (-1)^(di+dj+dk) * atan(x(di)*y(dj)/(z(dk)*sqrt(x(di)^2+y(dj)^2+z(dk)^2)));
                    end
                end
            end
            spatial_kernel(i,j,k) = -(1/(pi*4)) *sum(m(:));       
            
            
        end
    end
end

DC = (kx==0) & (ky==0) & (kz==0);
spatial_kernel(DC==1) = spatial_kernel(DC==1)+1/3;

kernel = fftn(fftshift(spatial_kernel));
kernel(1,1,1) = 0.0;
    
end


end

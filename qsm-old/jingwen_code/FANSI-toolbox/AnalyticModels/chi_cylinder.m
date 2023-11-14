function [chi, Bnuc] = chi_cylinder(FOV,N,p1,p2,radius, xin, xout)
% Calculates the susceptibility distribution and magnetization by a cylinder.
% This generates output field maps that show the sampled values at each voxel location.
% This function includes the Lorentz sphere approximation for the sensing proton.
%
% Parameters:
% FOV: field of view in x, y, and z directions
% N: no of samples in kx, ky, kz
% p1,p2: coordinates of two points that define the axis of the cylinder (px,py,px)
%        using the center of the volumen as origin.
% radius: radius of the cylinder in voxels.
% xin: susceptibility value inside the cylinder (full value, not ppm!)
% xout: susceptibility value outside the cylinder (full value, not ppm!)
% *** Susceptibility values should be in the full range, i.e. 1ppm = 1e-6 ***
%
% Output:
% chi: sampled susceptibility distribution
% Bnuc: analytic nuclear magnetic field
%
% Created by Carlos Milovic in 2019.05.26
% Last modified by Carlos Milovic in 2020.07.07


dp = p2-p1;
%na = norm(dp);
[alpha, theta, na] = cart2sph(dp(1),dp(2),dp(3));
st2 = sin(pi/2 - theta)^2;
ct2 = cos(pi/2 - theta)^2;
r2 = radius^2;

zu = dp/na;
yu = [zu(2) -zu(1) 0];
yu = yu/norm(yu);
xu = [-zu(1)*zu(3) -zu(2)*zu(3) zu(2)*zu(2)+zu(1)*zu(1)];
xu = xu/norm(xu);

kx = 1:N(1);
ky = 1:N(2);
kz = 1:N(3);

center = 1+floor(N/2);
kx = kx - center(1);
ky = ky - center(2);
kz = kz - center(3);

delta_kx = FOV(1)/N(1);
delta_ky = FOV(2)/N(2);
delta_kz = FOV(3)/N(3);


kx = single(kx) * delta_kx -p1(1);
ky = single(ky) * delta_ky-p1(2);
kz = single(kz) * delta_kz-p1(3);

kx = reshape(kx,[length(kx),1,1]);
ky = reshape(ky,[1,length(ky),1]);
kz = reshape(kz,[1,1,length(kz)]);

kx = repmat(kx,[1,N(2),N(3)]);
ky = repmat(ky,[N(1),1,N(3)]);
kz = repmat(kz,[N(1),N(2),1]);


k2 = kx.^2 + ky.^2 + kz.^2;

chi = zeros(N);

D2 = ( k2 - (zu(1)*kx + zu(2)*ky + zu(3)*kz).^2 );
PX = (xu(1)*kx + xu(2)*ky + xu(3)*kz) ;
PY = (yu(1)*kx + yu(2)*ky + yu(3)*kz) ;
C2phi = cos( 2*( atan2( (PY),(PX) ) ) );
C2phi(isnan(C2phi(:))) = 0.0;

chi(D2 > r2) = xout;

chi(D2 <= r2) = xin;

% Bnuc 
dX = xin-xout;
Bnuc = (1-chi*2/3); % Lorentz sphere approximation for sensing proton
Bmac = 0.5.*r2.*dX.*st2.*C2phi./(D2);
Bmac(D2 <= r2) = dX*(3*ct2 -1)/6; % Macroscopic field due to magnetization


Bnuc = Bmac.*Bnuc;

end

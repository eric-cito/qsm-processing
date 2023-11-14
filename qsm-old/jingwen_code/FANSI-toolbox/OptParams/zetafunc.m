function [zetaval] = zetafunc(e1,e2)
% Calculate the normalized discrepancy function between the energy of two
% regions defined by masks in the frequency domain.
% The normalization penalizes solutions with too attenuated energies.
%
% Paramaters:
% e1, e2: Mean mask amplitudes of two regions as defined by the mean squared 
%         values of the Fourier coefficients inside each of them, or the mean 
%         value of the magnitude of those coefficients.
%
% Output:
% zetaval: normalized discrepancy value.
%
% Last modified by Carlos Milovic in 2019.11.21

num = (e1-e2).^2;
den = (e1+e2).^2;

zetaval = num./den;

end


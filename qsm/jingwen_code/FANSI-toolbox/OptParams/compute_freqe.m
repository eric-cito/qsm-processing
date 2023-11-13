function [ e1, e2, e3 ] = compute_freqe( chi, m1, m2, m3, power )
% Calculates the Mean Amplitude of the energy inside masks
% m1, m2 and (optional) m3 in the Fourier domain of QSM reconstruction chi.
%
% The estimation of the Mean Amplitude may be the mean magnitude 
% of the coefficients (power = 1) or the squared mean of the 
% coefficients (power = 2). Both are robust estimators.
%
% Last modified by Carlos Milovic in 2020.07.08
%

if nargin < 5
    power = 1;
end
    
fchi = (abs(fftn(chi)));

n1 = sum(m1(:));
n2 = sum(m2(:));
if nargin > 3
n3 = sum(m3(:));
end

if power == 1
e1 = sum ( fchi(:).*m1(:) )/n1;
e2 = sum ( fchi(:).*m2(:) )/n2;
if nargin > 3
e3 = sum ( fchi(:).*m3(:) )/n3;
end

else
e1 = sqrt(sum( abs((fchi(:).*m1(:) )).^2 )/(n1));
e2 = sqrt(sum( abs((fchi(:).*m2(:) )).^2 )/(n2));
if nargin > 3
e3 = sqrt(sum( abs((fchi(:).*m3(:) )).^2 )/(n3));
end
end

if nargin < 4
    e3 = 0.0;
end
    


end

function [ cc ] = compute_cc( chi_recon, chi_true )
% Calculates the Correlation Coefficient index between two images.
%
% Last modified by Carlos Milovic in 2017.03.30
%
R = corrcoef(chi_recon, chi_true);
cc = R(2,1);

end


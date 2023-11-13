function [ mad ] = compute_mad( chi_recon, chi_true )
% Calculates the L1 norm of the difference between two images (Mean Absolute Error)
% 
% Last modified by Carlos Milovic in 2017.03.30
%


mad = 100 * sum( abs(chi_recon(:) - chi_true(:)) ) / sum(abs((chi_true(:))));


end


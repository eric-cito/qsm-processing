function [ metrics ] = compute_metrics( chi_recon, chi_true, mode )
% Calculates quality indices between two images.
%
% Parameters:
% mode: 0 - basic mode. Only standard metrics (RMSE, HFEN, SSIM and XSIM)
%       1 - include extended metrics (CC and MI)
%       2 - include extended and gradient domain metrics
% See README.txt for further information.
% Last modified by Carlos Milovic in 2020.07.07
%

if nargin < 3
mode = 0;
end

% Standard metrics
metrics.rmse = compute_rmse( chi_recon, chi_true );
metrics.hfen = compute_hfen( chi_recon, chi_true );
metrics.ssim = compute_ssim( chi_recon, chi_true ); % Not recommended for analysis. Use XSIM instead.
metrics.xsim = compute_xsim( chi_recon, chi_true );

 if mode > 0
% Extended metrics
metrics.cc = compute_cc( chi_recon, chi_true );
metrics.mi = compute_mi( chi_recon, chi_true );

  if mode > 1
% Gradient domain metrics
metrics.gxe = compute_rmse_GX( chi_recon, chi_true );
metrics.mad = compute_mad( chi_recon, chi_true );
metrics.madgx = compute_mad_GX( chi_recon, chi_true );

  end
 end


end


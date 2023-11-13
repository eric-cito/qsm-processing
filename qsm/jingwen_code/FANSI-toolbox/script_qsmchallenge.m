% This is a sample script to show how to use the functions in this toolbox, 
% and how to set the principal variables.
% This example uses the data from the 2016 QSM reconstruction challenge.
% 4th International Workshop on Phase Contrast and QSM:
% http://qsm.rocks
% Please cite Langkammer C, et al. MRM 2017.
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Modified by Carlos Milovic in 2017.12.27
% Last modified by Carlos Milovic in 2020.07.13

%%%-------------------------------------------------------------------------
%% Load data
%%-------------------------------------------------------------------------


set(0,'DefaultFigureWindowStyle','docked')

addpath(genpath(pwd)) % Please run this while at the FANSI root folder


load phs_tissue;            % tissue phase from transversal orientation (in ppm, normalized by gyro*TE*B0)
load spatial_res;           % voxel size
load msk;                   % brain mask => obtained by eroding the BET mask by 5 voxels (by setting peel=5 in LBV)

load magn;                  % magnitude from transversal orientation

N = size(msk);



imagesc3d2(msk, N/2, 1, [90,90,-90], [0,1], [], 'Mask')

imagesc3d2(phs_tissue, N/2, 2, [90,90,-90], [-0.05,0.05], [], 'Input Phase')
imagesc3d2(magn, N/2, 3, [90,90,-90], [0,0.5], [], 'Magnitude')


%%-------------------------------------------------------------------------
%% Create dipole kernel
%%-------------------------------------------------------------------------

kernel = dipole_kernel_fansi( N, spatial_res, 0 );

%%-------------------------------------------------------------------------
%% TKD recon
%%-------------------------------------------------------------------------

thre_tkd = 0.19;      % TKD threshold parameter
chi_tkd = tkd( phs_tissue, msk, kernel, thre_tkd, N );

imagesc3d2(chi_tkd, N/2, 11, [90,90,-90], [-0.10,0.14], [], 'TKD')


%%-------------------------------------------------------------------------
%% Closed-form L2 recon
%%-------------------------------------------------------------------------

l2beta = 2e-2;    % regularization parameter
chi_L2 = chiL2( phs_tissue, msk, kernel, l2beta, N );

imagesc3d2(chi_L2, N/2, 12, [90,90,-90], [-0.10,0.14], [], 'CF L2')


%%-------------------------------------------------------------------------
%% Linear vs Weighted Linear TV
%%-------------------------------------------------------------------------

num_iter = 50;
tol_update = 1;

params = [];

params.alpha1 = 2e-4;               % gradient L1 penalty

params.maxOuterIter = num_iter;
params.tol_update = tol_update;

params.K = kernel;
params.input = phs_tissue;
 

out = wTV(params); 

magn_use = magn .* msk;
magn_use = sum(msk(:))*magn_use / sum(abs(magn_use(:)));

params.weight = magn_use;
outw = wTV(params); 
chiw = outw.x;

imagesc3d2(out.x.*msk, N/2, 21, [90,90,-90], [-0.10,0.14], [], 'TV')
imagesc3d2(chiw.*msk, N/2, 22, [90,90,-90], [-0.10,0.14], [], 'wTV')
    
    
%%-------------------------------------------------------------------------
%% Nonlinear TV
%%-------------------------------------------------------------------------
    
TE = 25e-3;
B0 = 2.8936;
gyro = 2*pi*42.58;
phs_scale = TE * gyro * B0;
 
params = [];
params.input = phs_tissue * phs_scale;
magn_use = magn .* msk;
magn_use = magn_use / max(abs(magn_use(:)));
params.weight = magn_use;
params.K = kernel;

params.alpha1 = 3e-4;               % gradient L1 penalty
params.mu1 = 1e-2;                  % gradient consistency

outw2 = wTV(params);
chiw2 = outw2.x/phs_scale;

outnl = nlTV(params); 
chinl = outnl.x/phs_scale;

imagesc3d2(chiw2.*msk, N/2, 23, [90,90,-90], [-0.1,0.14], [], 'wTV2')
imagesc3d2(chinl.*msk, N/2, 24, [90,90,-90], [-0.1,0.14], [], 'nlTV')


%%-------------------------------------------------------------------------
%% Optimize reconstructions according to:
% L-curve 
% Frequency Analysis
% Challenge metrics (RMSE, HFEN and SSIM)
%%-------------------------------------------------------------------------
load chi_cosmos
load chi_sti

params = [];
params.input = phs_tissue * phs_scale;
params.K = kernel;
magn_use = magn .* msk;
magn_use = magn_use / max(abs(magn_use(:)));%norm(abs(magn_use(:)));
params.weight = magn_use;  

% Create the masks used in the Frequency Analysis
[m1, m2, m3 ] = create_freqmasks( spatial_res, kernel );

% For simplicity, here we demostrate how to do it with the linear TV method, weighted by the magnitude.
for ka = 1:30
    
    alpha(ka) = 10^(-1.5-ka/10); % Regularization weight
    params.alpha1 = alpha(ka);
    params.mu1 = 25*alpha(ka);                  % gradient consistency
    
    out = wTV(params); 
    
    % Calculate the costs associated to each term of the optimization functional
    [ data_cost, reg_cost ] = compute_costs( out.x.*msk, params.input.*msk, kernel );    
    dcw(ka) = data_cost;
    rcw(ka) = reg_cost;
    
    % Calculate the Mean Amplitudes for each mask
    [e1(ka), e2(ka), e3(ka)] = compute_freqe(out.x.*msk,m1,m2,m3);
    
    % Error evaluation with respect to different ground-truths
    r_cosmos(ka) = compute_rmse(out.x.*msk/phs_scale,chi_cosmos);
    r_x33(ka) = compute_rmse(out.x.*msk/phs_scale,chi_33);
    
    s_cosmos(ka) = compute_ssim(out.x.*msk/phs_scale,chi_cosmos);
    s_x33(ka) = compute_ssim(out.x.*msk/phs_scale,chi_33);% only for the wTV result, for simplicity
    
    % You may also use:
    % metrics_cosmos = compute_metrics( out.x.*msk/phs_scale,chi_cosmos);
    % metrics_33 = compute_metrics( out.x.*msk/phs_scale,chi_33);
    % This will calculate the RMSE, HFEN, SSIM and XSIM.
    % Build the the vector with RMSE:
    % r_cosmos(ka) = metrics_cosmos.rmse; % And so on...
    
    % To calculate the Correlation Coefficient and Mutual Information too you may use:
    % metrics_cosmos = compute_metrics( out.x.*msk/phs_scale,chi_cosmos,1);
    % To further include error metrics of the gradient domain, use:
    % metrics_cosmos = compute_metrics( out.x.*msk/phs_scale,chi_cosmos,2);
end

% L-curve analysis
[ Kappa ] = draw_lcurve( alpha, rcw, dcw, 31 );
[ KappaM ] = draw_lcurve_median( alpha, rcw, dcw, 32 ); % Use a median filter to smooth the curvature data.
% Optimal recommended reconstructions are those closer to the inflection point in the L-curve, i.e when
% the curvature is zero (changes its sign). Do the analysis from the right to the left to avoid unstabilities

% Frequency analysis
[opt23index, alpha_opt] = draw_freque(alpha,e1,e2,e3,33);
% The optimal reconstruction is found when the curves that describes the mask amplitudes for M2 and M3 
% (e2 and e3) intercept, or when zeta23 is minimized.
% Since this reconstruction tends to be slightly under-regularized, you may consider applying an
% up to 2x factor to this optimal predicted reconstruction to further reduce noise.

% Plot the RMSE scores
figure(28)
semilogx(alpha,r_cosmos,'b');
hold on;
semilogx(alpha,r_x33,'r');
hold off;

% Plot the SSIM scores
figure(51)
semilogx(alpha,s_cosmos,'b');
hold on;
semilogx(alpha,s_x33,'r');
hold off;





% This is a sample script to show how to use the functions in this toolbox, 
% and how to set the principal variables.
% This example uses a brain phantom based on a COSMOS reconstruction, part of the
% 2016 QSM Reconstruction Challenge dataset.
% Please cite Langkammer C, et al. MRM 2017.
%
% Based on the code by Bilgic Berkin at http://martinos.org/~berkin/software.html
% Last modified by Carlos Milovic in 2017.12.27
%

set(0,'DefaultFigureWindowStyle','docked')
addpath(genpath(pwd)) % Please run this while at the FANSI root folder

%% load COSMOS data

load spatial_res;           % voxel size
load msk;                   % brain mask => obtained by eroding the BET mask by 5 voxels (by setting peel=5 in LBV)


load magn;                  % magnitude from transversal orientation
load chi_cosmos;            % COSMOS from 12 orientations (in ppm)

N = size(chi_cosmos);

center = N/2 + 1;

TE = 25e-3;
B0 = 2.8936;
gyro = 2*pi*42.58;

phs_scale = TE * gyro * B0;

imagesc3d2(chi_cosmos, N/2, 1, [90,90,-90], [-0.10,0.14], [], 'COSMOS')


%% Create dipole kernel and susceptibility to field model
kernel = dipole_kernel_fansi( N, spatial_res, 0 ); % 0 for the continuous kernel by Salomir and Marques and Bowtell.


% Create raw phase
rphase = ifftn(fftn(chi_cosmos).*kernel);

signal = magn.*exp(1i*phs_scale*(rphase))+0.01*(randn(size(magn))+1i*randn(size(magn))); %SNR = 40
phase_use = angle(signal)/phs_scale;
mask_use = msk;
magn_use = magn.*mask_use;

clear magn msk signal rphase
imagesc3d2(magn_use, N/2, 2, [90,90,-90], [0,1], [], 'Magn')
imagesc3d2(mask_use, N/2, 3, [90,90,-90], [0,1], [], 'Mask')
imagesc3d2(phase_use, N/2, 4, [90,90,-90], [-0.10,0.14], [], 'Local Phase')

%% Reconstruction examples 
% Weighted Linear TV and TGV ADMM reconstructions


% Common parameters
params = [];

params.K = kernel;
params.input = mask_use.*phase_use*phs_scale; % Linear may work with PPM
% We use radians to mantain consistancy with nonlinear parameters


alpha = 8e-3;
params.alpha1 = alpha;
params.mu1 = 10*alpha;


% Linear, ROI weight
params.weight = mask_use; % Note that the input phase is masked.

% TV
out = wTV(params);
metrics_tv = compute_metrics(out.x.*mask_use/phs_scale,chi_cosmos); % Note that the solution is rescaled to ppm
imagesc3d2(mask_use.*out.x/phs_scale, N/2, 8, [90,90,-90], [-0.10,0.14], [], 'QSM: TV')

% TGV
out2 = wTGV(params);
metrics_tgv = compute_metrics(out2.x.*mask_use/phs_scale,chi_cosmos);
imagesc3d2(mask_use.*out2.x/phs_scale, N/2, 9, [90,90,-90], [-0.10,0.14], [], 'QSM: TGV')


% Linear, magnitude weighted
alpha = 4e-4;
params.alpha1 = alpha;
params.mu1 = 10*alpha;
params.weight = magn_use; % Magnitude data is in the [0,1] range

% TV
wout = wTV(params);
metrics_wtv = compute_metrics(out.x.*mask_use/phs_scale,chi_cosmos);
imagesc3d2(mask_use.*wout.x/phs_scale, N/2, 12, [90,90,-90], [-0.10,0.14], [], 'QSM: Weighted TV')

% TGV
wout2 = wTGV(params);
metrics_wtgv = compute_metrics(out2.x.*mask_use/phs_scale,chi_cosmos);
imagesc3d2(mask_use.*wout2.x/phs_scale, N/2, 13, [90,90,-90], [-0.10,0.14], [], 'QSM: Weighted TGV')


% Nonlinear reconstructions, with the same parameters

outnl = nlTV(params);
chi_nl = outnl.x/phs_scale; % Alternatively, another variable may store the rescaled result.
metrics_nltv = compute_metrics(chi_nl.*mask_use,chi_cosmos);
imagesc3d2(mask_use.*chi_nl, N/2, 14, [90,90,-90], [-0.10,0.14], [], 'QSM: nTV')

outnl2 = nlTGV(params);
chi_nl2 = outnl2.x/phs_scale;
metrics_nltgv = compute_metrics(chi_nl2.*mask_use,chi_cosmos);
imagesc3d2(mask_use.*chi_nl2, N/2, 15, [90,90,-90], [-0.10,0.14], [], 'QSM: nTGV')



%% Optimize for RMSE Example
% Fine tune the choise of the regularization weight, alpha.

params = [];
params.K = kernel;
params.input = mask_use.*phase_use*phs_scale; % Linear may work with PPM
%params.maxOuterIter = 50; % Lets do a fast example
%params.tol_update = 1;
% We'll find the optimal alpha for linear and nonlinear methods and different weighting, using TV.
for ka = 1:30
    disp(['Step number: ', num2str(ka)])
    % Set the regularization weight dependand parameters
    alpha(ka) = 10^(-1.5-ka/10);
    params.alpha1 = alpha(ka);
    params.mu1 = 10*alpha(ka);   
    
    
    params.weight = mask_use;
    out = wTV(params);
    chi = out.x/phs_scale;
    rmse_l2(ka) = compute_rmse(chi.*mask_use,chi_cosmos);
    
    params.weight = magn_use;
    out = wTV(params); % We are reusing variable names to save memory space
    chi = out.x/phs_scale;
    rmse_w(ka) = compute_rmse(chi.*mask_use,chi_cosmos);
    
    out = nlTV(params);
    chi = out.x/phs_scale;
    rmse_nl(ka) = compute_rmse(chi.*mask_use,chi_cosmos);
end


% Plot the RMSE curves
figure(20);
loglog(alpha,rmse_l2,'g');
hold on;
loglog(alpha,rmse_w,'r');
loglog(alpha,rmse_nl,'b');
hold off;

% Find the optimal value, kval, and its index ka:
[kval, ka] = min(rmse_l2);
% Optimal reconstruction
params.weight = ones(size(phase_use));
alphaOpt = 10^(-1.5-ka/10);
params.alpha1 = alphaOpt;
params.mu1 = 10*alphaOpt;                  % gradient consistency
out = wTV(params);
chi = out.x/phs_scale;

[kval, ka] = min(rmse_w);
params.weight = magn_use;
alphaOpt = 10^(-1.5-ka/10);
params.alpha1 = alphaOpt;
params.mu1 = 10*alphaOpt;                  % gradient consistency
outw = wTV(params);
chiw = outw.x/phs_scale;

[kval, ka] = min(rmse_nl);
params.weight = magn_use;
alphaOpt = 10^(-1.5-ka/10);
params.alpha1 = alphaOpt;
params.mu1 = 10*alphaOpt;         
outnl = nlTV(params);
chinl = outnl.x/phs_scale;


    
   
    

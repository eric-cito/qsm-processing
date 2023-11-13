% This is a sample script to show how to use the functions in this toolbox, 
% and how to set the principal variables.
% This example uses a data set from the Wellcome Trust Centre for
% NeuroImaging, University College London, UK. *** NOT PUBLICY AVAILABLE ***
% Use this as an example, to adapt to your own data.
%
% The original data was preprocessed using the MEDI Toolbox (Cornell).
% We kept the same nomenclature to facilitate the use of our toolbox.
%
% Last modified by Carlos Milovic in 2020.07.07

%%-------------------------------------------------------------------------
%% Load data
%%-------------------------------------------------------------------------


set(0,'DefaultFigureWindowStyle','docked')
addpath(genpath(pwd)) % Please run this while at the FANSI root folder

load RDF.mat
% This contained: 
% the local field map, RDF
% ROI mask, Mask
% magnitute, iMag
% vector with voxel dimensions in mm, voxel_size
% main field direction vector, B0_dir
% noise estimate, N_std

N = size(Mask);



imagesc3d2(Mask, N/2, 1, [90,90,-90], [0,1], [], 'Mask')
imagesc3d2(RDF, N/2, 2, [90,90,-90], [-1,1], [], 'Input Phase')
imagesc3d2(iMag/max(iMag(:)), N/2, 3, [90,90,-90], [0,0.5], [], 'Magnitude')



%%-------------------------------------------------------------------------
%% Create dipole kernel
%%-------------------------------------------------------------------------
spatial_res = voxel_size;

kernel = dipole_kernel_angulated( N, spatial_res, B0_dir );



%%-------------------------------------------------------------------------
%% Normalize and define variables
%%-------------------------------------------------------------------------
phs_use = RDF;
magn_use = min(2*Mask.*iMag/max(iMag(:)),1.0);
imagesc3d2(magn_use, N/2, 4, [90,90,-90], [0,1], [], 'rMagnitude')
imagesc3d2(N_std, N/2, 5, [90,90,-90], [0,2e-3], [], 'noise')

%% Phase to PPM scale

TE = 20e-3; % echo time in seconds
B0 = 3; % main field strength
gyro = 2*pi*42.58; % gyromagnetic ratio

phs_scale = TE * gyro * B0; % ppm-to-radian scaling factor



%%-------------------------------------------------------------------------
%% L-curve and Frequency analysis
%%-------------------------------------------------------------------------

% Create the masks used in the Frequency Analysis
[m1, m2, m3 ] = create_freqmasks( spatial_res, kernel );

% Create structure with parameters
params = [];

params.input = phs_use;
params.K = kernel;
params.maxOuterIter = 50; % To speed up the L-curve process, stop at an early convergence
params.tol_update = 1;

alpha = 10.^(-(1:17)/5-1);
for ka = 1:17
    
    params.alpha1 = alpha(ka);
    params.mu1 = 100*alpha(ka);    
    
    
    params.weight = magn_use;    
    out = wTV(params); 
    
    % Calculate the costs associated to each term of the optimization functional
    [ data_cost, reg_cost ] = compute_costs( out.x.*Mask, params.input.*Mask, kernel );    
    dcw(ka) = data_cost;
    rcw(ka) = reg_cost;
    
    % Calculate the Mean Amplitudes for each mask
    [e1(ka), e2(ka), e3(ka)] = compute_freqe(out.x.*Mask,m1,m2,m3);
end

% L-Curve analysis
[ Kappa ] = draw_lcurve( alpha, rcw, dcw, 11 );
[ KappaM ] = draw_lcurve_median( alpha, rcw, dcw, 12 );


% Frequency analysis
[opt23index, alpha_opt] = draw_freque(alpha,e1,e2,e3,13);

% Example optimal reconstruction
 
    params.alpha1 = 2.5*alpha(opt23index); % Slightly over-regularize to avoid noise amplification
    params.mu1 = 2.5*100*alpha(opt23index);    
    
    out = wTV(params); 
imagesc3d2(Mask.*out.x/phs_scale, N/2, 14, [90,90,-90], [-0.1,0.1], [], 'Frequency optimal')

%%-------------------------------------------------------------------------
% Example, FANSI wildcard function
%%-------------------------------------------------------------------------

options = [];
options.tgv = false;
options.nonlinear = true;

outf = FANSI( phs_use, magn_use, spatial_res, 0.0015, 0.0054, options, B0_dir );%alpha(9)


imagesc3d2(Mask.*outf.x, N/2, 21, [90,90,-90], [-1,1], [], 'Xfansi')

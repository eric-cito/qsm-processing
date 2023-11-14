function [param chi1_est chi_est chi_est_ppm fmap_est f1 obj_vals noise_thr fmap_lp mask1 mask2 mask3 nfm_final fmap_acq_corr, scales] = KOP_intsources_nfm_corr3_7T(subjnum_sph, subjnum_rois, ext_src_dir, dims) 
%
%
% subjnum_rois = '2020_2768'
%
% subjnum_ph = '0011_12192007'
%
%
%  ROI 1  %Right Pallidum/globus pallidus 
%  ROI 2  %Right Putamen 
%  ROI 3  %Right Caudate 
%  ROI 4  %Right Thalamus 
%  ROI 5  %Left Pallidum/globus pallidus 
%  ROI 6  %Left Putamen 
%  ROI 7  %Left Caudate 
%  ROI 8  %Left Thalamus
%
%%

%% SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% weights (lambdas)

LPWeight = 1e-8;          % 1e-8;     % term 1: penalty for disagreement between Laplacian of the estimated and true fieldmaps
OBJWeight = 1e-6;         % 1e-6;     % term 2: penalty for disagreement between estimated and 'true' (biasfield corrected) fieldmap 
EXTWeight = 1e10;         % 1e10;     % term 3: penalty for non-zero external susceptibility sources


% optimization parameters

lineSearchT0 = 0.1;      % step size to start with
Itnlim = 90; %30;             % number of CG iterations 
pNorm = 1;               % norm to use with term 1 (ie pNorm = 1 -> L1 norm)


% parameters for calculating the spatial dipole kernel from spatial (perturbation) formulation of the forward model

B0_b0calc = -1 * 2.675e8;   % rad/sec -> the estimated fieldmap will be in radians/second. (Do not change).
chi0 = 0 ; %0.4e-6;              % unitless   
delta = -9.1e-6 ; %-9.5e-6;            % unitless

% initial estimate of susceptibility (chi1) values

chi1_null = 0; %0.0421;                                 % unitless, chi1_null is the chi1 value that corresponds to chi = 0. (ie. solve chi = chi0 - delta * chi1 for chi1).
chi1_null_vol = chi1_null * ones(dims(1), dims(2), dims(3));     % set all voxels to chi1_null (the initial estimate of chi1 values)

% parameters for normalizing / converting the observed fieldmap to radians/second

B0 = 7.0;               % field strength of acquisition ( Tesla )
gyro = 2*pi*42.58;      % gyromagnetic ratio ( MHz/T * radians )
TE = 28e-3;             % echo time ( seconds )

%% READ IN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized field from sources outside the brain, unitless, in ppm
%filename = strcat(ext_src_dir,'/RORR_nfm_dipole'); 
filename = strcat(ext_src_dir,'/nfm_dipole'); 
nfm_dipole = read_avw(filename);
%nfm_dipole = imrotate(nfm_dipole,90);

% binary brain mask ( M_{0} )
%filename = 'RORR_mag_cmb_mask_ero_sphere20x2';  
filename = 'ROR_mag_cmb_mask_ero_sphere20x2'; 
[brain_mask, N, scales] = read_avw(filename);     
%brain_mask = imrotate(brain_mask,90);

% observed, unwrapped phase, radians
%filename = 'RORR_phase_cmb'; 
filename = 'ROR_phase_cmb'; 
uphase = read_avw(filename);
%uphase = imrotate(double(uphase.img),90);
uphase = uphase .* brain_mask;

% corresponding magnitude data
%filename = 'RORR_mag_cmb'; 
filename = 'ROR_mag_cmb'; 
mag = read_avw(filename);
%mag = imrotate(double(mag.img),90);

% binary mask of a region outside the head for computing noise
%noise_mask = read_avw('RORR_noise_mask');    
noise_mask = read_avw('ROR_noise_mask');  
%noise_mask = imrotate(double(noise_mask),90);

% ROI segmentation map
%filename = strcat(subjnum_rois,'_fdri_gp+put+caud+thal_roi_axp_2.img');
%[roi_labels] = load_nii(filename);
%roi_labels = imrotate(double(roi_labels.img),90);

%% GET ROI INDICES (RIGHT, LEFT, BILATERAL) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ind_pt1 = find(roi_labels == 2);
% ind_pt2 = find(roi_labels == 6);
% ind_pt = [ind_pt1;ind_pt2];
% ind_cd1 = find(roi_labels == 3);
% ind_cd2 = find(roi_labels == 7);
% ind_cd = [ind_cd1;ind_cd2];
% ind_pd1 = find(roi_labels == 1);
% ind_pd2 = find(roi_labels == 5);
% ind_pd = [ind_pd1;ind_pd2];
% ind_th1 = find(roi_labels == 4);
% ind_th2 = find(roi_labels == 8);
% ind_th = [ind_th1;ind_th2];

%% COMPUTE THE NORMALIZED OBSERVED FIELDMAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nfm = -uphase / (B0 * gyro * TE);                       % unitless, in ppm
nfm_pad = padarray(nfm,[(N(1)/2),(N(2)/2),(N(3)/2)]);


%% COMPUTE THE BIAS-CORRECTED OBSERVED FIELDMAP IN RAD/SEC ( THE FIELD DUE TO INTERNAL SOURCES ),  B_{i} =  B - B_{e}  

 
 nfm_final = nfm_pad - nfm_dipole;      % subtract the normalized field due to external sources from the normalized observed field
 
 %fmap_acq = nfm_final(257-128:256+128, 257-128:256+128, 65-32:64+32);
 fmap_acq = nfm_final((N(1)+1)-(N(1)/2):N(2)+(N(1)/2),(N(1)+1)-(N(1)/2):N(2)+(N(1)/2), (N(3)+1)-(N(3)/2):N(3)+(N(3)/2));
 
 %fmap_acq_corr = fmap_acq * -1 * 2.675e8;          %now in rad/sec * ppm
 
 fmap_acq_corr = fmap_acq * -1 * 2.675e8 * 1e-6;    %now in rad/sec    
 

%% COMPUTE THE SPATIAL DIPOLE KERNEL 
%
% The estimated fieldmap can then be computed by: K * chi1

[d2G_dz2_kernel] = d2Gdz2_kernel(N, scales);

K = KOP(d2G_dz2_kernel, B0_b0calc, chi0, delta);

%% COMPUTE MASKS FOR OBJECTIVE FUNCTION

% Term 1 mask ( W ) 

ind_noise = find(noise_mask == 1);
noise_vox = mag(ind_noise);
noise_thr = 2 * std(noise_vox)

ind_non_brain = find(brain_mask == 0);
ind_brain = find(brain_mask == 1);

% Compute Laplacian of observed fmap (in rad/sec)

fmap_lp = 0 * nfm;

% nfm_radsec = nfm * -1 * 2.675e8 * 1e-6;     % observed fmap in rad/sec
% 
% fmap_lp(2:end-1,2:end-1,2:end-1) =  -6 * nfm_radsec(2:end-1,2:end-1,2:end-1) + nfm_radsec(1:end-2,2:end-1,2:end-1) + nfm_radsec(3:end,2:end-1,2:end-1) ...
%                                          + nfm_radsec(2:end-1,1:end-2,2:end-1) + nfm_radsec(2:end-1,3:end,2:end-1)+ nfm_radsec(2:end-1,2:end-1,1:end-2) + nfm_radsec(2:end-1,2:end-1,3:end);                                 

% Note: Laplacian(B) = Laplacian(B_{i} + B_{e}) = Laplacian(B_{i}) assuming Laplacian(B_{e}) = 0 (ie. field from external sources is soln to Laplace Eq).

fmap_lp(2:end-1,2:end-1,2:end-1) =  -6 * fmap_acq_corr(2:end-1,2:end-1,2:end-1) + fmap_acq_corr(1:end-2,2:end-1,2:end-1) + fmap_acq_corr(3:end,2:end-1,2:end-1) ...
                                         + fmap_acq_corr(2:end-1,1:end-2,2:end-1) + fmap_acq_corr(2:end-1,3:end,2:end-1)+ fmap_acq_corr(2:end-1,2:end-1,1:end-2) + fmap_acq_corr(2:end-1,2:end-1,3:end);                               

                                     
fmap_lp = fmap_lp .* brain_mask;
                                     
mask1 = abs(fmap_lp);      

ind = find(mask1 < noise_thr);

mask1(ind) = 0;                      % TERM 1 MASK ( W )                                       

%---------------------------------------------------------------------

mask2 = brain_mask;                  % TERM 2 MASK ( M0 )

%---------------------------------------------------------------------

mask3 = zeros(N(1), N(2), N(3)); 
           
mask3(ind_non_brain) = ones(size(ind_non_brain),1);          % TERM 3 MASK ( M0C )    

%% INITIALIZE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param = init;                       % init.m returns a structure with default values
                            

%param.fmap_acq = nfm_radsec;                            
param.fmap_acq = fmap_acq_corr;    % param.fmap_acq = the observed field, B, in term 1.  Can set B = B_{i} if we assume Laplacian of B_{e}=0, then Laplacian(B)=Laplacian(B_{i})
param.data = fmap_acq_corr;         % param.data = the bias-corrected observed field, B_{i}, in term 2.  (B_{i} = B - B_{e})


param.x0 =  chi1_null_vol;                  % initial estimate of chi1 values    
param.chi1_null_vol = chi1_null_vol;        % used in term 3 to enforce chi outside brain = 0 (ie. chi1 = chi1_null outside brain)            


param.pNorm = pNorm;                    % type of norm to use (i.e. L1 L2 etc)
param.lineSearchT0 = lineSearchT0;      % step size to start with
param.Itnlim = Itnlim;                  % number of CG iterations 

param.K = K;  

param.LPWeight = LPWeight;           % weight for term 1
param.OBJWeight = OBJWeight;         % weight for term 2
param.EXTWeight = EXTWeight;         % weight for term 3

param.mask1 = mask1;                    % mask for term 1
param.mask2 = mask2;                    % mask for term 2
param.mask3 = mask3;                    % mask for term 3

param.scales = scales;
param.delta = delta;
param.chi0 = chi0;



%% OPTIMIZE: SEARCH FOR OPTIMAL CHI1 VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
[x, f1, obj_vals] = fnlCg_sc(param);
toc

chi1_est = x;                                               % Estimated chi1 values, unitless

chi_est = chi0 * ones(N(1),N(2),N(3)) + delta * chi1_est;   % Estimated chi values, unitless

fmap_est = K*chi1_est;                                      % Corresponding estimated fmap, radians/sec

%convert estimated chi values to ppm for easier visualization

chi_est_ppm = chi_est * 1e6;

%calc mean relative chi in each ROI
% susc_pt = chi_est(ind_pt);
% susc_cd = chi_est(ind_cd);
% susc_pd = chi_est(ind_pd);
% susc_th = chi_est(ind_th);
% 
% rois = ['th ','cd ','pt ','gp'];
% 
% mean_susc_rel = [mean(susc_th),  mean(susc_cd), mean(susc_pt), mean(susc_pd)]


% undo rotation of images
% chi1_est = imrotate(chi1_est,-90);
% chi_est = imrotate(chi_est,-90);
% chi_est_ppm = imrotate(chi_est_ppm,-90);
% fmap_est = imrotate(fmap_est,-90);
% fmap_lp = imrotate(fmap_lp,-90);
% mask1 = imrotate(mask1,-90);
% mask2 = imrotate(mask2,-90);
% mask3 = imrotate(mask3,-90);
% nfm_final = imrotate(nfm_final,-90);
% fmap_acq_corr = imrotate(fmap_acq_corr,-90);







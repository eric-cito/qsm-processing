function [param chi1_est chi_est chi_est_ppm fmap_est f1 obj_vals ...
    noise_thr fmap_lp mask1 mask2 mask3 nfm_final fmap_acq_corr, scales] = ...
    KOP_intsources_nfm_corr3_new(nfm_dipole, mask, magni, phase_uw, noise_mask, params) 

%% SET PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dims = size(phase_uw);
scales = params.voxelsize;


% weights (lambdas)

LPWeight = 1e-8;          % 1e-8;     % term 1: penalty for disagreement between Laplacian of the estimated and true fieldmaps
OBJWeight = 1e-6;         % 1e-6;     % term 2: penalty for disagreement between estimated and 'true' (biasfield corrected) fieldmap 
EXTWeight = 1e10;         % 1e10;     % term 3: penalty for non-zero external susceptibility sources


% optimization parameters

lineSearchT0 = 0.1;      % step size to start with
Itnlim = 50;             % number of CG iterations 
pNorm = 1;               % norm to use with term 1 (ie pNorm = 1 -> L1 norm)


% parameters for calculating the spatial dipole kernel from spatial (perturbation) formulation of the forward model

B0_b0calc = -1 * 2.675e8;   % rad/sec -> the estimated fieldmap will be in radians/second. (Do not change).
%chi0 = 0.4e-6;              % unitless   
%delta = -9.5e-6;            % unitless
chi0 = 0;              % unitless   
delta = -9.1e-6;            % unitless

% initial estimate of susceptibility (chi1) values

%chi1_null = 0.0421;                                 % unitless, chi1_null is the chi1 value that corresponds to chi = 0. (ie. solve chi = chi0 - delta * chi1 for chi1).
chi1_null = 0;
chi1_null_vol = chi1_null * ones(dims(1), dims(2), dims(3));     % set all voxels to chi1_null (the initial estimate of chi1 values)

% parameters for normalizing / converting the observed fieldmap to radians/second

B0 = params.bo;               % field strength of acquisition ( Tesla )
gyro = 2*pi*42.58;      % gyromagnetic ratio ( MHz/T * radians )
TE = params.te;             % echo time ( seconds )


%% READ IN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized field from sources outside the brain, unitless, in ppm
%filename = strcat(ext_src_dir,'/nfm_dipole');  
% filename = strcat(ext_src_dir,'/nfm_dipole'); 
% nfm_dipole = read_avw(filename);
%nfm_dipole = imrotate(nfm_dipole,90);

% filename = 'ROR_phase_bet_mask2_erode5';
% [brain_mask, dims, scales] = read_avw(filename);     
%brain_mask = imrotate(brain_mask,90);

% observed, unwrapped phase, radians
%filename = 'RO_phase_nan'; 
% filename = 'ROR_phase_cmb'; 
% uphase = read_avw(filename);
%uphase = imrotate(double(uphase.img),90);


% corresponding magnitude data
% filename = 'ROR_mag_cmb'; 
% mag = read_avw(filename);
%mag = imrotate(double(mag.img),90);

% binary mask of a region outside the head for computing noise
% noise_mask = read_avw('ROR_noise_mask');         

%% COMPUTE THE NORMALIZED OBSERVED FIELDMAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phase_uw = phase_uw .* mask;
nfm = -phase_uw / (B0 * gyro * TE);                       % unitless, in ppm
nfm_pad = padarray(nfm,[(dims(1)/2),(dims(2)/2),(dims(3)/2)]);


%% COMPUTE THE BIAS-CORRECTED OBSERVED FIELDMAP IN RAD/SEC ( THE FIELD DUE TO INTERNAL SOURCES ),  B_{i} =  B - B_{e}  

 
 nfm_final = nfm_pad - nfm_dipole;      % subtract the normalized field due to external sources from the normalized observed field
 fmap_acq = nfm_final((dims(1)+1)-(dims(1)/2):dims(1)+(dims(1)/2), (dims(2)+1)-(dims(2)/2):dims(2)+(dims(2)/2), (dims(3)+1)-(dims(3)/2):dims(3)+(dims(3)/2));
 fmap_acq_corr = fmap_acq * -1 * 2.675e8 * 1e-6;    %now in rad/sec    
 

%% COMPUTE THE SPATIAL DIPOLE KERNEL 
%
% The estimated fieldmap can then be computed by: K * chi1

[d2G_dz2_kernel] = d2Gdz2_kernel(dims, scales);

K = KOP(d2G_dz2_kernel, B0_b0calc, chi0, delta);

%% COMPUTE MASKS FOR OBJECTIVE FUNCTION

% Term 1 mask ( W ) 

ind_noise = find(noise_mask == 1);
noise_vox = magni(ind_noise);
noise_thr = 2 * std(noise_vox)

ind_non_brain = find(mask == 0);
ind_brain = find(mask == 1);

% Compute Laplacian of observed fmap (in rad/sec)

fmap_lp = 0 * nfm;
fmap_lp(2:end-1,2:end-1,2:end-1) =  -6 * fmap_acq_corr(2:end-1,2:end-1,2:end-1) + fmap_acq_corr(1:end-2,2:end-1,2:end-1) + fmap_acq_corr(3:end,2:end-1,2:end-1) ...
                                         + fmap_acq_corr(2:end-1,1:end-2,2:end-1) + fmap_acq_corr(2:end-1,3:end,2:end-1)+ fmap_acq_corr(2:end-1,2:end-1,1:end-2) + fmap_acq_corr(2:end-1,2:end-1,3:end);                               

                                     
fmap_lp = fmap_lp .* mask;                                    
mask1 = abs(fmap_lp);      
mask1(mask1 < noise_thr) = 0;                      % TERM 1 MASK ( W )                                       

%---------------------------------------------------------------------

mask2 = mask;                  % TERM 2 MASK ( M0 )

%---------------------------------------------------------------------

mask3 = zeros(dims(1), dims(2), dims(3)); 
           
mask3(ind_non_brain) = 1; % ones(size(ind_non_brain),1);          % TERM 3 MASK ( M0C )    

%% INITIALIZE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param = qsip_init;                       % init.m returns a structure with default values
                            

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

chi_est = chi0 * ones(dims(1),dims(2),dims(3)) + delta * chi1_est;   % Estimated chi values, unitless

fmap_est = K*chi1_est;                                      % Corresponding estimated fmap, radians/sec

%convert estimated chi values to ppm for easier visualization

chi_est_ppm = chi_est * 1e6;








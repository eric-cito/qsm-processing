%% navigate MATLAB to the 'eval' folder
% cd('~/deepMRI/xQSM/matlab/eval');

%% add NIFTI matlab toolbox for read and write NIFTI format
deepMRI_root = '/working/lupolab/jingwen/001_QSM/01_Code/deepMRI/'; % where deepMRI git repo is cloned to; change as yours; 
addpath(genpath(deepMRI_root));  %  add NIFTI saving and loading functions;

%% read in field map and COSMOS map (3rd dimension is z/B0 direction)

cd /working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308/output_HDBET/;
system('3dcopy localField_meanEcho_crop_ppm.nii.gz localField_meanEcho_crop_ppm.nii');

nii = load_nii('localField_meanEcho_crop_ppm.nii'); 
% replace the file name with yours. 
field = double(nii.img);
mask = field ~= 0; % brain tissue mask

% note the size of the field map input needs to be divisibel by 8
% otherwise 0 padding should be done first
imSize = size(field);
if mod(imSize, 8)
    [field, pos] = ZeroPadding(field, 8);
end

% illustration of one central axial slice of the input field 
figure;
imagesc(field(:,:,80)'); colormap gray; axis equal tight; colorbar; caxis([-0.05, 0.05]);
title('Slice 80 of the Input Field Map (ppm)');
drawnow;

%% read label (for evaluation purpose)
% nii = load_nii('../../cosmos_label.nii'); % replace the file name with yours. 
% label = double(nii.img);

%% label image normalization (mean of brain tissue region set to 0) for later comparison;
% label = label - sum(label(:)) / sum(mask(:));
% label = label .* mask;

% illustration of one central axial slice of the COSMOS label 
% figure; 
% imagesc(label(:,:,80)'); colormap gray; axis equal tight; colorbar; caxis([-0.1, 0.2]);
% title('Slice 80 of the COSMOS Label (ppm)');
% drawnow;

%% start recons
recon_methods_list = {'xQSM_invivo', 'xQSM_syn', 'xQSM_invivo_withNoiseLayer', 'Unet_invivo', 'Unet_syn'};

for k = 1:length(recon_methods_list)
    recon_method = recon_methods_list{k};
    fprintf('Reconstructing QSM using %s\n', recon_method);

    if canUseGPU()
        % (1) if your MATLAB is configured with CUDA GPU acceleration
        QSM_recon = Eval(field, recon_method, 'gpu');
    else
        % (2) otherwise if CUDA is not available, use CPU instead, this is much slower
        QSM_recon = Eval(field, recon_method, 'cpu');
    end

    % if zeropadding was performed, then do zero-removing before next step;
    if mod(imSize, 8)
        QSM_recon = ZeroRemoving(QSM_recon, pos);
    end

    %% image normalization (mean of brain tissue region set to 0)
    QSM_recon = QSM_recon - sum(QSM_recon(:)) / sum(mask(:));
    QSM_recon = -QSM_recon .* mask; 

    %% illustration of one central axial slice of the four different reconstructions; 
    figure,
    subplot(121), imagesc(QSM_recon(:,:,80)'); colormap gray; axis equal tight; colorbar; 
    caxis([-0.1, 0.15])
    title(['Slice 80 of the ', recon_method], 'Interpreter', 'none');
%     err  = QSM_recon - label;
%     subplot(122), imagesc(err(:,:,80)'); colormap gray; axis equal tight; colorbar; caxis([-0.1, 0.2])
%     title('Error');
    drawnow;
    
    %% use default pnsr and ssim 
%     PSNR_recon = psnr(QSM_recon, single(label));
%     fprintf('PSNR of %s is %f\n', recon_method, PSNR_recon);
%     SSIM_recon = ssim(QSM_recon, single(label));
%     fprintf('SSIM of %s is %f\n\n', recon_method, SSIM_recon);

    %% save the files for ROI measurements; 
    nii = make_nii(QSM_recon, [1, 1, 1]);
    save_nii(nii, ['Chi_', recon_method, '.nii']);
end

%% This demo shows the complete reconstruction pipeline for iQSM on Multi-echo MRI phase data
%% Assume your raw phase data is in NIFTI format


% (1) download or clone github repo for deepMRI: https://github.com/sunhongfu/deepMRI
% (2) download demo data and checkpoints here: https://www.dropbox.com/sh/9kmbytgf3jpj7bh/AACUZJ1KlJ1AFCPMIVyRFJi5a?dl=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data preparation guide: 

% 1. phase evolution type:
% The relationship between the phase data and filed pertubation (delta_B) 
% is assumed to satisfy the following equation: 
% "phase = -delta_B * gamma * TE" 
% Therefore, if your phase data is in the format of "phase = delta_B * gamma * TE;" 
% it will have to be preprocessed by multiplication by -1; 

% 2. For Ultra-high resolutin data:
% it is recommended that the phase data of ultra-high resolution (higher
% than 0.7 mm) should be interpoloated into 1 mm for better reconstruction results.  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set your own data paths and parameters
checkpoints_root  = [deepMRI_root '/iQSM/checkpoints'];

subj_path = '/working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308/output_HDBET/';

load([subj_path '/sepia_header.mat']);

PhasePath    = [subj_path 'Phase.nii.gz'];  % where raw phase data is (in NIFTI format)
ReconDir     = subj_path;  %% where to save reconstruction output
Eroded_voxel = 0;  %  set number of voxels for brain mask erosion; 0 means no erosion;
TE           = header.te; % set Echo Times (in second)
B0           = header.b0; % set B0 field (in Tesla)
vox          = header.voxelSize; % set voxel size a.k.a image resolution (in millimeter)
NetworkType  = 0; % network type: 0 for original iQSM, 1 for networks trained with data fidelity,
                  % 2 for networks trained with learnable Lap-Layer (15 learnable kernels) and data fidelity;

%% optional data paths to be set, simply comment out if not available
MaskPath = [subj_path 'brain_mask_HD.nii.gz']; %% brain mask; set to one will skip brain masking
MagPath = [subj_path 'Magni.nii.gz']; % magnitude image; set to one will skip magnitude weights in echo fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% add MATLAB paths
addpath(genpath([deepMRI_root,'/iQSM/iQSM_fcns/']));  % add necessary utility function for saving data and echo-fitting;
addpath(genpath([deepMRI_root,'/utils']));  %  add NIFTI saving and loading functions;

addpath('/home/jyao3/010_MATLAB_Utils/NIfTI');
addpath('/home/jyao3/030_QSM/01_Code');

%% 1. read in data
nii = load_nii(PhasePath);
phase = nii.img;

% interpolate the phase to isotropic
imsize = size(phase);
imsize2 = [round(imsize(1:3).*vox/min(vox)), imsize(4)];
vox2 = imsize(1:3).*vox/imsize2(1:3);
interp_flag = ~isequal(imsize,imsize2);

if interp_flag
    for echo_num = 1:imsize(4)
        phase2(:,:,:,echo_num) = angle(imresize3(exp(1j*phase(:,:,:,echo_num)),imsize2(1:3)));
    end
    phase = phase2;
    clear phase2
end

if ~ exist('MagPath','var') || isempty(MagPath)
    mag = ones(imsize2);
else
    nii = load_nii(MagPath);
    mag = nii.img;
    % interpolate the mag to isotropic
    if interp_flag
        for echo_num = 1:imsize(4)
            mag2(:,:,:,echo_num) = imresize3(mag(:,:,:,echo_num),imsize2(1:3));
        end
        mag = mag2;
        clear mag2
    end
end

if ~ exist('MaskPath','var') || isempty(MaskPath)
    mask = ones(imsize2(1:3));
else
    nii = load_nii(MaskPath);
    mask = nii.img;
    % interpolate the mask to isotropic
    if interp_flag
        mask = imresize3(mask,imsize2(1:3));
    end
end

%% mkdir for output folders
if ~exist(ReconDir, 'dir')
    mkdir(ReconDir)
end

[mask, pos] = ZeroPadding(mask, 16);

%% set inference.py path; 
switch NetworkType 
    case 0
        InferencePath = [deepMRI_root, '/iQSM/PythonCodes/Evaluation/Inference.py']; 
        checkpoints = [checkpoints_root, '/iQSM_and_iQFM'];
    case 1
        InferencePath = [deepMRI_root, '/iQSM/PythonCodes/Evaluation/DataFidelityVersion/Inference.py'];
        checkpoints = [checkpoints_root, '/iQSM_iQFM_DataFidelity']; 
    case 2
        InferencePath = [deepMRI_root, '/iQSM/PythonCodes/Evaluation/LearnableLapLayer/Inference.py'];
        checkpoints = [checkpoints_root, '/iQSM_learnableKernels']; 
end 

for echo_num = 1 : imsize(4)
    
    %% 2. save all information (B0, TE, phase) as .mat file for Network Reconstruction echo by echo
    tmp_TE = TE(echo_num);
    tmp_phase = phase(:,:,:,echo_num);
    
    tmp_phase = ZeroPadding(tmp_phase, 16);
    
    mask_eroded = Save_Input(tmp_phase, mask, tmp_TE, B0, Eroded_voxel, ReconDir);
    
    [exe_command] = iQSM_pyCmd(InferencePath, [ReconDir,'/Network_Input.mat'], ReconDir, checkpoints);
    
    cmd = ['source /netopt/rhel7/versions/python/Anaconda3-5.2.0/etc/profile.d/conda.sh; ' ...
        'conda activate /working/lupolab/jingwen/conda/envs/QSM_DL;' ...
        exe_command ';' ...
        'conda deactivate'];
    system(cmd);
    
    %% load reconstruction data and save as NIFTI
    load([ReconDir,'/iQSM.mat']);
    load([ReconDir,'/iQFM.mat']);
    
    pred_chi = ZeroRemoving(pred_chi, pos);
    pred_lfs = ZeroRemoving(pred_lfs, pos);
    
    chi(:,:,:,echo_num) = pred_chi;
    lfs(:,:,:,echo_num) = pred_lfs;
    
    clear tmp_phase;
end

%% save results of all echoes before echo fitting
nii = make_nii(chi, vox2);
save_nii(nii, [ReconDir, 'iQSM_all_echoes.nii']);

nii = make_nii(lfs, vox2);
save_nii(nii, [ReconDir, 'iQFM_all_echoes.nii']);

%% magnitude weighted echo-fitting and save as NIFTI

for echo_num = 1 : imsize(4)
    chi(:,:,:,echo_num) = TE(echo_num) .* chi(:,:,:,echo_num);
    lfs(:,:,:,echo_num) = TE(echo_num) .* lfs(:,:,:,echo_num);
end

chi_fitted = echofit(chi, mag, TE);
lfs_fitted = echofit(lfs, mag, TE);

if interp_flag
    
    nii = make_nii(chi_fitted, vox2);
    save_nii(nii, [ReconDir, 'iQSM_interp_echo_fitted.nii']);
    
    nii = make_nii(lfs_fitted, vox2);
    save_nii(nii, [ReconDir, 'iQFM_interp_echo_fitted.nii']);
    
    
    % back to original resolution if anisotropic
    chi_fitted = imresize3(chi_fitted,imsize(1:3));
    lfs_fitted = imresize3(lfs_fitted,imsize(1:3));
    
end

nii = make_nii(chi_fitted, vox);
save_nii(nii, [ReconDir,'/iQSM_echo_fitted.nii']);

nii = make_nii(lfs_fitted, vox);
save_nii(nii, [ReconDir,'/iQFM_echo_fitted.nii']);


delete([ReconDir,'/Network_Input.mat']);
delete([ReconDir,'/iQFM.mat']);
delete([ReconDir,'/iQSM.mat']);



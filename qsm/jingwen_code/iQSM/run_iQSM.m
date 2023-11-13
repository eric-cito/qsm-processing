function [] = run_iQSM(output_data_path, header, TEs)

deepMRI_root = '/working/lupolab/jingwen/001_QSM/01_Code/deepMRI/';
checkpoints_root  = [deepMRI_root '/iQSM/checkpoints'];

PhasePath    = [output_data_path '/Phase.nii.gz'];  % where raw phase data is (in NIFTI format)
ReconDir     = [output_data_path '/iQSM2/'];  %% where to save reconstruction output
Eroded_voxel = 1;  %  set number of voxels for brain mask erosion; 0 means no erosion;
TE           = TEs; % set Echo Times (in second)
B0           = header.b0; % set B0 field (in Tesla)
vox          = header.voxelSize; % set voxel size a.k.a image resolution (in millimeter)
NetworkType  = 2; % network type: 0 for original iQSM, 1 for networks trained with data fidelity,
                  % 2 for networks trained with learnable Lap-Layer (15 learnable kernels) and data fidelity;

MaskPath = [output_data_path '/brain_mask_HD.nii.gz']; %% brain mask; set to one will skip brain masking
MagPath = [output_data_path '/Magni.nii.gz']; % magnitude image; set to one will skip magnitude weights in echo fitting

%% add MATLAB paths

addpath(genpath([deepMRI_root,'/iQSM/iQSM_fcns/']));  % add necessary utility function for saving data and echo-fitting;
addpath(genpath([deepMRI_root,'/utils']));  %  add NIFTI saving and loading functions;

addpath('/home/jyao3/010_MATLAB_Utils/NIfTI');

%% 1. read in data

nii = load_nii(PhasePath);
phase = nii.img;

% interpolate the phase to 1mm isotropic
imsize = size(phase);
imsize(4) = length(TE);
imsize2 = [round(imsize(1:3).*vox./[1 1 1]), imsize(4)];
vox2 = imsize(1:3).*vox./imsize2(1:3);
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
        mask = imresize3(mask,imsize2(1:3)) > 0.5;
    end
end

%% mkdir for output folders

if ~exist(ReconDir, 'dir')
    mkdir(ReconDir)
end

[mask, pos] = ZeroPadding(mask, 16);

%% set inference.py path

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

nii = load_nii(MaskPath);

nii.img = chi_fitted;
save_nii(nii, [ReconDir,'/iQSM_echo_fitted.nii']);

nii.img = lfs_fitted;
save_nii(nii, [ReconDir,'/iQFM_echo_fitted.nii']);

delete([ReconDir,'/Network_Input.mat']);
delete([ReconDir,'/iQFM.mat']);
delete([ReconDir,'/iQSM.mat']);

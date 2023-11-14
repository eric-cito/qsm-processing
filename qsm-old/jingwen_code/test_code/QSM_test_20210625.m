clear; clc;
warning('off');

%% add path

addpath(genpath('/working/lupolab/jingwen/001_QSM/01_Code'));
addpath(genpath('/working/lupolab/jingwen/011_MATLAB_Utils'));

test_data_path = '/working/lupolab/jingwen/001_QSM/Test_data/b3025_t12518';
output_data_path = '/working/lupolab/jingwen/001_QSM/Test_data/output_crop_HDBET';
mkdir(output_data_path);

%% setup parameters

fprintf('## Setting up parameters \n');

algorParam.general.isInvert                 = false;
algorParam.general.isGPU                    = true;
algorParam.general.isBET                    = false;
algorParam.general.fractional_threshold   	= 0.5;
algorParam.general.gradient_threshold   	= 0;

% brain masking: BET or HD-BET
algorParam.brainextract.betMethod           = 'HD-BET';

% unwrapping
algorParam.unwrap.echoCombMethod            = 'Optimum weights';
algorParam.unwrap.unwrapMethod              = 'Laplacian (STI suite)';
algorParam.unwrap.isEddyCorrect             = 0;
algorParam.unwrap.excludeMaskThreshold      = Inf;
algorParam.unwrap.isSaveUnwrappedEcho       = 1;
algorParam.unwrap.excludeMethod             = 'Weighting map';
algorParam.unwrap.unit                      = 'Hz';

% background field removal
algorParam.bfr.method                       = 'VSHARP (STI suite)';
algorParam.bfr.erode_radius                 = 0;
algorParam.bfr.refine                       = 0;
algorParam.bfr.refine_method                = 'polyfit';
algorParam.bfr.refine_order                 = 4;

gyro = 42.57747892;
b0 = 7.0;
epsilon = 1e-15;

%% load data

fprintf('## Loading data from input directory \n');

% copy images to output -> float
if exist([output_data_path '/Magni.nii.gz'],'file') ~= 2
    system(['fslmaths ' test_data_path '/Magni.nii.gz ' ...
        output_data_path '/Magni.nii.gz -odt float']);
    system(['fslmaths ' test_data_path '/Phase.nii.gz ' ...
        output_data_path '/Phase.nii.gz -odt float']);
end

nii = load_untouch_nii([output_data_path '/Magni.nii.gz']);
iMag = double(nii.img);
matrix_size = nii.hdr.dime.dim(2:4);
Necho = nii.hdr.dime.dim(5);
voxel_size = nii.hdr.dime.pixdim(2:4); % mm

nii = load_untouch_nii([output_data_path '/Phase.nii.gz']);
iPha = double(nii.img);

load([test_data_path '/sepia_header.mat']);
TEs = header.te;

% split to subvolumes
cd(output_data_path);
system('fslsplit Magni.nii.gz Magni -t');
mkdir([output_data_path '/Magni_subvol']);
system('mv Magni0*.nii.gz Magni_subvol');

%% brain extraction

fprintf('## Brain Extraction \n');

if exist([output_data_path '/brain_mask.nii.gz'],'file') ~= 2
    fprintf(' - Brain extraction with BET \n');
    brain_mask = zeros(size(iMag));
    for ii = 1:Necho
        brain_mask(:,:,:,ii) = BET(iMag(:,:,:,ii),matrix_size,voxel_size,...
            algorParam.general.fractional_threshold,algorParam.general.gradient_threshold);
    end
    brain_mask = double(min(brain_mask,[],4));
    
    nii = load_untouch_nii([output_data_path '/Magni_subvol/Magni0000.nii.gz']);
    nii.img = brain_mask;
    save_untouch_nii(nii, [output_data_path '/brain_mask.nii.gz']);
end

if exist([output_data_path '/brain_mask_HD.nii.gz'],'file') ~= 2
    fprintf(' - Brain extraction with HD-BET \n');
    cd(output_data_path);
    if algorParam.general.isGPU
        status = system(['bash /working/lupolab/jingwen/001_QSM/01_Code/brain_extract_HD-BET.sh' ...
            ' Magni.nii.gz 1']);
    else
        status = system(['bash /working/lupolab/jingwen/001_QSM/01_Code/brain_extract_HD-BET.sh' ...
            ' Magni.nii.gz 0']);
    end
    
    if status
        fprintf(' - !! HD-BET failed, switching to BET \n');
        algorParam.brainextract.betMethod = 'BET';
    end
end

switch algorParam.brainextract.betMethod
    case 'BET'
        nii = load_untouch_nii([output_data_path '/brain_mask.nii.gz']);
        brain_mask = double(nii.img);
        
        nii_mask_path = [output_data_path '/brain_mask.nii.gz'];
    case 'HD-BET' 
        nii = load_untouch_nii([output_data_path '/brain_mask_HD.nii.gz']);
        brain_mask = double(nii.img);
        
        nii_mask_path = [output_data_path '/brain_mask_HD.nii.gz'];
    otherwise
        disp('! Brain extraction method not specified correctly.');
end

%% Crop the images for faster computation

fprintf('## Cropping images \n');

[x, y, z] = ind2sub(size(brain_mask), find(brain_mask > 0));
gap = 5; % set a magin of 5 pixels
x1 = max(1, min(x)-gap); x2 = min(size(brain_mask,1), max(x)+gap);
y1 = max(1, min(y)-gap); y2 = min(size(brain_mask,2), max(y)+gap);
z1 = max(1, min(z)-gap); z2 = min(size(brain_mask,3), max(z)+gap);
if mod(x2 - x1, 2) == 0
    x2 = x2 + 1;
end
if mod(y2 - y1, 2) == 0
    y2 = y2 + 1;
end
if mod(z2 - z1, 2) == 0
    z2 = z2 + 1;
end

crop.X = x1:x2;
crop.Y = y1:y2;
crop.Z = z1:z2;

% crop the images
brain_mask_ori = brain_mask;
brain_mask = brain_mask_ori(crop.X, crop.Y, crop.Z);

iMag_ori = iMag;
iMag = iMag_ori(crop.X, crop.Y, crop.Z, :);
iPha_ori = iPha;
iPha = iPha_ori(crop.X, crop.Y, crop.Z, :);

matrix_size_ori = matrix_size;
matrix_size = size(brain_mask);

header_ori = header;
header.matrixSize = matrix_size;

header.magn = iMag;
header.phase = iPha;

clear x y z x1 x2 y1 y2 z1 z2

%% echo combination and phase unwrapping + background field removal

% [totalField, noiseStd, fieldmapUnwrapAllEchoes] = ...
%     estimateTotalField(iPha, brain_mask, matrix_size, voxel_size, algorParam, header);
%
% % save maps
% nii = load_untouch_nii([output_data_path '/brain_mask.nii.gz']);
% nii.img = totalField;
% save_untouch_nii(nii, [output_data_path '/totalField.nii.gz']);
% nii.img = noiseStd;
% save_untouch_nii(nii, [output_data_path '/noiseStd.nii.gz']);
% nii = load_untouch_nii([output_data_path '/Magni.nii.gz']);
% nii.img = fieldmapUnwrapAllEchoes;
% save_untouch_nii(nii, [output_data_path '/Phase_unwrap.nii.gz']);
%
% localField = BackgroundRemovalMacro(totalField, brain_mask, matrix_size, voxel_size, ...
%              algorParam, header);
%
% % save maps
% nii = load_untouch_nii([output_data_path '/brain_mask.nii.gz']);
% nii.img = localField;
% save_untouch_nii(nii, [output_data_path '/localField.nii.gz']);

%% alternative of averaging accross all echoes

fprintf('## Phase unwrapping and background field removing \n');

if exist([output_data_path '/localField_meanEcho.nii.gz'],'file') ~= 2
    
    totalField_echo = zeros(size(iMag));
    localField_echo = zeros(size(iMag));
    totalField_echoNorm = zeros(size(iMag));
    localField_echoNorm = zeros(size(iMag));
    for ii = 1:Necho
        totalField_echo(:,:,:,ii) = ...
            UnwrapPhaseMacro(iPha(:,:,:,ii), brain_mask, matrix_size, ...
            voxel_size, algorParam, header);
        localField_echo(:,:,:,ii) = ...
            BackgroundRemovalMacro(totalField_echo(:,:,:,ii), brain_mask, matrix_size, ...
            voxel_size, algorParam, header);
        totalField_echoNorm(:,:,:,ii) = totalField_echo(:,:,:,ii)/header.te(ii)*TEs(end);
        localField_echoNorm(:,:,:,ii) = localField_echo(:,:,:,ii)/header.te(ii)*TEs(end);
    end
    
    totalField_mean = mean(totalField_echoNorm,4);
    localField_mean = mean(localField_echoNorm,4);
    
    totalField_output = zeros(matrix_size_ori);
    localField_output = zeros(matrix_size_ori);
    totalField_output(crop.X,crop.Y,crop.Z) = totalField_mean;
    localField_output(crop.X,crop.Y,crop.Z) = localField_mean;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = totalField_output;
    save_untouch_nii(nii, [output_data_path '/totalField_meanEcho.nii.gz']);
    nii.img = localField_output;
    save_untouch_nii(nii, [output_data_path '/localField_meanEcho.nii.gz']);
    
else
    
    nii = load_untouch_nii([output_data_path '/totalField_meanEcho.nii.gz']);
    totalField_output = double(nii.img);
    nii = load_untouch_nii([output_data_path '/localField_meanEcho.nii.gz']);
    localField_output = double(nii.img);
    
    totalField_mean = totalField_output(crop.X,crop.Y,crop.Z);
    localField_mean = localField_output(crop.X,crop.Y,crop.Z);
    
end

% calculate weight
tmp = sqrt(mean(iMag.^2,4));
wmap = (tmp./max(tmp(:))) .* brain_mask;
header.weights = wmap;

% one TE for all
header.te = TEs(end);

% save cropped images
if exist([output_data_path '/localField_meanEcho_crop.nii.gz'],'file') ~= 2
    cmd = sprintf('fslroi %s %s %i %i %i %i %i %i', ...
        [output_data_path '/localField_meanEcho.nii.gz'], ...
        [output_data_path '/localField_meanEcho_crop.nii.gz'], ...
        crop.X(1)-1, length(crop.X), ...
        crop.Y(1)-1, length(crop.Y), ...
        crop.Z(1)-1, length(crop.Z));
    system(cmd);
    
    % check dimension
    nii = load_untouch_nii([output_data_path '/localField_meanEcho_crop.nii.gz']);
    temp = double(nii.img);
    fprintf('Size: %i x %i x %i \n', size(temp));
end

%% QSM - iLSQR

if exist([output_data_path '/QSM_iLSQR_meanEcho.nii.gz'],'file') ~= 2
    
    algorParam.qsm.method                       = 'STI suite iLSQR';
    algorParam.qsm.maxiter                      = 100;
    algorParam.qsm.tol1                         = 0.01;
    algorParam.qsm.tol2                         = 0.001;
    algorParam.qsm.padsize                      = [0,0,6]; % Eason
    algorParam.qsm.reference_tissue             = 'None';
    
    QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_iLSQR_meanEcho.nii.gz']);
    
end

%% QSM - MEDI

if exist([output_data_path '/QSM_MEDI_meanEcho.nii.gz'],'file') ~= 2
    
    algorParam.qsm.method                       = 'MEDI';
    algorParam.qsm.isSMV                        = false;
    algorParam.qsm.merit                        = false;
    algorParam.qsm.lambda                       = 2000; % Eason
    
    QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_MEDI_meanEcho.nii.gz']);
    
end

%% QSM - FANSI

if exist([output_data_path '/QSM_FANSI_nonlinearTV_meanEcho.nii.gz'],'file') ~= 2
    
    algorParam.qsm.method                       = 'FANSI';
    algorParam.qsm.solver                       = 'nonlinear';
    algorParam.qsm.maxiter                      = 50;
    algorParam.qsm.lambda                       = 3e-5; % 3e-5; % 0.0015;
    algorParam.qsm.mu1                          = 5e-5; % 5e-5
    header.weights                              = iMag(:,:,:,1).*brain_mask;
    
    QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_FANSI_nonlinearTV_meanEcho.nii.gz']);
    
    header.weights = wmap;
    
end

%% QSM - Star-QSM

if exist([output_data_path '/QSM_STARQSM_meanEcho.nii.gz'],'file') ~= 2
    
    algorParam.qsm.method                       = 'Star-QSM';
    
    QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_STARQSM_meanEcho.nii.gz']);
    
end

%% QSM - HD-QSM

if exist([output_data_path '/QSM_HDQSM_meanEcho.nii.gz'],'file') ~= 2
    
    algorParam.qsm.method                       = 'HD-QSM';
    algorParam.qsm.tol_update       = 1.0;
    algorParam.qsm.maxOuterIterL2   = 280;
    algorParam.qsm.mu1L2            = 100 * 10^-4.785;
    algorParam.qsm.alphaL2          = 10^-4.785;
    algorParam.qsm.maxOuterIterL1   = 20;
    
    QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_HDQSM_meanEcho.nii.gz']);
    
end

%% QSM - QSMGAN

if exist([output_data_path '/QSM_QSMGAN_meanEcho.nii.gz'],'file') ~= 2
    
    % prepare local field map for QSMGAN
    
    nii = load_untouch_nii([output_data_path '/localField_meanEcho_crop.nii.gz']);
    temp = double(nii.img)/(b0*gyro*header.te);
    temp(temp == 0) = epsilon;
    temp(isnan(temp)) = epsilon;
    temp(isinf(temp)) = epsilon;
    nii.img = temp;
    save_untouch_nii(nii, [output_data_path '/localField_meanEcho_crop_ppm.nii.gz']);
    
    cd('/working/lupolab/jingwen/001_QSM/01_Code/QSMGAN/Code');
    
    cmd = ['source /netopt/rhel7/versions/python/Anaconda3-5.2.0/etc/profile.d/conda.sh; ' ...
        'conda activate QSM_DL; '...
        'python make_swan_qsm_DL_JY.py ' output_data_path '/localField_meanEcho_crop_ppm.nii.gz; '...
        'conda deactivate; '];
    
    system(cmd);
    
    cd('/working/lupolab/jingwen');
    
    nii = load_untouch_nii([output_data_path '/localField_meanEcho_crop_ppm_susc_DL.nii.gz']);
    QSM = double(nii.img);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_QSMGAN_meanEcho.nii.gz']);
    
end

%% QSM - QSMnet+ 

if exist([output_data_path '/QSM_QSMnet_meanEcho.nii.gz'],'file') ~= 2
    
    % create isotropic tissue phase and mask
    FOV = header.voxelSize .* header.matrixSize;
    matrixSize_iso = ceil(FOV);
    
    % convert to even dimensions though it will not be precisely isotropic
    if mod(matrixSize_iso(1), 2) == 1
        matrixSize_iso(1) = matrixSize_iso(1) + 1;
    end
    if mod(matrixSize_iso(2), 2) == 1
        matrixSize_iso(2) = matrixSize_iso(2) + 1;
    end
    if mod(matrixSize_iso(3), 2) == 1
        matrixSize_iso(3) = matrixSize_iso(3) + 1;
    end
    
    temp = localField_mean/(b0*gyro*header.te); % Hz -> ppm
    localField_iso = imresize3(temp,matrixSize_iso);
    
    mask_iso = imresize3(brain_mask,matrixSize_iso) > 0.5;
    
    localField_iso = -flip(localField_iso,2); % not sure why but needed the negative sign
    mask_iso = flip(mask_iso,2);
    
    % check orientation
    % load('/working/lupolab/jingwen/001_QSM/01_Code/QSMnet-master/Data/Test/Input/test_input1.mat');
    %
    % figure('Position',[100,100,900,400]);
    % subplot(1,2,1); plot3D(phs_tissue);
    % subplot(1,2,2); plot3D(localField_iso);
    
    localField_iso(localField_iso == 0) = epsilon;
    localField_iso(isnan(localField_iso)) = epsilon;
    localField_iso(isinf(localField_iso)) = epsilon;
    
    phs_tissue = localField_iso;
    mask = mask_iso;
    save([output_data_path '/QSMnet_input.mat'], 'phs_tissue', 'mask');
    
    cd('/working/lupolab/jingwen/001_QSM/01_Code/QSMnet-master/Code');
    
    cmd = ['source /netopt/rhel7/versions/python/Anaconda3-5.2.0/etc/profile.d/conda.sh; ' ...
        'conda activate QSM_DL; '...
        'python inference.py ' output_data_path '; '...
        'conda deactivate; '];
    
    system(cmd);
    
    % check results
    
    % figure('Position',[100,100,900,400]);
    % load('Data/Test/Prediction/subject1_QSMnet+_64_25.mat');
    % subplot(1,2,1); plot3D(sus);
    % load('Data/Test/Prediction/subject99_QSMnet+_64_25.mat');
    % subplot(1,2,2); plot3D(sus); caxis([-0.03 0.03])
    
    load([output_data_path '/QSMnet_QSMnet+_64_25.mat']);
    QSM = flip(imresize3(sus,header.matrixSize),2);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_QSMnet_meanEcho.nii.gz']);
    
    cd('/working/lupolab/jingwen');
    
end

%% QSM - SS_TGV

if exist([output_data_path '/QSM_SSTGV_meanEcho.nii.gz'],'file') ~= 2
    
    algorParam.qsm.method                       = 'SS-TGV';
    algorParam.qsm.SS.method                    = 'SS-TGV';
    algorParam.qsm.SSTGV.alpha1                 = 7e-3;
    algorParam.qsm.SSTGV.mu0                    = 1e-1;
    
    QSM = -QSMMacro(totalField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_SSTGV_meanEcho.nii.gz']);
    
end

%% QSM - SS_TV

if exist([output_data_path '/QSM_SSTV_meanEcho.nii.gz'],'file') ~= 2
    
    algorParam.qsm.method                       = 'SS-TGV';
    algorParam.qsm.SS.method                    = 'SS-TV';
    algorParam.qsm.SSTV.alpha                   = 7e-3;
    algorParam.qsm.SSTV.mu1                     = 1e-1;
    
    QSM = -QSMMacro(totalField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_SSTV_meanEcho.nii.gz']);
    
end

%% QSM - QSIP

if exist([output_data_path '/QSM_QSIP_meanEcho.nii.gz'],'file') ~= 2
    
    tmp = sqrt(mean(iMag.^2,4));
    header.magn = tmp;
    
    algorParam.qsm.method                       = 'QSIP';
    algorParam.qsm.QSIP.atlas_thr               = 0.96;
    algorParam.qsm.QSIP.atlas_flag              = 1;
    algorParam.qsm.QSIP.num_iter                = 300;
    
    QSM = QSMMacro(totalField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_QSIP_meanEcho.nii.gz']);
    
    header.magn = iMag;
    
end

%%



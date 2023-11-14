clc; clear;
warning('off');

%% add path

addpath(genpath('/working/lupolab/jingwen/001_QSM/01_Code'));
addpath(genpath('/working/lupolab/jingwen/011_MATLAB_Utils'));

combine_opt = 'combQSM';

%% set up parameters

input_data_path = '/working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308';
output_data_path = [input_data_path '/output_HDBET_' combine_opt];
mkdir(output_data_path);

copyfile([input_data_path '/output_HDBET/Magni.nii.gz'], output_data_path);
copyfile([input_data_path '/output_HDBET/Phase.nii.gz'], output_data_path);
copyfile([input_data_path '/output_HDBET/brain_mask_HD.nii.gz'], output_data_path);
copyfile([input_data_path '/output_HDBET/brain_mask.nii.gz'], output_data_path);

input.magnitudeFile     = [input_data_path '/Magni.nii.gz'];
input.phaseFile         = [input_data_path '/Phase.nii.gz'];
input.headerFile        = [input_data_path '/sepia_header.mat'];

opts.writeLog           = 1;
opts.isGPU              = 1;
opts.BETmethod          = 'HD-BET';
opts.iLSQR              = 1;
opts.QSMGAN             = 1;
opts.All                = 0;

fprintf('## Setting up parameters \n');

algorParam.general.isInvert                 = false;
algorParam.general.isGPU                    = opts.isGPU;
algorParam.general.isBET                    = false;
algorParam.general.fractional_threshold   	= 0.5;
algorParam.general.gradient_threshold   	= 0;

% brain masking: BET or HD-BET
algorParam.brainextract.betMethod           = opts.BETmethod;

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

%% load data

fprintf('## Loading data from input directory \n');

nii = load_untouch_nii([output_data_path '/Magni.nii.gz']);
iMag = double(nii.img);
matrix_size = nii.hdr.dime.dim(2:4);
Necho = nii.hdr.dime.dim(5);
voxel_size = nii.hdr.dime.pixdim(2:4); % mm

nii = load_untouch_nii([output_data_path '/Phase.nii.gz']);
iPha = double(nii.img);

load(input.headerFile);
TEs = header.te;
deltaTE = header.delta_TE;
gyro = 42.57747892;
b0 = header.b0;
epsilon = 1e-15;

cd(output_data_path);

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

if exist([output_data_path '/localField.nii.gz'],'file') ~= 2
    
    switch combine_opt
        case 'combPhase'
            [totalField, noiseStd, fieldmapUnwrapAllEchoes] = ...
                estimateTotalField(iPha, brain_mask, matrix_size, voxel_size, algorParam, header);
            % save maps
            nii = load_untouch_nii([output_data_path '/brain_mask.nii.gz']);
            IMG_output = zeros(matrix_size_ori);
            IMG_output(crop.X,crop.Y,crop.Z) = totalField;
            nii.img = IMG_output;
            save_untouch_nii(nii, [output_data_path '/totalField.nii.gz']);
            IMG_output = zeros(matrix_size_ori);
            IMG_output(crop.X,crop.Y,crop.Z) = noiseStd;
            nii.img = IMG_output;
            save_untouch_nii(nii, [output_data_path '/noiseStd.nii.gz']);
            
            localField = BackgroundRemovalMacro(totalField, brain_mask, matrix_size, voxel_size, ...
                algorParam, header);
            % save maps
            IMG_output = zeros(matrix_size_ori);
            IMG_output(crop.X,crop.Y,crop.Z) = localField;
            nii.img = IMG_output;
            save_untouch_nii(nii, [output_data_path '/localField.nii.gz']);
        case 'combQSM'
            for ii = 1:length(TEs)
                header.te = TEs(ii);
                header.delta_TE = TEs(ii);
                header.magn = iMag(:,:,:,ii);
                header.phase = iPha(:,:,:,ii);
                [totalField(:,:,:,ii), noiseStd(:,:,:,ii), ~] = ...
                    estimateTotalField(iPha(:,:,:,ii), brain_mask, matrix_size, voxel_size, algorParam, header);
                localField(:,:,:,ii) = ...
                    BackgroundRemovalMacro(totalField(:,:,:,ii), brain_mask, matrix_size, voxel_size, ...
                    algorParam, header);
            end
            % save maps
            nii = load_untouch_nii([output_data_path '/Magni.nii.gz']);
            temp_size = size(nii.img);
            IMG_output = zeros(temp_size);
            IMG_output(crop.X,crop.Y,crop.Z,:) = totalField;
            nii.img = IMG_output;
            save_untouch_nii(nii, [output_data_path '/totalField.nii.gz']);
            IMG_output = zeros(temp_size);
            IMG_output(crop.X,crop.Y,crop.Z,:) = noiseStd;
            nii.img = IMG_output;
            save_untouch_nii(nii, [output_data_path '/noiseStd.nii.gz']);
            IMG_output = zeros(temp_size);
            IMG_output(crop.X,crop.Y,crop.Z,:) = localField;
            nii.img = IMG_output;
            save_untouch_nii(nii, [output_data_path '/localField.nii.gz']);
    end

else
    
    nii = load_untouch_nii([output_data_path '/totalField.nii.gz']);
    totalField_output = double(nii.img);
    nii = load_untouch_nii([output_data_path '/localField.nii.gz']);
    localField_output = double(nii.img);
    
    totalField = totalField_output(crop.X,crop.Y,crop.Z,:);
    localField = localField_output(crop.X,crop.Y,crop.Z,:);
    
end

switch combine_opt
    case 'combPhase'
        % one TE for all
        header.te = 1; % s
        header.delta_TE = 1; % s
end

% calculate weight
tmp = sqrt(mean(iMag.^2,4));
wmap = (tmp./max(tmp(:))) .* brain_mask;
header.weights = wmap;

% save cropped images
if exist([output_data_path '/localField_crop.nii.gz'],'file') ~= 2
    cmd = sprintf('fslroi %s %s %i %i %i %i %i %i', ...
        [output_data_path '/localField.nii.gz'], ...
        [output_data_path '/localField_crop.nii.gz'], ...
        crop.X(1)-1, length(crop.X), ...
        crop.Y(1)-1, length(crop.Y), ...
        crop.Z(1)-1, length(crop.Z));
    system(cmd);
    
    % check dimension
    nii = load_untouch_nii([output_data_path '/localField_crop.nii.gz']);
    temp = double(nii.img);
    switch combine_opt
        case 'combPhase'
            fprintf('Size: %i x %i x %i \n', size(temp));
        case 'combQSM'
            fprintf('Size: %i x %i x %i x %i \n', size(temp));
    end
    
end

%% QSM - iLSQR

if opts.iLSQR && exist([output_data_path '/QSM_iLSQR_' combine_opt '.nii.gz'],'file') ~= 2
    
    tic;
    
    algorParam.qsm.method                       = 'STI suite iLSQR';
    algorParam.qsm.maxiter                      = 100;
    algorParam.qsm.tol1                         = 0.01;
    algorParam.qsm.tol2                         = 0.001;
    algorParam.qsm.padsize                      = [0,0,6]; % Eason
    algorParam.qsm.reference_tissue             = 'None';
    
    switch combine_opt
        case 'combPhase'
            QSM = -QSMMacro(localField, brain_mask, matrix_size, voxel_size, algorParam, header);
        case 'combQSM'
            for ii = 1:length(TEs)
                header.te = TEs(ii);
                header.delta_TE = TEs(ii);
                header.magn = iMag(:,:,:,ii);
                header.phase = iPha(:,:,:,ii);
                QSM(:,:,:,ii) = -QSMMacro(localField(:,:,:,ii), brain_mask, matrix_size, voxel_size, algorParam, header);
            end
            QSM = mean(QSM,4);
    end
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM;
    
    Tcomp.iLSQR = toc;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_iLSQR_' combine_opt '.nii.gz']);
    
end

%% QSM - QSMGAN

if opts.QSMGAN && exist([output_data_path '/QSM_QSMGAN_' combine_opt '.nii.gz'],'file') ~= 2
    
    tic;
    
    % prepare local field map for QSMGAN
    
    nii = load_untouch_nii([output_data_path '/localField_crop.nii.gz']);
    switch combine_opt
        case 'combPhase'
            temp = double(nii.img)/(b0*gyro*header.te);
        case 'combQSM'
            temp = double(nii.img);
            for ii = 1:length(TEs)
                temp(:,:,:,ii) = double(nii.img(:,:,:,ii))/(b0*gyro*TEs(ii));
            end
    end
    temp(temp == 0) = epsilon;
    temp(isnan(temp)) = epsilon;
    temp(isinf(temp)) = epsilon;
    nii.img = temp;
    save_untouch_nii(nii, [output_data_path '/localField_crop_ppm.nii.gz']);
    
    switch combine_opt
        case 'combPhase'
            cd('/working/lupolab/jingwen/001_QSM/01_Code/QSMGAN/Code');
            cmd = ['source /netopt/rhel7/versions/python/Anaconda3-5.2.0/etc/profile.d/conda.sh; ' ...
                'conda activate QSM_DL; '...
                'python make_swan_qsm_DL_JY.py ' output_data_path '/localField_crop_ppm.nii.gz; '...
                'conda deactivate; '];
            system(cmd);
            cd('/working/lupolab/jingwen');
            nii = load_untouch_nii([output_data_path '/localField_crop_ppm_susc_DL.nii.gz']);
            QSM = double(nii.img);
        case 'combQSM'
            cd(output_data_path);
            if exist([output_data_path '/ppm_subvol/ppm0000.nii.gz'],'file') ~= 2
                system('fslsplit localField_crop_ppm.nii.gz ppm -t');
                mkdir([output_data_path '/ppm_subvol']);
                system('mv ppm000*.nii.gz ppm_subvol');
            end
            for ii = 1:length(TEs)
                cd('/working/lupolab/jingwen/001_QSM/01_Code/QSMGAN/Code');
                cmd = ['source /netopt/rhel7/versions/python/Anaconda3-5.2.0/etc/profile.d/conda.sh; ' ...
                    'conda activate QSM_DL; '...
                    'python make_swan_qsm_DL_JY.py ' output_data_path '/ppm_subvol/ppm000' num2str(ii-1) '.nii.gz; '...
                    'conda deactivate; '];
                system(cmd);
                cd('/working/lupolab/jingwen');
                nii = load_untouch_nii([output_data_path '/ppm_subvol/ppm000' num2str(ii-1) '_susc_DL.nii.gz']);
                QSM(:,:,:,ii) = double(nii.img);
            end
            QSM = mean(QSM,4);
    end
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM;
    
    Tcomp.QSMGAN = toc;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_QSMGAN_' combine_opt '.nii.gz']);
    
end
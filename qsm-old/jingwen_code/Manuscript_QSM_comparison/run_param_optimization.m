function [] = run_param_optimization(outputPath)

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

%% subject for optimization

% outputPath = '/working/lupolab/jingwen/001_QSM/Test_data/b4469_t12309/';
% '/working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308/';
% b4469_t12309
output_data_path = outputPath;

%% load data

fprintf('## Loading data from input directory \n');

nii = load_untouch_nii([output_data_path '/Magni.nii.gz']);
iMag = double(nii.img);
matrix_size = nii.hdr.dime.dim(2:4);
Necho = nii.hdr.dime.dim(5);
voxel_size = nii.hdr.dime.pixdim(2:4); % mm

nii = load_untouch_nii([output_data_path '/Phase.nii.gz']);
iPha = double(nii.img);

load([output_data_path '/../sepia_header.mat']);
TEs = header.te;
deltaTE = header.delta_TE;
gyro = 42.57747892;
b0 = header.b0;
epsilon = 1e-15;

cd(output_data_path);

nii = load_untouch_nii([output_data_path '/brain_mask_HD.nii.gz']);
brain_mask = double(nii.img);

nii_mask_path = [output_data_path '/brain_mask_HD.nii.gz'];

%% Crop the images for faster computation

fprintf('## Cropping images \n');

[x, y, z] = ind2sub(size(brain_mask), find(brain_mask > 0));
gap = 5; % set a magin of 5 pixels
x1 = max(1, min(x)-gap); x2 = min(size(brain_mask,1), max(x)+gap);
y1 = max(1, min(y)-gap); y2 = min(size(brain_mask,2), max(y)+gap);
z1 = max(1, min(z)-gap); z2 = min(size(brain_mask,3), max(z)+gap);
if mod(x2 - x1, 2) == 0
    x2 = x2 - 1;
end
if mod(y2 - y1, 2) == 0
    y2 = y2 - 1;
end
if mod(z2 - z1, 2) == 0
    z2 = z2 - 1;
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

%% local field

nii = load_untouch_nii([output_data_path '/totalField_meanEcho.nii.gz']);
totalField_output = double(nii.img);
nii = load_untouch_nii([output_data_path '/localField_meanEcho.nii.gz']);
localField_output = double(nii.img);

totalField_mean = totalField_output(crop.X,crop.Y,crop.Z);
localField_mean = localField_output(crop.X,crop.Y,crop.Z);

% calculate weight
tmp = sqrt(mean(iMag.^2,4));
wmap = (tmp./max(tmp(:))) .* brain_mask;
header.weights = wmap;

% one TE for all
header.te = TEs(end); % s
header.delta_TE = TEs(end); % s

%% QSM - MEDI

output_data_path = [outputPath '/MEDI'];
mkdir(output_data_path);

algorParam.qsm.method                       = 'MEDI';
algorParam.qsm.isSMV                        = false;
algorParam.qsm.merit                        = false;
algorParam.qsm.lambda                       = 1000; % Eason 2000

for lambda  = 10.^[1:0.25:3 3.1:0.1:5]
    
    paraStr = sprintf('l10E%.2f', log10(lambda));
    paraStr = strrep(paraStr,'.','p');
    
    if exist([output_data_path '/QSM_MEDI_' paraStr '.nii.gz'],'file') == 0
        
        algorParam.qsm.lambda                       = lambda;
        
        QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
        
        QSM_output = zeros(matrix_size_ori);
        QSM_output(crop.X,crop.Y,crop.Z) = QSM;
        
        % save maps
        nii = load_untouch_nii(nii_mask_path);
        nii.img = QSM_output;
        
        save_untouch_nii(nii, [output_data_path '/QSM_MEDI_' paraStr '.nii.gz']);
    end
    
end

%% QSM - FANSI

output_data_path = [outputPath '/FANSI'];
mkdir(output_data_path);

algorParam.qsm.method                       = 'FANSI';
algorParam.qsm.solver                       = 'nonlinear';
algorParam.qsm.maxiter                      = 150;
header.weights                              = iMag(:,:,:,1).*brain_mask;

for lambda  = 10.^[-5:0.25:-4.5 -4.3:0.1:-3.5 -3.25:0.25:-2]
    
    paraStr = sprintf('l10E%.2f', log10(lambda));
    paraStr = strrep(paraStr,'.','p');
    
    if exist([output_data_path '/QSM_FANSI_' paraStr '.nii.gz'],'file') == 0
        
        mu = 100*lambda;
        algorParam.qsm.lambda                   = lambda;
        algorParam.qsm.mu1                      = mu;
        
        QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
        
        QSM_output = zeros(matrix_size_ori);
        QSM_output(crop.X,crop.Y,crop.Z) = QSM;
        
        % save maps
        nii = load_untouch_nii(nii_mask_path);
        nii.img = QSM_output;
        
        save_untouch_nii(nii, [output_data_path '/QSM_FANSI_' paraStr '.nii.gz']);
    end
    
end

header.weights = wmap;

%% QSM - HD-QSM

output_data_path = [outputPath '/HDQSM'];
mkdir(output_data_path);

algorParam.qsm.method           = 'HD-QSM';
algorParam.qsm.tol_update       = 1.0;
algorParam.qsm.maxOuterIterL2   = 280;
algorParam.qsm.maxOuterIterL1   = 20;

for lambda  = 10.^[-6:0.25:-5.25 -5.2:0.1:-3.8 -3.75:0.25:-2]
    
    paraStr = sprintf('l10E%.2f', log10(lambda));
    paraStr = strrep(paraStr,'.','p');
    
    if exist([output_data_path '/QSM_HDQSM_' paraStr '.nii.gz'],'file') == 0
        
        algorParam.qsm.alphaL2          = lambda;
        algorParam.qsm.mu1L2            = 100 * lambda;
        
        QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
        
        QSM_output = zeros(matrix_size_ori);
        QSM_output(crop.X,crop.Y,crop.Z) = QSM;
        
        % save maps
        nii = load_untouch_nii(nii_mask_path);
        nii.img = QSM_output;
        
        save_untouch_nii(nii, [output_data_path '/QSM_HDQSM_' paraStr '.nii.gz']);
    end
    
end

%% QSM - SS_TGV

output_data_path = [outputPath '/SSTGV'];
mkdir(output_data_path);

algorParam.qsm.method                       = 'SS-TGV';
algorParam.qsm.SS.method                    = 'SS-TGV';
algorParam.qsm.SSTGV.mu0                    = 1e-1;

for lambda  = 10.^[-3:0.1:-1.3]
    
    paraStr = sprintf('l10E%.2f', log10(lambda));
    paraStr = strrep(paraStr,'.','p');
    
    if exist([output_data_path '/QSM_SSTGV_' paraStr '.nii.gz'],'file') == 0
        
        algorParam.qsm.SSTGV.alpha0                 = lambda;
        
        QSM = -QSMMacro(totalField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
        
        QSM_output = zeros(matrix_size_ori);
        QSM_output(crop.X,crop.Y,crop.Z) = QSM;
        
        % save maps
        nii = load_untouch_nii(nii_mask_path);
        nii.img = QSM_output;
        
        save_untouch_nii(nii, [output_data_path '/QSM_SSTGV_' paraStr '.nii.gz']);
    end
    
end

%% QSM - SS_TV

output_data_path = [outputPath '/SSTV'];
mkdir(output_data_path);

algorParam.qsm.method                       = 'SS-TGV';
algorParam.qsm.SS.method                    = 'SS-TV';
algorParam.qsm.SSTV.mu1                     = 1e-1;

for lambda  = 10.^[-3:0.1:-1.3]
    
    paraStr = sprintf('l10E%.2f', log10(lambda));
    paraStr = strrep(paraStr,'.','p');
    
    if exist([output_data_path '/QSM_SSTV_' paraStr '.nii.gz'],'file') == 0
        
        algorParam.qsm.SSTV.alpha                 = lambda;
        
        QSM = -QSMMacro(totalField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
        
        QSM_output = zeros(matrix_size_ori);
        QSM_output(crop.X,crop.Y,crop.Z) = QSM;
        
        % save maps
        nii = load_untouch_nii(nii_mask_path);
        nii.img = QSM_output;
        
        save_untouch_nii(nii, [output_data_path '/QSM_SSTV_' paraStr '.nii.gz']);
    end
    
end

end

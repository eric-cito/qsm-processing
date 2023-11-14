function [Tcomp] = QSM_processing_CPU(input, output_data_path, opts)

% =========================================================================
% Input
% --------------
% input_data        : structure containing input file info
%       - input_data.magnitudeFile  
%               file path to magnitude nifti, default: pwd/Magni.nii.gz
%       - input_data.phaseFile  
%               file path to phase nifti, default: pwd/Phase.nii.gz
%       - input_data.headerFile  
%               file path to header mat file, default: pwd/sepia_header.mat
%               - structure header with the following fields
%               - header.voxelSize      % mm 3D
%               - header.matrixSize     % matrix size 3D
%               - header.b0             % Tesla
%               - header.b0dir          % [0;0;1]
%               - header.CF             % central frequency Hz
%               - header.te             % s
%               - header.delta_TE       % s
% output_data_path  : file path to output data
%     	- default: pwd/output
% opts              : structure containing fields with options
%       - opts.writeLog         % creates log file in the output path, default 1
%       - opts.isGPU            % availability of GPU, default 1
%       - opts.Necho            % number of echoes, default 4
%       - opts.BETmethod        % default 'BET'
%                               % recommend 'HD-BET' with GPU
%       - opts.iLSQR            % flag to run, default 1
%       - opts.QSMGAN           % flag to run, default 1
%       - opts.All              % flag to run all algorithms, default 0
%
% Output
% --------------
% Tcomp             : structure containing computation times
%
% Description: This is a wrapper function to run QSM processing
%
% Jingwen Yao
% jingwen.yao@ucsf.ucsf.edu
% Date created: 1 July 2021
% Date modified: 1 July 2021
% =========================================================================

%% add path

warning('off');
addpath(genpath('/data/morrison/scripts/qsm/jingwen_code'));

%% set default arguments

currentPath = pwd;

% if nargin < 3 || isempty(opts)
%     opts = struct;
% end
% 
% if nargin < 2 || isempty(output_data_path)
%     output_data_path = sprintf('%s/output', currentPath);
% end
% 
% if nargin < 1 || isempty(input)
%     input = struct;
% end

input = check_and_set_default_input(input, currentPath);

opts = check_and_set_default_opts(opts);

if exist(input.headerFile, 'file') ~= 2
    set_SEPIA_header;
    input.headerFile = [pwd '/sepia_header.mat'];
end
load(input.headerFile);

if opts.writeLog
    mkdir(sprintf('%s/Log/', output_data_path));
    logfile = sprintf('%s/Log/Log_%s', output_data_path, datestr(now,'yyyy-mm-dd-HHMM'));
    diary(logfile);
end

Tcomp = struct;

%% add path

warning('off');

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils'));

mkdir(output_data_path);
mkdir([output_data_path '/Log']);

%% setup parameters

fprintf('## Setting up parameters \n');

algorParam.general.isInvert                 = false;
algorParam.general.isGPU                    = opts.isGPU;
algorParam.general.isBET                    = true;
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

gyro = 42.57747892;
b0 = header.b0;
epsilon = 1e-15;

%% load data

fprintf('## Loading data from input directory \n');

%copy images to output -> float
if exist([output_data_path '/Magni.nii.gz'],'file') ~= 2
    system(sprintf('fslmaths %s %s/Magni.nii.gz -odt float', ...
        input.magnitudeFile, output_data_path));
    system(sprintf('fslmaths %s %s/Phase.nii.gz -odt float', ...
        input.phaseFile, output_data_path));
end

% system(['cp ' input.magnitudeFile ' ' output_data_path '/Magni.nii.gz']);
% system(['cp ' input.phaseFile ' ' output_data_path '/Phase.nii.gz']);


nii = load_untouch_nii([output_data_path 'Magni.nii.gz']);
iMag = double(nii.img);
matrix_size = nii.hdr.dime.dim(2:4);
Necho = nii.hdr.dime.dim(5);
voxel_size = nii.hdr.dime.pixdim(2:4); % mm

nii = load_untouch_nii([output_data_path 'Phase.nii.gz']);
iPha = double(nii.img);

% rescale phase to -pi to pi
maxPha = max(iPha(:));
minPha = min(iPha(:));
iPha = (iPha + minPha)/(maxPha - minPha)*2*pi - pi;

TEs = header.te;
deltaTE = header.delta_TE;

% split to subvolumes
cd(output_data_path);
if exist([output_data_path 'Magni_subvol/Magni0000.nii.gz'],'file') ~= 2
    system('fslsplit Magni.nii.gz Magni -t');
    mkdir([output_data_path 'Magni_subvol']);
    system('mv Magni0*.nii.gz Magni_subvol');
end


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
% 
% if exist([output_data_path '/brain_mask_HD.nii.gz'],'file') ~= 2
%     fprintf(' - Brain extraction with HD-BET \n');
%     cd(output_data_path);
%     if algorParam.general.isGPU
%         status = system(['bash /home/jyao3/030_QSM/01_Code/brain_extract_HD-BET.sh' ...
%             ' Magni.nii.gz 1']);
%     else
%         status = system(['bash /home/jyao3/030_QSM/01_Code/brain_extract_HD-BET.sh' ...
%             ' Magni.nii.gz 0']);
%     end
%     
%     if status
%         fprintf(' - !! HD-BET failed, switching to BET \n');
%         algorParam.brainextract.betMethod = 'BET';
%     end
% end

% switch algorParam.brainextract.betMethod
%     case 'BET'
        nii = load_untouch_nii([output_data_path '/brain_mask.nii.gz']);
        brain_mask = double(nii.img);
        
        nii_mask_path = [output_data_path '/brain_mask.nii.gz'];
%     case 'HD-BET' 
%         nii = load_untouch_nii([output_data_path '/brain_mask_HD.nii.gz']);
%         brain_mask = double(nii.img);
%         
%         nii_mask_path = [output_data_path '/brain_mask_HD.nii.gz'];
%     otherwise
%         disp('! Brain extraction method not specified correctly.');
% end

%% Crop the images for faster computation

fprintf('## Cropping images \n');

[x, y, z] = ind2sub(size(brain_mask), find(brain_mask > 0));
gap = 5; % set a magin of 5 pixels
x1 = max(1, min(x)-gap); x2 = min(size(brain_mask,1), max(x)+gap);
y1 = max(1, min(y)-gap); y2 = min(size(brain_mask,2), max(y)+gap);
z1 = max(1, min(z)-gap); z2 = min(size(brain_mask,3), max(z)+gap);
if mod(x2 - x1, 2) == 0
    if x2 == size(brain_mask,1)
        x1 = x1 - 1;
    else
        x2 = x2 + 1;
    end
end
if mod(y2 - y1, 2) == 0
    if y2 == size(brain_mask,2)
        y1 = y1 - 1;
    else
        y2 = y2 + 1;
    end
end
if mod(z2 - z1, 2) == 0 && z2 == size(brain_mask,3)
    z1 = z1 - 1;
elseif mod(z2 - z1, 2) == 0
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
    
    tic;
    
    totalField_echo = zeros(size(iMag));
    localField_echo = zeros(size(iMag));
    totalField_echoNorm = zeros(size(iMag));
    localField_echoNorm = zeros(size(iMag));
    for ii = 1:Necho
        totalField_echo(:,:,:,ii) = ...
            UnwrapPhaseMacro(iPha(:,:,:,ii), brain_mask, matrix_size, ...
            voxel_size, algorParam, header); % rad
        localField_echo(:,:,:,ii) = ...
            BackgroundRemovalMacro(totalField_echo(:,:,:,ii), brain_mask, matrix_size, ...
            voxel_size, algorParam, header); % rad
        totalField_echoNorm(:,:,:,ii) = totalField_echo(:,:,:,ii)/header.te(ii)*TEs(end); % rad
        localField_echoNorm(:,:,:,ii) = localField_echo(:,:,:,ii)/header.te(ii)*TEs(end); % rad
    end
    
    totalField_mean = mean(totalField_echoNorm,4);
    localField_mean = mean(localField_echoNorm,4)/(2*pi); % rad -> Hz at echo(0.0165ms)
    
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
    
    Tcomp.unwrap_bfr = toc;
    
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
header.te = TEs(end); % s
header.delta_TE = TEs(end); % s

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

if opts.iLSQR && exist([output_data_path '/QSM_iLSQR_meanEcho.nii.gz'],'file') ~= 2
    
    tic;
    
    algorParam.qsm.method                       = 'STI suite iLSQR';
    algorParam.qsm.maxiter                      = 100;
    algorParam.qsm.tol1                         = 0.01;
    algorParam.qsm.tol2                         = 0.001;
    algorParam.qsm.padsize                      = [0,0,6]; % Eason
    algorParam.qsm.reference_tissue             = 'None';
    
    QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM.*brain_mask;
        
    Tcomp.iLSQR = toc;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_iLSQR_meanEcho.nii.gz']);
    
end

%% QSM - MEDI

if opts.All && exist([output_data_path '/QSM_MEDI_meanEcho.nii.gz'],'file') ~= 2
    
    tic;
    
    algorParam.qsm.method                       = 'MEDI';
    algorParam.qsm.isSMV                        = false;
    algorParam.qsm.merit                        = false;
    algorParam.qsm.lambda                       = 1.78e2; % Eason 2000, paper 1000
    
    QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM.*brain_mask;
        
    Tcomp.MEDI = toc;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_MEDI_meanEcho.nii.gz']);
    
end

%% QSM - FANSI

if opts.All && exist([output_data_path '/QSM_FANSI_nonlinearTV_meanEcho.nii.gz'],'file') ~= 2
    
    tic;
    
    algorParam.qsm.method                       = 'FANSI';
    algorParam.qsm.solver                       = 'nonlinear';
    algorParam.qsm.maxiter                      = 150;
    algorParam.qsm.lambda                       = 1e-4; % 3e-5; % 0.0015;
    algorParam.qsm.mu1                          = 100*algorParam.qsm.lambda; % 5e-5
    header.weights                              = iMag(:,:,:,1).*brain_mask;
    
    QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM.*brain_mask;
        
    Tcomp.FANSI = toc;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_FANSI_nonlinearTV_meanEcho.nii.gz']);
    
    header.weights = wmap;
    
end

%% QSM - Star-QSM

if opts.All && exist([output_data_path '/QSM_STARQSM_meanEcho.nii.gz'],'file') ~= 2
    
    tic;
    
    algorParam.qsm.method                       = 'Star-QSM';
    
    QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM.*brain_mask;
        
    Tcomp.STARQSM = toc;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_STARQSM_meanEcho.nii.gz']);
    
end

%% QSM - HD-QSM

if opts.All && exist([output_data_path '/QSM_HDQSM_meanEcho.nii.gz'],'file') ~= 2
    
    tic;
    
    algorParam.qsm.method           = 'HD-QSM';
    algorParam.qsm.tol_update       = 1.0;
    algorParam.qsm.maxOuterIterL2   = 280;
    algorParam.qsm.alphaL2          = 6.3e-5; % default 10^-4.785
    algorParam.qsm.mu1L2            = 100 * algorParam.qsm.alphaL2;
    algorParam.qsm.maxOuterIterL1   = 20;
    
    QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM.*brain_mask;
        
    Tcomp.HDQSM = toc;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_HDQSM_meanEcho.nii.gz']);
    
end

%% QSM - SS_TGV

if opts.All && exist([output_data_path '/QSM_SSTGV_meanEcho.nii.gz'],'file') ~= 2
    
    tic;
    
    algorParam.qsm.method                       = 'SS-TGV';
    algorParam.qsm.SS.method                    = 'SS-TGV';
    algorParam.qsm.SSTGV.alpha0                 = 7.9e-3; % default 7e-3
    algorParam.qsm.SSTGV.mu0                    = 1e-1;
    % algorParam.qsm.SSTGV.tol_soln               = 0.1; % default 1
    % algorParam.qsm.SSTGV.maxOuterIter           = 300; % default 100
    
    QSM = -QSMMacro(totalField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM.*brain_mask;
    
    Tcomp.SSTGV = toc;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_SSTGV_meanEcho.nii.gz']);
    
end

%% QSM - SS_TV

if opts.All && exist([output_data_path '/QSM_SSTV_meanEcho.nii.gz'],'file') ~= 2
    
    tic;
    
    algorParam.qsm.method                       = 'SS-TGV';
    algorParam.qsm.SS.method                    = 'SS-TV';
    algorParam.qsm.SSTV.alpha                   = 4.0e-3; % default 7e-3
    algorParam.qsm.SSTV.mu1                     = 1e-1;
    % algorParam.qsm.SSTV.tol_soln                = 0.1;
    % algorParam.qsm.SSTV.maxOuterIter            = 300; % default 100
    
    QSM = -QSMMacro(totalField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM.*brain_mask;
    
    Tcomp.SSTV = toc;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_SSTV_meanEcho.nii.gz']);
    
end

%% QSM - QSIP

if opts.All && exist([output_data_path '/QSM_QSIP_meanEcho.nii.gz'],'file') ~= 2
    
    tic;
    
    tmp = sqrt(mean(iMag.^2,4));
    header.magn = tmp;
    
    algorParam.qsm.method                       = 'QSIP';
    algorParam.qsm.QSIP.atlas_thr               = 0.96;
    algorParam.qsm.QSIP.atlas_flag              = 1;
    algorParam.qsm.QSIP.num_iter                = 300;
    
    QSM = QSMMacro(totalField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM.*brain_mask;
    
    Tcomp.QSIP = toc;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_QSIP_meanEcho.nii.gz']);
    
    header.magn = iMag;
    
end

%% QSM - QSMGAN

if opts.QSMGAN && exist([output_data_path '/QSM_QSMGAN_meanEcho.nii.gz'],'file') ~= 2
    
    tic;
    
    % prepare local field map for QSMGAN
    
    nii = load_untouch_nii([output_data_path '/localField_meanEcho_crop.nii.gz']);
    temp = double(nii.img)/(b0*gyro*header.te);
    temp(temp == 0) = epsilon;
    temp(isnan(temp)) = epsilon;
    temp(isinf(temp)) = epsilon;
    nii.img = temp;
    save_untouch_nii(nii, [output_data_path '/localField_meanEcho_crop_ppm.nii.gz']);
    
    try
        cd('/working/lupolab/jingwen/001_QSM/01_Code/QSMGAN/Code');
        
        cmd = ['source /netopt/rhel7/versions/python/Anaconda3-5.2.0/etc/profile.d/conda.sh; ' ...
            'conda activate /working/lupolab/jingwen/conda/envs/QSM_DL; '...
            'python make_swan_qsm_DL_JY.py ' output_data_path '/localField_meanEcho_crop_ppm.nii.gz; '...
            'conda deactivate; '];
        
        system(cmd);
    catch
        disp('Python DL failed.');
    end
    
    cd('/working/lupolab/jingwen');
    
    nii = load_untouch_nii([output_data_path '/localField_meanEcho_crop_ppm_susc_DL.nii.gz']);
    QSM = double(nii.img);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM.*brain_mask;
        
    Tcomp.QSMGAN = toc;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_QSMGAN_meanEcho.nii.gz']);
    
end

%% QSM - QSMnet+ 

if opts.All && exist([output_data_path '/QSM_QSMnet_meanEcho.nii.gz'],'file') ~= 2
    
    tic;
    
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
    
    try
        cd('/working/lupolab/jingwen/001_QSM/01_Code/QSMnet-master/Code');
        
        cmd = ['source /netopt/rhel7/versions/python/Anaconda3-5.2.0/etc/profile.d/conda.sh; ' ...
            'conda activate /working/lupolab/jingwen/conda/envs/QSM_DL; '...
            'python inference.py ' output_data_path '; '...
            'conda deactivate; '];
        
        system(cmd);
    catch
        disp('Python DL failed.');
    end
    
    % check results
    
    % figure('Position',[100,100,900,400]);
    % load('Data/Test/Prediction/subject1_QSMnet+_64_25.mat');
    % subplot(1,2,1); plot3D(sus);
    % load('Data/Test/Prediction/subject99_QSMnet+_64_25.mat');
    % subplot(1,2,2); plot3D(sus); caxis([-0.03 0.03])
    
    load([output_data_path '/QSMnet_QSMnet+_64_25.mat']);
    QSM = flip(imresize3(sus,header.matrixSize),2);
    
    QSM_output = zeros(matrix_size_ori);
    QSM_output(crop.X,crop.Y,crop.Z) = QSM.*brain_mask;
    
    Tcomp.QSMnet = toc;
    
    % save maps
    nii = load_untouch_nii(nii_mask_path);
    nii.img = QSM_output;
    save_untouch_nii(nii, [output_data_path '/QSM_QSMnet_meanEcho.nii.gz']);
    
    cd('/working/lupolab/jingwen');
    
end

%% save output

fprintf('## Done!! \n');

outfile = sprintf('%s/Log/Output_%s.mat', output_data_path, datestr(now,'yyyy-mm-dd-HHMM'));
save(outfile, 'Tcomp');

if opts.writeLog
    diary off
end

end

%% set default parameter for unspecific input

function input = check_and_set_default_input(input_data, currentPath)

input = input_data;

try     input.magnitudeFile     = input_data.magnitudeFile;       
catch;  input.magnitudeFile 	= sprintf('%s/Magni.nii.gz', currentPath); end
try     input.phaseFile         = input_data.phaseFile;       
catch;  input.phaseFile         = sprintf('%s/Phase.nii.gz', currentPath); end
try     input.headerFile        = input_data.headerFile;       
catch;  input.headerFile        = sprintf('%s/sepia_header.mat', currentPath); end

end

function opts2 = check_and_set_default_opts(opts)

opts2 = opts;

try     opts2.writeLog  = opts.writeLog;    catch;  opts2.writeLog  = 1; end
try     opts2.isGPU     = opts.isGPU;       catch;  opts2.isGPU     = 1; end
try     opts2.BETmethod = opts.BETmethod;   catch;  opts2.BETmethod = 'BET'; end
try     opts2.iLSQR     = opts.iLSQR;       catch;  opts2.iLSQR     = 1; end
try     opts2.QSMGAN    = opts.QSMGAN;      catch;  opts2.QSMGAN    = 1; end
try     opts2.All       = opts.All;         catch;  opts2.All       = 0; end
try     opts2.Necho     = opts.Necho;       catch;  opts2.Necho     = 4; end

end
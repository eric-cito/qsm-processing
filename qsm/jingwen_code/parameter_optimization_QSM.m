clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

%% subject for optimization

outputPath = '/working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308/';
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

load([output_data_path '/sepia_header.mat']);
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

% output_data_path = [outputPath '/MEDI'];
% mkdir(output_data_path);
% 
% algorParam.qsm.method                       = 'MEDI';
% algorParam.qsm.isSMV                        = false;
% algorParam.qsm.merit                        = false;
% algorParam.qsm.lambda                       = 1000; % Eason 2000
% 
% for lambda  = 10.^[3:0.1:5]
%     
%     algorParam.qsm.lambda                       = lambda;
%     
%     QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
%     
%     QSM_output = zeros(matrix_size_ori);
%     QSM_output(crop.X,crop.Y,crop.Z) = QSM;
%     
%     % save maps
%     nii = load_untouch_nii(nii_mask_path);
%     nii.img = QSM_output;
%     
%     paraStr = sprintf('l10E%.2f', log10(lambda));
%     paraStr = strrep(paraStr,'.','p');
%     save_untouch_nii(nii, [output_data_path '/QSM_MEDI_' paraStr '.nii.gz']);
%     
% end

%% QSM - FANSI

% output_data_path = [outputPath '/output_HDBET/FANSI'];
% mkdir(output_data_path);
% 
% algorParam.qsm.method                       = 'FANSI';
% algorParam.qsm.solver                       = 'nonlinear';
% algorParam.qsm.maxiter                      = 150;
% header.weights                              = iMag(:,:,:,1).*brain_mask;
% 
% for lambda  = 10.^[-3.6:-0.1:-4.4]
%     
% %     for mu = 10.^[log10(lambda):0.5:log10(100*lambda)]
% 
%         mu = 100*lambda;
%         algorParam.qsm.lambda                   = lambda;
%         algorParam.qsm.mu1                      = mu;
%         
%         QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
%         
%         QSM_output = zeros(matrix_size_ori);
%         QSM_output(crop.X,crop.Y,crop.Z) = QSM;
%         
%         % save maps
%         nii = load_untouch_nii(nii_mask_path);
%         nii.img = QSM_output;
%         
%         paraStr = sprintf('l10E%.2f', log10(lambda));
%         paraStr = strrep(paraStr,'.','p');
%         save_untouch_nii(nii, [output_data_path '/QSM_FANSI_' paraStr '.nii.gz']);
%         
% %     end
%     
% end
% 
% header.weights = wmap;

%% QSM - HD-QSM

% output_data_path = [outputPath '/HDQSM'];
% mkdir(output_data_path);
% 
% algorParam.qsm.method           = 'HD-QSM';
% algorParam.qsm.tol_update       = 1.0;
% algorParam.qsm.maxOuterIterL2   = 280;
% algorParam.qsm.maxOuterIterL1   = 20;
% 
% for lambda  = 10.^[-5.4:0.1:-3.8]
%     
%     algorParam.qsm.alphaL2          = lambda;
%     algorParam.qsm.mu1L2            = 100 * lambda;
%     
%     QSM = -QSMMacro(localField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
%     
%     QSM_output = zeros(matrix_size_ori);
%     QSM_output(crop.X,crop.Y,crop.Z) = QSM;
%     
%     % save maps
%     nii = load_untouch_nii(nii_mask_path);
%     nii.img = QSM_output;
%     
%     paraStr = sprintf('l10E%.2f', log10(lambda));
%     paraStr = strrep(paraStr,'.','p');
%     save_untouch_nii(nii, [output_data_path '/QSM_HDQSM_' paraStr '.nii.gz']);
%     
% end

%% QSM - SS_TGV

% output_data_path = [outputPath '/SSTGV'];
% mkdir(output_data_path);
% 
% algorParam.qsm.method                       = 'SS-TGV';
% algorParam.qsm.SS.method                    = 'SS-TGV';
% algorParam.qsm.SSTGV.mu0                    = 1e-1;
% 
% for lambda  = 10.^[-1.9:0.1:-1]
%     
%     algorParam.qsm.SSTGV.alpha0                 = lambda;
%     
%     QSM = -QSMMacro(totalField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
%     
%     QSM_output = zeros(matrix_size_ori);
%     QSM_output(crop.X,crop.Y,crop.Z) = QSM;
%     
%     % save maps
%     nii = load_untouch_nii(nii_mask_path);
%     nii.img = QSM_output;
%     
%     paraStr = sprintf('l10E%.2f', log10(lambda));
%     paraStr = strrep(paraStr,'.','p');
%     save_untouch_nii(nii, [output_data_path '/QSM_SSTGV_' paraStr '.nii.gz']);
%     
% end

%% QSM - SS_TV

% output_data_path = [outputPath '/SSTV'];
% mkdir(output_data_path);
% 
% algorParam.qsm.method                       = 'SS-TGV';
% algorParam.qsm.SS.method                    = 'SS-TV';
% algorParam.qsm.SSTV.mu1                     = 1e-1;
% 
% for lambda  = 10.^[-3.5:0.1:-1]
%     
%     algorParam.qsm.SSTV.alpha                 = lambda;
%     
%     QSM = -QSMMacro(totalField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
%     
%     QSM_output = zeros(matrix_size_ori);
%     QSM_output(crop.X,crop.Y,crop.Z) = QSM;
%     
%     % save maps
%     nii = load_untouch_nii(nii_mask_path);
%     nii.img = QSM_output;
%     
%     paraStr = sprintf('l10E%.2f', log10(lambda));
%     paraStr = strrep(paraStr,'.','p');
%     save_untouch_nii(nii, [output_data_path '/QSM_SSTV_' paraStr '.nii.gz']);
%     
% end

%% QSM - QSIP

output_data_path = [outputPath '/QSIP'];
mkdir(output_data_path);

algorParam.qsm.method                       = 'QSIP';
algorParam.qsm.QSIP.atlas_thr               = 0.96;
algorParam.qsm.QSIP.atlas_flag              = 1;
algorParam.qsm.QSIP.num_iter                = 300;
algorParam.qsm.QSIP.EXTWeight               = 1e10;

for lambda1 = 10.^[-10:0.5:-5] % default 10-8
    
    for lambda2 = 10.^[-8:0.5:-3] % default 10-6
        
        paraStr = sprintf('l10E%.2fl10E%.2f', log10(lambda1), log10(lambda2));
        paraStr = strrep(paraStr,'.','p');
        if exist([output_data_path '/QSM_QSIP_' paraStr '.nii.gz'],'file') == 2
            continue;
        end
        
        algorParam.qsm.QSIP.LPWeight                 = lambda1;
        algorParam.qsm.QSIP.OBJWeight                = lambda2;
        
        QSM = -QSMMacro(totalField_mean, brain_mask, matrix_size, voxel_size, algorParam, header);
        
        QSM_output = zeros(matrix_size_ori);
        QSM_output(crop.X,crop.Y,crop.Z) = QSM;
        
        % save maps
        nii = load_untouch_nii(nii_mask_path);
        nii.img = QSM_output;

        save_untouch_nii(nii, [output_data_path '/QSM_QSIP_' paraStr '.nii.gz']);
        
    end
    
end

%% Frequency masking

kernel = dipole_kernel_angulated(header.matrixSize, header.voxelSize, header.b0dir);

% Create the masks used in the Frequency Analysis
[m1, m2, m3] = create_freqmasks(header.voxelSize, kernel);

method = 'QSIP';
lambda1 = 10.^[-10:0.5:-5];
lambda2 = 10.^[-8:0.5:-3];
% FANSI: -5:0.25:-4.5 -4.3:0.1:-3.5 -3.25:0.25:-2
% MEDI: 1:0.25:4
% SSTGV/SSTV: -3:0.1:-1.3
% HDQSM: -6:0.25:-5.25 -5.2:0.1:-3.8 -3.75:0.25:-2
e1 = zeros(length(lambda1), length(lambda2));
e2 = zeros(length(lambda1), length(lambda2));
e3 = zeros(length(lambda1), length(lambda2));
for ii = 1:length(lambda1)
    for jj = 1:length(lambda2)
        
        % load image
        paraStr = sprintf('l10E%.2fl10E%.2f', log10(lambda1(ii)), log10(lambda2(jj)));
        paraStr = strrep(paraStr,'.','p');
        fprintf(' >> Loading data %s \n', paraStr);
        
        try
            nii = load_untouch_nii([method '/QSM_' method '_' paraStr '.nii.gz']);
            QSM = -nii.img(crop.X,crop.Y,crop.Z);
            % Calculate the Mean Amplitudes for each mask
            [e1(ii,jj), e2(ii,jj), e3(ii,jj)] = compute_freqe(QSM.*brain_mask,m1,m2,m3);
        catch
            e1(ii,jj) = nan; e2(ii,jj) = nan; e3(ii,jj) = nan; 
        end
        
    end
end

%% Frequency analysis
% [opt23index, alpha_opt] = draw_freque(lambda,e1,e2,e3,13);
% zeta23 = zetafunc(e2,e3);
% opt23  = find(zeta23 == min(zeta23));
% 
% set(0,'defaultAxesFontSize',14);
% set(0,'defaultLineLineWidth',1);
% 
% figure('position',[100 100 400 400]); 
% plot(lambda, zeta23,'-o','Color','k','LineWidth',1.0); hold on;
% set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
% xlabel('Regularization weight'); xlim([lambda(1) lambda(end)]);
% ylabel('\zeta_{23} function')
% scatter(lambda(opt23), zeta23(opt23),'filled','kd');
% legend({'\zeta_{23}', 'Optimal \lambda'},'location','best','box','off');
% title('Frequency Equalization Plot');
% 
% fprintf('Minimum Zeta23: %.2e, 10^%.2f \n', lambda(opt23), log10(lambda(opt23)));
% export_fig(['optimization_result/' method], '-png','-transparent'); % close;

%% Frequency analysis - 2D
zeta23 = zetafunc(e2,e3);
[i,j]  = find(zeta23 == min(zeta23(:)));
fprintf('Minimum Zeta23 lambda1: 10^%.2f, lambda2: 10^%.2f \n', log10(lambda1(i)), log10(lambda2(j)));

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

[Y,X] = meshgrid(log10(lambda1), log10(lambda2));
figure('position',[100 100 400 400]); 
% scatter3(X(:), Y(:), zeta23(:)); hold on;
surface(X, Y, zeta23); hold on;
view(3);
xlabel('lambda 1'); ylabel('lambda 2');
zlabel('\zeta_{23} function');

% zeta23 = zetafunc(e2,e3);
% opt23  = find(zeta23 == min(zeta23));
% 
% set(0,'defaultAxesFontSize',14);
% set(0,'defaultLineLineWidth',1);
% 
% figure('position',[100 100 400 400]); 
% plot(lambda, zeta23,'-o','Color','k','LineWidth',1.0); hold on;
% set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log')
% xlabel('Regularization weight'); xlim([lambda(1) lambda(end)]);
% ylabel('\zeta_{23} function')
% scatter(lambda(opt23), zeta23(opt23),'filled','kd');
% legend({'\zeta_{23}', 'Optimal \lambda'},'location','best','box','off');
% title('Frequency Equalization Plot');
% 
% fprintf('Minimum Zeta23: %.2e, 10^%.2f \n', lambda(opt23), log10(lambda(opt23)));

%% plot Frequency ROI

% kernel = dipole_kernel_angulated(header.matrixSize, header.voxelSize, header.b0dir);
% 
% % Create the masks used in the Frequency Analysis
% [m1, m2, m3] = create_freqmasks(header.voxelSize, kernel);
% 
% nii = load_untouch_nii([outputPath '/output_HDBET/HDQSM/QSM_HDQSM_l10E-3p50.nii.gz']);
% QSM = double(nii.img(crop.X,crop.Y,crop.Z));
% fchi = (abs(fftn(QSM)));
% QSMkplot = fftshift(fchi).^0.2;
% maskkplot = fftshift(m1+2*m2+3*m3);
% 
% indX = floor(size(maskkplot,1)/2);
% indY = floor(size(maskkplot,2)/2);
% indZ = floor(size(maskkplot,3)/2);
% 
% X = rot90(squeeze(QSMkplot(:,indY,:)),1);
% Y = rot90(squeeze(maskkplot(:,indY,:)),1);
% 
% figure('position',[100 100 700 700]);
% ax1 = axes; imagesc(X); colormap(ax1,'gray');
% axis tight; axis off; axis equal;
% ax2 = axes; imagesc(Y,'alphadata',Y>0); colormap(ax2,flip(parula(3)));
% caxis(ax2,[min(nonzeros(Y)) max(nonzeros(Y))]);
% ax2.Visible = 'off';
% axis tight; axis off; axis equal;
% linkprop([ax1 ax2],'Position');
% export_fig(['optimization_result/HDQSM-350'], '-png'); % close;

%% save images

% foldername = 'b4468_t12308';
% savename = 'HDQSM_meanEcho';
% QSMname = [foldername '/QSM_' savename '.nii.gz'];
% targetSize = [220 220 148]; 
% 
% nii = load_untouch_nii([outputPath '/output_HDBET/' QSMname]);
% QSM = nii.img(crop.X,crop.Y,crop.Z).*brain_mask;
% currSize = size(QSM);
% padSize = (targetSize - currSize)/2;
% 
% QSM_Ax = padarray(flip(rot90(QSM(:,:,66),1)),padSize([2 1]));
% QSM_Cor = padarray(flip(rot90(squeeze(QSM(87,:,:)),1)),padSize([3 2]));
% QSM_Sag = padarray(flip(rot90(squeeze(QSM(:,83,:)),1)),padSize([3 1]));
% 
% fig = figure('position',[100 100 500 500]); 
% imshow([QSM_Ax; QSM_Cor; QSM_Sag],[-0.15 0.15]); colormap gray;
% pause;
% export_fig(['optimization_result/' savename], '-png'); close;

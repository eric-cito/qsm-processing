clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

%% subject for optimization

subjPath_list = {'/working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308/output_HDBET/', ...
    '/working/lupolab/jingwen/001_QSM/Test_data/b4469_t12309/', ...
    '/working/lupolab/jingwen/001_QSM/Test_data/b4473_t12317/', ...
    '/working/lupolab/jingwen/001_QSM/Test_data/b4489_t12374/', ...
    '/working/lupolab/jingwen/001_QSM/Test_data/b4536_t12517/'};

output_data_path = subjPath_list{1};

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

%% load header

load([output_data_path '/../sepia_header.mat']);
TEs = header.te;
deltaTE = header.delta_TE;
gyro = 42.57747892;
b0 = header.b0;
epsilon = 1e-15;

header_ori = header;

%% optimization methods

method_list = {'FANSI','MEDI','SSTGV','SSTV','HDQSM'};
lambda_list = {[-5:0.25:-4.5 -4.3:0.1:-3.5 -3.25:0.25:-2], ...
    [1:0.25:3 3.1:0.1:5], ...
    [-3:0.1:-1.5], [-3:0.1:-1.5], ...
    [-6:0.25:-5.25 -5.2:0.1:-3.8 -3.75:0.25:-2]};

%% loop through subjects

for mm = 5 % 1:length(method_list)
    
    method = method_list{mm};
    lambda = 10.^lambda_list{mm};
    
    e1 = zeros(length(lambda), length(subjPath_list));
    e2 = zeros(length(lambda), length(subjPath_list));
    e3 = zeros(length(lambda), length(subjPath_list));
    
    figure('position',[100 100 400 400]);
    
    for ii = 1:length(subjPath_list)
        
        output_data_path = subjPath_list{ii};
        
        % load data
        cd(output_data_path);
        nii = load_untouch_nii([output_data_path '/brain_mask_HD.nii.gz']);
        brain_mask = double(nii.img);
        
        % Crop the images
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
        
        matrix_size_ori = size(brain_mask_ori);
        matrix_size = size(brain_mask);
        
        header.matrixSize = matrix_size;
        
        % Frequency masking
        kernel = dipole_kernel_angulated(header.matrixSize, header.voxelSize, header.b0dir);
        
        % Create the masks used in the Frequency Analysis
        [m1, m2, m3] = create_freqmasks(header.voxelSize, kernel);
        
        % Calculate frequency components
        for ll = 1:length(lambda)
            % load image
            paraStr = sprintf('l10E%.2f', log10(lambda(ll)));
            paraStr = strrep(paraStr,'.','p');
            
            % fprintf('Loading %s %s\n', method, paraStr);
            try
                nii = load_untouch_nii([method '/QSM_' method '_' paraStr '.nii.gz']);
                QSM = nii.img(crop.X,crop.Y,crop.Z);
                % Calculate the Mean Amplitudes for each mask
                [e1(ll,ii), e2(ll,ii), e3(ll,ii)] = compute_freqe(QSM.*brain_mask,m1,m2,m3);
            catch
                e1(ll,ii) = nan; e2(ll,ii) = nan; e3(ll,ii) = nan;
            end
        end
        
        % Frequency analysis
        zeta23 = zetafunc(e2(:,ii),e3(:,ii));
        opt23  = find(zeta23 == min(zeta23));
        
        plot(lambda, zeta23,'-o','LineWidth',1.0); hold on;
        fprintf('%s Minimum Zeta23: %.2e, 10^%.2f \n', method, lambda(opt23), log10(lambda(opt23)));
        
    end
    
    % set(gca, 'YScale', 'log'); 
    ylim([0 0.01]);
    set(gca, 'XScale', 'log')
    xlabel('Regularization weight'); xlim([lambda(1) lambda(end)]);
    ylabel('\zeta_{23} function')
    title([method ' - Frequency Equalization Plot']);
    
    legend({'Subj1', 'Subj2', 'Subj3', 'SUbj4', 'SUbj5'},'location','best','box','off');
    
end

%% save images

img_root = '/working/lupolab/jingwen/001_QSM/temp';

method = 'FANSI';
title(method);
export_fig([img_root '/optimization_' method], '-png','-transparent'); % close;

%% plot Frequency ROI

kernel = dipole_kernel_angulated(header.matrixSize, header.voxelSize, header.b0dir);

% Create the masks used in the Frequency Analysis
[m1, m2, m3] = create_freqmasks(header.voxelSize, kernel);

nii = load_untouch_nii([output_data_path '/HDQSM/QSM_HDQSM_l10E-6p00.nii.gz']);
QSM = double(nii.img(crop.X,crop.Y,crop.Z));
fchi = (abs(fftn(QSM)));
QSMkplot = fftshift(fchi).^0.5;
maskkplot = fftshift(m2+2*m3);

indX = floor(size(maskkplot,1)/2)+1;
indY = floor(size(maskkplot,2)/2)+1;
indZ = floor(size(maskkplot,3)/2)+1;

X = rot90(squeeze(QSMkplot(:,indY,:)),1);
Y = rot90(squeeze(maskkplot(:,indY,:)),1);

figure('position',[100 100 700 700]);
ax1 = axes; imagesc(X, prctile(X(:),[2 99])); colormap(ax1,'gray');
axis tight; axis off; axis equal;
ax2 = axes; imagesc(Y,'alphadata',Y>0); colormap(ax2,flip(parula(3)));
caxis(ax2,[0 3]);
ax2.Visible = 'off';
axis tight; axis off; axis equal;
linkprop([ax1 ax2],'Position');
export_fig([img_root '/optim_kexample'], '-png'); % close;

%% plot all k & images

paraList = {'l10E-6p00','l10E-4p20','l10E-3p25'};

targetSize = [220 220 148];

for ii = 1:length(paraList)
    
    nii = load_untouch_nii([output_data_path '/HDQSM/QSM_HDQSM_' paraList{ii} '.nii.gz']);
    QSM = double(nii.img(crop.X,crop.Y,crop.Z));
    fchi = (abs(fftn(QSM)));
    QSMkplot = fftshift(fchi).^0.5;
    
    X = rot90(squeeze(QSMkplot(:,indY,:)),1);
    Y = rot90(squeeze(maskkplot(:,indY,:)),1);
    
    figure('position',[100 100 700 700]);
    ax1 = axes; imagesc(X, prctile(X(:),[2 99])); colormap(ax1,'gray');
    axis tight; axis off; axis equal;
    pause(1);
    export_fig([img_root '/optim_' paraList{ii}], '-png'); % close;
    
%     currSize = size(QSM);
%     padSize = (targetSize - currSize)/2;
%     
%     QSM_Ax = padarray(flip(rot90(QSM(:,:,66),1)),padSize([2 1]));
%     QSM_Cor = padarray(flip(rot90(squeeze(QSM(87,:,:)),1)),padSize([3 2]));
%     QSM_Sag = padarray(flip(rot90(squeeze(QSM(:,83,:)),1)),padSize([3 1]));
%     
%     figure('position',[100 100 500 500]);
%     imshow([QSM_Ax; QSM_Cor; QSM_Sag],[-0.15 0.15]); colormap gray;
    
end

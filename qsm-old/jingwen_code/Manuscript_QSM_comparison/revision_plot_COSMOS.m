clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath('/home/jyao3/010_MATLAB_Utils/');
addpath('/home/jyao3/010_MATLAB_Utils/violinplot');
addpath('/home/jyao3/010_MATLAB_Utils/export_fig');
addpath('/home/jyao3/010_MATLAB_Utils/NIfTI');

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

%% subject list from COSMOS

mainPath = '/working/lupolab/eason/DL_QSM/data/7T_cosmos';
subjFolders = dir(mainPath);
subjFolders(~[subjFolders.isdir]) = [];
subjFolders(~cellfun(@isempty,strfind({subjFolders.name},'.'))) = [];
subjFolders(cellfun(@isempty,strfind({subjFolders.name},'volunteer'))) = [];

%% Loop through subjects

matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data';
load([matout_root '/COSMOS_stats.mat'], 'NMSE', 'dNMSE', 'SSIM', 'XSIM', 'HFEN');

QSMfile_list = {
    'QSM_iLSQR_meanEcho' ...
    'QSM_QSIP_meanEcho' ...
    'QSM_SSTGV_meanEcho' ...
    'QSM_SSTV_meanEcho' ...
    'QSM_STARQSM_meanEcho' ...
    'QSM_FANSI_nonlinearTV_meanEcho' ...
    'QSM_HDQSM_meanEcho' ...
    'QSM_MEDI_meanEcho' ...
    'QSM_QSMGAN_meanEcho' ...
    'QSM_QSMnet_meanEcho' ...
    'QSM_xQSM2_meanEcho' ...
    'QSM_iQSM2_meanEcho'};

% NMSE = zeros(12,length(subjFolders));
% SSIM = zeros(12,length(subjFolders));
% XSIM = zeros(12,length(subjFolders));
% HFEN = zeros(12,length(subjFolders));
% dNMSE = zeros(12,length(subjFolders));
for ii = 1:length(subjFolders)
    
    fprintf('Processing %s \n', subjFolders(ii).name);
    
    dataPath = [mainPath '/' subjFolders(ii).name '/scan1/swan_qsm/COSMOSmask_allQSM'];
    
    if exist([dataPath '/QSM_COSMOS.nii.gz'],'file') ~= 2
        cmd = sprintf('3dcopy %s %s', [mainPath '/' subjFolders(ii).name '/scan1_cosmos_VSHARP.nii'], ...
            [dataPath '/temp.nii.gz']);
        system(cmd);
        cmd = sprintf('3drefit -orient RAS -duporigin %s %s', [dataPath '/QSM_iLSQR_meanEcho.nii.gz'], ...
            [dataPath '/temp.nii.gz']);
        system(cmd);
        cmd = sprintf('3dresample -master %s -prefix %s -input %s', [dataPath '/QSM_iLSQR_meanEcho.nii.gz'], ...
            [dataPath '/QSM_COSMOS.nii.gz'], [dataPath '/temp.nii.gz']);
        system(cmd);
    end
    
    if exist([dataPath '/QSM_iLSQR_meanEcho_reg.nii.gz'],'file') == 0
        cmd = sprintf('flirt -v -in %s -ref %s -out %s -omat %s', ...
            [dataPath '/QSM_iLSQR_meanEcho.nii.gz'], ...
            [dataPath '/QSM_COSMOS.nii.gz'], ...
            [dataPath '/QSM_iLSQR_meanEcho_reg.nii.gz'], ...
            [dataPath '/qsm2cosmos.mat']);
        system(cmd);
    end
    
    % delete([dataPath '/QSM_iQSM2_meanEcho.nii.gz']);
    % delete([dataPath '/QSM_iQSM2_meanEcho_reg.nii.gz']);
    if exist([dataPath '/QSM_iQSM2_meanEcho.nii.gz'],'file') == 0
        cmd = sprintf('3dresample -master %s -prefix %s -input %s', [dataPath '/QSM_COSMOS.nii.gz'], ...
            [dataPath '/QSM_iQSM2_meanEcho.nii.gz'], [dataPath '/iQSM2/iQSM_echo_fitted.nii']);
        system(cmd);
    end
    
    nii = load_nii([dataPath '/QSM_COSMOS.nii.gz']);
    cosmos = double(nii.img);
    
    nii = load_nii([dataPath '/brain_mask_HD.nii.gz']);
    brain_mask = double(nii.img);
    
    SE = strel('sphere',2);
    brain_mask = imerode(brain_mask,SE);
    
    for qq = 1:length(QSMfile_list)
        
        fprintf('  # %s \n', QSMfile_list{qq});
        
        if exist([dataPath '/' QSMfile_list{qq} '_reg.nii.gz'],'file') == 0
            cmd = sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s', ...
                [dataPath '/' QSMfile_list{qq} '.nii.gz'], ...
                [dataPath '/QSM_COSMOS.nii.gz'], ...
                [dataPath '/' QSMfile_list{qq} '_reg.nii.gz'], ...
                [dataPath '/qsm2cosmos.mat']);
            system(cmd);
        end
        
        nii = load_nii([dataPath '/' QSMfile_list{qq} '_reg.nii.gz']);
        qsm = double(nii.img);
        
        if qq == 9
            qsm = qsm/0.5684;
        end
        
        ind = brain_mask > 0; % ~isnan(cosmos) & (brain_mask > 0) & ~isnan(qsm);
        cosmos_mask = cosmos.*ind;
        qsm_mask = qsm.*ind;
        
        NMSE(qq,ii) = compute_rmse(qsm_mask, cosmos_mask)/100;
        SSIM(qq,ii) = ssim(qsm_mask, cosmos_mask);
        XSIM(qq,ii) = compute_xsim(qsm_mask, cosmos_mask);
        HFEN(qq,ii) = compute_hfen(qsm_mask, cosmos_mask)/100;
        
        P1 = polyfit(cosmos_mask, qsm_mask, 1);
        P(1) = 1 / P1(1);
        P(2) = -P1(2) / P1(1);
        chi_detrended = polyval(P, qsm_mask);

        dNMSE(qq,ii) = compute_rmse(chi_detrended, cosmos_mask)/100;
    end
    
end

%% save results

matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data';
save([matout_root '/COSMOS_stats.mat'], 'NMSE', 'dNMSE', 'SSIM', 'XSIM', 'HFEN');

%% plot metric

QSM_name = {'iLSQR' 'QSIP' 'SSTGV' 'SSTV' 'STARQSM' 'FANSI' 'HDQSM' 'MEDI' ...
    'QSMGAN' 'QSMnet+' 'xQSM' 'iQSM'};

figure('position', [100 100 1200 900]);
subplot(321);
[SSIMmean] = plotMet(SSIM, QSM_name);
ylabel('SSIM with COSMOS'); ylim([0.88 1]);
% export_fig('NMSE', '-png','-transparent'); % close;

subplot(322);
[XSIMmean] = plotMet(XSIM, QSM_name);
ylabel('XSIM with COSMOS'); ylim([0 0.8]);
% export_fig('SSIM', '-png','-transparent'); % close;

subplot(323);
[NMSEmean] = plotMet(NMSE, QSM_name);
ylabel('NRMSE with COSMOS'); ylim([0 1.5]);

subplot(324);
[dNMSEmean] = plotMet(dNMSE, QSM_name);
ylabel('dNRMSE with COSMOS'); ylim([0 2.2]);

subplot(325);
[HFENmean] = plotMet(HFEN, QSM_name);
ylabel('HFEN with COSMOS'); ylim([0 2]);

subplot(326);
scatter(dNMSE, XSIM);

T_COSMOS = table(QSM_name', NMSEmean, dNMSEmean, XSIMmean);

save([matout_root '/QSMcomp_all.mat'], 'T_COSMOS', '-append');

%% save plot - box

img_path = '/working/lupolab/jingwen/001_QSM/temp';

figure('position', [100 100 800 300]);
e = errorbar(mean(XSIM,2), std(XSIM,[],2), 'k'); hold on;
e.LineStyle = 'none';
b = bar(mean(XSIM,2), 0.5); hold on;
b.FaceColor = 'none';
ax = gca;
ax.XTick = 1:12;
ax.XTickLabel = QSM_name;
ax.Box = 'off';
xlim([0 13]); ylim([0 0.8]);
ylabel('XSIM with COSMOS'); pause(1)
export_fig([img_path '/XSIM_bar'], '-png','-transparent');

figure('position', [100 100 800 300]);
e = errorbar(mean(NMSE,2), std(NMSE,[],2), 'k'); hold on;
e.LineStyle = 'none';
b = bar(mean(NMSE,2), 0.5); hold on;
b.FaceColor = 'none';
ax = gca;
ax.XTick = 1:12;
ax.XTickLabel = QSM_name;
ax.Box = 'off';
xlim([0 13]); ylim([0 1.2]);
ylabel('NRMSE with COSMOS'); pause(1)
export_fig([img_path '/NRMSE_bar'], '-png','-transparent');

%% functions

function [SSIMmean] = plotMet(SSIM, QSMname)

SSIMmean = mean(SSIM,2);
SSIMstd = std(SSIM,[],2);

plot(1:size(SSIM,1), SSIMmean, 'o-', 'Color', [0 0 0], ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5]/2); hold on; %  + [0 0.4470 0.7410]
x = [1:12];
fill([x fliplr(x)], [SSIMmean'-SSIMstd' fliplr(SSIMmean'+SSIMstd')], ...
    [1 1 1]*0.2, 'FaceAlpha',0.1, 'EdgeColor','none');
ax = gca;
ax.XTick = 1:12;
ax.XTickLabel = QSMname(:);
xlim([0 13]);
ax.Box = 'off';
% ax.XTickLabelRotation = 30; 

end

function nmse = calNMSE(orgSig,recSig)

% mse=norm(orgSig(:)-recSig(:),2)^2/length(orgSig(:));
% sigEner=norm(orgSig(:))^2;
% nmse=(mse/sigEner);

nmse = sum((orgSig - recSig).^2)/sum(orgSig.^2);

end
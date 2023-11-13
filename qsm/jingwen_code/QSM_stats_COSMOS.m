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

QSMfile_list = {'QSM_FANSI_nonlinearTV_meanEcho.nii.gz' ...
    'QSM_HDQSM_meanEcho.nii.gz' ...
    'QSM_iLSQR_meanEcho.nii.gz' ...
    'QSM_MEDI_meanEcho.nii.gz' ...
    'QSM_QSIP_meanEcho.nii.gz' ...
    'QSM_QSMGAN_meanEcho.nii.gz' ...
    'QSM_QSMnet_meanEcho.nii.gz' ...
    'QSM_SSTGV_meanEcho.nii.gz' ...
    'QSM_SSTV_meanEcho.nii.gz' ...
    'QSM_STARQSM_meanEcho.nii.gz'};

NMSE = zeros(10,length(subjFolders));
SSIM = zeros(10,length(subjFolders));
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
    
    nii = load_untouch_nii([dataPath '/QSM_COSMOS.nii.gz']);
    cosmos = double(nii.img);
    
    nii = load_untouch_nii([dataPath '/brain_mask_HD.nii.gz']);
    brain_mask = double(nii.img);
    
    SE = strel('sphere',2);
    brain_mask = imerode(brain_mask,SE);
    
    for qq = 1:length(QSMfile_list)
        nii = load_untouch_nii([dataPath '/' QSMfile_list{qq}]);
        qsm = double(nii.img);
        
        if qq == 6
            qsm = qsm/0.5684;
        end
        
        ind = brain_mask > 0; % ~isnan(cosmos) & (brain_mask > 0) & ~isnan(qsm);
        NMSE(qq,ii) = calNMSE(cosmos(ind),qsm(ind));
        
        cosmos_mask = cosmos.*ind;
        qsm_mask = qsm.*ind;
        SSIM(qq,ii) = ssim(qsm_mask, cosmos_mask);
    end
    
end

% reorder
QSMorder = [3 10 1 2 4 5 8 9 6 7];
NMSE = NMSE(QSMorder,:);
SSIM = SSIM(QSMorder,:);

%% plot metric

QSMfile_list = {'QSM_FANSI' ...
    'QSM_HDQSM' ...
    'QSM_iLSQR' ...
    'QSM_MEDI' ...
    'QSM_QSIP' ...
    'QSM_QSMGAN' ...
    'QSM_QSMnet+' ...
    'QSM_SSTGV' ...
    'QSM_SSTV' ...
    'QSM_STARQSM'};
QSMfile_list = QSMfile_list(QSMorder);

QSMname = strrep(QSMfile_list,'_','-');
QSMname = strrep(QSMname,'QSM-','');
QSMname = strrep(QSMname,'-meanEcho','');

RMSEmean = mean(NMSE,2);
SSIMmean = mean(SSIM,2);
RMSEstd = std(NMSE,[],2);
SSIMstd = std(SSIM,[],2);

T = table(QSMfile_list', RMSEmean, SSIMmean);

figure('position', [100 100 900 300]);
plot(1:10, RMSEmean, 'o-', 'Color', [0 0 0], ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5]/2); hold on; %  + [0 0.4470 0.7410]
x = [1:10];
fill([x fliplr(x)], [RMSEmean'-RMSEstd' fliplr(RMSEmean'+RMSEstd')], ...
    [1 1 1]*0.2, 'FaceAlpha',0.1, 'EdgeColor','none');
ylabel('NMSE with COSMOS');
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30; xlim([0 11]); ylim([0.4 2]);
ax.Box = 'off';
export_fig('NMSE', '-png','-transparent'); % close;

figure('position', [100 100 900 300]);
plot(1:10, SSIMmean, 'o-', 'Color', [0 0 0], ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5]/2); hold on;
x = [1:10];
fill([x fliplr(x)], [SSIMmean'-SSIMstd' fliplr(SSIMmean'+SSIMstd')], ...
    [1 1 1]*0.2, 'FaceAlpha',0.1, 'EdgeColor','none');
ylabel('SSIM with COSMOS');
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30; xlim([0 11]); ylim([0.89 0.98]);
ax.Box = 'off';
export_fig('SSIM', '-png','-transparent'); % close;

%% functions

function nmse = calNMSE(orgSig,recSig)

% mse=norm(orgSig(:)-recSig(:),2)^2/length(orgSig(:));
% sigEner=norm(orgSig(:))^2;
% nmse=(mse/sigEner);

nmse = sum((orgSig - recSig).^2)/sum(orgSig.^2);

end
clc;clear;close;

addpath('/home/jyao3/010_MATLAB_Utils/export_fig');
addpath('/home/jyao3/010_MATLAB_Utils/NIfTI');
addpath('/home/jyao3/010_MATLAB_Utils/imoverlay');
addpath('/home/jyao3/010_MATLAB_Utils/');

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

display_root = ('/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/temp');

%% load MNI images
MNI_root = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/NIFTI_HC';
nii = load_nii([MNI_root '/MNI_T1w.nii.gz']);
Anat = double(nii.img);

nii = load_nii([MNI_root '/Whole_brain_segmentation_DGM.nii']);
brainROI = double(nii.img);

nii = load_nii([MNI_root '/MNI_GMWM.nii.gz']);
brainMask = double(nii.img);
% brainMask(brainMask == 0) = nan;

QSMNii_list = {'QSM_iLSQR_meanEcho' ...
    'QSM_STARQSM_meanEcho' ...
    'QSM_FANSI_nonlinearTV_meanEcho' ...
    'QSM_HDQSM_meanEcho' ...
    'QSM_MEDI_meanEcho' ...
    'QSM_QSIP_meanEcho' ...
    'QSM_SSTGV_meanEcho' ...
    'QSM_SSTV_meanEcho' ...
    'QSM_QSMGAN_meanEcho' ...
    'QSM_QSMnet_meanEcho'};

%% specify subjects

% blist = {'b4595' 'temp_110320_keep'};
% tlist = {'t12703' 'for_archiving'};
% agelist = [50.64 52.23];

blist = {'b4604' 'temp_110620_keep'};
tlist = {'t12773' 'for_archiving'};
agelist = [35.62 37.18];

QSMroot1 = ['/data/7T_hunt/' blist{1} '/QSM_longitudinal/TP1'];
QSMroot2 = ['/data/7T_hunt/' blist{1} '/QSM_longitudinal/TP2'];

%% register tp2 to tp1

for ii = 1:length(QSMNii_list)
    newName = strrep(QSMNii_list{ii},'_meanEcho','');
    cmd = sprintf(['flirt -in %s/%s.nii.gz -ref %s/%s.nii.gz -out %s/%s_tp2.nii.gz '...
        '-applyxfm -init %s/tp2to1.mat'],...
        QSMroot2, QSMNii_list{ii}, QSMroot1, QSMNii_list{ii}, QSMroot1, newName, QSMroot1);
    system(cmd);
end

%% load images

for qq = 1:length(QSMNii_list)
    newName = strrep(QSMNii_list{qq},'_meanEcho','');
    nii = load_nii([QSMroot1 '/' QSMNii_list{qq} '.nii.gz']);
    QSM_TP1 = double(nii.img);
    if qq == 9; QSM_TP1 = QSM_TP1/0.5684; end
    QSMmap(:,:,:,qq,1) = QSM_TP1;
    nii = load_nii([QSMroot1 '/' newName '_tp2.nii.gz']);
    QSM_TP2 = double(nii.img);
    if qq == 9; QSM_TP2 = QSM_TP2/0.5684; end
    QSMmap(:,:,:,qq,2) = QSM_TP2;
    
%     QSM_diff = (QSM_TP2 - QSM_TP1).*brainMask;
%     nii.img = QSM_diff;
%     save_nii(nii,[MNI_root '/' QSMNii_list{qq} '_' blist{1} '_diff.nii.gz']);
%     QSM_SD = std(cat(4,QSM_TP1,QSM_TP2),0,4).*brainMask;
%     nii.img = QSM_SD;
%     save_nii(nii,[MNI_root '/' QSMNii_list{qq} '_' blist{1} '_SD.nii.gz']);
end

%% read in ROI index

ROItxt = fileread('/working/lupolab/jingwen/090_Brain_ROIs_QSM_Atlas/ROIfromWord.txt');
ROItxt = strsplit(ROItxt, '\n');
ROItxt(cellfun('isempty', ROItxt)) = [];
ROItxt(strcmp(ROItxt, ' ')) = [];

ROIlist = cell(length(ROItxt),2);

for ii = 1:length(ROItxt)
    temp = strsplit(ROItxt{ii},' ');
    temp(cellfun('isempty', temp)) = [];
    if length(temp) > 1
        ROIlist{ii,1} = [temp{1:end-1}];
        ROIlist{ii,2} = str2double(temp{end});
    end
end

nii = load_nii([QSMroot1 '/QSM_atlas_ROI.nii.gz']);
brainROI = double(nii.img);

nii = load_nii(['/data/7T_hunt/' blist{1} '/' tlist{1} '/swan_qsm/HDBET_allQSM/brain_mask_HD.nii.gz']);
brainMask = double(nii.img);

%% plot maps

QSM_diff = QSMmap(:,:,:,:,2) - QSMmap(:,:,:,:,1);
QSM_SD = std(QSMmap,0,5);

sliceN = linspace(50,110,7);

dispMap(brainMask, QSM_diff(:,:,:,1:5), brainMask, sliceN, [-1 1]*0.05, jet);
export_fig([blist{1} '_diffQSM'], '-png','-transparent'); close;
pause(2);
dispMap(brainMask, QSM_diff(:,:,:,6:10), brainMask, sliceN, [-1 1]*0.05, jet);
export_fig([blist{1} '_diffQSM2'], '-png','-transparent'); close;

%% compute ROI medians

QSMstats_TP1 = QSMmedian(QSMmap(:,:,:,:,1), brainROI, brainMask, ROIlist, 1);
QSMstats_TP2 = QSMmedian(QSMmap(:,:,:,:,2), brainROI, brainMask, ROIlist, 1);

%% plot heatmap

QSMmedian_TP1 = reshape([QSMstats_TP1.median],[],10);
QSMmedian_TP2 = reshape([QSMstats_TP2.median],[],10);

QSM_diff_ROI = QSMmedian_TP2 - QSMmedian_TP1;
QSM_diff_SD = std(cat(3,QSMmedian_TP1,QSMmedian_TP2),0,3);

set(0,'defaultAxesFontSize',12);

QSMname = strrep(QSMNii_list,'_','-');
QSMname = strrep(QSMname,'QSM-','');
QSMname = strrep(QSMname,'-meanEcho','');
QSMname = strrep(QSMname,'-nonlinearTV','');

cmap1 = flip(colormap_RWG([-1 1],100));
r = ones(100, 1); g = linspace(1, 0, 100)'; 
cmap2 = [r g g];

figure('position', [100 100 800 400]);
ax(1) = subplot(211);
plotHeatMap(QSM_diff_ROI, cmap1, QSMname, [-0.05 0.05]);
title('QSM difference');
ax(2) = subplot(212);
plotHeatMap(QSM_diff_SD, cmap2, QSMname, [0 0.05]);
title('QSM SD accross time points');

colormap(ax(1), cmap1); colormap(ax(2), cmap2);
export_fig([blist{1} '_diff_ROI'], '-png','-transparent');

%% functions

function [crop] = cropMask(brain_mask)

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

end

function [] = plotHeatMap(data, cmap, QSMname, colorrange)

imagesc(data(1:209,:)',colorrange); colormap(cmap); colorbar; hold on;
plot([109.5 109.5],[0 11],'k','LineWidth',2);
plot([176.5 176.5],[0 11],'k','LineWidth',2);
plot([190.5 190.5],[0 11],'k','LineWidth',2);
ax = gca;
ax.YTick = 1:10;
ax.YTickLabel = QSMname(:);
ax.YTickLabelRotation = 0;
ax.XTick = floor([(1+109)/2 (109+176)/2 (176+191)/2 (191+210)/2]);
ax.XTickLabelRotation = 0;
ax.XTickLabel = {'Cortical GM','WM','BG','TH'};
title('QSM-Age Correlation Coefficient');

end

function [] = dispMap(underlay, overlay, mask, sliceN, colorrange, colormap)

Nqsm = size(overlay,4);

[crop] = cropMask(mask);
underlay = underlay(crop.X, crop.Y, crop.Z, :);
overlay = overlay(crop.X, crop.Y, crop.Z, :);
mask = mask(crop.X, crop.Y, crop.Z, :);

imA = underlay.*mask; imA = imA(:,:,sliceN-crop.Z(1)); imA = rot90(imA);
imA = repmat(imA(:,:),[Nqsm 1]);

imB = overlay.*repmat(mask,[1 1 1 Nqsm]); 
imB = imB(:,:,sliceN-crop.Z(1),:); imB = rot90(imB);
imB = reshape(permute(imB,[1 4 2 3]), size(imA));
imB(imB == 0) = nan;

imoverlay(imA(:,:),imB(:,:),colorrange,[],colormap,1);

end

function [QSMstats] = ...
    QSMmedian(QSMmaps, QSM_ROI, brain_mask, ROIlist, flag_erode)

if nargin < 5
    flag_erode = 0;
end
se = strel('disk',1);

QSM_ROI(brain_mask == 0) = 0;

% loop through QSM files
QSMstats(size(QSMmaps,4)) = struct;
for nn = 1:size(QSMmaps,4)
    
    QSMmap = QSMmaps(:,:,:,nn);

    for rr = 1:size(ROIlist,1)
        ROIindex = ROIlist{rr,2};
        ROIname = ROIlist(rr,1);
        
        ROImask = QSM_ROI == ROIindex;
        if flag_erode; ROImask = imerode(ROImask,se); end
        ROIdata = nonzeros(QSMmap(ROImask));
        ROIdata(isnan(ROIdata)) = [];
        QSMstats(nn).median(rr) = median(ROIdata);
        QSMstats(nn).ROIindex(rr) = ROIindex;
        QSMstats(nn).ROIname(rr) = ROIname;
    end
    
    ROIname = {'Thalamus_L'};
    ROImask = (QSM_ROI == 205) | (QSM_ROI == 207) | (QSM_ROI == 209) | (QSM_ROI == 211) ...
        | (QSM_ROI == 213) | (QSM_ROI == 215);
    if flag_erode; ROImask = imerode(ROImask,se); end
    ROIdata = nonzeros(QSMmap(ROImask));
    ROIdata(isnan(ROIdata)) = [];
    QSMstats(nn).median(size(ROIlist,1)+1) = median(ROIdata);
    QSMstats(nn).ROIindex(size(ROIlist,1)+1) = ROIindex;
    QSMstats(nn).ROIname(size(ROIlist,1)+1) = ROIname;
    
    ROIname = {'Thalamus_R'};
    ROImask = (QSM_ROI == 206) | (QSM_ROI == 208) | (QSM_ROI == 210) | (QSM_ROI == 212) ...
        | (QSM_ROI == 214) | (QSM_ROI == 216);
    if flag_erode; ROImask = imerode(ROImask,se); end
    ROIdata = nonzeros(QSMmap(ROImask));
    ROIdata(isnan(ROIdata)) = [];
    QSMstats(nn).median(size(ROIlist,1)+2) = median(ROIdata);
    QSMstats(nn).ROIindex(size(ROIlist,1)+2) = ROIindex;
    QSMstats(nn).ROIname(size(ROIlist,1)+2) = ROIname;
    
    ROIname = {'SubstantiaNigra_L'};
    ROImask = (QSM_ROI == 197) | (QSM_ROI == 199);
    if flag_erode; ROImask = imerode(ROImask,se); end
    ROIdata = nonzeros(QSMmap(ROImask));
    ROIdata(isnan(ROIdata)) = [];
    QSMstats(nn).median(size(ROIlist,1)+3) = median(ROIdata);
    QSMstats(nn).ROIindex(size(ROIlist,1)+3) = ROIindex;
    QSMstats(nn).ROIname(size(ROIlist,1)+3) = ROIname;
    
    ROIname = {'SubstantiaNigra_R'};
    ROImask = (QSM_ROI == 198) | (QSM_ROI == 200);
    if flag_erode; ROImask = imerode(ROImask,se); end
    ROIdata = nonzeros(QSMmap(ROImask));
    ROIdata(isnan(ROIdata)) = [];
    QSMstats(nn).median(size(ROIlist,1)+4) = median(ROIdata);
    QSMstats(nn).ROIindex(size(ROIlist,1)+4) = ROIindex;
    QSMstats(nn).ROIname(size(ROIlist,1)+4) = ROIname;
    
    ROIname = {'LatVen'};
    ROIdata = 1e-16;
    QSMstats(nn).median(size(ROIlist,1)+5) = median(ROIdata);
    QSMstats(nn).ROIindex(size(ROIlist,1)+5) = ROIindex;
    QSMstats(nn).ROIname(size(ROIlist,1)+5) = ROIname;
    
    ROIname = {'WholeBrain'};
    ROIdata = QSMmap(brain_mask > 0);
    ROIdata(isnan(ROIdata)) = [];
    QSMstats(nn).median(size(ROIlist,1)+6) = median(ROIdata);
    QSMstats(nn).ROIindex(size(ROIlist,1)+6) = ROIindex;
    QSMstats(nn).ROIname(size(ROIlist,1)+6) = ROIname;
    
end

end
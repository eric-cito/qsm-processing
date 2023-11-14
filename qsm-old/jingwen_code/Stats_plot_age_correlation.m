clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath('/home/jyao3/010_MATLAB_Utils/');
addpath('/home/jyao3/010_MATLAB_Utils/violinplot');
addpath('/home/jyao3/010_MATLAB_Utils/export_fig');
addpath('/home/jyao3/010_MATLAB_Utils/NIfTI');
addpath('/home/jyao3/010_MATLAB_Utils/imoverlay');
addpath('/home/jyao3/010_MATLAB_Utils/shadedErrorBar');

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220103.xlsx','Sheet','NoRep');
statusList = [T.status];

subjList = T.b_num; subjList(~strcmp(statusList,'HC')) = [];
examList = T.t_num; examList(~strcmp(statusList,'HC')) = [];
sexList = T.sex; sexList(~strcmp(statusList,'HC')) = [];
ageList = T.age; ageList(~strcmp(statusList,'HC')) = [];

matout_root = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/Data';

Mind = strcmp(sexList,'M');
Find = strcmp(sexList,'F');

fprintf('Age mean %.3f SD %.3f \n', mean(ageList), std(ageList));
fprintf('M %i F %i \n', sum(Mind), sum(Find));

%% load stats on ROI of interest

ROIlist = 1:213;
QSMfile_list = {'FANSI' ...
    'HDQSM' ...
    'iLSQR' ...
    'MEDI' ...
    'QSIP' ...
    'QSMGAN' ...
    'QSMnet+' ...
    'SSTGV' ...
    'SSTV' ...
    'STARQSM'};

meanROI = zeros(length(subjList), length(ROIlist), length(QSMfile_list));
stdROI = zeros(length(subjList), length(ROIlist), length(QSMfile_list));
medianROI = zeros(length(subjList), length(ROIlist), length(QSMfile_list));
madROI = zeros(length(subjList), length(ROIlist), length(QSMfile_list));
for ii = 1:length(subjList)
    
    exam_id = [subjList{ii} '_' examList{ii}];
    fprintf('# Reading subj %s \n', exam_id);
    dataPath = [matout_root '/' exam_id '_erode.mat'];
    
    load(dataPath, 'QSMstats');
    for QSMnum = 1:length(QSMfile_list)
        meanROI(ii,:,QSMnum) = table2array(QSMstats(QSMnum).QSMtable(ROIlist,'ROImean'));
        stdROI(ii,:,QSMnum) = table2array(QSMstats(QSMnum).QSMtable(ROIlist,'ROIstd'));
        medianROI(ii,:,QSMnum) = table2array(QSMstats(QSMnum).QSMtable(ROIlist,'ROImedian'));
        madROI(ii,:,QSMnum) = table2array(QSMstats(QSMnum).QSMtable(ROIlist,'ROImad'));
    end
    
end

ROIname = table2cell(QSMstats(1).QSMtable(ROIlist,'ROIname'));

% scale the QSMGAN
meanROI(:,:,6) = meanROI(:,:,6)/0.5684;
stdROI(:,:,6) = stdROI(:,:,6)/0.5684;
medianROI(:,:,6) = medianROI(:,:,6)/0.5684;
madROI(:,:,6) = madROI(:,:,6)/0.5684;

% reorder
QSMorder = [3 10 1 2 4 5 8 9 6 7];
QSMfile_list = QSMfile_list(QSMorder);
meanROI = meanROI(:,:,QSMorder);
stdROI = stdROI(:,:,QSMorder);
medianROI = medianROI(:,:,QSMorder);
madROI = madROI(:,:,QSMorder);

%% output age correlation stats

meanROI_R = zeros(length(ROIlist), length(QSMfile_list));
meanROI_P = zeros(length(ROIlist), length(QSMfile_list));
meanROI_slope = zeros(length(ROIlist), length(QSMfile_list));
medianROI_R = zeros(length(ROIlist), length(QSMfile_list));
medianROI_P = zeros(length(ROIlist), length(QSMfile_list));
medianROI_slope = zeros(length(ROIlist), length(QSMfile_list));
stdROI_R = zeros(length(ROIlist), length(QSMfile_list));
stdROI_P = zeros(length(ROIlist), length(QSMfile_list));
stdROI_slope = zeros(length(ROIlist), length(QSMfile_list));
for rr = 1:length(ROIlist)
    for ii = 1:10
        % mean
        data = squeeze(meanROI(:,rr,ii)) - squeeze(meanROI(:,212,ii)); 
        [R,P] = corr(ageList(:),data(:),'Type','Pearson');
        meanROI_R(rr,ii) = R;
        meanROI_P(rr,ii) = P;
        try
            F = fit(ageList(~isnan(data)),data(~isnan(data)),'poly1');
            meanROI_slope(rr,ii) = F.p1;
        catch
            meanROI_slope(rr,ii) = nan;
        end
        % median
        data = squeeze(medianROI(:,rr,ii)) - squeeze(medianROI(:,212,ii)); 
        [R,P] = corr(ageList(:),data(:),'Type','Pearson');
        medianROI_R(rr,ii) = R;
        medianROI_P(rr,ii) = P;
        try
            F = fit(ageList(~isnan(data)),data(~isnan(data)),'poly1');
            medianROI_slope(rr,ii) = F.p1;
        catch
            medianROI_slope(rr,ii) = nan;
        end
        % std
        data = squeeze(stdROI(:,rr,ii)); 
        [R,P] = corr(ageList(:),data(:),'Type','Pearson');
        stdROI_R(rr,ii) = R;
        stdROI_P(rr,ii) = P;
        try
            F = fit(ageList(~isnan(data)),data(~isnan(data)),'poly1');
            stdROI_slope(rr,ii) = F.p1;
        catch
            stdROI_slope(rr,ii) = nan;
        end
    end
end

%% sort correlation

set(0,'defaultAxesFontSize',12);
set(0,'defaultLineLineWidth',1);

median_R = mean(medianROI_R,2);
mad_R = std(medianROI_R,0,2);

invalidInd = isnan(median_R);

figure('position', [100 100 800 300]); 
plot([82.5 82.5], [-1 1], ':', 'Color', [1 1 1]*0.8); hold on;
plot([108.5 108.5], [-1 1], 'Color', [1 1 1]*0.8);
plot([176.5 176.5], [-1 1], 'Color', [1 1 1]*0.8);
plot([190.5 190.5], [-1 1], 'Color', [1 1 1]*0.8);
H = shadedErrorBar(1:length(median_R(~invalidInd)), median_R(~invalidInd), mad_R(~invalidInd),...
    'lineProps','-b'); 
H.mainLine.Color = [0 0.4470 0.7410];
H.edge(1).Color = 'none'; H.edge(2).Color = 'none';
H.patch.FaceColor = [0 0.4470 0.7410];
plot([0 210],[0 0],'k:');
xlim([0 206]); ylim([-0.6 0.8]);
ylabel('Correlation with age R');
xlabel('Brain segmentations');
export_fig('Age-Corr-Rmean', '-png','-transparent');

[median_Rsort,i] = sort(median_R);
ROIname_sort = ROIname(i);

T = table(ROIname_sort,median_Rsort);
% disp(T);

save('age_median_R.mat','medianROI_R');

%% plot heatmap

set(0,'defaultAxesFontSize',12);

QSMname = strrep(QSMfile_list,'_','-');
QSMname = strrep(QSMname,'QSM-','');
QSMname = strrep(QSMname,'-meanEcho','');

cmap1 = flip(colormap_RWG([-1 1],100));
r = ones(100, 1); g = linspace(1, 0, 100)'; 
cmap2 = [r g g];

% median
data = [medianROI_R]; data(isnan(data)) = 0;
dataP = [medianROI_P]; dataP(isnan(dataP)) = 1;
dataS = [medianROI_slope]; dataS(isnan(dataS)) = 0;

figure('position', [100 100 800 750]);
ax(1) = subplot(311);
plotHeatMap(data, cmap1, QSMname, [-1 1]);
title('Median QSM-Age Correlation Coefficient');
ax(2) = subplot(312);
plotHeatMap(dataP, flip(cmap2), QSMname, [0 0.005]);
title('Median QSM-Age Correlation P-value');
ax(3) = subplot(313);
plotHeatMap(dataS, flip(cmap2), QSMname, [-0.001 0.001]);
title('Median QSM-Age Linear Regression Slope');

colormap(ax(1), cmap1); colormap(ax(2), flip(cmap2)); colormap(ax(3), cmap1);
export_fig('Age-Corr-Median', '-png','-transparent');

% mean
data = [meanROI_R]; data(isnan(data)) = 0;
dataP = [meanROI_P]; dataP(isnan(dataP)) = 1;
dataS = [meanROI_slope]; dataS(isnan(dataS)) = 0;

figure('position', [100 100 800 750]);
ax(1) = subplot(311);
plotHeatMap(data, cmap1, QSMname, [-1 1]);
title('Mean QSM-Age Correlation Coefficient');
ax(2) = subplot(312);
plotHeatMap(dataP, flip(cmap2), QSMname, [0 0.005]);
title('Mean QSM-Age Correlation P-value');
ax(3) = subplot(313);
plotHeatMap(dataS, flip(cmap2), QSMname, [-0.001 0.001]);
title('Mean QSM-Age Linear Regression Slope');

colormap(ax(1), cmap1); colormap(ax(2), flip(cmap2)); colormap(ax(3), cmap1);
export_fig('Age-Corr-Mean', '-png','-transparent');

% std
data = [stdROI_R]; data(isnan(data)) = 0;
dataP = [stdROI_P]; dataP(isnan(dataP)) = 1;
dataS = [stdROI_slope]; dataS(isnan(dataS)) = 0;

figure('position', [100 100 800 750]);
ax(1) = subplot(311);
plotHeatMap(data, cmap1, QSMname, [-1 1]);
title('SD QSM-Age Correlation Coefficient');
ax(2) = subplot(312);
plotHeatMap(dataP, flip(cmap2), QSMname, [0 0.005]);
title('SD QSM-Age Correlation P-value');
ax(3) = subplot(313);
plotHeatMap(dataS, flip(cmap2), QSMname, [-0.001 0.001]);
title('SD QSM-Age Linear Regression Slope');

colormap(ax(1), cmap1); colormap(ax(2), flip(cmap2)); colormap(ax(3), cmap1);
export_fig('Age-Corr-SD', '-png','-transparent');

%% load MNI images
MNI_root = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/NIFTI_HC';
nii = load_nii([MNI_root '/MNI_T1w.nii.gz']);
Anat = double(nii.img);

nii = load_nii([MNI_root '/Whole_brain_segmentation_DGM.nii']);
brainROI = double(nii.img);

nii = load_nii([MNI_root '/MNI_GMWM.nii.gz']);
brainMask = double(nii.img);
% brainMask(brainMask == 0) = nan;

%% load ROI index

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

%% convert to map

sliceN = linspace(50,110,7);

Rmap = zeros([numel(Anat) length(QSMfile_list)]);
Pmap = zeros([numel(Anat) length(QSMfile_list)]);
Smap = zeros([numel(Anat) length(QSMfile_list)]);
for rr = 1:length(ROIlist)
    ROImask = brainROI == ROIlist{rr,2};
    for qq = 1:length(QSMfile_list)
        Rmap(ROImask,qq) = medianROI_R(rr,qq);
        Pmap(ROImask,qq) = medianROI_P(rr,qq);
        Smap(ROImask,qq) = medianROI_slope(rr,qq);
    end
end
Rmap = reshape(Rmap,[size(Anat) length(QSMfile_list)]);
Pmap = reshape(Pmap,[size(Anat) length(QSMfile_list)]);
Smap = reshape(Smap,[size(Anat) length(QSMfile_list)]);

%% display maps - R

dispMap(Anat, Rmap(:,:,:,1:5), brainMask, sliceN, [-1 1], jet);
export_fig('Age-Corr-ROI-R1', '-png','-transparent');
pause(1);
dispMap(Anat, Rmap(:,:,:,6:10), brainMask, sliceN, [-1 1], jet);
export_fig('Age-Corr-ROI-R2', '-png','-transparent');

%% display maps - P

dispMap(Anat, Pmap(:,:,:,1:5), brainMask, sliceN, [0 0.05], autumn);
export_fig('Age-Corr-ROI-P1', '-png','-transparent');
pause(1);
dispMap(Anat, Pmap(:,:,:,6:10), brainMask, sliceN, [0 0.05], autumn);
export_fig('Age-Corr-ROI-P2', '-png','-transparent');

%% display maps - slppe

dispMap(Anat, Smap(:,:,:,1:5), brainMask, sliceN, [-0.001 0.001], jet);
export_fig('Age-Corr-ROI-S1', '-png','-transparent');
pause(1);
dispMap(Anat, Smap(:,:,:,6:10), brainMask, sliceN, [-0.001 0.001], jet);
export_fig('Age-Corr-ROI-S2', '-png','-transparent');

%% voxelwise maps

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

Rmap = zeros([size(Anat) length(QSMfile_list)]);
Pmap = zeros([size(Anat) length(QSMfile_list)]);
Smap = zeros([size(Anat) length(QSMfile_list)]);
for qq = 1:length(QSMfile_list)
    nii = load_nii([MNI_root '/' QSMNii_list{qq} '_R.nii.gz']);
    tmp = double(nii.img);
    Rmap(:,:,:,qq) = tmp;
    nii = load_nii([MNI_root '/' QSMNii_list{qq} '_P.nii.gz']);
    tmp = double(nii.img);
    Pmap(:,:,:,qq) = tmp;
    nii = load_nii([MNI_root '/' QSMNii_list{qq} '_Slope.nii.gz']);
    tmp = double(nii.img);
    Smap(:,:,:,qq) = tmp;
end

%% display maps - R

dispMap(Anat, Rmap(:,:,:,1:5), brainMask, sliceN, [-1 1], jet);
export_fig('Age-Corr-Vox-R1', '-png','-transparent');
pause(1);
dispMap(Anat, Rmap(:,:,:,6:10), brainMask, sliceN, [-1 1], jet);
export_fig('Age-Corr-Vox-R2', '-png','-transparent');

%% display maps - P

dispMap(Anat, Pmap(:,:,:,1:5), brainMask, sliceN, [0 0.05], autumn);
export_fig('Age-Corr-Vox-P1', '-png','-transparent');
pause(1);
dispMap(Anat, Pmap(:,:,:,6:10), brainMask, sliceN, [0 0.05], autumn);
export_fig('Age-Corr-Vox-P2', '-png','-transparent');

%% display maps - slope

dispMap(Anat, Smap(:,:,:,1:5), brainMask, sliceN, [-0.001 0.001], jet);
export_fig('Age-Corr-Vox-S1', '-png','-transparent');
pause(1);
dispMap(Anat, Smap(:,:,:,6:10), brainMask, sliceN, [-0.001 0.001], jet);
export_fig('Age-Corr-Vox-S2', '-png','-transparent');

%% functions

function [] = dispMap(underlay, overlay, mask, sliceN, colorrange, colormap)

Nqsm = size(overlay,4);

imA = underlay.*mask; imA = imA(:,:,sliceN); imA = rot90(imA);
imA = repmat(imA(:,:),[Nqsm 1]);

imB = overlay.*repmat(mask,[1 1 1 Nqsm]); 
imB = imB(:,:,sliceN,:); imB = rot90(imB);
imB = reshape(permute(imB,[1 4 2 3]), size(imA));
imB(imB == 0) = nan;

imoverlay(imA(:,:),imB(:,:),colorrange,[],colormap,1);

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
clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/010_MATLAB_Utils/'));

dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/SubjImgDataBG_ANTS_iLSQR_1004.mat';
img_root = '/working/lupolab/jingwen/001_QSM/temp';

%% Load data

load(dataPath);

indHC = strcmp({HD_BGanalysis.group},'HC');
indPM = strcmp({HD_BGanalysis.group},'PM');
indEM = strcmp({HD_BGanalysis.group},'EM');
indMan = strcmp({HD_BGanalysis.group},'Manifest');

statusList = nan(1,length(HD_BGanalysis));

age = [HD_BGanalysis.age];
sex = {HD_BGanalysis.sex};
CAG = [HD_BGanalysis.CAG];
volList = zeros(3,length(HD_BGanalysis));

statusList(indHC) = 1;
statusList(indPM) = 2;
statusList(indEM) = 3;
statusList(indMan) = 4;
statusList(CAG < 36) = nan;

%% demographics

CAPS = age.*(CAG-33.66)/432.3326;
CAPS(isnan(statusList)) = nan;

CAG = [HD_BGanalysis.CAG];
TMS = [HD_BGanalysis.TMS];
DCL = [HD_BGanalysis.DCL];
TFC = [HD_BGanalysis.TFC];

DART = [HD_BGanalysis.DART];
DARTtiming = [HD_BGanalysis.DARTtiming];
Flanker = [HD_BGanalysis.Flanker];
Match = [HD_BGanalysis.Match];
SetShift = [HD_BGanalysis.SetShift];

ageOnset = zeros(1,length(HD_BGanalysis));
for ii = 1:length(HD_BGanalysis)
    if isnan(CAG(ii))
        ageOnset(ii) = nan;
    else
        ageOnset(ii) = medianYearOnset(CAG(ii));
    end
end

AOO = ageOnset;
YTO = AOO - age;

% YTO = exp(4.4196 - 2.8102*CAPS);

indPM = statusList == 2;
indEM = statusList == 3;
indMan = statusList == 4;
indPMfar = indPM & YTO > 15;
indPMnear = indPM & YTO <= 15;

statusList5 = statusList;
statusList5(indHC) = 0;
statusList5(indPMfar) = 1;
statusList5(indPMnear) = 2;
statusList5(indEM | indMan) = 3;

statusList3 = statusList;
statusList3(indHC) = 0;
statusList3(indPM) = 1;
statusList3(indEM | indMan) = 2;

%% load QSM data

nii = load_nii(['/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/allQSM_ANTS.nii.gz']);
allQSM = double(nii.img);

nii = load_nii(['/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/allFA_ANTS.nii.gz']);
allFA = double(nii.img);
allFA(:,:,38) = nan;

nii = load_nii(['/working/lupolab/jingwen/001_QSM/03_QSM_HD/Atlas/Seg_Subcortical.nii.gz']);
ROImask = double(nii.img);

nii = load_nii(['/working/lupolab/jingwen/001_QSM/03_QSM_HD/Atlas/MNI_T1_brain.nii.gz']);
Brain_img = double(nii.img);
Brain_mask = double(nii.img > 0);

%% bivariate histogram

Nsubj = length(HD_BGanalysis);
allQSM_data = reshape(allQSM,[],Nsubj);
allFA_data = reshape(allFA,[],Nsubj);
indInclude = statusList3 > -1; % 0 HC 1 PMfar 2 PMnear 3 EM+M

set(0,'DefaultAxesFontSize', 14);

figure('position', [100 100 1200 800]);

ROItarget = ROImask == 1 | ROImask == 2;
nBins = 64;
subplot(3,4,1);
[r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, statusList3 == 0, nBins);
xlim([0 0.4]); ylim([-0.07 0.1]);
ylabel('Susceptibility (ppm)'); % xlabel('FA'); 
title('Healthy control');
text(0.01,-0.06,sprintf('r = %.3f',r),'Color','w');
subplot(3,4,2);
[r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, statusList3 == 1, nBins);
xlim([0 0.4]); ylim([-0.07 0.1]);
title('Premanifest HD');
% xlabel('FA'); % ylabel('Susceptibility (ppm)');
text(0.01,-0.06,sprintf('r = %.3f',r),'Color','w');
subplot(3,4,3);
[r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, statusList3 == 2, nBins);
xlim([0 0.4]); ylim([-0.07 0.1]);
title('Manifest HD');
% xlabel('FA'); % ylabel('Susceptibility (ppm)');
text(0.01,-0.06,sprintf('r = %.3f',r),'Color','w');

subplot(3,4,4);
nBins = 20;
plotContour(allQSM_data, allFA_data, ROItarget, indInclude, nBins, statusList3);
xlim([0 0.4]); ylim([-0.07 0.1]);
legend({'','','','Healthy control','Premanifest HD','Manifest HD'},'location','best','FontSize',10);
% xlabel('FA'); % ylabel('Susceptibility (ppm)');

ROItarget = ROImask == 3 | ROImask == 4;
nBins = 64;
subplot(3,4,5);
[r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, statusList3 == 0, nBins);
xlim([0.05 0.45]); ylim([-0.07 0.12]);
ylabel('Susceptibility (ppm)'); % xlabel('FA'); 
text(0.06,-0.06,sprintf('r = %.3f',r),'Color','w');
subplot(3,4,6);
[r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, statusList3 == 1, nBins);
xlim([0.05 0.45]); ylim([-0.07 0.12]);
% xlabel('FA'); % ylabel('Susceptibility (ppm)');
text(0.06,-0.06,sprintf('r = %.3f',r),'Color','w');
subplot(3,4,7);
[r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, statusList3 == 2, nBins);
xlim([0.05 0.45]); ylim([-0.07 0.12]);
% xlabel('FA'); % ylabel('Susceptibility (ppm)');
text(0.06,-0.06,sprintf('r = %.3f',r),'Color','w');

subplot(3,4,8);
nBins = 20;
plotContour(allQSM_data, allFA_data, ROItarget, indInclude, nBins, statusList3);
xlim([0.05 0.45]); ylim([-0.07 0.12]);
% xlabel('FA');

ROItarget = ROImask == 5 | ROImask == 6 | ROImask == 7 | ROImask == 8;
nBins = 40;
subplot(3,4,9);
[r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, statusList3 == 0, nBins);
xlim([0.05 0.75]); ylim([-0.07 0.25]);
ylabel('Susceptibility (ppm)'); xlabel('FA'); 
text(0.06,-0.05,sprintf('r = %.3f',r),'Color','w');
subplot(3,4,10);
[r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, statusList3 == 1, nBins);
xlim([0.05 0.75]); ylim([-0.07 0.25]);
xlabel('FA'); % ylabel('Susceptibility (ppm)');
text(0.06,-0.05,sprintf('r = %.3f',r),'Color','w');
subplot(3,4,11);
[r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, statusList3 == 2, nBins);
xlim([0.05 0.75]); ylim([-0.07 0.25]);
xlabel('FA'); % ylabel('Susceptibility (ppm)');
text(0.06,-0.05,sprintf('r = %.3f',r),'Color','w');

subplot(3,4,12);
nBins = 20;
plotContour(allQSM_data, allFA_data, ROItarget, indInclude, nBins, statusList3);
xlim([0.05 0.75]); ylim([-0.07 0.25]);
xlabel('FA');

% pause(1); export_fig([img_root '/BiHist'], '-png','-transparent');

%% DN

figure('position', [100 100 1200 800]);

ROItarget = ROImask == 17 | ROImask == 18;
nBins = 64;
subplot(3,4,1);
[r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, statusList3 == 0, nBins);
xlim([0 0.6]); ylim([-0.07 0.1]);
ylabel('Susceptibility (ppm)'); % xlabel('FA'); 
title('Healthy control');
text(0.01,-0.06,sprintf('r = %.3f',r),'Color','w');
subplot(3,4,2);
[r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, statusList3 == 1, nBins);
xlim([0 0.6]); ylim([-0.07 0.1]);
title('Premanifest HD');
% xlabel('FA'); % ylabel('Susceptibility (ppm)');
text(0.01,-0.06,sprintf('r = %.3f',r),'Color','w');
subplot(3,4,3);
[r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, statusList3 == 2, nBins);
xlim([0 0.6]); ylim([-0.07 0.1]);
title('Manifest HD');
% xlabel('FA'); % ylabel('Susceptibility (ppm)');
text(0.01,-0.06,sprintf('r = %.3f',r),'Color','w');

subplot(3,4,4);
nBins = 20;
plotContour(allQSM_data, allFA_data, ROItarget, indInclude, nBins, statusList3);
xlim([0 0.6]); ylim([-0.07 0.1]);
% xlabel('FA'); % ylabel('Susceptibility (ppm)');

%% 3D rendering of ROIs

figure('position', [100 100 500 500]);
render3d(Brain_img, ROImask == 1 | ROImask == 2);
pause(1); export_fig([img_root '/CNmask'], '-png','-transparent');

figure('position', [100 100 500 500]);
render3d(Brain_img, ROImask == 3 | ROImask == 4);
pause(1); export_fig([img_root '/PUmask'], '-png','-transparent');

figure('position', [100 100 500 500]);
render3d(Brain_img, ROImask == 5 | ROImask == 6 | ROImask == 7 | ROImask == 8);
pause(1); export_fig([img_root '/GPmask'], '-png','-transparent');

%% DN

figure('position', [100 100 500 500]);
render3d(Brain_img, ROImask == 17 | ROImask == 18);
pause(1); export_fig([img_root '/DNmask'], '-png','-transparent');

%% Define line

set(0,'DefaultAxesFontSize', 14);

yind = 100:130;
xind = 120;

figure('position', [100 100 1200 600]);
ax = subplot(2,4,1);
imagesc(rot90(Brain_img(:,:,74))); axis equal tight off; hold on; colormap(ax,'gray');
plot(yind,(size(Brain_img,2)-xind)*ones(size(yind)),'k-');

subplot(2,4,2);
QSMdata = squeeze(allQSM(yind,xind,74,:));
ROIdata = squeeze(ROImask(yind,xind,74));
plotLineProfile(ROIdata, QSMdata, statusList3, [-0.15 0.28]);

subplot(2,4,3);
FAdata = squeeze(allFA(yind,xind,74,:));
ROIdata = squeeze(ROImask(yind,xind,74));
plotLineProfile(ROIdata, FAdata, statusList3, [0 1]);
ylabel('FA'); legend('off');

subplot(2,4,4);
MDdata = squeeze(allMD(yind,xind,74,:))*1e3;
ROIdata = squeeze(ROImask(yind,xind,74));
plotLineProfile(ROIdata, MDdata, statusList3, [0 0.8]);
ylabel('MD (x10^{-3} mm^2/s)'); legend('off');

yind = 100:130;
xind = 140;

ax = subplot(2,4,5);
imagesc(rot90(Brain_img(:,:,74))); axis equal tight off; hold on; colormap(ax,'gray');
plot(yind,(size(Brain_img,2)-xind)*ones(size(yind)),'k-');

subplot(2,4,6);
QSMdata = squeeze(allQSM(yind,xind,74,:));
ROIdata = squeeze(ROImask(yind,xind,74));
plotLineProfile(ROIdata, QSMdata, statusList3, [-0.1 0.1]);
legend('off');

subplot(2,4,7);
FAdata = squeeze(allFA(yind,xind,74,:));
ROIdata = squeeze(ROImask(yind,xind,74));
plotLineProfile(ROIdata, FAdata, statusList3, [0 0.6]);
ylabel('FA'); legend('off');

subplot(2,4,8);
MDdata = squeeze(allMD(yind,xind,74,:))*1e3;
ROIdata = squeeze(ROImask(yind,xind,74));
plotLineProfile(ROIdata, MDdata, statusList3, [0.2 1.2]);
ylabel('MD (x10^{-3} mm^2/s)'); legend('off');

% imagesc(rot90(ROImask(:,:,74))); axis equal tight off; hold on;
% plot(yind,(size(Brain_img,2)-xind)*ones(size(yind)),'k-');

% imagesc(allQSM(:,:,74,1)); axis equal tight off; hold on;
% plot(xind*ones(size(yind)),yind,'k-');

pause(1); export_fig([img_root '/LineProfile'], '-png','-transparent');

%% helper function

function [] = plotLineProfile(ROIdata, QSMdata, statusList3, yrange)

ind0 = find(ROIdata == 0,1)-0.5;
ind1 = find(ROIdata == 8,1)-0.5;
ind2 = find(ROIdata == 6,1)-0.5;
ind3 = find(ROIdata == 4,1)-0.5;
ind4 = find(ROIdata == 4,1,'last')+0.5;

meanQSM_HC = mean(QSMdata(:,statusList3 == 0),2);
sdQSM_HC = std(QSMdata(:,statusList3 == 0),0,2);
meanQSM_PM = mean(QSMdata(:,statusList3 == 1),2);
sdQSM_PM = std(QSMdata(:,statusList3 == 1),0,2);
meanQSM_EM = mean(QSMdata(:,statusList3 == 2),2);
sdQSM_EM = std(QSMdata(:,statusList3 == 2),0,2);

% plot(QSMdata(:,statusList3 == 0),'-','Color',[1 1 1]*0.9); hold on;
% ([1 1 1]-[0 0.4470 0.7410])*0.75 + [0 0.4470 0.7410]
% plot(QSMdata(:,statusList3 == 1),'-','Color',[1 1 1]*0.9);
% ([1 1 1]-[0.8500 0.3250 0.0980])*0.75 + [0.8500 0.3250 0.0980]
% plot(QSMdata(:,statusList3 == 2),'-','Color',[1 1 1]*0.9);
% ([1 1 1]-[0.9290 0.6940 0.1250])*0.75 + [0.9290 0.6940 0.1250]
p = fill([1:length(meanQSM_HC) length(meanQSM_HC):-1:1],...
    [meanQSM_HC-sdQSM_HC; meanQSM_HC(end:-1:1)+sdQSM_HC(end:-1:1)],'k'); hold on;
p.FaceColor = [0 0.4470 0.7410]; 
p.FaceAlpha = 0.1;
p.EdgeColor = 'none'; 

p = fill([1:length(meanQSM_PM) length(meanQSM_PM):-1:1],...
    [meanQSM_PM-sdQSM_PM; meanQSM_PM(end:-1:1)+sdQSM_PM(end:-1:1)],'k');
p.FaceColor = [0.8500 0.3250 0.0980]; 
p.FaceAlpha = 0.1;
p.EdgeColor = 'none';

p = fill([1:length(meanQSM_EM) length(meanQSM_EM):-1:1],...
    [meanQSM_EM-sdQSM_EM; meanQSM_EM(end:-1:1)+sdQSM_EM(end:-1:1)],'k');
p.FaceColor = [0.9290 0.6940 0.1250]; 
p.FaceAlpha = 0.1;
p.EdgeColor = 'none';

l(1) = plot(meanQSM_HC,'-','Color',[0 0.4470 0.7410], 'LineWidth',1); hold on;
l(2) = plot(meanQSM_PM,'-','Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
l(3) = plot(meanQSM_EM,'-','Color',[0.9290 0.6940 0.1250], 'LineWidth',1);

if ~isempty(ind0); plot([1 1]*ind0,yrange,':','Color',[1 1 1]*0); end
if ~isempty(ind1); plot([1 1]*ind1,yrange,':','Color',[1 1 1]*0); end
if ~isempty(ind2); plot([1 1]*ind2,yrange,':','Color',[1 1 1]*0); end
if ~isempty(ind3); plot([1 1]*ind3,yrange,':','Color',[1 1 1]*0); end
if ~isempty(ind4); plot([1 1]*ind4,yrange,':','Color',[1 1 1]*0); end
axis tight;
legend(l,{'Healthy control','Premanifest HD','Manifest HD'},'box','off','location','best','FontSize',10);
ylabel('Susceptibility (ppm)');
xlabel('Line profile (mm)');

end

function [] = render3d(Brain_img, ROImask)

plotMask = smooth3(Brain_img > 0,'gaussian',[5 5 5],2);
hisoBG = patch(isosurface(plotMask,0.1),...
    'FaceColor',[0.9 0.9 0.9], ...
    'FaceAlpha',0.2, ...
    'EdgeColor','none');
isonormals(plotMask,hisoBG); hold on;

plotCNMask = smooth3(ROImask,'gaussian',[5 5 5],2);
hisoCN = patch(isosurface(plotCNMask,0.1),...
    'FaceColor', [0.3010 0.7450 0.9330], ...
    'EdgeColor','none');
isonormals(plotCNMask,hisoCN);

view(30,0);
axis tight equal off

lightangle(45,45);
lighting gouraud

end

function [] = plotContour(allQSM_data, allFA_data, ROItarget, indInclude, nBins, statusList3)

allQSM_PU = allQSM_data(ROItarget, indInclude);
allFA_PU = allFA_data(ROItarget, indInclude);

indInclude = allQSM_PU ~=0 & allFA_PU ~=0;

allQSM_PU(allQSM_PU == 0) = nan;
allFA_PU(allFA_PU == 0) = nan;

[contValues, FAedge, QSMedge] = ...
    bihistContourValues(allFA_PU, allQSM_PU, nBins, statusList3);

medianQSM_PU = nanmedian(allQSM_PU,1);
medianFA_PU = nanmedian(allFA_PU,1);

[QSMcont, FAcont] = ndgrid(QSMedge(1:end-1)+(QSMedge(2)-QSMedge(1))/2, ...
    FAedge(1:end-1)+(FAedge(2)-FAedge(1))/2);

nContour = 6;

contour(FAcont, QSMcont, contValues(:,:,1)',nContour,'-','Color',[0 0.4470 0.7410]); hold on
contour(FAcont, QSMcont, contValues(:,:,2)',nContour,'-','Color',[0.8500 0.3250 0.0980]);
contour(FAcont, QSMcont, contValues(:,:,3)',nContour,'-','Color',[0.9290 0.6940 0.1250]);
scatter(medianFA_PU(statusList3 == 0), medianQSM_PU(statusList3 == 0), 50, ...
    'o', 'filled', 'MarkerFaceAlpha', 0.75);
scatter(medianFA_PU(statusList3 == 1), medianQSM_PU(statusList3 == 1), 50, ...
    'o', 'filled', 'MarkerFaceAlpha', 0.75);
scatter(medianFA_PU(statusList3 == 2), medianQSM_PU(statusList3 == 2), 50, ...
    'o', 'filled', 'MarkerFaceAlpha', 0.75);

end

function [r,p] = plotBiHist(allQSM_data, allFA_data, ROItarget, indInclude, nBins)

allQSM_PU = allQSM_data(ROItarget, indInclude);
allFA_PU = allFA_data(ROItarget, indInclude);

indInclude = allQSM_PU ~=0 & allFA_PU ~=0;

[r,p] = corr(allFA_PU(indInclude),allQSM_PU(indInclude),'Type','Pearson');

QSMedge = linspace(-0.15, 0.25,nBins);
% linspace(prctile(allQSM_PU(:),0),prctile(allQSM_PU(:),99.5),nBins);
FAedge = linspace(0,1,nBins);
% linspace(prctile(allFA_PU(:),0),prctile(allFA_PU(:),99),nBins);

h = histogram2(allFA_PU(indInclude),allQSM_PU(indInclude),FAedge,QSMedge,...
    'DisplayStyle','tile','ShowEmptyBins','on','Normalization','probability',...
    'EdgeAlpha',0,'FaceAlpha',1); hold on;

end

function [contValues, FAedge, QSMedge] = bihistContourValues(allFA_PU, allQSM_PU, nBins, statusList)

status = unique(statusList);
status = sort(status);

QSMedge = linspace(prctile(allQSM_PU(:),0),prctile(allQSM_PU(:),99.5),nBins);
FAedge = linspace(0,prctile(allFA_PU(:),99),nBins);
% linspace(prctile(allFA_PU(:),0),prctile(allFA_PU(:),99),nBins);

contValues = zeros(length(FAedge)-1, length(QSMedge)-1, length(status));

for ii = 1:length(status)
    h1 = histogram2(allFA_PU(:,statusList == status(ii)),allQSM_PU(:,statusList == status(ii)),...
        FAedge,QSMedge,'Normalization','probability');
    contValues(:,:,ii) = h1.Values;
end

end

function [ageOnset] = medianYearOnset(CAG)

syms f(x)
f(x) = (1+exp(pi/sqrt(3)*(-21.54-exp(9.56-0.146*CAG)+x)./(sqrt(35.55+exp(17.72-0.327*CAG))))).^-1 - 0.5;

tmp = vpasolve(f);
ageOnset = double(tmp);

end

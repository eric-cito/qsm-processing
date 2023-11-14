clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/'));

% dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/SubjImgDataBG_MNI_QSM.mat';
% dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/SubjImgDataBG_ANTS_iLSQR_1004.mat';
% dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/SubjImgDataBG_ANTS_iQSM_1228_erode.mat';
dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/SubjImgDataBG_ANTS_iQSM_0213_erode.mat';
img_root = '/working/lupolab/jingwen/001_QSM/temp';

%% Load data

load(dataPath);

% HD_BGanalysis([HD_BGanalysis.CAG] < 39) = [];

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

fprintf('=== HC === \n');
printDemo(indHC, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

indPM = statusList == 2;

fprintf('=== PM === \n');
printDemo(indPM, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

indEM = statusList == 3;

fprintf('=== EM === \n');
printDemo(indEM, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

indMan = statusList == 4;

fprintf('=== Manifest === \n');
printDemo(indMan, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

fprintf('=== EM+Manifest === \n');
printDemo(indEM | indMan, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

indPMfar = indPM & YTO > 15;
indPMnear = indPM & YTO <= 15;

fprintf('=== PM far === \n');
printDemo(indPMfar, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

fprintf('=== PM near === \n');
printDemo(indPMnear, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

statusList5 = statusList;
statusList5(indHC) = 0;
statusList5(indPMfar) = 1;
statusList5(indPMnear) = 2;
statusList5(indEM | indMan) = 3;

statusList3 = statusList;
statusList3(indHC) = 0;
statusList3(indPM) = 1;
statusList3(indEM | indMan) = 2;

%% histogram of YTO

figure;
histogram(YTO(indEM | indMan), [-15:5:30]); hold on;
histogram(YTO(indPMfar), [-15:5:30]);
histogram(YTO(indPMnear), [-15:5:30]);
legend({'Manifest', 'PM far', 'PM near'}, 'location', 'best');
xlabel('YTO (year)');

%% remove bad data

HD_BGanalysis(68).imData.FA(:) = nan;
HD_BGanalysis(68).imData.MD(:) = nan;
HD_BGanalysis(68).imData.RD(:) = nan;

%% Imaging metrics

Nroi = size(HD_BGanalysis(1).imData,1)/2;

fieldname = 'Volume';
[volMat, volMatL, volMatR] = extractMetric(HD_BGanalysis, fieldname);

fieldname = 'QSMmedian';
[qsmMat, qsmMatL, qsmMatR] = extractMetric(HD_BGanalysis, fieldname);

fieldname = 'QSMmad';
[qsmSDMat, qsmSDMatL, qsmSDMatR] = extractMetric(HD_BGanalysis, fieldname);

fieldname = 'FA';
[faMat, faMatL, faMatR] = extractMetric(HD_BGanalysis, fieldname);

fieldname = 'MD';
[mdMat, mdMatL, mdMatR] = extractMetric(HD_BGanalysis, fieldname);

fieldname = 'RD';
[rdMat, rdMatL, rdMatR] = extractMetric(HD_BGanalysis, fieldname);

ROI_name = {'CN','PU','GPe','GPi','TH','RN','SN','STN','DN'};

%% Age correction

qsmMatAC = qsmMat;
volMatAC = volMat;
faMatAC = faMat;
mdMatAC = mdMat;
rdMatAC = rdMat;

slope = zeros(Nroi,5);
rAge = zeros(Nroi,5);
pAge = zeros(Nroi,5);
for ii = 1:Nroi
    [qsmMatAC(:,ii), slope(ii,1), rAge(ii,1), pAge(ii,1)] = ageCorr(age, qsmMat(:,ii), indHC);
    [volMatAC(:,ii), slope(ii,2), rAge(ii,2), pAge(ii,2)] = ageCorr(age, volMat(:,ii), indHC);
    [faMatAC(:,ii), slope(ii,3), rAge(ii,3), pAge(ii,3)] = ageCorr(age, faMat(:,ii), indHC);
    [mdMatAC(:,ii), slope(ii,4), rAge(ii,4), pAge(ii,4)] = ageCorr(age, mdMat(:,ii), indHC);
    [rdMatAC(:,ii), slope(ii,5), rAge(ii,5), pAge(ii,5)] = ageCorr(age, rdMat(:,ii), indHC);
end

%% group comparison

set(0,'DefaultAxesFontSize', 14);

statusListPlot = statusList3;
indROI = [1:7 9];

data = qsmMatAC(:,indROI);
figure('position', [100 0 1200 900]);
subaxis(2,2,1,'SpacingHoriz', 0.06);
[~, pCorrQSM] = plotANOVA(data, ROI_name(indROI), statusListPlot, [-0.07 0.33]);
ylabel('Susceptibility (ppm)'); ylim([-0.07 0.33]);
legend({'Healthy control','Premanifest HD', 'Manifest HD'});

data = volMatAC(:,indROI);
subaxis(2,2,2);
[~, pCorrVol] = plotANOVA(data, ROI_name(indROI), statusListPlot, [0 7.2]);
ylabel('Corrected Volume (ml)'); ylim([0 7.2]);

axes('pos', [0.78 0.68 0.128 0.215]);
data = volMatAC(:,[3 4 6 7]);
plotANOVA_noCorr(data, ROI_name([3 4 6 7]), statusListPlot, [0 1.5], 1);
set(gca,'FontSize', 10);

data = faMat(:,indROI);
subaxis(2,2,3);
[~, pCorrFA] = plotANOVA(data, ROI_name(indROI), statusListPlot, [0 0.75]);
ylabel('FA'); ylim([0 0.75]);

data = mdMat(:,indROI)*1e3;
subaxis(2,2,4);
[~, pCorrMD] = plotANOVA(data, ROI_name(indROI), statusListPlot, [0 0.95]);
ylabel('MD (x10^{-3} mm^2/s)'); ylim([0 0.95]);

pause(1); export_fig([img_root '/All_groupComp1'], '-png','-transparent');

%% group comparison - finer groups

set(0,'DefaultAxesFontSize', 14);

statusListPlot = statusList5;

data = qsmMatAC(:,indROI);
figure('position', [100 0 1200 900]);
subaxis(2,2,1,'SpacingHoriz', 0.06);
[~, pCorrQSM, chiANOVAQSM] = plotANOVA(data, ROI_name(indROI), statusListPlot, [-0.07 03325]);
ylabel('Susceptibility (ppm)'); ylim([-0.07 0.33]);
legend({'Healthy control','Premanifest HD far', 'Premanifest HD near', 'Manifest HD'});

data = volMatAC(:,indROI);
subaxis(2,2,2);
[~, pCorrVol, chiANOVAVol] = plotANOVA(data, ROI_name(indROI), statusListPlot, [0 7.2]);
ylabel('Corrected Volume (ml)'); ylim([0 7.2]);

axes('pos', [0.78 0.68 0.128 0.215]);
data = volMatAC(:,[3 4 6 7]);
plotANOVA_noCorr(data, ROI_name([3 4 6 7]), statusListPlot, [0 1.5], 1);
set(gca,'FontSize', 10);

data = faMat(:,indROI);
subaxis(2,2,3);
[~, pCorrFA, chiANOVAFA] = plotANOVA(data, ROI_name(indROI), statusListPlot, [0 0.75]);
ylabel('FA'); ylim([0 0.75]);

data = mdMat(:,indROI)*1e3;
subaxis(2,2,4);
[~, pCorrMD, chiANOVAMD] = plotANOVA(data, ROI_name(indROI), statusListPlot, [0 0.95]);
ylabel('MD (x10^{-3} mm^2/s)'); ylim([0 0.95]);

pause(1); export_fig([img_root '/All_groupComp2'], '-png','-transparent');

%% create table

T = table(chiANOVAQSM', pCorrQSM', chiANOVAVol', pCorrVol', ...
    chiANOVAFA', pCorrFA', chiANOVAMD', pCorrMD');

%% clustering

indvalid = ~isnan(faMat(:,9));

ROIind = [1:9];
VarMRI = [qsmMatAC(:,ROIind) volMatAC(:,ROIind) faMat(:,ROIind) mdMat(:,ROIind)];

[coeff,score,latent,tsquared,explained] = pca(VarMRI);

figure;
scatter(score(indHC,1), score(indHC,2), 'o'); hold on;
scatter(score(indPMfar,1), score(indPMfar,2), 'o');
scatter(score(indPMnear,1), score(indPMnear,2), 'o');
scatter(score(indEM,1), score(indEM,2), 'o');
scatter(score(indMan,1), score(indMan,2), 'o');

figure;
boxchart(statusList5,score(:,1));

%% TSNE

Y = tsne(VarMRI);
gscatter(Y(:,1),Y(:,2),statusList5(indvalid));

%% test

% nii = load_untouch_nii('/data/7T_hunt/temp_052022_keep/for_archiving/swan_qsm/HDBET_allQSM/FSseg/MD_regT1.nii.gz');
% mdMap = nii.img;
% nii = load_untouch_nii('/data/7T_hunt/temp_052022_keep/for_archiving/swan_qsm/HDBET_allQSM/FSseg/Seg_ANTS_manual.nii.gz');
% mask = nii.img;
% 
% histogram(mdMap(mask > 1.5 & mask < 2.5));

%% helper function

function [] = printDemo(indEM, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC)

fprintf('Age mean %.3f SD %.3f \n', mean(age(indEM)), std(age(indEM)));
fprintf('M %i F %i \n', sum(strcmp(sex(indEM),'M')), sum(strcmp(sex(indEM),'F')));
fprintf('CAG mean %.3f SD %.3f \n', mean(CAG(indEM)), std(CAG(indEM)));
fprintf('CAG range %.3f - %.3f \n', min(CAG(indEM)), max(CAG(indEM)));
fprintf('CAPS median %.3f range %.3f - %.3f \n', median(CAPS(indEM)), min(CAPS(indEM)), max(CAPS(indEM)));
fprintf('AOO median %.3f range %.3f - %.3f \n', median(AOO(indEM)), min(AOO(indEM)), max(AOO(indEM)));
fprintf('YTO median %.3f range %.3f - %.3f \n', median(YTO(indEM)), min(YTO(indEM)), max(YTO(indEM)));
fprintf('TMS median %.3f range %.3f - %.3f \n', median(TMS(indEM)), min(TMS(indEM)), max(TMS(indEM)));
fprintf('DCL median %.3f range %.3f - %.3f \n', median(DCL(indEM)), min(DCL(indEM)), max(DCL(indEM)));
fprintf('TFC median %.3f range %.3f - %.3f \n', median(TFC(indEM)), min(TFC(indEM)), max(TFC(indEM)));

end

function [ageOnset] = medianYearOnset(CAG)

syms f(x)
f(x) = (1+exp(pi/sqrt(3)*(-21.54-exp(9.56-0.146*CAG)+x)./(sqrt(35.55+exp(17.72-0.327*CAG))))).^-1 - 0.5;

tmp = vpasolve(f);
ageOnset = double(tmp);

end

function [dataAC, slope, r, p] = ageCorr(age, data, indHC)

indValid = ~isnan(data) & indHC';

F = fit(age(indValid)',data(indValid),'poly1');
dataAC = data - F.p1*(age'-median(age));

slope = F.p1;

[r,p] = corr(age(indValid)',data(indValid));

end

function [matAll, matL, matR] = extractMetric(HD_BGanalysis, fieldname)

matL = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2);
matR = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2);
for ii = 1:length(HD_BGanalysis)
    matL(ii,:) = HD_BGanalysis(ii).imData.(fieldname)(1:2:end);
    matR(ii,:) = HD_BGanalysis(ii).imData.(fieldname)(2:2:end);
end
matAll = (matL + matR)/2;

end

function [pCorr, pANOVACorr, chiANOVA] = plotANOVA(data, ROI_name, statusList3, ylines, statsFlag)

if nargin < 5
    statsFlag = 1;
end

ymax = max(data(:)); ymin = min(data(:)); yrange = ymax - ymin;
colorList = repmat(statusList3',[1 length(ROI_name)]);
roiList = repmat(1:length(ROI_name), [size(data,1) 1]);

b = boxchart(roiList(:), data(:), 'GroupByColor', colorList(:), ...
    'MarkerStyle', '.', 'Notch', 'off'); hold on;

% stats
if statsFlag
    pAll = [];
    for ii = 1:length(ROI_name)
        statusUnique = unique(statusList3);
        for ss = 1:length(statusUnique)
            hh(ss) = swtest(data(statusList3 == statusUnique(ss),ii));
        end
        if length(statusUnique) == 3
            cmpGroups = {[1 2],[1 3],[2 3]};
        elseif length(statusUnique) == 4
            cmpGroups = {[1 2],[1 3],[1 4],[2 3],[2 4],[3 4]};
        end
        for gg = 1:length(cmpGroups)
            cmpInd = cmpGroups{gg};
            if sum(hh(cmpInd)) == 0
                [~,pAll(gg,ii)] = ttest2(data(statusList3 == cmpInd(1)-1,ii),...
                    data(statusList3 == cmpInd(2)-1,ii));
            else
                pAll(gg,ii) = ranksum(data(statusList3 == cmpInd(1)-1,ii),...
                    data(statusList3 == cmpInd(2)-1,ii));
            end
        end
        c = pAll(:,ii);
        [pTmp, tbl] = kruskalwallis(data(:,ii), statusList3, 'off');
        pANOVA(ii) = pTmp;
        chiANOVA(ii) = tbl{2,5};
        % c = multcompare(stats,'Display','off','CType','dunn-sidak');
        % pAll(:,ii) = c(:,6);
    end
    % multiple comparison correction
    pANOVACorr = reshape(fdr_BH(pANOVA, 0.05, false),size(pANOVA));
    for ii = 1:length(ROI_name)
        if pANOVACorr(ii) > 0.05
            pAll(:,ii) = nan;
        end
    end
    pCorr = reshape(fdr_BH(pAll, 0.05, false),size(pAll));

end

if statsFlag
    for ii = 1:length(ROI_name)
%         pStat = pAll(:,ii);
%         pStatCorr = pCorr(:,ii);
        pStat = pCorr(:,ii);
        pStatCorr = nan*pCorr(:,ii);
        if size(c,1) == 3
            sigstar({ii+[-0.33,0],ii+[-0.33,0.33],ii+[0, 0.33]},pStat,0,max(data(:,ii)),0,pStatCorr);
        elseif size(c,1) == 6
            sigstar({ii+[-0.375, -0.125],ii+[-0.375, 0.125],ii+[-0.375, 0.375], ...
                ii+[-0.125,0.125],ii+[-0.125,0.375],ii+[0.125,0.375]},...
                pStat,0,max(data(:,ii)),0,pStatCorr);
        elseif size(c,1) == 10
            sigstar({ii+[-0.4,-0.2],ii+[-0.4,0],ii+[-0.4,0.2],ii+[-0.4,0.4], ...
                ii+[-0.2,0],ii+[-0.2,0.2],ii+[-0.2,0.4],ii+[0,0.2],ii+[0,0.4],ii+[0.2,0.4]},...
                pStat,0,max(data(:,ii)),0,pStatCorr);
        end
    end

end

% separation lines
for ii = 1:length(ROI_name)+1
    plot([ii-0.5 ii-0.5],ylines,...
        ':','Color', [0.5 0.5 0.5]);
end
axis tight;

% set(gca, 'YScale', 'log');
xticks(1:length(ROI_name));
xticklabels(ROI_name);

end

function [pCorr, pAll] = plotANOVA_noCorr(data, ROI_name, statusList3, ylines, statsFlag)

if nargin < 5
    statsFlag = 1;
end

ymax = max(data(:)); ymin = min(data(:)); yrange = ymax - ymin;
colorList = repmat(statusList3',[1 length(ROI_name)]);
roiList = repmat(1:length(ROI_name), [size(data,1) 1]);

b = boxchart(roiList(:), data(:), 'GroupByColor', colorList(:), ...
    'MarkerStyle', '.', 'Notch', 'off'); hold on;

% stats
if statsFlag
    pAll = [];
    for ii = 1:length(ROI_name)
        [~,~,stats] = kruskalwallis(data(:,ii), statusList3, 'off');
        c = multcompare(stats,'Display','off','CType','dunn-sidak');
        pAll(:,ii) = c(:,6);
    end
    % multiple comparison correction
    pCorr = reshape(fdr_BH(pAll, 0.05, false),size(pAll));

end

if statsFlag
    for ii = 1:length(ROI_name)
        pStat = pAll(:,ii);
        if size(c,1) == 3
            sigstar({ii+[-0.33,0],ii+[-0.33,0.33],ii+[0, 0.33]},pStat,0,max(data(:,ii)),0);
        elseif size(c,1) == 6
            sigstar({ii+[-0.375, -0.125],ii+[-0.375, 0.125],ii+[-0.375, 0.375], ...
                ii+[-0.125,0.125],ii+[-0.125,0.375],ii+[0.125,0.375]},...
                pStat,0,max(data(:,ii)),0);
        elseif size(c,1) == 10
            sigstar({ii+[-0.4,-0.2],ii+[-0.4,0],ii+[-0.4,0.2],ii+[-0.4,0.4], ...
                ii+[-0.2,0],ii+[-0.2,0.2],ii+[-0.2,0.4],ii+[0,0.2],ii+[0,0.4],ii+[0.2,0.4]},...
                pStat,0,max(data(:,ii)),0);
        end
    end

end

% separation lines
for ii = 1:length(ROI_name)+1
    plot([ii-0.5 ii-0.5],ylines,...
        ':','Color', [0.5 0.5 0.5]);
end
axis tight;

% set(gca, 'YScale', 'log');
xticks(1:length(ROI_name));
xticklabels(ROI_name);

end
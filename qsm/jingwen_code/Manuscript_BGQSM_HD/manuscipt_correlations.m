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

%% read table

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20230213.xlsx','Sheet','NoRep');
T(T.CAG < 36,:) = [];
T(strcmp(T.status_reclass,'MISSING'),:) = [];

%% demographics

CAPS = age.*(CAG-33.66)/432.3326;
CAPS(isnan(statusList)) = nan;

CAG = [HD_BGanalysis.CAG];
TMS = [HD_BGanalysis.TMS];
DCL = [HD_BGanalysis.DCL];
TFC = [HD_BGanalysis.TFC];

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

%% Clinical correlations

DART = [T.DART_ToC];
Flanker = [T.Flanker_Score];
Match = [T.Match_CorrectTotal];
SetShift = [T.SetShift_Shift_Score];
DARTtiming = [T.DART_Total_AvgReaction];
FavRec = [T.FavRec_totalCorrect];
FavDelay = [T.FavDelay_totalCorrect];
LOtime = [T.LO_AvgTimePerTrial];
LOscore = [T.LO_ThresholdScore_3_10];

% multiple correlation correction
% Nyholt Dale R.. A Simple Correction for Multiple Testing for Single-Nucleotide Polymorphisms
% in Linkage Disequilibrium with Each Other. The American Journal of Human Genetics.
% 2004;74:765–769. (not used) -> BH

set(0,'DefaultAxesFontSize', 13);

indInclude = statusList < 6;
indROI = [1:7 9];

% CliData = [CAPS' YTO' TMS' DART DARTtiming Flanker Match SetShift ...
%     FavRec FavDelay LOscore LOtime];
% Cli_name = {'CAPS','YTO','TMS','DART', 'DARTtime', 'Flanker','Match','Set Shift',...
%     'FavRec','FavDelay','LOscore','LOtime'};

CliData = [CAPS' YTO' TMS' DART Flanker Match SetShift];
Cli_name = {'CAPS','YTO','TMS','DART', 'Flanker','Match','Set Shift'};

figure('position', [100 100 700 800]);
subaxis(2,2,1, 'SpacingHoriz', 0.07);
ImgData = [volMatAC(:,indROI)];
[rho, pCorr] = correlograph(ImgData(indInclude,:), CliData(indInclude,:), ROI_name(indROI), Cli_name); 
axis tight; xticks([]); colorbar off; 
title('Corrected Volume');

subaxis(2,2,2);
ImgData = [qsmMatAC(:,indROI)];
[rho, pCorr] = correlograph(ImgData(indInclude,:), CliData(indInclude,:), ROI_name(indROI), Cli_name);
axis tight; xticks([]); colorbar off; 
title('Susceptibility');

subaxis(2,2,3);
ImgData = [faMat(:,indROI)];
[rho, pCorr] = correlograph(ImgData(indInclude,:), CliData(indInclude,:), ROI_name(indROI), Cli_name); 
axis tight; colorbar off;
title('FA');

subaxis(2,2,4);
ImgData = [mdMat(:,indROI)];
[rho, pCorr] = correlograph(ImgData(indInclude,:), CliData(indInclude,:), ROI_name(indROI), Cli_name); 
axis tight; colorbar off; 
title('MD');

pause(1); export_fig([img_root '/All_corrClinical'], '-png','-transparent');

%% Clinical correlations

indInclude = statusList5 < 6;

DemData = [grp2idx(sex) age'];

figure('position', [100 100 700 800]);
subaxis(2,2,1, 'SpacingHoriz', 0.07);
ImgData = [volMatAC(:,indROI)];
[rho, pCorr] = correlograph_PC(ImgData(indInclude,:), CliData(indInclude,:), ROI_name(indROI), Cli_name, DemData(indInclude,:)); 
axis tight; xticks([]); colorbar off; 
title('Corrected Volume');

subaxis(2,2,2);
ImgData = [qsmMatAC(:,indROI)];
[rho, pCorr] = correlograph_PC(ImgData(indInclude,:), CliData(indInclude,:), ROI_name(indROI), Cli_name, DemData(indInclude,:)); 
axis tight; xticks([]); colorbar off; 
title('Susceptibility');

subaxis(2,2,3);
ImgData = [faMat(:,indROI)];
[rho, pCorr] = correlograph_PC(ImgData(indInclude,:), CliData(indInclude,:), ROI_name(indROI), Cli_name, DemData(indInclude,:)); 
axis tight; colorbar off;
title('FA');

subaxis(2,2,4);
ImgData = [mdMat(:,indROI)];
[rho, pCorr] = correlograph_PC(ImgData(indInclude,:), CliData(indInclude,:), ROI_name(indROI), Cli_name, DemData(indInclude,:)); 
axis tight; colorbar off; 
title('MD');

pause(1); export_fig([img_root '/All_corrClinical_PC'], '-png','-transparent');

%% plot Correlations

set(0,'DefaultAxesFontSize', 14);

indInclude = statusList < 5;

figure('position', [100 100 700 800]);

subplot(3,2,1);
score = CAPS; xStr = 'CAPS';
ROIselect = 1;
data = volMatAC(:,ROIselect,1);
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr);
xlabel([ROI_name{ROIselect} ' corrected volume (ml)']);
if p < 0.0001
    text(2.8,1.4,sprintf('r = %.2f, p < 0.0001',r));
else
    text(2.8,1.4,sprintf('r = %.2f, p = %.4f',r, p));
end
title(''); axis tight

subplot(3,2,2);
score = TMS; xStr = 'TMS';
ROIselect = 1;
data = volMatAC(:,ROIselect,1);
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr);
xlabel([ROI_name{ROIselect} ' corrected volume (ml)']);
if p < 0.0001
    text(1.8,-4,sprintf('r = %.2f, p < 0.0001',r));
else
    text(1.8,-4,sprintf('r = %.2f, p = %.4f',r, p));
end
title(''); axis tight; ylim([-8 48]);
legend({'','','HC','PM far','PM near','Manifest'},'location','best','box','off','FontSize',10);

subplot(3,2,3);
score = Match'; xStr = 'Match Score';
ROIselect = 1;
data = volMatAC(:,ROIselect,1);
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr);
xlabel([ROI_name{ROIselect} ' corrected volume (ml)']);
if p < 0.0001
    text(2.8,30,sprintf('r = %.2f, p < 0.0001',r));
else
    text(2.8,30,sprintf('r = %.2f, p = %.4f',r, p));
end
title(''); axis tight

subplot(3,2,4);
score = Match'; xStr = 'Match Score';
ROIselect = 2;
data = volMatAC(:,ROIselect,1);
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr);
xlabel([ROI_name{ROIselect} ' corrected volume (ml)']);
if p < 0.0001
    text(3.5,30,sprintf('r = %.2f, p < 0.0001',r));
else
    text(3.5,30,sprintf('r = %.2f, p = %.4f',r, p));
end
title(''); axis tight
legend('off');

subplot(3,2,5);
score = Match'; xStr = 'Match Score';
ROIselect = 2;
data = qsmMatAC(:,ROIselect,1);
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr);
xlabel([ROI_name{ROIselect} ' susceptibility (ppm)']);
if p < 0.0001
    text(0.01,70,sprintf('r = %.2f, p < 0.0001',r));
else
    text(0.01,70,sprintf('r = %.2f, p = %.4f',r, p));
end
title(''); axis tight
legend('off');

subplot(3,2,6);
score = Match'; xStr = 'Match Score';
ROIselect = 2;
data = faMat(:,ROIselect,1);
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr);
xlabel([ROI_name{ROIselect} ' FA']);
if p < 0.0001
    text(0.15,70,sprintf('r = %.2f, p < 0.0001',r));
else
    text(0.15,70,sprintf('r = %.2f, p = %.4f',r, p));
end
title(''); axis tight
legend('off');

pause(0.1); export_fig([img_root '/Example_corrClinical'], '-png','-transparent');

%% plot Correlations

set(0,'DefaultAxesFontSize', 14);

indInclude = statusList < 3;

figure('position', [100 100 400 800])
subplot(3,1,1);
score = CAPS; xStr = 'CAPS';
ROIselect = 2;
data = volMatAC(:,ROIselect,1);
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr);
xlabel([ROI_name{ROIselect} ' corrected volume (ml)']);
text(4.5,0.96,sprintf('r = %.2f, p = %.3f',r,p));
title(''); axis tight
legend({'','','','PM far','PM near'},'location','best','box','off','FontSize',10);

subplot(3,1,2);
score = CAPS; xStr = 'CAPS';
ROIselect = 2;
data = qsmMatAC(:,ROIselect,1);
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr); 
xlabel([ROI_name{ROIselect} ' susceptibility (ppm)']); % xlim([0.03 0.15]);
text(-0.01,1,sprintf('r = %.2f, p = %.3f',r,p));
title(''); axis tight
legend('off');

subplot(3,1,3);
score = TMS; xStr = 'TMS';
ROIselect = 7;
data = volMatAC(:,ROIselect,1);
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr); ylim([-1 6]);
xlabel([ROI_name{ROIselect} ' corrected volume (ml)']); xlim([0.4 0.55]);
text(0.46,5.5,sprintf('r = %.2f, p = %.3f',r,p));
title('');
legend('off');

pause(0.1); export_fig([img_root '/Example_corrClinical2'], '-png','-transparent');

%% correlation between metrics

indInclude = statusList > 1;
metName = {'Vol.','Susc.','FA','MD'};

figure('position', [100 100 800 800]);
for ii = 1:length(ROI_name)
    subaxis(3,3,ii);
    indROI = ii;
    ImgData = [volMatAC(:,indROI) qsmMatAC(:,indROI) faMat(:,indROI) mdMat(:,indROI)];
    [rho, pCorr] = correlographMRI(ImgData, ROI_name{indROI}, indInclude, metName);
    axis equal tight; colorbar off;
    if ii < 7; xticks([]); end
    if mod(ii,3) ~= 1; yticks([]); end
end

pause(0.1); export_fig([img_root '/All_corrImg'], '-png','-transparent');

%% correlation between metrics - multiple correction for all

indInclude = statusList > 1;
metName = {'Vol.','Susc.','FA','MD'};

indROI = [1:7 9];

ImgData = [volMatAC(:,indROI) qsmMatAC(:,indROI) faMat(:,indROI) mdMat(:,indROI)];

figure('position', [100 100 800 800]);
[rho_all, pCorr] = correlographMRIall(ImgData, ROI_name(indROI), indInclude, metName);

pause(0.1); export_fig([img_root '/All_corrImg2'], '-png','-transparent');

%% plot correlation

indInclude = statusList3 > 0;

figure('position', [100 100 800 800]);

met1 = volMatAC;
met2 = qsmMatAC;

subplot(221);
ROIind = [1 2];
[h, legend_str] = plotCorrMetricsAll(met1, met2, ROIind, ROI_name, [1 2]);
xlim([1.7 6]); ylim([-0.01 0.1]);
legend(h, legend_str,'location','best','FontSize',10);
xlabel('Corrected Volume (mL)');
ylabel('Susceptibility (ppm)');

subplot(222);
ROIind = [3 4 9];
[h, legend_str] = plotCorrMetricsAll(met1, met2, ROIind, ROI_name, [3 6 4]);
xlim([0.3 1.4]); ylim([-0.04 0.21]);
legend(h, legend_str,'location','best','FontSize',10);
xlabel('Corrected Volume (mL)');
ylabel('Susceptibility (ppm)');

met1 = faMat;
met2 = qsmMatAC;

subplot(223);
ROIind = [1 2 3 9];
[h, legend_str] = plotCorrMetricsAll(met1, met2, ROIind, ROI_name, [1 2 3 4]);
xlim([0.1 0.5]); ylim([-0.02 0.25]);
legend(h, legend_str,'location','best','FontSize',10);
xlabel('FA');
ylabel('Susceptibility (ppm)');

met1 = faMat;
met2 = mdMat*1e3;

subplot(224);
ROIind = [2 3 5];
[h, legend_str] = plotCorrMetricsAll(met1, met2, ROIind, ROI_name, [2 3 5]);
xlim([0.1 0.5]); ylim([0.07 0.65]);
xlabel('FA');
ylabel('MD (x10^{-3} mm^2/s)');
legend(h, legend_str,'location','best','FontSize',10);

pause(0.1); export_fig([img_root '/Example_corrImg'], '-png','-transparent');

%% plot correlation - all

% met1 = mdMat*1e3; 
% met2 = faMat; 
% indHD = ~indHC;
% 
% xStr = 'MD (x10^{-3} mm^2/s)'; % MD (x10^{-3} mm^2/s) % FA 
% yStr = 'FA'; % Corrected Volume (mL) % Susceptibility (ppm)
% 
% figure('position', [100 100 1200 800]);
% plotCorrMetricsSubplots(met1, met2, ROI_name, indHD, xStr, yStr);
% 
% pause(0.1); export_fig([img_root '/PlotAll_corrImg1'], '-png','-transparent');

%% TSNE

indvalid = ~isnan(faMat(:,9));

ROIind = [1 2 3 4 9];
VarMRI = [qsmMatAC(:,ROIind) volMatAC(:,ROIind) faMat(:,ROIind) mdMat(:,ROIind)];

Y = tsne(VarMRI);
gscatter(Y(:,1),Y(:,2),statusList5(indvalid));

%% helper function

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
matAll = matR; % (matL + matR)/2;

end

function [R,P, F,gof] = plotCorr(data, score, ROI_name, indHC, indPMfar, indPMnear, indM, xStr, plotType)

if nargin < 9
    plotType = 'poly1';
end

xrange = max(score) - min(score);
xplot = linspace(min(score)-0.1*xrange, max(score)+0.1*xrange, 50); hold on;
[F,gof] = fit(score(~isnan(data) & ~isnan(score')),data(~isnan(data) & ~isnan(score'))',plotType);
ci = predint(F, xplot, 0.95, 'functional','off');
h = plot(xplot, F(xplot)); h.Color = [0.5 0.5 0.5];
p = fill([xplot xplot(end:-1:1)],[ci(:,1); ci(end:-1:1,2)],'k');
p.FaceColor = 'k'; 
p.FaceAlpha = 0.1;
p.EdgeColor = 'none';  

set(gca,'ColorOrderIndex',1);

scatter(score(indHC), data(indHC), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
scatter(score(indPMfar), data(indPMfar), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
scatter(score(indPMnear), data(indPMnear), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
scatter(score(indM), data(indM), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
% scatter(score(indEM2), data(indEM2,pp), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
xlabel(ROI_name);
ylabel(xStr);

[R,P] = corr(score(~isnan(data)),data(~isnan(data))','Type','Pearson','Rows','complete');
title(sprintf('r %.2f p %.4f', R, P));

% legend({'HC','PM far','PM near','Manifest'},'location','best','FontSize',10);

end

function [rho_all, pCorr] = correlographMRIall(ImgData, ROI_name, indInclude, metName)

mask = [0 1 1 1; 0 0 1 1; 0 0 0 1; 0 0 0 0] > 0;
mask2 = flip(flip(mask,1),2);

J_rho = customcolormap_preset('red-white-blue');

rho_all = zeros(length(metName), length(metName), length(ROI_name));
p_all = zeros(length(metName), length(metName), length(ROI_name));
for ii = 1:length(ROI_name)
    ImgData_ROI = ImgData(:,ii:length(ROI_name):end);
    [rho, pval] = corr(ImgData_ROI, 'Type', 'Pearson', 'Rows', 'Pairwise');
    [rho2, pval2] = corr(ImgData_ROI(indInclude,:), 'Type', 'Pearson', 'Rows', 'Pairwise');
    
    rho_all(:,:,ii) = rho.*mask + rho2.*~mask;
    p_all(:,:,ii) = pval.*mask + pval2.*~mask;
end

ptmp = reshape(p_all,[], length(ROI_name));
ptmp = reshape(fdr_BH(ptmp(mask | mask2, :), 0.05, false), size(ptmp(mask | mask2, :)));
pCorr = reshape(p_all,[], length(ROI_name));
pCorr(mask | mask2, :) = ptmp;
pCorr = reshape(pCorr, size(p_all));

for nn = 1:length(ROI_name)
    subaxis(3,3,nn);
    imagesc(rho_all(:,:,nn), [-1 1]); colorbar; colormap(J_rho); hold on;
    for ii = 1:size(pCorr,1)
        for jj = 1:size(pCorr,2)
            if p_all(ii,jj,nn) < 0.0001
                t = text(jj,ii,'****','HorizontalAlignment','center');
            elseif p_all(ii,jj,nn) < 0.001
                t = text(jj,ii,'***','HorizontalAlignment','center');
            elseif p_all(ii,jj,nn) < 0.01
                t = text(jj,ii,'**','HorizontalAlignment','center');
            elseif p_all(ii,jj,nn) < 0.05
                t = text(jj,ii,'*','HorizontalAlignment','center');
            end
            
            if pCorr(ii,jj,nn) < 0.05
                t.Color = 'w';
            end
        end
    end
    axis equal tight; colorbar off;
    title(ROI_name{nn});
    yticks(1:4);
    yticklabels(metName);
    xticks(1:4);
    xticklabels(metName);
    if nn < 7; xticks([]); end
    if mod(nn,3) ~= 1; yticks([]); end
end

end

function [rho_combined, pCorr] = correlographMRI(ImgData, ROI_name, indInclude, metName)

[rho, pval] = corr(ImgData, 'Type', 'Pearson', 'Rows', 'Pairwise');
[rho2, pval2] = corr(ImgData(indInclude,:), 'Type', 'Pearson', 'Rows', 'Pairwise');

J_rho = customcolormap_preset('red-white-blue');

mask = [0 1 1 1; 0 0 1 1; 0 0 0 1; 0 0 0 0] > 0;
mask2 = flip(flip(mask,1),2);
mask_diag = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] > 0;
rho_combined = rho.*mask + rho2.*~mask;
p_combined = pval.*mask + pval2.*mask2;
p_combined(mask_diag) = nan;
pCorr_tmp = fdr_BH(p_combined(mask | mask2), 0.05, false);
pCorr = pval;
pCorr(mask | mask2) = pCorr_tmp;

imagesc(rho_combined, [-1 1]); colorbar; colormap(J_rho); hold on;

for ii = 1:size(p_combined,1)
    for jj = 1:size(p_combined,2)
        if p_combined(ii,jj) < 0.0001
            t = text(jj,ii,'****','HorizontalAlignment','center');
        elseif p_combined(ii,jj) < 0.001
            t = text(jj,ii,'***','HorizontalAlignment','center');
        elseif p_combined(ii,jj) < 0.01
            t = text(jj,ii,'**','HorizontalAlignment','center');
        elseif p_combined(ii,jj) < 0.05
            t = text(jj,ii,'*','HorizontalAlignment','center');
        end
        
        if pCorr(ii,jj) < 0.05
            t.Color = 'w';
        end
    end
end

title(ROI_name);

yticks(1:4); 
yticklabels(metName);
xticks(1:4); 
xticklabels(metName);

end

function [rho, pCorr] = correlograph_PC(ImgData, CliData, ROI_name, Cli_name, DemData)

[rho, pval] = partialcorr(ImgData, CliData, DemData, 'Type', 'Pearson', 'Rows', 'pairwise');
% pCorr = pval*MeffMRI*MeffClin; 
pCorr = reshape(fdr_BH(pval, 0.05, false),size(pval));

J_rho = customcolormap_preset('red-white-blue');

imagesc(rho, [-1 1]); colorbar; colormap(J_rho); hold on;
yticks([1:length(ROI_name)]); 
yticklabels(ROI_name);
xticks([1:size(CliData,2)]); 
xticklabels(Cli_name);

for ii = 1:size(pCorr,1)
    for jj = 1:size(pCorr,2)
        if pval(ii,jj) < 0.0001
            t = text(jj,ii,'****','HorizontalAlignment','center');
        elseif pval(ii,jj) < 0.001
            t = text(jj,ii,'***','HorizontalAlignment','center');
        elseif pval(ii,jj) < 0.01
            t = text(jj,ii,'**','HorizontalAlignment','center');
        elseif pval(ii,jj) < 0.05
            t = text(jj,ii,'*','HorizontalAlignment','center');
        end
        
        if pCorr(ii,jj) < 0.05
            t.Color = 'w';
        end
    end
end

end

function [rho, pCorr] = correlograph(ImgData, CliData, ROI_name, Cli_name)

[rMRI,~] = corr(ImgData,ImgData,'type','Pearson','Rows','pairwise');
[~,S,~] = svd(rMRI);
MeffMRI = 1 + (size(ImgData,2)-1)*(1-var(diag(S))/size(ImgData,2));

[rClin,~] = corr(CliData,CliData,'type','Pearson','Rows','pairwise');
[~,S,~] = svd(rClin);
MeffClin = 1 + (size(CliData,2)-1)*(1-var(diag(S))/size(CliData,2));

alphaEffective = 0.05/MeffMRI/MeffClin;
disp(['corrected alpha = ' num2str(alphaEffective)]);

[rho, pval] = corr(ImgData, CliData, 'Type', 'Pearson', 'Rows', 'Pairwise');
% pCorr = pval*MeffMRI*MeffClin; 
pCorr = reshape(fdr_BH(pval, 0.05, false),size(pval));

J_rho = customcolormap_preset('red-white-blue');

imagesc(rho, [-1 1]); colorbar; colormap(J_rho); hold on;
yticks([1:length(ROI_name)]); 
yticklabels(ROI_name);
xticks([1:size(CliData,2)]); 
xticklabels(Cli_name);

for ii = 1:size(pCorr,1)
    for jj = 1:size(pCorr,2)
        if pval(ii,jj) < 0.0001
            t = text(jj,ii,'****','HorizontalAlignment','center');
        elseif pval(ii,jj) < 0.001
            t = text(jj,ii,'***','HorizontalAlignment','center');
        elseif pval(ii,jj) < 0.01
            t = text(jj,ii,'**','HorizontalAlignment','center');
        elseif pval(ii,jj) < 0.05
            t = text(jj,ii,'*','HorizontalAlignment','center');
        end
        
        if pCorr(ii,jj) < 0.05
            t.Color = 'w';
        end
    end
end

end

function [h, legend_str] = plotCorrMetricsAll(met1, met2, ROIind, ROI_name, ColorInd)

clear h;
for ii = 1:length(ROIind)
    hold on;
    set(gca,'ColorOrderIndex',ColorInd(ii));
    h(ii) = scatter(met1(:,ROIind(ii)),met2(:,ROIind(ii)), 50, 'o', 'filled', 'MarkerFaceAlpha', 0.75);
    
    indValid = ~isnan(met1(:,ROIind(ii))') & ~isnan(met2(:,ROIind(ii))');
    F = fit(met1(indValid,ROIind(ii)),met2(indValid,ROIind(ii)),'poly1');
    
    set(gca,'ColorOrderIndex',ColorInd(ii));
    xrange = [min(met1(indValid,ROIind(ii))) max(met1(indValid,ROIind(ii)))];
    xrange = xrange + [-0.2 0.2]*(xrange(2) - xrange(1));
    plot(xrange, F(xrange), '-');
    
    [rho, pval] = corr(met1(indValid,ROIind(ii)),met2(indValid,ROIind(ii)), ...
        'Type', 'Pearson', 'Rows', 'Pairwise');
    
    if pval < 0.0001
        legend_str{ii} = sprintf('%s r = %.3f, p < 0.0001', ROI_name{ROIind(ii)}, rho);
    elseif pval < 0.001
        legend_str{ii} = sprintf('%s r = %.3f, p < 0.001', ROI_name{ROIind(ii)}, rho);
    else
        legend_str{ii} = sprintf('%s r = %.3f, p = %.3f', ROI_name{ROIind(ii)}, rho, pval);
    end
end

end

function [] = plotCorrMetricsSubplots(met1, met2, ROI_name, indHD, xStr, yStr)

for ii = 1:length(ROI_name)
    subplot(3,3,ii)
    scatter(met1(~indHD,ii),met2(~indHD,ii), 50, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
    scatter(met1(indHD,ii),met2(indHD,ii), 50, 'o', 'filled', 'MarkerFaceAlpha', 0.75);
    
    indValid = ~isnan(met1(:,ii)') & ~isnan(met2(:,ii)');
    F = fit(met1(indValid,ii),met2(indValid,ii),'poly1');
    plot([min(met1(indValid,ii)) max(met1(indValid,ii))], ...
        F([min(met1(indValid,ii)) max(met1(indValid,ii))]), 'k-');
    
    [r,p] = corr(met1(indValid,ii),met2(indValid,ii));
    
    title(ROI_name{ii});
    
    if ischar(yStr)
        if mod(ii,3) == 1
            ylabel(yStr);
        end
    else
        ylabel(yStr{ii});
    end
    if ii >= 7
        xlabel(xStr);
    end
    
    [xloc, yloc] = calcTextLocation(xlim, ylim);
    text(xloc, yloc, sprintf('r = %.2f p = %.3f', r, p));
end

legend({'HC','HD'},'location','best','FontSize',10);

end

function [xloc, yloc] = calcTextLocation(xlim, ylim)

xrange = xlim(2) - xlim(1);
yrange = ylim(2) - ylim(1);

xloc = xlim(1)+0.05*xrange;
yloc = ylim(1)+0.05*yrange;

end
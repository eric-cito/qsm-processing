clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/'));

dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/SubjImgDataBG_SSTGV.mat';
img_root = '/home/jyao3/030_QSM/img_temp';

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

% CAPS = [HD_BGanalysis.CAPS];
CAPS = age.*(CAG-35.5);
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

indPM = statusList == 2;

fprintf('=== PM === \n');
printDemo(indPM, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

indEM = statusList == 3;

fprintf('=== EM === \n');
printDemo(indEM, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

indPMfar = indPM & YTO > 15;
indPMnear = indPM & YTO <= 15;

indEM1 = (indEM | indMan) & DCL < 4;
indEM2 = (indEM | indMan) & DCL >= 4;

fprintf('=== PM far === \n');
printDemo(indPMfar, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

fprintf('=== PM near === \n');
printDemo(indPMnear, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

fprintf('=== EM DCL < 4 === \n');
printDemo(indEM1, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

fprintf('=== EM DCL >= 4 === \n');
printDemo(indEM2, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

%% Imaging metrics

Nroi = size(HD_BGanalysis(1).imData,1)/2;

DTI_flag = 1;

volMatL = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2);
qsmMatL = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2);
qsmSDMatL = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2);
qsmMatIIL = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2);
dtiMatL = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2,4); % FA MD RD ODI
volMatR = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2);
qsmMatR = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2);
qsmSDMatR = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2);
qsmMatIIR = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2);
dtiMatR = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1)/2,4);
for ii = 1:length(HD_BGanalysis)
    volMatL(ii,:) = HD_BGanalysis(ii).imData.NormVolume(1:2:end);
    qsmMatL(ii,:) = HD_BGanalysis(ii).imData.QSMmedian(1:2:end);
    qsmSDMatL(ii,:) = HD_BGanalysis(ii).imData.QSMmad(1:2:end);
    qsmMatIIL(ii,:) = HD_BGanalysis(ii).imData.QSMIImedian(1:2:end);
    if DTI_flag
        dtiMatL(ii,:,1) = HD_BGanalysis(ii).imData.FAmedian(1:2:end);
        dtiMatL(ii,:,2) = HD_BGanalysis(ii).imData.MDmedian(1:2:end);
        dtiMatL(ii,:,3) = HD_BGanalysis(ii).imData.RDmedian(1:2:end);
        dtiMatL(ii,:,4) = HD_BGanalysis(ii).imData.ODImedian(1:2:end);
    end
    
    volMatR(ii,:) = HD_BGanalysis(ii).imData.NormVolume(2:2:end);
    qsmMatR(ii,:) = HD_BGanalysis(ii).imData.QSMmedian(2:2:end);
    qsmSDMatR(ii,:) = HD_BGanalysis(ii).imData.QSMmad(2:2:end);
    qsmMatIIR(ii,:) = HD_BGanalysis(ii).imData.QSMIImedian(2:2:end);
    if DTI_flag
        dtiMatR(ii,:,1) = HD_BGanalysis(ii).imData.FAmedian(2:2:end);
        dtiMatR(ii,:,2) = HD_BGanalysis(ii).imData.MDmedian(2:2:end);
        dtiMatR(ii,:,3) = HD_BGanalysis(ii).imData.RDmedian(2:2:end);
        dtiMatR(ii,:,4) = HD_BGanalysis(ii).imData.ODImedian(2:2:end);
    end
    
end
volMat = (volMatL + volMatR);
qsmMat = (qsmMatL + qsmMatR)/2;
qsmSDMat = (qsmSDMatL + qsmSDMatR)/2;
qsmIIMat = (qsmMatIIL + qsmMatIIR)/2;
dtiMat = (dtiMatL + dtiMatR)/2;

ROI_name = {'CN','PU','GP','SNpr','SNpc','RN','DN','SN','FS-CN','FS-PU','FS-GP'};

%% Age correction

volMatAC = volMat;
qsmMatAC = qsmMat;
qsmIIMatAC = qsmIIMat;
dtiMatAC = dtiMat;
for ii = 1:Nroi
    volMatAC(:,ii) = ageCorr(age, volMat(:,ii), indHC);
    qsmMatAC(:,ii) = ageCorr(age, qsmMat(:,ii), indHC);
    qsmIIMatAC(:,ii) = ageCorr(age, qsmIIMat(:,ii), indHC);
    if DTI_flag
        dtiMatAC(:,ii,1) = ageCorr(age, dtiMat(:,ii,1), indHC);
        dtiMatAC(:,ii,2) = ageCorr(age, dtiMat(:,ii,2), indHC);
        dtiMatAC(:,ii,3) = ageCorr(age, dtiMat(:,ii,3), indHC);
        dtiMatAC(:,ii,4) = ageCorr(age, dtiMat(:,ii,4), indHC);
    end
end

%% group comparison

statusListEx = statusList;
statusListEx(indHC) = 0;
statusListEx(indPMfar) = 1;
statusListEx(indPMnear) = 2;
statusListEx(indEM1) = 3;
statusListEx(indEM2) = 4;

figure('position', [100 100 1500 900]);
data = qsmSDMat(:,:,1);
for pp = 1:size(data,2)
    subplot(3,4,pp);
    [p,pPMfar,pPMnear,pEM] = plotROI_ANOVAex(data(:,pp), statusListEx, ROI_name{pp});
    legend('off');
end
pause(1); export_fig([img_root '/vol_group1'], '-png','-transparent'); % close;

figure('position', [100 100 1500 900]);
for pp = 1:size(data,2)
    subplot(3,4,pp);
    [p,pPM,pEM] = plotROI_ANOVA(data(:,pp), statusList, ROI_name{pp});
    legend('off');
end
pause(1); export_fig([img_root '/vol_group2'], '-png','-transparent'); % close;

%% CAPS correlation

score = CAPS;

figure('position', [100 100 1500 900]);
data = qsmSDMat(:,:,1);
for pp = 1:size(data,2)
    subplot(3,4,pp);
    
    scatter(score(indPMfar), data(indPMfar,pp), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
    scatter(score(indPMnear), data(indPMnear,pp), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
    scatter(score(indEM1), data(indEM1,pp), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
    scatter(score(indEM2), data(indEM2,pp), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
    xlabel('CAPS');
    
    [R,P] = corr(score(~isnan(score))',data(~isnan(score),pp),'Type','Spearman');
    title(sprintf('r %.2f p %.4f', R, P));
    
end
legend({'PM1','PM2','EM1','EM2'},'location','best');

pause(0.1); export_fig([img_root '/vol_caps'], '-png','-transparent');

%% correlation between metrics

indHD = statusList > 1;
met1 = dtiMat(:,:,1);
met2 = qsmMat; % dtiMat(:,:,1)

ROIind = [9 10 11];
legend_str = {};

clear h

figure('position', [100 100 800 500]);
for ii = 1:length(ROIind)
    
    indValid = indHD & ~isnan(met1(:,ROIind(ii))') & ~isnan(met2(:,ROIind(ii))');
    set(gca,'ColorOrderIndex',ii);
    h(ii) = scatter(met1(indHD,ROIind(ii)),met2(indHD,ROIind(ii)), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
    set(gca,'ColorOrderIndex',ii);
    F = fit(met1(indValid,ROIind(ii)),met2(indValid,ROIind(ii)),'poly1');
    plot([min(met1(indValid,ROIind(ii))) max(met1(indValid,ROIind(ii)))], ...
        F([min(met1(indValid,ROIind(ii))) max(met1(indValid,ROIind(ii)))]), '-');
    
    [R(ii),P(ii)] = corr(met1(indValid,ROIind(ii)),met2(indValid,ROIind(ii)),'Type','Spearman');
    legend_str{ii} = sprintf('%s (r %.2f p %.3f)', ROI_name{ROIind(ii)}, R(ii), P(ii));
    
end

% Corrected Volume (mL)
% Susceptibility (ppm)
% RII Susceptibility (ppm)
% FA

legend(h,legend_str,'box','off','location','best'); % 
xlabel('FA'); 
ylabel('Susceptibility (ppm)');
pause(0.1); export_fig([img_root '/vol_qsm5'], '-png','-transparent');

%% clustering

ROIind = [1 2 3 6 7 8];
VarMRI = [qsmMatAC(:,ROIind) volMatAC(:,ROIind) dtiMatAC(:,ROIind,1) dtiMatAC(:,ROIind,2)];

n_kcluster = 3;

rng(1);
[idx,C] = kmeans(VarMRI,n_kcluster);

figure;
for ii = 1:n_kcluster
    plot3(VarMRI(idx==ii,1),VarMRI(idx==ii,7),VarMRI(idx==ii,13),'.','MarkerSize',12); hold on
end
plot3(VarMRI(indHC,1),VarMRI(indHC,7),VarMRI(indHC,13),'ko','MarkerSize',12)
plot3(VarMRI(indPM,1),VarMRI(indPM,7),VarMRI(indPM,13),'co','MarkerSize',12)
plot3(VarMRI(indEM,1),VarMRI(indEM,7),VarMRI(indEM,13),'mo','MarkerSize',12)
plot3(VarMRI(indMan,1),VarMRI(indMan,7),VarMRI(indMan,13),'ro','MarkerSize',12)
xlabel('CN susc.');
ylabel('CN vol.');
zlabel('CN FA');
legend({'','','','HC','PM','EM','M'});

figure;
subplot(241); boxchart(idx,TMS);
subplot(242); boxchart(statusList,TMS);
subplot(243); boxchart(idx,DCL);
subplot(244); boxchart(statusList,DCL);
subplot(245); boxchart(idx,statusList);
subplot(246); boxchart(statusList,statusList);
subplot(247); boxchart(idx,YTO);
subplot(248); boxchart(statusList,YTO);

%% TSNE

ROIind = [1 2 3 6 7 8];
VarMRI = [qsmMatAC(:,ROIind) volMatAC(:,ROIind)];

Y = tsne(VarMRI);
gscatter(Y(:,1),Y(:,2),statusList);

%% 3dplot

ROIind = 2;

QSM = qsmIIMat(:,ROIind);
FA = dtiMat(:,ROIind,1);
Vol = volMat(:,ROIind);

figure;
scatter3(FA(indHD),Vol(indHD),QSM(indHD),'filled'); hold on;
scatter3(FA(~indHD),Vol(~indHD),QSM(~indHD),'filled');

legend('HD','HC');
xlabel('FA'); ylabel('Volume'); zlabel('Susceptibility');
pause(0.1); export_fig([img_root '/corr3d2'], '-png','-transparent');

mdl = fitglm([FA,Vol,age'],QSM)

%% TMS GLM

ROIind = 1;

QSM = qsmIIMat(:,ROIind);
FA = dtiMat(:,ROIind,1);
Vol = volMat(:,ROIind);

sexList = strcmp(sex,'M');

mdl = fitglm([age', sexList', QSM], TMS)

%% plot all

plotRoot = [img_root '/volAC'];
Tvol = plotAll(age, volMatAC(:,:,1), ROI_name, statusList, CAG, CAPS, AOO, TMS, DCL, TFC, plotRoot);

% plotRoot = [img_root '/qsmAC'];
% Tqsm = plotAll(age, qsmMatAC(:,:,1), ROI_name, statusList, CAG, CAPS, AOO, TMS, DCL, TFC, plotRoot);
% 
% plotRoot = [img_root '/qsmIIAC'];
% TqsmII = plotAll(age, qsmIIMatAC(:,:,1), ROI_name, statusList, CAG, CAPS, AOO, TMS, DCL, TFC, plotRoot);

if DTI_flag
%     plotRoot = [img_root '/faAC'];
%     Tfa = plotAll(age, dtiMatAC(:,:,1), ROI_name, statusList, CAG, CAPS, AOO, TMS, DCL, TFC, plotRoot);
%     
%     plotRoot = [img_root '/mdAC'];
%     Tmd = plotAll(age, dtiMatAC(:,:,2), ROI_name, statusList, CAG, CAPS, AOO, TMS, DCL, TFC, plotRoot);
    
%     plotRoot = [img_root '/rdAC'];
%     Trd = plotAll(age, dtiMatAC(:,:,3), ROI_name, statusList, CAG, CAPS, AOO, TMS, DCL, TFC, plotRoot);
%     
%     plotRoot = [img_root '/odiAC'];
%     Todi = plotAll(age, dtiMatAC(:,:,4), ROI_name, statusList, CAG, CAPS, AOO, TMS, DCL, TFC, plotRoot);
end

%% test non-monotonic relationship

CAPS = [HD_BGanalysis.CAPS];
CAPS(isnan(statusList)) = nan;

indValid = ~isnan(CAPS);
% CAPS(isnan(CAPS)) = 0;

data = qsmMatAC(:,1,1);

figure;
plot(CAPS, data,'o'); hold on;
[r,p] = corr(CAPS', data, 'Rows', 'complete', 'Type', 'Spearman');
[F, gof] = fit(CAPS(indValid)', data(indValid),'smoothingspline')
plot(F,CAPS',data);

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

function [dataAC] = ageCorr(age, data, indHC)

indValid = ~isnan(data) & indHC';

F = fit(age(indValid)',data(indValid),'poly1');
dataAC = data - F.p1*(age'-median(age));

end

function [T] = plotAll(age, metricList, ROI_name, statusList, CAG, CAPS, AOO, TMS, DCL, TFC, plotRoot)

indHC = find(statusList == 1);
indPM = find(statusList == 2);
indEM = find(statusList == 3);
indMan = find(statusList == 4);

figure('position', [100 100 1500 900]);
for pp = 1:size(metricList,2)
    subplot(3,4,pp);
    plotAge(age, metricList(:,pp), indHC, indPM, indEM, indMan);
    ylabel(ROI_name{pp});
    if pp > 1; legend('off'); end
end
pause(1); export_fig([plotRoot '_age'], '-png','-transparent'); close;

figure('position', [100 100 1500 900]);
for pp = 1:size(metricList,2)
    subplot(3,4,pp);
    [r,p] = plotCAG(CAG, metricList(:,pp), indPM, indEM, indMan);
    ylabel(ROI_name{pp}); xlabel('CAG repeats');
    legend('off');
    Rcag(pp) = r; Pcag(pp) = p;
end
pause(1); export_fig([plotRoot '_CAG'], '-png','-transparent'); close;

figure('position', [100 100 1500 900]);
for pp = 1:size(metricList,2)
    subplot(3,4,pp);
    [r,p] = plotCAG(CAPS, metricList(:,pp), indPM, indEM, indMan);
    ylabel(ROI_name{pp}); xlabel('CAPS');
    legend('off');
    Rcaps(pp) = r; Pcaps(pp) = p;
end
pause(1); export_fig([plotRoot '_CAPS'], '-png','-transparent'); close;

YTO = AOO - age;
figure('position', [100 100 1500 900]);
for pp = 1:size(metricList,2)
    subplot(3,4,pp);
    [r,p] = plotCAG(YTO, metricList(:,pp), indPM, indEM, indMan);
    ylabel(ROI_name{pp}); xlabel('YTO');
    legend('off');
    Ryto(pp) = r; Pyto(pp) = p;
end
pause(1); export_fig([plotRoot '_YTO'], '-png','-transparent'); close;

figure('position', [100 100 1500 900]);
for pp = 1:size(metricList,2)
    subplot(3,4,pp);
    [r,p] = plotCAG(TMS, metricList(:,pp), indPM, indEM, indMan);
    ylabel(ROI_name{pp}); xlabel('TMS');
    legend('off');
    Rtms(pp) = r; Ptms(pp) = p;
end
pause(1); export_fig([plotRoot '_TMS'], '-png','-transparent'); close;

figure('position', [100 100 1500 900]);
for pp = 1:size(metricList,2)
    subplot(3,4,pp);
    [r,p] = plotCAG(DCL, metricList(:,pp), indPM, indEM, indMan);
    ylabel(ROI_name{pp}); xlabel('DCL');
    legend('off');
    Rdcl(pp) = r; Pdcl(pp) = p;
end
pause(1); export_fig([plotRoot '_DCL'], '-png','-transparent'); close;

figure('position', [100 100 1500 900]);
for pp = 1:size(metricList,2)
    subplot(3,4,pp);
    [r,p] = plotCAG(TFC, metricList(:,pp), indPM, indEM, indMan);
    ylabel(ROI_name{pp}); xlabel('TFC');
    legend('off');
    Rtfc(pp) = r; Ptfc(pp) = p;
end
pause(1); export_fig([plotRoot '_TFC'], '-png','-transparent'); close;

figure('position', [100 100 1500 900]);
for pp = 1:size(metricList,2)
    subplot(3,4,pp);
    [p,pPM,pEM] = plotROI_ANOVA(metricList(:,pp), statusList, ROI_name{pp});
    legend('off');
    Panova(pp) = p;
    PanovaPM(pp) = pPM;
    PanovaEM(pp) = pEM;
end
pause(1); export_fig([plotRoot '_ANOVA'], '-png','-transparent'); close;

T = table(ROI_name', Rcag', Rcaps', Ryto', Rtms', Rtfc', Rdcl', ...
    Pcag', Pcaps', Pyto', Ptms', Ptfc', Pdcl', Panova', PanovaPM', PanovaEM');

end

function [] = plotAge(ageList, volList, indHC, indPM, indEM, indMan)

plot(ageList(indHC), volList(indHC), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2); hold on;
plot(ageList(indPM), volList(indPM), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.8500 0.3250 0.0980]/2);
plot(ageList(indEM), volList(indEM), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.9290 0.6940 0.1250]/2);
plot(ageList(indMan), volList(indMan), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.4940 0.1840 0.5560]/2);

meanHC = mean(volList(indHC));
stdHC = std(volList(indHC));

plot([20 100], [1 1]*meanHC, '-', 'Color', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2);
plot([20 100], [1 1]*(meanHC+1.96*stdHC), ':', 'Color', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2);
plot([20 100], [1 1]*(meanHC-1.96*stdHC), ':', 'Color', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2);

legend({'HC','PM','EM','M'},'location','best');
xlabel('Age (year)'); ylabel('Susc.');
xlim([20 90]);

ylimL = min(volList) - 0.1*(max(volList) - min(volList));
ylimU = max(volList) + 0.1*(max(volList) - min(volList));
ylim([ylimL ylimU]);

end

function [R,P] = plotCAG(cagList, volList, indPM, indEM, indMan)

indHD = [indPM indEM indMan];

plot(cagList(indPM), volList(indPM), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.8500 0.3250 0.0980]/2, ...
    'MarkerEdgeColor', [0.8500 0.3250 0.0980]); hold on;
plot(cagList(indEM), volList(indEM), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.9290 0.6940 0.1250]/2, ...
    'MarkerEdgeColor', [0.9290 0.6940 0.1250]);
plot(cagList(indMan), volList(indMan), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.4940 0.1840 0.5560]/2, ...
    'MarkerEdgeColor', [0.4940 0.1840 0.5560]);

legend({'PM','EM', 'Manifest'},'location','best');
% xlabel('age x (CAG-30)/6.27'); ylabel('Susc.');

ylimL = min(cagList(indHD)) - 0.1*(max(cagList(indHD)) - min(cagList(indHD)));
ylimU = max(cagList(indHD)) + 0.1*(max(cagList(indHD)) - min(cagList(indHD)));
xlim([ylimL ylimU]);

ylimL = min(volList) - 0.1*(max(volList) - min(volList));
ylimU = max(volList) + 0.1*(max(volList) - min(volList));
ylim([ylimL ylimU]);

A = cagList(indHD)'; B = volList(indHD);
invalidInd = isnan(A) | isnan(B);

[R,P] = corr(A(~invalidInd),B(~invalidInd),'Type','Spearman');
title(sprintf('r %.2f p %.4f', R, P));

% Mdl = fitrgam(A(~invalidInd),B(~invalidInd))

end

function [p,pPM,pEM] = plotROI_ANOVA(data, statusList, ROIname)

statusList(statusList == 4) = 3;
violinplot(data, statusList, 'showMean', true);

ymax = max(data);
ymin = min(data);
yrange = ymax - ymin;

[p,~,stats] = kruskalwallis(data, statusList, 'off');
if p < 0.0001
    text(0.6,ymin+0.1*yrange,sprintf('p < 0.0001'),'FontSize',12); % -0.03
elseif p < 0.001
    text(0.6,ymin+0.1*yrange,sprintf('p < 0.001'),'FontSize',12);
else
    text(0.6,ymin+0.1*yrange,sprintf('p = %.3f',p),'FontSize',12);
end

c = multcompare(stats,'Display','off','CType','dunn-sidak');

ylim([ymin-0.1*yrange ymax+0.2*yrange]);
ylabel([ROIname]);

xlim([0.5 3.5]);
ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel = {'HC','PM','EM+M'};
sigstar({[1,2],[1,3],[2,3]},c(:,6));

pPM = c(1,6); pEM = c(2,6);

end

function [p,pPMfar,pPMnear,pEM] = plotROI_ANOVAex(data, statusListEx, ROIname)

% 0 HC 1 PMfar 2 PMnear 3 EM1 4 M2

% statusListEx(statusListEx == 4) = 3;
violinplot(data, statusListEx, 'showMean', true);

ymax = max(data);
ymin = min(data);
yrange = ymax - ymin;

[p,~,stats] = kruskalwallis(data, statusListEx, 'off');
if p < 0.0001
    text(0.6,ymin+0.1*yrange,sprintf('p < 0.0001'),'FontSize',12); % -0.03
elseif p < 0.001
    text(0.6,ymin+0.1*yrange,sprintf('p < 0.001'),'FontSize',12);
else
    text(0.6,ymin+0.1*yrange,sprintf('p = %.3f',p),'FontSize',12);
end

c = multcompare(stats,'Display','off','CType','dunn-sidak');

ylim([ymin-0.1*yrange ymax+0.2*yrange]);
ylabel([ROIname]);

xlim([0.5 5.5]);
ax = gca;
ax.XTick = [1 2 3 4 5];
ax.XTickLabel = {'HC','PM1','PM2','EM1','EM2'};
sigstar({[1,2],[1,3],[1,4],[1,5]},c(1:4,6));

pPMfar = c(1,6); pPMnear = c(2,6); pEM = c(3:4,6);

end
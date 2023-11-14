clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/'));

% dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/SubjImgDataBG_MNI_QSM.mat';
dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/SubjImgDataBG_ANTS_iLSQR_1004.mat';
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

HD_BGanalysis(38).imData.FA(:) = nan;
HD_BGanalysis(38).imData.MD(:) = nan;
HD_BGanalysis(38).imData.RD(:) = nan;

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
indROI = [1:4];

data = qsmMatAC(:,indROI);
figure('position', [100 0 800 800]);
subaxis(2,2,1,'SpacingHoriz', 0.11);
pCorrQSM = plotANOVA(data, ROI_name(indROI), statusListPlot, [-0.07 0.25]);
ylabel('Susceptibility (ppm)'); ylim([-0.07 0.25]);

data = volMatAC(:,indROI);
subaxis(2,2,2);
pCorrVol = plotANOVA(data, ROI_name(indROI), statusListPlot, [0 7.2]);
ylabel('Corrected Volume (ml)'); ylim([0 7.2]);

data = faMat(:,indROI);
subaxis(2,2,3);
pCorrFA = plotANOVA(data, ROI_name(indROI), statusListPlot, [0 0.75]);
ylabel('FA'); ylim([0 0.75]);
legend({'Healthy control','Premanifest HD', 'Manifest HD'},'location','best');

data = mdMat(:,indROI)*1e3;
subaxis(2,2,4);
pCorrMD = plotANOVA(data, ROI_name(indROI), statusListPlot, [0 0.95]);
ylabel('MD (x10^{-3} mm^2/s)'); ylim([0 0.95]);

pause(1); export_fig([img_root '/All_groupComp1'], '-png','-transparent');

%% Clinical correlations

set(0,'DefaultAxesFontSize', 13);

indInclude = statusList < 5;
indROI = [1:4];

CliData = [CAPS' YTO' TMS' DART' DARTtiming' Flanker' Match' SetShift'];
Cli_name = {'CAPS','YTO','TMS','DART', 'DARTtime', 'Flanker','Match','Set Shift'};

figure('position', [100 100 800 400]);
subaxis(2,2,1, 'PaddingBottom',0.05);
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
axis tight; colorbar off; xtickangle(30);
title('FA');

subaxis(2,2,4);
ImgData = [mdMat(:,indROI)];
[rho, pCorr] = correlograph(ImgData(indInclude,:), CliData(indInclude,:), ROI_name(indROI), Cli_name); 
axis tight; colorbar off; xtickangle(30);
title('MD');

pause(1); export_fig([img_root '/All_corrClinical'], '-png','-transparent');

%% plot Correlations

set(0,'DefaultAxesFontSize', 14);

indInclude = statusList < 5;

figure('position', [100 100 800 1000]);

subplot(3,2,1);
score = CAPS'; xStr = 'CAPS';
ROIselect = 1;
data = volMatAC(:,ROIselect,1);
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr);
xlabel([ROI_name{ROIselect} ' corrected volume (ml)']);
text(2.8,1.4,sprintf('r = %.2f, p < 0.0001',r), 'FontSize',12);
title(''); axis tight

subplot(3,2,2);
score = TMS'; xStr = 'TMS';
ROIselect = 1;
data = volMatAC(:,ROIselect,1);
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr);
xlabel([ROI_name{ROIselect} ' corrected volume (ml)']);
text(1.9,-5,sprintf('r = %.2f, p < 0.0001',r), 'FontSize',12);
legend({'','','HC','PM','Manifest'},'location','best','box','off','FontSize',11);
title(''); axis tight

subplot(3,2,3);
score = TMS'; xStr = 'TMS';
ROIselect = 1;
data = mdMat(:,ROIselect,1)*1e3;
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr);
xlabel([ROI_name{ROIselect} ' MD (x10^{-3} mm^2/s)']);
text(0.5,32,sprintf('r = %.2f, p = %.4f',r,p), 'FontSize',12);
title(''); axis tight

subplot(3,2,5);
score = Match'; xStr = 'Match Score';
ROIselect = 2;
data = qsmMatAC(:,ROIselect,1);
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr);
xlabel([ROI_name{ROIselect} ' susceptibility (ppm)']);
text(0.01,70,sprintf('r = %.2f, p = %.4f',r,p), 'FontSize',12);
title(''); axis tight
legend('off');

subplot(3,2,4);
score = Match'; xStr = 'Match Score';
ROIselect = 1;
data = mdMat(:,ROIselect,1)*1e3;
[r,p] = plotCorr(score(indInclude), data(indInclude), ROI_name(ROIselect), ...
    indHC(indInclude), indPMfar(indInclude), indPMnear(indInclude), indEM(indInclude) | indMan(indInclude), xStr);
ylabel(xStr);
xlabel([ROI_name{ROIselect} ' MD (x10^{-3} mm^2/s)']);
text(0.5,30,sprintf('r = %.2f, p < 0.0001',r), 'FontSize',12);
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
text(0.12,30,sprintf('r = %.2f, p < 0.0001',r), 'FontSize',12);
title(''); axis tight
legend('off');

pause(0.1); export_fig([img_root '/Example_corrClinical'], '-png','-transparent');

%% correlation between metrics - multiple correction for all

indInclude = statusList > 1;
metName = {'Vol.','Susc.','FA','MD'};

indROI = 1:4;

ImgData = [volMatAC(:,indROI) qsmMatAC(:,indROI) faMat(:,indROI) mdMat(:,indROI)];

figure('position', [100 100 800 200]);
[rho_all, pCorr] = correlographMRIall(ImgData, ROI_name(indROI), indInclude, metName);
colorbar;

pause(0.1); export_fig([img_root '/All_corrImg2'], '-png','-transparent');

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

function [pCorr, pAll] = plotANOVA(data, ROI_name, statusList3, ylines, statsFlag)

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
        pStatCorr = pCorr(:,ii);
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

function [R,P, F,gof] = plotCorr(data, score, ROI_name, indHC, indPMfar, indPMnear, indM, xStr, plotType)

if nargin < 9
    plotType = 'poly1';
end

xrange = max(score) - min(score);
xplot = linspace(min(score)-0.1*xrange, max(score)+0.1*xrange, 50); hold on;
[F,gof] = fit(score(~isnan(data(:)) & ~isnan(score(:))),data(~isnan(data(:)) & ~isnan(score(:))),plotType);
ci = predint(F, xplot, 0.95, 'functional','off');
h = plot(xplot, F(xplot)); h.Color = [0.5 0.5 0.5];
p = fill([xplot xplot(end:-1:1)],[ci(:,1); ci(end:-1:1,2)],'k');
p.FaceColor = 'k'; 
p.FaceAlpha = 0.1;
p.EdgeColor = 'none';  

set(gca,'ColorOrderIndex',1);

scatter(score(indHC), data(indHC), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
scatter(score(indPMfar | indPMnear), data(indPMfar | indPMnear), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
scatter(score(indM), data(indM), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
% scatter(score(indEM2), data(indEM2,pp), 75, 'o', 'filled', 'MarkerFaceAlpha', 0.75); hold on;
xlabel(ROI_name);
ylabel(xStr);

[R,P] = corr(score(~isnan(data)),data(~isnan(data)),'Type','Pearson','Rows','complete');
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
    subplot(1,4,nn);
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
    if mod(nn,4) ~= 1; yticks([]); end
end

end
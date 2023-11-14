clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/'));

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

img_root = '/home/jyao3/030_QSM/img_temp/';

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220209.xlsx','Sheet','NoRep');
statusList = nan(1,length(T.status_reclass));
statusList(strcmp(T.status_reclass,'HC')) = 1;
statusList(strcmp(T.status_reclass,'PM')) = 2;
statusList(strcmp(T.status_reclass,'EM')) = 3;
% statusList(strcmp(T.status_reclass,'Manifest')) = 4;

subjList = T.b_num; subjList(isnan(statusList)) = [];
examList = T.t_num; examList(isnan(statusList)) = [];
typeList = T.status_reclass; typeList(isnan(statusList)) = [];
ageList = T.age; ageList(isnan(statusList)) = [];
sexList = T.sex; sexList(isnan(statusList)) = [];
cagList = T.CAG; cagList(isnan(statusList)) = [];
pairList = T.PairID; pairList(isnan(statusList)) = [];
DTIvalidList = T.DTIReg>0; DTIvalidList(isnan(statusList)) = [];

statusList(isnan(statusList)) = [];
[pairList,i] = sort(pairList);

subjList = subjList(i);
examList = examList(i);
typeList = typeList(i);
ageList = ageList(i);
sexList = sexList(i);
cagList = cagList(i);
statusList = statusList(i);
DTIvalidList = DTIvalidList(i);

matout_root = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/Data';

indHC = find(strcmp(typeList,'HC'));
indPM = find(strcmp(typeList,'PM'));
indEM = find(strcmp(typeList,'EM'));

fprintf('PM: Age mean %.3f SD %.3f \n', mean(ageList(indPM)), std(ageList(indPM)));
fprintf('M %i F %i \n', sum(strcmp(sexList(indPM),'M')), sum(strcmp(sexList(indPM),'F')));
fprintf('CAG mean %.3f SD %.3f \n', mean(cagList(indPM)), std(cagList(indPM)));
fprintf('CAG range %.3f - %.3f \n', min(cagList(indPM)), max(cagList(indPM)));

fprintf('EM: Age mean %.3f SD %.3f \n', mean(ageList(indEM)), std(ageList(indEM)));
fprintf('M %i F %i \n', sum(strcmp(sexList(indEM),'M')), sum(strcmp(sexList(indEM),'F')));
fprintf('CAG mean %.3f SD %.3f \n', mean(cagList(indEM)), std(cagList(indEM)));
fprintf('CAG range %.3f - %.3f \n', min(cagList(indEM)), max(cagList(indEM)));

%% load stats on ROI of interest

ROIlist = 1:213;

meanROI = nan(length(subjList), length(ROIlist), 11);
stdROI = nan(length(subjList), length(ROIlist), 11);
medianROI = nan(length(subjList), length(ROIlist), 11);
madROI = nan(length(subjList), length(ROIlist), 11);
for ii = 1:length(subjList)
    
    exam_id = [subjList{ii} '_' examList{ii}];
    fprintf('# Reading subj %s \n', exam_id);
    
    dataPath = [matout_root '/' exam_id '_erode.mat'];
    load(dataPath, 'QSMstats');
    for qq = 1:10
        medianROI(ii,:,qq) = table2array(QSMstats(qq).QSMtable(ROIlist,'ROImedian'));
    end
    
    if DTIvalidList(ii)
        dataPath = [matout_root '/' exam_id '_FA_erode.mat'];
        load(dataPath, 'DTIstats');
        medianROI(ii,:,11) = table2array(DTIstats.QSMtable(ROIlist,'ROImedian'));
    end
    
end

ROIname = table2cell(QSMstats(1).QSMtable(ROIlist,'ROIname'));

%% QSM list

QSMname_list =  {'FANSI','HDQSM','iLSQR','MEDI','QSIP',...
    'QSMGAN','QSMnet+','SSTGV','SSTV','STARQSM'};
medianROI(:,:,6) = medianROI(:,:,6)/0.5684;

%% plot results - all

% [177 179 181 187 189 210 208 115]; % 208 115 183 127
% 153 155 111 192 194
% {'CN' 'PU' 'GP' 'RN' 'DN' 'SN' 'TH' 'CST'}; %  'TH' 'CST' 'SNpr' 'PLIC'
% 'IFO' 'UNC' 'GCC' 'TLA' 'TLM'
indList = [177 179 127 153 111];
indname = {'CN' 'PU' 'PLIC' 'IFO' 'GCC'};

ValidInd = true(size(pairList));
indQSM = 7;

figure('position', [100 100 1200 500]);
for rr = 1:length(indList)
    ROIind = indList(rr); % [177 179 181 187 189 210 208 115]
    ROInameT = indname{rr};
    
    % QSM
    data_all = mean(squeeze(medianROI(ValidInd,ROIind:ROIind+1,indQSM)),2) - squeeze(medianROI(ValidInd,212,indQSM));
    ymin = min(data_all(:)); ymax = max(data_all(:)); yrange = ymax - ymin;
    ymin = ymin - 0.05*yrange; ymax = ymax + 0.05*yrange;
    
    data = mean(squeeze(medianROI(:,ROIind:ROIind+1,indQSM)),2) - squeeze(medianROI(:,212,indQSM));
    F = fit(ageList(indHC),data(indHC),'poly1');
    
    data_ageCorr = data - F.p1*(ageList-median(ageList));
    
    subplot(2,5,rr);
    plotROI_ANOVA(data_ageCorr(ValidInd), statusList(ValidInd), ymin, ymax, yrange, {}, ROInameT, 1);
    title(['QSM-' indname{rr}]);
    
    % DTI
    data_all = mean(squeeze(medianROI(ValidInd,ROIind:ROIind+1,11)),2);
    ymin = min(data_all(:)); ymax = max(data_all(:)); yrange = ymax - ymin;
    ymin = ymin - 0.05*yrange; ymax = ymax + 0.05*yrange;
    
    data = mean(squeeze(medianROI(:,ROIind:ROIind+1,11)),2);
    F = fit(ageList(statusList == 1 & DTIvalidList'),data(statusList == 1 & DTIvalidList'),'poly1');
    
    data_ageCorr = data - F.p1*(ageList-median(ageList));
    
    subplot(2,5,rr+5);
    plotROI_ANOVA(data_ageCorr(ValidInd), statusList(ValidInd), ymin, ymax, yrange, {}, ROInameT, 1);
    ylabel([ROInameT ' FA']);
    title(['FA-' indname{rr}]);
end

export_fig([img_root 'QSM-FA-ROIcomp_' QSMname_list{indQSM}], '-png','-transparent'); close;

%% output ANOVA test stats

% group HD together
statusListHD = statusList > 1; % 1 HD 0 HC

% corrected alpha
alpha = 1-(1-0.05)^(1/3);
z = @(p) -sqrt(2) * erfcinv(p*2);
zcritical = z(1-alpha);

medianROI_P = zeros(length(ROIlist),11,5,2);
medianROI_T = zeros(length(ROIlist),11,5,2);
medianROIcorr = medianROI;
for rr = 1:length(ROIlist)
    % median
    try
        % QSM
        for qq = 1:10
            data = squeeze(medianROI(:,rr,qq)) - squeeze(medianROI(:,212,qq));
            F = fit(ageList(indHC),data(indHC),'poly1');
            
            % no correction
            [p,tbl,stats] = kruskalwallis(data, statusList, 'off');
            c = multcompare(stats,'Display','off','CType','dunn-sidak');
            medianROI_P(rr,qq,1,1) = p;
            medianROI_T(rr,qq,1,1) = tbl{2,5};
            medianROI_P(rr,qq,2:4,1) = c(:,6);
            [~,~,stats] = ranksum(data(indPM),data(indHC));
            medianROI_T(rr,qq,2,1) = stats.zval;
            [~,~,stats] = ranksum(data(indEM),data(indHC));
            medianROI_T(rr,qq,3,1) = stats.zval;
            [p,~,stats] = ranksum(data(statusListHD),data(~statusListHD));
            medianROI_T(rr,qq,5,1) = stats.zval;
            medianROI_P(rr,qq,5,1) = p;
            
            % linear correction
            data_ageCorr = data - F.p1*(ageList-median(ageList));
            medianROIcorr(:,rr,qq) = data_ageCorr;
            [p,tbl,stats] = kruskalwallis(data_ageCorr, statusList, 'off');
            c = multcompare(stats,'Display','off','CType','dunn-sidak');
            medianROI_P(rr,qq,1,2) = p;
            medianROI_T(rr,qq,1,2) = tbl{2,5};
            medianROI_P(rr,qq,2:4,2) = c(:,6);
            [~,~,stats] = ranksum(data_ageCorr(indPM),data_ageCorr(indHC));
            medianROI_T(rr,qq,2,2) = stats.zval;
            [~,~,stats] = ranksum(data_ageCorr(indEM),data_ageCorr(indHC));
            medianROI_T(rr,qq,3,2) = stats.zval;
            [p,~,stats] = ranksum(data_ageCorr(statusListHD),data_ageCorr(~statusListHD));
            medianROI_T(rr,qq,5,2) = stats.zval;
            medianROI_P(rr,qq,5,2) = p;
        end
        
        % DTI
        data = squeeze(medianROI(:,rr,11));
        F = fit(ageList(statusList == 1 & DTIvalidList'),data(statusList == 1 & DTIvalidList'),'poly1');
        
        % no correction
        [p,tbl,stats] = kruskalwallis(data, statusList, 'off');
        c = multcompare(stats,'Display','off','CType','dunn-sidak');
        medianROI_P(rr,11,1,1) = p;
        medianROI_T(rr,11,1,1) = tbl{2,5};
        medianROI_P(rr,11,2:4,1) = c(:,6);
        [~,~,stats] = ranksum(data(indPM),data(indHC));
        medianROI_T(rr,11,2,1) = stats.zval;
        [~,~,stats] = ranksum(data(indEM),data(indHC));
        medianROI_T(rr,11,3,1) = stats.zval;
        [p,~,stats] = ranksum(data(statusListHD),data(~statusListHD));
        medianROI_T(rr,11,5,1) = stats.zval;
        medianROI_P(rr,11,5,1) = p;
        
        % linear correction
        data_ageCorr = data - F.p1*(ageList-median(ageList));
        medianROIcorr(:,rr,11) = data_ageCorr;
        [p,tbl,stats] = kruskalwallis(data_ageCorr, statusList, 'off');
        c = multcompare(stats,'Display','off','CType','dunn-sidak');
        medianROI_P(rr,11,1,2) = p;
        medianROI_T(rr,11,1,2) = tbl{2,5};
        medianROI_P(rr,11,2:4,2) = c(:,6);
        [~,~,stats] = ranksum(data_ageCorr(indPM),data_ageCorr(indHC));
        medianROI_T(rr,11,2,2) = stats.zval;
        [~,~,stats] = ranksum(data_ageCorr(indEM),data_ageCorr(indHC));
        medianROI_T(rr,11,3,2) = stats.zval;
        [p,~,stats] = ranksum(data_ageCorr(statusListHD),data_ageCorr(~statusListHD));
        medianROI_T(rr,11,5,2) = stats.zval;
        medianROI_P(rr,11,5,2) = p;
        
    catch
        medianROI_P(rr,:,:,:) = 1;
        medianROI_T(rr,:,:,:) = 0;
    end
end

%% plot heatmap

set(0,'defaultAxesFontSize',12);

cmap1 = flip(colormap_RWG([-1 1],100));
r = ones(100, 1); g = linspace(1, 0, 100)';
cmap2 = [r g g];

dataT = medianROI_T(:,:,:,1); dataT(isnan(dataT)) = 0;
dataP = medianROI_P(:,:,:,1); dataP(isnan(dataP)) = 1;

figure('position', [100 100 600 900]);
plotHeatMap_HD(dataT, dataP, cmap1, cmap2)
export_fig([img_root 'QSMFA-noCorr'], '-png','-transparent');

dataT = medianROI_T(:,:,:,2); dataT(isnan(dataT)) = 0;
dataP = medianROI_P(:,:,:,2); dataP(isnan(dataP)) = 1;

figure('position', [100 100 600 900]);
plotHeatMap_HD(dataT, dataP, cmap1, cmap2)
export_fig([img_root 'QSMFA-linCorr'], '-png','-transparent');

%% sort HC-PM

indQSM = 10;

set(0,'defaultAxesFontSize',12);
set(0,'defaultLineLineWidth',1);

dataT = medianROI_T(:,:,:,2); dataT(isnan(dataT)) = 0;
dataP = medianROI_P(:,:,:,2); dataP(isnan(dataP)) = 1;

indQSM = 10;
medianPM_T = mean(dataT(:,indQSM,2),2);
madPM_R = std(dataT(:,indQSM,2),0,2);
medianEM_T = mean(dataT(:,indQSM,3),2);
madEM_R = std(dataT(:,indQSM,3),0,2);
medianHD_T = mean(dataT(:,indQSM,5),2);
madHD_R = std(dataT(:,indQSM,5),0,2);

medianPM_FA = mean(dataT(:,11,2),2);
madPM_FA = std(dataT(:,11,2),0,2);
medianEM_FA = mean(dataT(:,11,3),2);
madEM_FA = std(dataT(:,11,3),0,2);
medianHD_FA = mean(dataT(:,11,5),2);
madHD_FA = std(dataT(:,11,5),0,2);

[median_Rsort,i] = sort(medianPM_T);
ROIname_sort = ROIname(i);

THD = table(ROIname_sort,median_Rsort);

invalidInd = isnan(medianPM_T);

figure('position', [100 100 1000 500]);
subplot(211);
H1 = shadedErrorBar(1:length(medianPM_T(~invalidInd)), medianPM_T(~invalidInd), ...
    madPM_R(~invalidInd),'lineProps','-b'); hold on;
H2 = shadedErrorBar(1:length(medianEM_T(~invalidInd)), medianEM_T(~invalidInd), ...
    madEM_R(~invalidInd),'lineProps','-b');
H1.mainLine.Color = [0 0.4470 0.7410];
H1.edge(1).Color = 'none'; H1.edge(2).Color = 'none';
H1.patch.FaceColor = [0 0.4470 0.7410];
H2.mainLine.Color = [0.8500 0.3250 0.0980];
H2.edge(1).Color = 'none'; H2.edge(2).Color = 'none';
H2.patch.FaceColor = [0.8500 0.3250 0.0980];
gridPlot(zcritical);
ylabel('HC-PM-EM z-value');
title('QSM');
legend([H1.mainLine,H2.mainLine],{'HC-PM','HC-EM'},'box','off','location','southwest');

subplot(212);
H1 = shadedErrorBar(1:length(medianPM_FA(~invalidInd)), medianPM_FA(~invalidInd), ...
    madPM_FA(~invalidInd),'lineProps','-b'); hold on;
H2 = shadedErrorBar(1:length(medianEM_FA(~invalidInd)), medianEM_FA(~invalidInd), ...
    madEM_FA(~invalidInd),'lineProps','-b');
H1.mainLine.Color = [0 0.4470 0.7410];
H1.edge(1).Color = 'none'; H1.edge(2).Color = 'none';
H1.patch.FaceColor = [0 0.4470 0.7410];
H2.mainLine.Color = [0.8500 0.3250 0.0980];
H2.edge(1).Color = 'none'; H2.edge(2).Color = 'none';
H2.patch.FaceColor = [0.8500 0.3250 0.0980];
gridPlot(zcritical);
ylabel('HC-PM-EM z-value');
title('FA');
legend([H1.mainLine,H2.mainLine],{'HC-PM','HC-EM'},'box','off','location','southwest');

export_fig([img_root 'QSM-FA-z'], '-png','-transparent');

%% QSM vs. FA

figure('position', [100 100 800 800]);

subplot(221);
scatter(medianPM_T(1:108),medianPM_FA(1:108),50,'o','filled','MarkerFaceAlpha',0.5); hold on;
scatter(medianPM_T(109:176),medianPM_FA(109:176),50,'o','filled','MarkerFaceAlpha',0.5);
scatter(medianPM_T(177:190),medianPM_FA(177:190),50,'o','filled','MarkerFaceAlpha',0.5);
scatter(medianPM_T(191:206),medianPM_FA(191:206),50,'o','filled','MarkerFaceAlpha',0.5);
xlabel('QSM HC-PM z-stat'); xlim([-5 5]);
ylabel('FA HC-PM z-stat'); ylim([-5 5]);
plot([0 0],[-5 5],'k:');
plot([-5 5],[0 0],'k:');

subplot(222);
scatter(medianEM_T(1:108),medianEM_FA(1:108),50,'o','filled','MarkerFaceAlpha',0.5); hold on;
scatter(medianEM_T(109:176),medianEM_FA(109:176),50,'o','filled','MarkerFaceAlpha',0.5);
scatter(medianEM_T(177:190),medianEM_FA(177:190),50,'o','filled','MarkerFaceAlpha',0.5);
scatter(medianEM_T(191:206),medianEM_FA(191:206),50,'o','filled','MarkerFaceAlpha',0.5);
xlabel('QSM HC-EM z-stat'); xlim([-5 5]);
ylabel('FA HC-EM z-stat'); ylim([-5 5]);
plot([0 0],[-5 5],'k:');
plot([-5 5],[0 0],'k:');

legend({'Cortical GM','WM','BG','TH'},'box','off','location','best');

subplot(223);
scatter(medianPM_T(109:176),medianPM_FA(109:176),50,'o','filled','MarkerFaceAlpha',0.5); hold on;
scatter(medianEM_T(109:176),medianEM_FA(109:176),50,'o','filled','MarkerFaceAlpha',0.5);
xlabel('QSM WM z-stat'); xlim([-5 5]);
ylabel('FA WM z-stat'); ylim([-5 5]);

xFit = linspace(-5,5,100);
[F,gof] = fit([medianPM_T(109:176); medianEM_T(109:176)],...
    [medianPM_FA(109:176); medianEM_FA(109:176)],'poly1');
plot(xFit,F(xFit),'-','Color',[1 1 1]*0.5);
plot([0 0],[-5 5],'k:');
plot([-5 5],[0 0],'k:');

[R,P] = corr([medianPM_T(109:176); medianEM_T(109:176)],...
    [medianPM_FA(109:176); medianEM_FA(109:176)]);

legend({'HC-PM','HC-EM'},'box','off','location','best');
title(sprintf('Pearson R %.2f P %.4f',R,P));

% QSM vs FA patient-wise
subplot(224);
scatter(reshape(medianROIcorr(:,1:108,indQSM),1,[]), reshape(medianROIcorr(:,1:108,11),1,[]),'.'); hold on;
scatter(reshape(medianROIcorr(:,109:176,indQSM),1,[]), reshape(medianROIcorr(:,109:176,11),1,[]),'.');
scatter(reshape(medianROIcorr(:,177:190,indQSM),1,[]), reshape(medianROIcorr(:,177:190,11),1,[]),'.');
scatter(reshape(medianROIcorr(:,191:206,indQSM),1,[]), reshape(medianROIcorr(:,191:206,11),1,[]),'.');
legend({'Cortical GM','WM','BG','TH'},'box','off','location','best');
xlabel('QSM (ppm)'); ylabel('FA');

% export_fig([img_root 'QSM-FA-corr'], '-png','-transparent');

% subplot(222);
% scatter(mean(medianROIcorr(indHC,153:154,1),2), ...
%     mean(medianROIcorr(indHC,153:154,2),2),50,'o','filled','MarkerFaceAlpha',0.5); hold on;
% scatter(mean(medianROIcorr(indPM,153:154,1),2), ...
%     mean(medianROIcorr(indPM,153:154,2),2),50,'o','filled','MarkerFaceAlpha',0.5);
% scatter(mean(medianROIcorr(indEM,153:154,1),2), ...
%     mean(medianROIcorr(indEM,153:154,2),2),50,'o','filled','MarkerFaceAlpha',0.5);
% legend({'HC','PM','EM'},'box','off','location','best');
% xlabel('QSM (ppm)'); ylabel('FA');

%% FA vs. QSM in WM in HC

indQSM = 10;

figure('position', [100 100 800 800]);
subplot(221);
dataQSM_HC_WM = medianROIcorr(:, 109:176, indQSM);
dataQSM_HC_WM = averageLR_WM(dataQSM_HC_WM);
dataFA_HC_WM = medianROIcorr(:, 109:176, 11);
dataFA_HC_WM = averageLR_WM(dataFA_HC_WM);
scatter(dataQSM_HC_WM(:), dataFA_HC_WM(:), '.'); hold on;
xlim([-0.1 0.05]);
xlabel('QSM susc.'); ylabel('FA');
title('All subjects');

indValid = ~isnan(dataQSM_HC_WM(:)) & ~isnan(dataFA_HC_WM(:));
xFit = linspace(-0.1,0.05,100);
[F,gof] = fit(dataQSM_HC_WM(indValid), dataFA_HC_WM(indValid),'poly1');
h = plot(xFit,F(xFit),':','Color',[1 1 1]*0.5);
[R,P] = corr(dataQSM_HC_WM(indValid), dataFA_HC_WM(indValid));
legend(h,{sprintf('Fitting R %.2f P %.2e \n',R,P)},'box','off','location','best');

subplot(223);
errorbarxy(nanmean(dataQSM_HC_WM,1), nanmean(dataFA_HC_WM,1), ...
    nanstd(dataQSM_HC_WM,1), nanstd(dataFA_HC_WM,1),{'ko','k','k'}); hold on;

xFit = linspace(-0.1,0.05,100);
[F,gof] = fit(nanmean(dataQSM_HC_WM,1)', nanmean(dataFA_HC_WM,1)','poly1');
h = plot(xFit,F(xFit),':','Color',[1 1 1]*0.5);
[R,P] = corr(nanmean(dataQSM_HC_WM,1)', nanmean(dataFA_HC_WM,1)');
legend(h,{sprintf('Fitting R %.2f P %.2e \n',R,P)},'box','off','location','best');
xlabel('QSM susc.'); ylabel('FA');
title('Mean±SD of all subjects');

subplot(222);
dataQSM_HC_WM = medianROIcorr(indHC, 109:176, indQSM);
dataQSM_HC_WM = averageLR_WM(dataQSM_HC_WM);
dataFA_HC_WM = medianROIcorr(indHC, 109:176, 11);
dataFA_HC_WM = averageLR_WM(dataFA_HC_WM);
scatter(dataQSM_HC_WM(:), dataFA_HC_WM(:), '.'); hold on;
xlim([-0.1 0.05]);
xlabel('QSM susc.'); ylabel('FA');
title('HC subjects');

indValid = ~isnan(dataQSM_HC_WM(:)) & ~isnan(dataFA_HC_WM(:));
xFit = linspace(-0.1,0.05,100);
[F,gof] = fit(dataQSM_HC_WM(indValid), dataFA_HC_WM(indValid),'poly1');
h = plot(xFit,F(xFit),':','Color',[1 1 1]*0.5);
[R,P] = corr(dataQSM_HC_WM(indValid), dataFA_HC_WM(indValid));
legend(h,{sprintf('Fitting R %.2f P %.2e \n',R,P)},'box','off','location','best');

subplot(224);
errorbarxy(nanmean(dataQSM_HC_WM,1), nanmean(dataFA_HC_WM,1), ...
    nanstd(dataQSM_HC_WM,1), nanstd(dataFA_HC_WM,1),{'ko','k','k'}); hold on;

xFit = linspace(-0.1,0.05,100);
[F,gof] = fit(nanmean(dataQSM_HC_WM,1)', nanmean(dataFA_HC_WM,1)','poly1');
h = plot(xFit,F(xFit),':','Color',[1 1 1]*0.5);
[R,P] = corr(nanmean(dataQSM_HC_WM,1)', nanmean(dataFA_HC_WM,1)');
legend(h,{sprintf('Fitting R %.2f P %.2e \n',R,P)},'box','off','location','best');
xlabel('QSM susc.'); ylabel('FA');
title('Mean±SD of HC subjects');

export_fig([img_root 'QSM-FA-corr_' QSMname_list{indQSM}], '-png','-transparent');

%% QSM vs. FA for all methods

dataFA_HC_WM = medianROIcorr(indHC, 109:176, 11);
dataFA_HC_WM = averageLR_WM(dataFA_HC_WM);

figure('position', [100 100 1600 800]);
for indQSM = 1:10
    dataQSM_HC_WM = medianROIcorr(indHC, 109:176, indQSM);
    dataQSM_HC_WM = averageLR_WM(dataQSM_HC_WM);
    
    subplot(2,5,indQSM);
    errorbarxy(nanmean(dataQSM_HC_WM,1), nanmean(dataFA_HC_WM,1), ...
        nanstd(dataQSM_HC_WM,1), nanstd(dataFA_HC_WM,1),{'ko','k','k'}); hold on;
    
    xFit = linspace(-0.1,0.05,100);
    [F,gof] = fit(nanmean(dataQSM_HC_WM,1)', nanmean(dataFA_HC_WM,1)','poly1');
    h = plot(xFit,F(xFit),':','Color',[1 1 1]*0.5);
    [R,P] = corr(nanmean(dataQSM_HC_WM,1)', nanmean(dataFA_HC_WM,1)');
    legend(h,{sprintf('R %.2f P %.2e \n',R,P)},'box','off','location','north');
    xlabel([QSMname_list{indQSM} ' QSM']); ylabel('FA');
    ylim([0 0.6]);
    title('Mean±SD of HC subjects');
end

export_fig([img_root 'QSM-FA-corr_mean_' QSMname_list{indQSM}], '-png','-transparent');

figure('position', [100 100 1600 800]);
for indQSM = 1:10
    dataQSM_HC_WM = medianROIcorr(indHC, 109:176, indQSM);
    dataQSM_HC_WM = averageLR_WM(dataQSM_HC_WM);
    
    subplot(2,5,indQSM);
    scatter(dataQSM_HC_WM(:), dataFA_HC_WM(:), '.'); hold on;
    
    indValid = ~isnan(dataQSM_HC_WM(:)) & ~isnan(dataFA_HC_WM(:));
    xFit = linspace(-0.1,0.05,100);
    [F,gof] = fit(dataQSM_HC_WM(indValid), dataFA_HC_WM(indValid),'poly1');
    h = plot(xFit,F(xFit),':','Color',[1 1 1]*0.5);
    [R,P] = corr(dataQSM_HC_WM(indValid), dataFA_HC_WM(indValid));
    legend(h,{sprintf('R %.2f P %.2e \n',R,P)},'box','off','location','north');
    xlabel([QSMname_list{indQSM} ' QSM']); ylabel('FA');
    ylim([0 0.6]);
    title('All HC subjects');
end

export_fig([img_root 'QSM-FA-corr_data_' QSMname_list{indQSM}], '-png','-transparent');

%% functions

function [dataAve] = averageLR_WM(data)

dataAve = cat(2,data(:,1:6),0.5*(data(:,7:2:end)+data(:,8:2:end)));

end

function [] = gridPlot(zcritical)

plot([82.5 82.5], [-4 4], ':', 'Color', [1 1 1]*0.8);
plot([108.5 108.5], [-4 4], 'Color', [1 1 1]*0.8);
plot([176.5 176.5], [-4 4], 'Color', [1 1 1]*0.8);
plot([190.5 190.5], [-4 4], 'Color', [1 1 1]*0.8);
plot([0 210],[0 0],'k-');
plot([0 210],[zcritical zcritical],'k:');
plot([0 210],[-zcritical -zcritical],'k:');
xlim([0 206]);
ylim([-4 4]);
ax = gca;
ax.XTick = floor([(1+109)/2 (109+176)/2 (176+191)/2 (191+206)/2]);
ax.XTickLabelRotation = 0;
ax.XTickLabel = {'Cortical GM','WM','BG','TH'};

end

function [] = plotROI_ANOVA(data, statusList, ymin, ymax, yrange, QSMname, ROIname, ii)

violinplot(data, statusList, 'showMean', true);

[p,~,stats] = kruskalwallis(data, statusList, 'off');
if p < 0.0001
    text(0.6,ymin+0.1*yrange,sprintf('p < 0.0001'),'FontSize',12); % -0.03
elseif p < 0.001
    text(0.6,ymin+0.1*yrange,sprintf('p < 0.001'),'FontSize',12);
else
    text(0.6,ymin+0.1*yrange,sprintf('p = %.3f',p),'FontSize',12);
end

c = multcompare(stats,'Display','off','CType','dunn-sidak');
title([strrep(strrep(QSMname,'_','-'),'-meanEcho','')]);

ylim([ymin ymax]); % [-0.05 0.15] [0 0.2]
if mod(ii,10) == 1
    ylabel([ROIname ' Susc. (ppm)']);
else
    ylabel('')
end
xlim([0.5 3.5]);
ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel = {'HC','PM','EM'};
sigstar({[1,2],[1,3],[2,3]},c(:,6));

end

function [] = plotHeatMap(data, cmap, colorrange)

imagesc(data(1:209,:)',colorrange); colormap(cmap); colorbar; hold on;
plot([109.5 109.5],[0 11],'k','LineWidth',2);
plot([176.5 176.5],[0 11],'k','LineWidth',2);
plot([190.5 190.5],[0 11],'k','LineWidth',2);
ax = gca;
ax.YTick = 1:11;
ax.YTickLabel = {'QSM-FANSI','QSM-HDQSM','QSM-iLSQR','QSM-MEDI','QSM-QSIP',...
    'QSM-QSMGAN','QSM-QSMnet+','QSM-SSTGV','QSM-SSTV','QSM-STARQSM','FA'};
ax.YTickLabelRotation = 0;
ax.XTick = floor([(1+109)/2 (109+176)/2 (176+191)/2 (191+210)/2]);
ax.XTickLabelRotation = 0;
ax.XTickLabel = {'Cortical GM','WM','BG','TH'};

end

function [] = plotHeatMap_HD(dataT, dataP, cmap1, cmap2)

ax(1) = subplot(411);
plotHeatMap(dataT(:,:,2), cmap1, [-4 4]);
title('Median QSM/FA HC-PM Z-value');
ax(2) = subplot(412);
plotHeatMap(dataP(:,:,2), flip(cmap2), [0 0.05]);
title('Median QSM/FA HC-PM P-value');
ax(3) = subplot(413);
plotHeatMap(dataT(:,:,3), cmap1, [-4 4]);
title('Median QSM/FA HC-EM Z-value');
ax(4) = subplot(414);
plotHeatMap(dataP(:,:,3), flip(cmap2), [0 0.05]);
title('Median QSM/FA HC-EM P-value');

colormap(ax(1), cmap1); colormap(ax(2), flip(cmap2));
colormap(ax(3), cmap1); colormap(ax(4), flip(cmap2));

end
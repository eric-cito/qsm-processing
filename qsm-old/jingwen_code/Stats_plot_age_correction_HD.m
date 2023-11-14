clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/'));

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220103.xlsx','Sheet','NoRep');
statusList = nan(1,length(T.status_reclass));
statusList(strcmp(T.status_reclass,'HC')) = 1;
statusList(strcmp(T.status_reclass,'PM')) = 2;
statusList(strcmp(T.status_reclass,'EM')) = 3;

subjList = T.b_num; subjList(isnan(statusList)) = [];
examList = T.t_num; examList(isnan(statusList)) = [];
typeList = T.status_reclass; typeList(isnan(statusList)) = [];
ageList = T.age; ageList(isnan(statusList)) = [];
sexList = T.sex; sexList(isnan(statusList)) = [];
cagList = T.CAG; cagList(isnan(statusList)) = [];
pairList = T.PairID; pairList(isnan(statusList)) = [];

statusList(isnan(statusList)) = [];
[pairList,i] = sort(pairList);

subjList = subjList(i);
examList = examList(i);
typeList = typeList(i);
ageList = ageList(i);
sexList = sexList(i);
cagList = cagList(i);
statusList = statusList(i);

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

% load
% load('QSM_age_fit.mat');

%% load stats on ROI of interest - voxCorrected

medianROI_voxCorr = zeros(length(subjList), length(ROIlist), length(QSMfile_list));
for ii = 1:length(subjList)
    
    exam_id = [subjList{ii} '_' examList{ii}];
    fprintf('# Reading subj %s \n', exam_id);
    dataPath = [matout_root '/' exam_id '_voxCorr.mat'];
    
    load(dataPath, 'QSMstats');
    for QSMnum = 1:length(QSMfile_list)
        medianROI_voxCorr(ii,:,QSMnum) = table2array(QSMstats(QSMnum).QSMtable(ROIlist,'ROImedian'));
    end
    
end

% scale the QSMGAN
medianROI_voxCorr(:,:,6) = medianROI_voxCorr(:,:,6)/0.5684;

% reorder
medianROI_voxCorr = medianROI_voxCorr(:,:,QSMorder);

%% plot results - all

% [177 179 181 187 189 210 208 115]; % 208 115 183 127 
% 153 155 111 192 194
% {'CN' 'PU' 'GP' 'RN' 'DN' 'SN' 'TH' 'CST'}; %  'TH' 'CST' 'SNpr' 'PLIC'
% 'IFO' UNC GCC TLA TLM
indList = [192 194]; 
indname = {'TLA' 'TLM'}; 

ValidInd = true(size(pairList));

for rr = 1:length(indList)
    ROIind = indList(rr); % [177 179 181 187 189 210 208 115]
    ROInameT = indname{rr};
    
    % calcultae y range
    data_all = mean(squeeze(medianROI(ValidInd,ROIind:ROIind+1,:)),2) - squeeze(medianROI(ValidInd,212,:));
    ymin = min(data_all(:)); ymax = max(data_all(:)); yrange = ymax - ymin;
    ymin = ymin - 0.05*yrange; ymax = ymax + 0.05*yrange;
    
    figure('position', [100 100 1800 1000]);
    for ii = 1:10
        
        data = mean(squeeze(medianROI(:,ROIind:ROIind+1,ii)),2) - squeeze(medianROI(:,212,ii));
        F = fit(ageList(indHC),data(indHC),'poly1');
        
        data_ageCorr = data - F.p1*(ageList-median(ageList));
        data_ageDiv = data./ageList*median(ageList);
        data_voxCorr = mean(squeeze(medianROI_voxCorr(:,ROIind:ROIind+1,ii)),2);
        
        QSMname = QSMfile_list{ii};
        
        subplot(4,10,ii);
        plotROI_ANOVA(data(ValidInd), statusList(ValidInd), ymin, ymax, yrange, QSMname, ROInameT, ii);
        subplot(4,10,ii+10);
        plotROI_ANOVA(data_ageCorr(ValidInd), statusList(ValidInd), ymin, ymax, yrange, QSMname, ROInameT, ii);
        title('');
        subplot(4,10,ii+20);
        plotROI_ANOVA(data_ageDiv(ValidInd), statusList(ValidInd), ymin, ymax, yrange, QSMname, ROInameT, ii);
        title('');
        subplot(4,10,ii+30);
        plotROI_ANOVA(data_voxCorr(ValidInd), statusList(ValidInd), ymin, ymax, yrange, QSMname, ROInameT, ii);
        title('');

    end
    export_fig(['HD-' indname{rr} '_compCorr'], '-png','-transparent'); % close;
end

%% plot results - paired

% indList = [183]; % [177 179 181 187 189 210 208 115]; % 208 115
% indname = {'SNpr'}; % {'CN' 'PU' 'GP' 'RN' 'DN' 'SN' 'TH' 'CST'}; %  'TH' 'CST'

ValidInd = ~isnan(pairList); % ones(size(pairList));

for rr = 1:length(indList)
    ROIind = indList(rr); % [177 179 181 187 189 210 208 115]
    ROInameT = indname{rr};
    
    % calcultae y range
    data_all = mean(squeeze(medianROI(ValidInd,ROIind:ROIind+1,:)),2) - squeeze(medianROI(ValidInd,212,:));
    ymin = min(data_all(:)); ymax = max(data_all(:)); yrange = ymax - ymin;
    ymin = ymin - 0.05*yrange; ymax = ymax + 0.05*yrange;
    
    figure('position', [100 100 1800 1000]);
    for ii = 1:10
        
        data = mean(squeeze(medianROI(:,ROIind:ROIind+1,ii)),2) - squeeze(medianROI(:,212,ii));
        F = fit(ageList(indHC),data(indHC),'poly1');
        
        data_ageCorr = data - F.p1*(ageList-median(ageList));
        data_ageDiv = data./ageList*median(ageList);
        data_voxCorr = mean(squeeze(medianROI_voxCorr(:,ROIind:ROIind+1,ii)),2);
        
        QSMname = QSMfile_list{ii};
        
        subplot(5,10,ii);
        plotROI_ANOVA(data(ValidInd), statusList(ValidInd), ymin, ymax, yrange, QSMname, ROInameT, ii);
        subplot(5,10,ii+10);
        plotROI_ANOVA(data_ageCorr(ValidInd), statusList(ValidInd), ymin, ymax, yrange, QSMname, ROInameT, ii);
        title('');
        subplot(5,10,ii+20);
        plotROI_ANOVA(data_ageDiv(ValidInd), statusList(ValidInd), ymin, ymax, yrange, QSMname, ROInameT, ii);
        title('');
        subplot(5,10,ii+30);
        plotROI_ANOVA(data_voxCorr(ValidInd), statusList(ValidInd), ymin, ymax, yrange, QSMname, ROInameT, ii);
        title('');
        subplot(5,10,ii+40);
        plotROI_pair(data, typeList, pairList, ymin, ymax, QSMname, ROInameT, ii);
        title('');

    end
    export_fig(['HD-' indname{rr} '_compCorr_pair'], '-png','-transparent'); % close;
end

%% output ANOVA test stats

% group HD together
statusListHD = statusList > 1; % 1 HD 0 HC

% corrected alpha
alpha = 1-(1-0.05)^(1/3);
z = @(p) -sqrt(2) * erfcinv(p*2);                                    
zcritical = z(1-alpha);

medianROI_P = zeros(length(ROIlist), length(QSMfile_list),5,4);
medianROI_T = zeros(length(ROIlist), length(QSMfile_list),5,4);
for rr = 1:length(ROIlist)
    for ii = 1:10
        % median
        try
            data = squeeze(medianROI(:,rr,ii)) - squeeze(medianROI(:,212,ii));
            F = fit(ageList(indHC),data(indHC),'poly1');
            
            % no correction
            [p,tbl,stats] = kruskalwallis(data, statusList, 'off');
            c = multcompare(stats,'Display','off','CType','dunn-sidak');
            medianROI_P(rr,ii,1,1) = p;
            medianROI_T(rr,ii,1,1) = tbl{2,5};
            medianROI_P(rr,ii,2:4,1) = c(:,6);
            [~,~,stats] = ranksum(data(indPM),data(indHC));
            medianROI_T(rr,ii,2,1) = stats.zval;
            [~,~,stats] = ranksum(data(indEM),data(indHC));
            medianROI_T(rr,ii,3,1) = stats.zval;
            [p,~,stats] = ranksum(data(statusListHD),data(~statusListHD));
            medianROI_T(rr,ii,5,1) = stats.zval;
            medianROI_P(rr,ii,5,1) = p;
            
            % linear correction
            data_ageCorr = data - F.p1*(ageList-median(ageList));
            [p,tbl,stats] = kruskalwallis(data_ageCorr, statusList, 'off');
            c = multcompare(stats,'Display','off','CType','dunn-sidak');
            medianROI_P(rr,ii,1,2) = p;
            medianROI_T(rr,ii,1,2) = tbl{2,5};
            medianROI_P(rr,ii,2:4,2) = c(:,6);
            [~,~,stats] = ranksum(data_ageCorr(indPM),data_ageCorr(indHC));
            medianROI_T(rr,ii,2,2) = stats.zval;
            [~,~,stats] = ranksum(data_ageCorr(indEM),data_ageCorr(indHC));
            medianROI_T(rr,ii,3,2) = stats.zval;
            [p,~,stats] = ranksum(data_ageCorr(statusListHD),data_ageCorr(~statusListHD));
            medianROI_T(rr,ii,5,2) = stats.zval;
            medianROI_P(rr,ii,5,2) = p;
            
            % division correction
            data_ageCorr = data./ageList*median(ageList);
            [p,tbl,stats] = kruskalwallis(data_ageCorr, statusList, 'off');
            c = multcompare(stats,'Display','off','CType','dunn-sidak');
            medianROI_P(rr,ii,1,3) = p;
            medianROI_T(rr,ii,1,3) = tbl{2,5};
            medianROI_P(rr,ii,2:4,3) = c(:,6);
            [~,~,stats] = ranksum(data_ageCorr(indPM),data_ageCorr(indHC));
            medianROI_T(rr,ii,2,3) = stats.zval;
            [~,~,stats] = ranksum(data_ageCorr(indEM),data_ageCorr(indHC));
            medianROI_T(rr,ii,3,3) = stats.zval;
            [p,~,stats] = ranksum(data_ageCorr(statusListHD),data_ageCorr(~statusListHD));
            medianROI_T(rr,ii,5,3) = stats.zval;
            medianROI_P(rr,ii,5,3) = p;
            
            % voxelwise linear correction
            data_voxCorr = squeeze(medianROI_voxCorr(:,rr,ii));
            [p,tbl,stats] = kruskalwallis(data_voxCorr, statusList, 'off');
            c = multcompare(stats,'Display','off','CType','dunn-sidak');
            medianROI_P(rr,ii,1,4) = p;
            medianROI_T(rr,ii,1,4) = tbl{2,5};
            medianROI_P(rr,ii,2:4,4) = c(:,6);
            [~,~,stats] = ranksum(data_voxCorr(indPM),data_voxCorr(indHC));
            medianROI_T(rr,ii,2,4) = stats.zval;
            [~,~,stats] = ranksum(data_voxCorr(indEM),data_voxCorr(indHC));
            medianROI_T(rr,ii,3,4) = stats.zval;
            [p,~,stats] = ranksum(data_voxCorr(statusListHD),data_voxCorr(~statusListHD));
            medianROI_T(rr,ii,5,4) = stats.zval;
            medianROI_P(rr,ii,5,4) = p;
            
        catch
            medianROI_P(rr,ii,:) = 1;
            medianROI_T(rr,ii,:) = 0;
        end
    end
end

%% plot heatmap

set(0,'defaultAxesFontSize',12);

QSMname = strrep(QSMfile_list,'_','-');
QSMname = strrep(QSMname,'QSM-','');
QSMname = strrep(QSMname,'-meanEcho','');

cmap1 = flip(colormap_RWG([-1 1],100));
r = ones(100, 1); g = linspace(1, 0, 100)'; 
cmap2 = [r g g];

dataT = medianROI_T(:,:,:,1); dataT(isnan(dataT)) = 0;
dataP = medianROI_P(:,:,:,1); dataP(isnan(dataP)) = 1;

figure('position', [100 100 600 900]);
plotHeatMap_HD(dataT, dataP, cmap1, cmap2, QSMname)
export_fig('HC-HD-noCorr', '-png','-transparent');

dataT = medianROI_T(:,:,:,2); dataT(isnan(dataT)) = 0;
dataP = medianROI_P(:,:,:,2); dataP(isnan(dataP)) = 1;

figure('position', [100 100 600 900]);
plotHeatMap_HD(dataT, dataP, cmap1, cmap2, QSMname)
export_fig('HC-HD-linCorr', '-png','-transparent');

dataT = medianROI_T(:,:,:,3); dataT(isnan(dataT)) = 0;
dataP = medianROI_P(:,:,:,3); dataP(isnan(dataP)) = 1;

figure('position', [100 100 600 900]);
plotHeatMap_HD(dataT, dataP, cmap1, cmap2, QSMname)
export_fig('HC-HD-divCorr', '-png','-transparent');

dataT = medianROI_T(:,:,:,4); dataT(isnan(dataT)) = 0;
dataP = medianROI_P(:,:,:,4); dataP(isnan(dataP)) = 1;

figure('position', [100 100 600 900]);
plotHeatMap_HD(dataT, dataP, cmap1, cmap2, QSMname)
export_fig('HC-HD-voxCorr', '-png','-transparent');

%% sort HC-PM

set(0,'defaultAxesFontSize',12);
set(0,'defaultLineLineWidth',1);

dataT = medianROI_T(:,:,:,2); dataT(isnan(dataT)) = 0;
dataP = medianROI_P(:,:,:,2); dataP(isnan(dataP)) = 1;

% figure('position', [100 100 800 900]);
% plotHeatMap_HD(dataT, dataP, cmap1, cmap2, QSMname)
% export_fig('temp', '-png','-transparent');

medianPM_T = mean(dataT(:,:,2),2);
madPM_R = std(dataT(:,:,2),0,2);
medianEM_T = mean(dataT(:,:,3),2);
madEM_R = std(dataT(:,:,3),0,2);
medianHD_T = mean(dataT(:,:,5),2);
madHD_R = std(dataT(:,:,5),0,2);

[median_Rsort,i] = sort(medianPM_T);
ROIname_sort = ROIname(i);

THD = table(ROIname_sort,median_Rsort);

invalidInd = isnan(medianPM_T);

figure('position', [100 100 800 300]); 
plot([82.5 82.5], [-4 4], ':', 'Color', [1 1 1]*0.8); hold on;
plot([108.5 108.5], [-4 4], 'Color', [1 1 1]*0.8);
plot([176.5 176.5], [-4 4], 'Color', [1 1 1]*0.8);
plot([190.5 190.5], [-4 4], 'Color', [1 1 1]*0.8);
H1 = shadedErrorBar(1:length(medianPM_T(~invalidInd)), medianPM_T(~invalidInd), ...
    madPM_R(~invalidInd),'lineProps','-b'); 
H2 = shadedErrorBar(1:length(medianEM_T(~invalidInd)), medianEM_T(~invalidInd), ...
    madEM_R(~invalidInd),'lineProps','-r'); 
H3 = shadedErrorBar(1:length(medianHD_T(~invalidInd)), medianHD_T(~invalidInd), ...
    madHD_R(~invalidInd),'lineProps','-g'); 
H1.mainLine.Color = [0 0.4470 0.7410];
H1.edge(1).Color = 'none'; H1.edge(2).Color = 'none';
H1.patch.FaceColor = [0 0.4470 0.7410];
H2.mainLine.Color = [0.8500 0.3250 0.0980];
H2.edge(1).Color = 'none'; H2.edge(2).Color = 'none';
H2.patch.FaceColor = [0.8500 0.3250 0.0980];
H3.mainLine.Color = [0.4940 0.1840 0.5560]; 
H3.edge(1).Color = 'none'; H3.edge(2).Color = 'none';
H3.patch.FaceColor = [0.4940 0.1840 0.5560];
plot([0 210],[0 0],'k-');
plot([0 210],[zcritical zcritical],'k:');
plot([0 210],[-zcritical -zcritical],'k:');
xlim([0 206]); 
ylim([-4 4]);
ylabel('HC-PM-EM z-value');
xlabel('Brain segmentations');
legend([H1.mainLine,H2.mainLine],{'HC-PM','HC-EM'},'box','off','location','southwest');
export_fig('HC-PM-EM-z', '-png','-transparent');

%% change in age vs. change in PM

load('age_median_R.mat','medianROI_R');

median_R = mean(medianROI_R,2);
mad_R = std(medianROI_R,0,2);


figure;
scatter(median_R, median_T, 'o');

%% functions

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

function [] = plotROI_pair(data, typeList, pairList, ymin, ymax, QSMname, ROIname, ii)

indPMmatch = find(strcmp(typeList,'PM') & ~isnan(pairList));
indPMmatchHC = indPMmatch+1;
indEMmatch = find(strcmp(typeList,'EM') & ~isnan(pairList));
indEMmatchHC = indEMmatch+1;

plot([1 2], [data(indPMmatchHC) data(indPMmatch)], '-', 'Color', [0 0.4470 0.7410]); hold on;
scatter([1 2], [data(indPMmatchHC) data(indPMmatch)], 'o', 'filled', ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 0.4470 0.7410], ...
    'MarkerFaceAlpha', 0.5);
plot([3 4], [data(indEMmatchHC) data(indEMmatch)], '-', 'Color', [0.8500 0.3250 0.0980]);
scatter([3 4], [data(indEMmatchHC) data(indEMmatch)], 'o', 'filled', ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
    'MarkerFaceAlpha', 0.5);

[p1,~,stats] = signrank(data(indPMmatchHC), data(indPMmatch));
[p2,~,stats] = signrank(data(indEMmatchHC), data(indEMmatch));

title([strrep(strrep(QSMname,'_','-'),'-meanEcho','')]);

ylim([ymin ymax]); % [-0.05 0.15] [0 0.2]
if mod(ii,10) == 1
    ylabel([ROIname ' Susc. (ppm)']);
else
    ylabel('')
end
xlim([0.5 4.5]);
ax = gca;
ax.XTick = [1 2 3 4];
ax.XTickLabel = {'HC','PM','HC','EM'};
sigstar({[1,2],[3,4]},[p1 p2]);

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

function [] = plotHeatMap_HD(dataT, dataP, cmap1, cmap2, QSMname)

ax(1) = subplot(411);
plotHeatMap(dataT(:,:,2), cmap1, QSMname, [-4 4]);
title('Median QSM HC-PM Z-value');
ax(2) = subplot(412);
plotHeatMap(dataP(:,:,2), flip(cmap2), QSMname, [0 0.05]);
title('Median QSM HC-PM P-value');
ax(3) = subplot(413);
plotHeatMap(dataT(:,:,3), cmap1, QSMname, [-4 4]);
title('Median QSM HC-EM Z-value');
ax(4) = subplot(414);
plotHeatMap(dataP(:,:,3), flip(cmap2), QSMname, [0 0.05]);
title('Median QSM HC-EM P-value');

colormap(ax(1), cmap1); colormap(ax(2), flip(cmap2)); 
colormap(ax(3), cmap1); colormap(ax(4), flip(cmap2)); 

end
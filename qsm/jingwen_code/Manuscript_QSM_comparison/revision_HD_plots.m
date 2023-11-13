clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/'));

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

img_path = '/working/lupolab/jingwen/001_QSM/temp';

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','NoRep');
statusList = nan(1,length(T.status_reclass));
statusList(strcmp(T.status_reclass,'HC')) = 1;
statusList(strcmp(T.status_reclass,'PM')) = 2;
statusList(strcmp(T.status_reclass,'EM')) = 3;
% statusList(65) = nan;
% statusList(strcmp(T.status_reclass,'Manifest')) = 3;

subjList = T.b_num; subjList(isnan(statusList)) = [];
examList = T.t_num; examList(isnan(statusList)) = [];
ageList = T.age; ageList(isnan(statusList)) = [];
sexList = T.sex; sexList(isnan(statusList)) = [];
cagList = T.CAG; cagList(isnan(statusList)) = [];

statusList(isnan(statusList)) = [];
[statusList,i] = sort(statusList);

subjList = subjList(i);
examList = examList(i);
ageList = ageList(i);
sexList = sexList(i);
cagList = cagList(i);

indHC = statusList == 1;
indPM = statusList == 2; 
indEM = statusList == 3;

fprintf('PM: Age mean %.3f SD %.3f \n', mean(ageList(indPM)), std(ageList(indPM)));
fprintf('M %i F %i \n', sum(strcmp(sexList(indPM),'M')), sum(strcmp(sexList(indPM),'F')));
fprintf('CAG mean %.3f SD %.3f \n', mean(cagList(indPM)), std(cagList(indPM)));
fprintf('CAG range %.3f - %.3f \n', min(cagList(indPM)), max(cagList(indPM)));

fprintf('EM: Age mean %.3f SD %.3f \n', mean(ageList(indEM)), std(ageList(indEM)));
fprintf('M %i F %i \n', sum(strcmp(sexList(indEM),'M')), sum(strcmp(sexList(indEM),'F')));
fprintf('CAG mean %.3f SD %.3f \n', mean(cagList(indEM)), std(cagList(indEM)));
fprintf('CAG range %.3f - %.3f \n', min(cagList(indEM)), max(cagList(indEM)));

matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data';
load([matout_root '/SubjImgData_QSMcomp_erode.mat'],'DataStruct');

%% load stats on ROI of interest

ROIlist = 1:233;
QSM_name = {'iLSQR' 'QSIP' 'SSTGV' 'SSTV' 'STARQSM' 'FANSI' 'HDQSM' 'MEDI' ...
    'QSMGAN' 'QSMnet+' 'xQSM' 'iQSM'};

QSMfile_list = {'QSM_iLSQR' 'QSM_QSIP' 'QSM_SSTGV' 'QSM_SSTV' ...
    'QSM_STARQSM' 'QSM_FANSI_nonlinearTV' 'QSM_HDQSM' 'QSM_MEDI' ...
    'QSM_QSMGAN' 'QSM_QSMnet' 'QSM_xQSM2' 'QSM_iQSM2'};

meanROI = zeros(length(subjList), length(ROIlist), length(QSMfile_list));
stdROI = zeros(length(subjList), length(ROIlist), length(QSMfile_list));
medianROI = zeros(length(subjList), length(ROIlist), length(QSMfile_list));
madROI = zeros(length(subjList), length(ROIlist), length(QSMfile_list));
for ii = 1:length(subjList)
    
    exam_id = [subjList{ii} '_' examList{ii}];
    fprintf('# Reading subj %s \n', exam_id);
    
    indSubj = find(strcmp({DataStruct.subjID}, subjList{ii}));
    
    for QSMnum = 1:length(QSMfile_list)
        meanROI(ii,:,QSMnum) = DataStruct(indSubj).(QSMfile_list{QSMnum}).ROImean;
        stdROI(ii,:,QSMnum) = DataStruct(indSubj).(QSMfile_list{QSMnum}).ROIstd;
        medianROI(ii,:,QSMnum) = DataStruct(indSubj).(QSMfile_list{QSMnum}).ROImedian;
        madROI(ii,:,QSMnum) = DataStruct(indSubj).(QSMfile_list{QSMnum}).ROImad;
    end
    
end

% scale the QSMGAN
meanROI(:,:,9) = meanROI(:,:,9)/0.5684;
stdROI(:,:,9) = stdROI(:,:,9)/0.5684;
medianROI(:,:,9) = medianROI(:,:,9)/0.5684;
madROI(:,:,9) = madROI(:,:,9)/0.5684;

%% plot results

meanStat = medianROI;
REFind = 212;

indList = [177 179 181 187 189 210 208 115]; % 208 115
indname = {'CN' 'PU' 'GP' 'RN' 'DN' 'SN' 'TH' 'CST'}; %  'TH' 'CST'

for rr = 1:3 % :length(indList)
    ROIind = indList(rr); % [177 179 181 187 189 210 208 115]
    
    % calcultae y range
    data_all = mean(squeeze(meanStat(:,ROIind:ROIind+1,:)),2) - squeeze(meanStat(:,REFind,:));
    ymin = min(data_all(:)); ymax = max(data_all(:)); yrange = ymax - ymin;
    ymin = ymin - 0.05*yrange; ymax = ymax + 0.05*yrange;
    
    figure('position', [100 0 1200 900]);
    for ii = 1:12
        subaxis(3,4,ii,'SpacingVert',0.11)
        data = mean(squeeze(meanStat(:,ROIind:ROIind+1,ii)),2) - squeeze(meanStat(:,REFind,ii));
        F = fit(ageList(indHC),data(indHC),'poly1');
        data_ageCorr = data - F.p1*(ageList-40);
        % data_ageCorr = data;
        % data_ageCorr = data./ageList*median(ageList);
        violinplot(data_ageCorr, statusList, 'showMean', true);

        [p,~,stats] = kruskalwallis(data_ageCorr, statusList, 'off');
        if p < 0.0001
            text(0.6,ymin+0.1*yrange,sprintf('p < 0.0001'),'FontSize',12); % -0.03
        elseif p < 0.001
            text(0.6,ymin+0.1*yrange,sprintf('p < 0.001'),'FontSize',12);
        else
            text(0.6,ymin+0.1*yrange,sprintf('p = %.3f',p),'FontSize',12);
        end
        
        c = multcompare(stats,'Display','off','CType','dunn-sidak');
        
        title(QSM_name{ii});
        % ylim([ymin ymax]); % [-0.05 0.15] [0 0.2]
        
        if mod(ii,4) == 1
            ylabel([indname{rr} ' Susc. (ppm)']);
        else
            ylabel('')
        end
        xlim([0.5 3.5]);
        ax = gca;
        ax.XTick = [1 2 3];
        ax.XTickLabel = {'HC','PM','EM'};
        
        ylim([ymin ymax]);
%         if ii == 12
%             ylim([0.01 0.47]);
%         elseif ii == 3 || ii == 4
%             ylim([0.01 0.22]);
%         else
%             ylim([0.01 0.16]);
%         end
        
        sigstar({[1,2],[1,3],[2,3]},c(:,6));
    end
    export_fig([img_path '/HD-' indname{rr} '_AgeCorr'], '-png','-transparent'); % close;
end

%% Comparison p values

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

% GP, RN, SN, PU, DN, CN
% [181 187 210 179 177 189]
% [232 224 226 216 214 230]
% BGind = [181 182 187 188 210 211 179 180 189 190 177 178];
% GP, PU, CN
BGind = [181 182 179 180 177 178];
% BGind = [232 233 216 217 214 215];
QSM_L = meanStat(:,BGind(1:2:end),:);
QSM_L = cat(2,QSM_L,mean(QSM_L,2));
QSM_R = meanStat(:,BGind(2:2:end),:);
QSM_R = cat(2,QSM_R,mean(QSM_R,2));
QSM_LR = (QSM_L + QSM_R)/2 - repmat(meanStat(:,REFind,:),[1 length(BGind)/2+1 1]);

% corrected alpha
alpha = 1-(1-0.05)^(1/3);
z = @(p) -sqrt(2) * erfcinv(p*2);
zcritical = z(1-alpha);

medianROI_P_BG = zeros(length(BGind)/2+1, length(QSMfile_list),4);
medianROI_T_BG = zeros(length(BGind)/2+1, length(QSMfile_list),4);
medianROI_ES_BG = zeros(length(BGind)/2+1, length(QSMfile_list),4);
medianROI_SS_BG = zeros(length(BGind)/2+1, length(QSMfile_list),4);
for rr = 1:length(BGind)/2+1
    for ii = 1:12
        data = QSM_LR(:,rr,ii);
        F = fit(ageList(indHC),data(indHC),'poly1');
        data_ageCorr = data - F.p1*(ageList-40);
        [p,tbl,stats] = kruskalwallis(data_ageCorr, statusList, 'off');
        c = multcompare(stats,'Display','off','CType','dunn-sidak');
        medianROI_P_BG(rr,ii,1) = p;
        medianROI_T_BG(rr,ii,1) = tbl{2,5};
        medianROI_P_BG(rr,ii,2:4) = c(:,6);
        [~,~,stats] = ranksum(data_ageCorr(indPM),data_ageCorr(indHC));
        medianROI_T_BG(rr,ii,2) = stats.zval;
%         medianROI_ES_BG(rr,ii,2) = stats.zval/sqrt(length(indPM)+length(indHC));
        statsES = mes(data_ageCorr(indPM),data_ageCorr(indHC),{'hedgesg','cles'});
        medianROI_ES_BG(rr,ii,2) = statsES.hedgesg;
%         medianROI_ES_BG(rr,ii,2) = computeCohen_d(data_ageCorr(indPM),data_ageCorr(indHC));
        medianROI_SS_BG(rr,ii,2) = ...
            sampsizepwr('t2',[mean(data_ageCorr(indHC)) std(data_ageCorr(indHC))],...
            median(data_ageCorr(indPM)),0.9,[],'Ratio',2);
        
        [~,~,stats] = ranksum(data_ageCorr(indEM),data_ageCorr(indHC));
        medianROI_T_BG(rr,ii,3) = stats.zval;
%         medianROI_ES_BG(rr,ii,3) = stats.zval/sqrt(length(indEM)+length(indHC));
        statsES = mes(data_ageCorr(indEM),data_ageCorr(indHC),{'hedgesg','cles'});
        medianROI_ES_BG(rr,ii,3) = statsES.hedgesg;
%         medianROI_ES_BG(rr,ii,3) = computeCohen_d(data_ageCorr(indEM),data_ageCorr(indHC));
        medianROI_SS_BG(rr,ii,3) = ...
            sampsizepwr('t2',[mean(data_ageCorr(indHC)) std(data_ageCorr(indHC))],...
            median(data_ageCorr(indEM)),0.9,[],'Ratio',2);
    end
end

% effect size

figure('position', [100 100 800 275]);
data = medianROI_ES_BG(:,:,2);
plot_scatter_Z(data, QSM_name);
hold on;
plot([0 length(QSM_name)+1], 0.5*[1 1], 'k:');
plot([0 length(QSM_name)+1], 0.8*[1 1], 'k:');
% legend({'GP','PU','CN','Mean of ROIs'},'box','off','location','best');
ylabel('Effect Size (HC-PM)'); % ylim([0.6 0.9]);
export_fig([img_path '/HD-BG_ES_compPM'], '-png','-transparent'); % close;

figure('position', [100 100 800 275]);
data = medianROI_ES_BG(:,:,3);
plot_scatter_Z(data, QSM_name);
hold on; 
plot([0 length(QSM_name)+1], 0.5*[1 1], 'k:');
plot([0 length(QSM_name)+1], 0.8*[1 1], 'k:');
% legend({'GP','PU','CN','Mean of ROIs'},'box','off','location','best');
ylabel('Effect Size (HC-EM)'); % ylim([0.6 0.9]);
export_fig([img_path '/HD-BG_ES_compEM'], '-png','-transparent'); % close;

%% save

save_stat = medianROI_ES_BG;
ES_GP_PM = save_stat(1,:,2)';
ES_PU_PM = save_stat(2,:,2)';
ES_CN_PM = save_stat(3,:,2)';
ES_all_PM = save_stat(4,:,2)';
ES_GP_EM = save_stat(1,:,3)';
ES_PU_EM = save_stat(2,:,3)';
ES_CN_EM = save_stat(3,:,3)';
ES_all_EM = save_stat(4,:,3)';

T_HD = table(QSM_name', ES_GP_PM, ES_PU_PM, ES_CN_PM, ES_all_PM, ES_GP_EM, ES_PU_EM, ES_CN_EM, ES_all_EM);

% save([matout_root '/QSMcomp_all.mat'], 'T_HD', '-append');

%% Plot

% figure('position', [100 100 750 300]);
% data = medianROI_P_BG(:,:,1);
% plot_violin_P(data, QSM_name);
% plot([0 13], log10([0.05 0.05]), 'k--');
% ylabel('log(P) in BG (ANOVA)'); % ylim([-2.5 0]);
% export_fig([img_path '/HD-BG_ANOVAPc'], '-png','-transparent'); % close;

% figure('position', [100 100 750 300]);
% data = medianROI_P_BG(:,:,2);
% plot_violin_P(data, QSM_name);
% plot([0 13], log10([alpha alpha]), 'k--');
% ylabel('log(P) in BG (HC-PM)'); % ylim([-2.5 0]);
% export_fig([img_path '/HD-BG_ANOVAPc_compPM'], '-png','-transparent'); % close;

% figure('position', [100 100 750 300]);
% data = medianROI_P_BG(:,:,3);
% plot_violin_P(data, QSM_name);
% plot([0 13], log10([alpha alpha]), 'k--');
% ylabel('log(P) in BG (HC-EM)'); % ylim([-5 1]);
% export_fig([img_path '/HD-BG_ANOVAPc_compEM'], '-png','-transparent'); % close;

% figure('position', [100 100 750 300]);
% data = medianROI_T_BG(:,:,1);
% plot_scatter_Z(data, QSM_name);
% legend({'GP','PU','CN','Mean of ROIs'},'box','off','location','best');
% ylabel('chi square (ANOVA)'); % ylim([-2.5 0]);
% export_fig([img_path '/HD-BG_ANOVAZc'], '-png','-transparent'); % close;

figure('position', [100 100 800 275]);
data = medianROI_T_BG(:,:,2);
plot_scatter_Z(data, QSM_name);
hold on; plot([0 13], zcritical*[1 1], 'k:');
legend({'GP','PU','CN'},'box','off','location','best','FontSize',10);
ylabel('z-value (HC-PM)'); ylim([1 4]);
export_fig([img_path '/HD-BG_ANOVAZc_compPM'], '-png','-transparent'); % close;

figure('position', [100 100 800 275]);
data = medianROI_T_BG(:,:,3);
plot_scatter_Z(data, QSM_name);
hold on; plot([0 13], zcritical*[1 1], 'k:');
% legend({'GP','PU','CN','Mean of ROIs'},'box','off','location','best');
ylabel('z-value (HC-EM)'); ylim([1.5 5.5]);
export_fig([img_path '/HD-BG_ANOVAZc_compEM'], '-png','-transparent'); % close;

%% sample size

figure('position', [100 100 750 300]);
data = medianROI_SS_BG(:,:,2);
plot_scatter_Z(data, QSM_name);
ylabel('Sample Size (HC-PM)'); set(gca, 'YScale', 'log')
export_fig([img_path '/HD-BG_SS_compPM'], '-png','-transparent'); % close;

figure('position', [100 100 750 300]);
data = medianROI_SS_BG(:,:,3);
plot_scatter_Z(data, QSM_name);
ylabel('Sample Size (HC-EM)'); set(gca, 'YScale', 'log')
export_fig([img_path '/HD-BG_SS_compEM'], '-png','-transparent'); % close;

%% output ANOVA test stats

% medianROI_P = zeros(length(ROIlist), length(QSMfile_list),4);
% for rr = 1:length(ROIlist)
%     for ii = 1:10
%         % median
%         try
%             data = squeeze(medianROI(:,rr,ii)) - squeeze(medianROI(:,212,ii));
%             F = fit(ageList(indHC),data(indHC),'poly1');
%             data_ageCorr = data - F.p1*(ageList-40);
%             [p,~,stats] = kruskalwallis(data_ageCorr, statusList, 'off');
%             c = multcompare(stats,'Display','off','CType','dunn-sidak');
%             medianROI_P(rr,ii,1) = p;
%             medianROI_P(rr,ii,2:4) = c(:,6);
%         catch
%             medianROI_P(rr,ii,:) = 1;
%         end
%     end
% end

%% plot heatmap

r = ones(100, 1);
g = linspace(1, 0, 100)'; 
cmap = flip([r g g]);

alpha = 0.05;

data = [medianROI_P(:,:,1)];
data(isnan(data)) = 1;
plotHeatMap(data, cmap, QSM_name, alpha);
title('HC-PM-EM Comparison P-value');
% export_fig('ANOVA1-Median', '-png','-transparent'); % close;

data = [medianROI_P(:,:,2)];
data(isnan(data)) = 1;
plotHeatMap(data, cmap, QSM_name, alpha);
title('HC-PM Comparison P-value');
% export_fig('ANOVA1-PM-Median', '-png','-transparent'); % close;

data = [medianROI_P(:,:,3)];
data(isnan(data)) = 1;
plotHeatMap(data, cmap, QSM_name, alpha);
title('HC-EM Comparison P-value');
% export_fig('ANOVA1-EM-Median', '-png','-transparent'); % close;

%% functions

function [] = plotHeatMap(data, cmap, QSMname, alpha)

figure('position', [100 100 500 700]);
imagesc(data(1:211,:), [0 alpha]); colormap(cmap); colorbar; hold on;
plot([0 lenghth(QSMname)+1],[109.5 109.5],'k');
plot([0 lenghth(QSMname)+1],[177.5 177.5],'k');
ax = gca;
ax.XTick = 1:lenghth(QSMname);
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30;
ax.YTick = floor([(1+109)/2 (109+177)/2 (177+212)/2]);
ax.YTickLabelRotation = 90;
ax.YTickLabel = {'Cortical GM','WM','BG'};

end

function [] = plot_violin_P(data, QSMname)

v = violinplot(log10(data));
v(1).ViolinColor = [0.4940    0.1840    0.5560];
for ii = 2:4
    v(ii).ViolinColor = [0    0.4470    0.7410];
end
for ii = 5:8
    v(ii).ViolinColor = [0.8500    0.3250    0.0980];
end
for ii = 9:12
    v(ii).ViolinColor = [0.4660    0.6740    0.1880];
end
ax = gca;
ax.XTick = 1:length(QSMname);
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30; xlim([0 length(QSMname)+1]);
hold on;

end

function [] = plot_scatter_Z(data, QSMname)

plot(1:length(QSMname), data(1,:), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2); hold on;
plot(1:length(QSMname), data(2,:), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.8500 0.3250 0.0980]/2);
plot(1:length(QSMname), data(3,:), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.9290 0.6940 0.1250]/2);
plot(1:length(QSMname), data(4,:), 'o-', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.4940 0.1840 0.5560]/2);
% hold on; plot([0 11], 1.96*[1 1], 'k--');
ax = gca;
ax.XTick = 1:length(QSMname);
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30; xlim([0 length(QSMname)+1]);
ax.Box = 'off';

end
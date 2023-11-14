clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath('/home/jyao3/010_MATLAB_Utils/');

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList.xlsx');
statusList = [T.PairID];

subjList = T.b_num; subjList(isnan(statusList)) = [];
examList = T.t_num; examList(isnan(statusList)) = [];
typeList = T.status; typeList(isnan(statusList)) = [];

statusList(isnan(statusList)) = [];
[statusList,i] = sort(statusList);

subjList = subjList(i);
examList = examList(i);
typeList = typeList(i);

matout_root = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/Data';

%% load stats on ROI of interest

ROIlist = 1:213;
QSMfile_list = {'QSM_FANSI' ...
    'QSM_HDQSM' ...
    'QSM_iLSQR' ...
    'QSM_MEDI' ...
    'QSM_QSIP' ...
    'QSM_QSMGAN' ...
    'QSM_QSMnet' ...
    'QSM_SSTGV' ...
    'QSM_SSTV' ...
    'QSM_STARQSM'};

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

%% plot results

% [177 179 181 187 189 210 208 115]
% CN PU GP RN DN SN TH CST

indPM = find(strcmp(typeList,'PM'));
indPMHC = indPM + 1;
indEM = find(strcmp(typeList,'EM'));
indEMHC = indEM + 1;

for ROIind = 177 % [177 179 181 187 189 210 208 115]
    figure('position', [100 100 1500 500]);
    for ii = 1:10
        subplot(2,5,ii)
        data = mean(squeeze(medianROI(:,ROIind:ROIind+1,ii)),2) - squeeze(medianROI(:,212,ii));
        for nn = 2:2:30
            if strcmp(typeList(nn-1),'PM')
                a = plot([1 2],data([nn nn-1]),'o-', 'MarkerSize', 7, 'Color', [0 0.4470 0.7410],...
                    'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2); hold on;
            else
                plot([3 4],data([nn nn-1]),'o-', 'MarkerSize', 7, 'Color', [0.8500 0.3250 0.0980],...
                    'MarkerFaceColor', [0.5 0.5 0.5] + [0.8500 0.3250 0.0980]/2); hold on;
            end
        end
        pPM = signrank(data(indPM), data(indPMHC));
        pEM = signrank(data(indEM), data(indEMHC));
        title([strrep(strrep(QSMfile_list(ii),'_','-'),'-meanEcho','')]);
        ylim([-0.05 0.2]);
        ylabel('Susceptibility (ppm)');
        xlim([0.5 4.5]);
        ax = gca;
        ax.XTick = [1 2 3 4];
        ax.XTickLabel = {'HC','PM','HC','EM'};
        sigstar({[1,2],[3,4]},[pPM pEM]);
    end
    export_fig([ROIname{ROIind} '-compHD'], '-png','-transparent'); % close;
end

%% output paired test stats

medianROI_P = zeros(length(ROIlist), length(QSMfile_list),2);
stdROI_P = zeros(length(ROIlist), length(QSMfile_list),2);
for rr = 1:length(ROIlist)
    for ii = 1:10
        % median
        data = squeeze(medianROI(:,rr,ii)) - squeeze(medianROI(:,212,ii));
        try
            P = signrank(data(indPM), data(indPMHC));
            medianROI_P(rr,ii,1) = P;
        catch
            medianROI_P(rr,ii,1) = 1;
        end
        try
            P = signrank(data(indEM), data(indEMHC));
            medianROI_P(rr,ii,2) = P;
        catch
            medianROI_P(rr,ii,2) = 1;
        end
        
        % std
        data = squeeze(stdROI(:,rr,ii));
        try
            P = signrank(data(indPM), data(indPMHC));
            stdROI_P(rr,ii,1) = P;
        catch
            stdROI_P(rr,ii,1) = 1;
        end
        try
            P = signrank(data(indEM), data(indEMHC));
            stdROI_P(rr,ii,2) = P;
        catch
            stdROI_P(rr,ii,2) = 1;
        end
    end
end

%% correction for multiple comparison

for ii = 1:10
    test = medianROI_P(1:211,ii,1);
    [c_pvalues, ~, ~] = fdr_BH(test, 0.05, false); % , 'corr+'
    medianROI_P(1:211,ii,1) = c_pvalues;
    test = medianROI_P(1:211,ii,2);
    [c_pvalues, ~, ~] = fdr_BH(test, 0.05, false);
    medianROI_P(1:211,ii,2) = c_pvalues;
    
    test = stdROI_P(1:211,ii,1);
    [c_pvalues, ~, ~] = fdr_BH(test, 0.05, false);
    stdROI_P(1:211,ii,1) = c_pvalues;
    test = stdROI_P(1:211,ii,2);
    [c_pvalues, ~, ~] = fdr_BH(test, 0.05, false);
    stdROI_P(1:211,ii,2) = c_pvalues;
end

%% plot heatmap

QSMname = strrep(QSMfile_list,'_','-');
QSMname = strrep(QSMname,'QSM-','');
QSMname = strrep(QSMname,'-meanEcho','');

r = ones(100, 1);
g = linspace(1, 0, 100)'; 
cmap = flip([r g g]);

alpha = 0.05;

data = [medianROI_P(:,:,1)];
data(isnan(data)) = 1;
plotHeatMap(data, cmap, QSMname, alpha);
title('HC-PM Comparison P-value');
export_fig('PMcomp-Median', '-png','-transparent'); % close;

data = [medianROI_P(:,:,2)];
data(isnan(data)) = 1;
plotHeatMap(data, cmap, QSMname, alpha);
title('HC-EM Comparison P-value');
export_fig('EMcomp-Median', '-png','-transparent'); % close;

data = [stdROI_P(:,:,1)];
data(isnan(data)) = 1;
plotHeatMap(data, cmap, QSMname, alpha);
title('HC-PM Comparison P-value');
export_fig('PMcomp-STD', '-png','-transparent'); % close;

data = [stdROI_P(:,:,2)];
data(isnan(data)) = 1;
plotHeatMap(data, cmap, QSMname, alpha);
title('HC-EM Comparison P-value');
export_fig('EMcomp-STD', '-png','-transparent'); % close;

%% Comparison p values

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

% GP, RN, SN, PU, DN, CN
% BGind = [181 182 187 188 210 211 179 180 189 190 177 178];
% GP, PU, CN
BGind = [181 182 179 180 177 178];
QSM_L = medianROI(:,BGind(1:2:end),:);
QSM_L = cat(2,QSM_L,mean(QSM_L,2));
QSM_R = medianROI(:,BGind(2:2:end),:);
QSM_R = cat(2,QSM_R,mean(QSM_R,2));
QSM_LR = (QSM_L + QSM_R)/2 - repmat(medianROI(:,212,:),[1 length(BGind)/2+1 1]);

medianROI_P_BG = zeros(length(BGind)/2+1, length(QSMfile_list),2);
medianROI_Z_BG = zeros(length(BGind)/2+1, length(QSMfile_list),2);
for rr = 1:length(BGind)/2+1
    for ii = 1:10
        data = QSM_LR(:,rr,ii); 
        [P,~,stats] = signrank(data(indPM), data(indPMHC),'method','approximate');
        medianROI_P_BG(rr,ii,1) = P;
        medianROI_Z_BG(rr,ii,1) = stats.zval;
        [P,~,stats] = signrank(data(indEM), data(indEMHC),'method','approximate');
        medianROI_P_BG(rr,ii,2) = P;
        medianROI_Z_BG(rr,ii,2) = stats.zval;
    end
end

figure('position', [100 100 900 300]);
data = medianROI_P_BG(:,:,1);
plot_violin_P(data, QSMname);
ylabel('log(P) in BG (HC-PM)'); ylim([-2.5 0]);
export_fig('BG_P_compPM', '-png','-transparent'); % close;

figure('position', [100 100 900 300]);
data = medianROI_P_BG(:,:,2);
plot_violin_P(data, QSMname);
ylabel('log(P) in BG (HC-EM)'); ylim([-5 1]);
export_fig('BG_P_compEM', '-png','-transparent'); % close;

figure('position', [100 100 900 300]);
data = medianROI_Z_BG(:,:,1);
plot_scatter_Z(data, QSMname);
ylabel('z-value (HC-PM)'); % ylim([-2.5 0]);
% legend({'GP','PU','CN','Mean of ROIs'},'box','off','location','best');
export_fig('BG_Z_compPM', '-png','-transparent'); % close;

figure('position', [100 100 900 300]);
data = medianROI_Z_BG(:,:,2);
plot_scatter_Z(data, QSMname);
legend({'GP','PU','CN','Mean of ROIs'},'box','off','location','best');
ylabel('z-value (HC-EM)'); % ylim([-2.5 0]);
export_fig('BG_Z_compEM', '-png','-transparent'); % close;

%% functions

function [] = plotHeatMap(data, cmap, QSMname, alpha)

figure('position', [100 100 500 700]);
imagesc(data(1:211,:), [0 alpha]); colormap(cmap); colorbar; hold on;
plot([0 11],[109.5 109.5],'k');
plot([0 11],[177.5 177.5],'k');
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30;
ax.YTick = floor([(1+109)/2 (109+177)/2 (177+212)/2]);
ax.YTickLabelRotation = 90;
ax.YTickLabel = {'Cortical GM','WM','BG'};

end

function [] = plot_violin_P(data, QSMname)

v = violinplot(log10(data));
v(1).ViolinColor = [0.4940    0.1840    0.5560];
for ii = 2:5
    v(ii).ViolinColor = [0    0.4470    0.7410];
end
for ii = 6:8
    v(ii).ViolinColor = [0.8500    0.3250    0.0980];
end
for ii = 9:10
    v(ii).ViolinColor = [0.4660    0.6740    0.1880];
end
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30; xlim([0 11]);
hold on; plot([0 11], log10([0.05 0.05]), 'k--');

end

function [] = plot_scatter_Z(data, QSMname)

plot(1:10, data(1,:), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2); hold on;
plot(1:10, data(2,:), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.8500 0.3250 0.0980]/2);
plot(1:10, data(3,:), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.9290 0.6940 0.1250]/2);
plot(1:10, data(4,:), 'o-', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.4940 0.1840 0.5560]/2);
hold on; plot([0 11], 1.96*[1 1], 'k--');
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30; xlim([0 11]);
ax.Box = 'off';

end
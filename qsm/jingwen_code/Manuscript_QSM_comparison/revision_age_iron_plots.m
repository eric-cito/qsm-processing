clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath('/home/jyao3/010_MATLAB_Utils/');
addpath('/home/jyao3/010_MATLAB_Utils/violinplot');
addpath('/home/jyao3/010_MATLAB_Utils/export_fig');

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

img_path = '/working/lupolab/jingwen/001_QSM/temp';

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','NoRep');
statusList = [T.status];

subjList = T.b_num; subjList(~strcmp(statusList,'HC')) = [];
examList = T.t_num; examList(~strcmp(statusList,'HC')) = [];
sexList = T.sex; sexList(~strcmp(statusList,'HC')) = [];
ageList = T.age; ageList(~strcmp(statusList,'HC')) = [];

Mind = strcmp(sexList,'M');
Find = strcmp(sexList,'F');

fprintf('Age mean %.3f SD %.3f \n', mean(ageList), std(ageList));
fprintf('M %i F %i \n', sum(Mind), sum(Find));

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

ROIname = DataStruct(1).(QSMfile_list{1}).ROIname;

%% plot ROI example

meanStat = medianROI;
REFind = 212;

set(0,'defaultAxesFontSize',14);

% GP, RN, SN, PU, DN, CN
% [181 187 210 179 189 177];

tmpList = [177 179 181 187 210 189];
tmpName = {'CN','PU','GP','RN','SN','DN'};

tt = 5;

ROIind = tmpList(tt);
ROIname = tmpName{tt};

figure('position', [100 0 1200 900]);
% data range
data = mean(squeeze(meanStat(:,ROIind:ROIind+1,:)),2) - squeeze(meanStat(:,REFind,:));
ymin = min(data(:)); ymax = max(data(:)); yrange = ymax - ymin;
ymin = ymin - 0.01*yrange; ymax = ymax + 0.3*yrange;
for ii = 1:12
    
    subaxis(3,4,ii,'SpacingVert',0.11)
    QSMcomp_age_subplot(3, 4, ii, meanStat, ageList, sexList, ROIname, ROIind, REFind, [ymin ymax]);
    title(QSM_name(ii));
    
end

export_fig([img_path '/Age_' tmpName{tt}], '-png','-transparent'); 

%% output age correlation stats

meanROI_R = zeros(length(ROIlist), length(QSMfile_list));
meanROI_P = zeros(length(ROIlist), length(QSMfile_list));
meanROI_slope = zeros(length(ROIlist), length(QSMfile_list));
for rr = 1:length(ROIlist)
    for ii = 1:length(QSMfile_list)
        data = squeeze(meanStat(:,rr,ii)) - squeeze(meanStat(:,REFind,ii)); 
        [R,P] = corr(ageList(:),data(:),'Type','Pearson');
        meanROI_R(rr,ii) = R;
        meanROI_P(rr,ii) = P;
        try
            F = fit(ageList(~isnan(data)),data(~isnan(data)),'poly1');
            meanROI_slope(rr,ii) = F.p1;
        catch
            meanROI_slope(rr,ii) = nan;
        end
    end
end

% save('QSM_age_fit.mat', 'medianROI_slope', 'medianROI_intercept');

%% plot heatmap - correlation

data = [meanROI_R];
data(isnan(data)) = 0;
cmap = flip(colormap_RWG([-1 1],100));

figure('position', [100 100 750 450]);
imagesc(data(1:209,:)',[-1 1]); colormap(cmap); colorbar; hold on;
plot([109.5 109.5],[0 13],'k','LineWidth',2);
plot([176.5 176.5],[0 13],'k','LineWidth',2);
plot([190.5 190.5],[0 13],'k','LineWidth',2);
ax = gca;
ax.YTick = 1:length(QSM_name);
ax.YTickLabel = QSM_name(:);
ax.YTickLabelRotation = 0;
ax.XTick = floor([(1+109)/2 (109+176)/2 (176+191)/2 (191+210)/2]);
ax.XTickLabelRotation = 0;
ax.XTickLabel = {'Cortical GM','WM','BG','TH'};
title('QSM-Age Correlation Coefficient');
export_fig([img_path '/Age-Corr-Median-R'], '-png','-transparent'); % close;

r = ones(100, 1);
g = linspace(1, 0, 100)'; 
cmap = [r g g];

data = [meanROI_P];
data(isnan(data)) = 1;

figure('position', [100 100 750 300]);
imagesc(data(1:209,:)', [0 0.05]); colormap(flip(cmap)); colorbar; hold on;
plot([109.5 109.5],[0 13],'k','LineWidth',2);
plot([176.5 176.5],[0 13],'k','LineWidth',2);
plot([190.5 190.5],[0 13],'k','LineWidth',2);
ax = gca;
ax.YTick = 1:length(QSM_name);
ax.YTickLabel = QSM_name(:);
ax.YTickLabelRotation = 0;
ax.XTick = floor([(1+109)/2 (109+176)/2 (176+191)/2 (191+210)/2]);
ax.XTickLabelRotation = 0;
ax.XTickLabel = {'Cortical GM','WM','BG','TH'};
title('QSM-Age Correlation P-value');
% export_fig([img_path '/Age-Corr-Median-P'], '-png','-transparent'); % close;

%% BG correlation coefficients

meanStat = medianROI;
REFind = 212;

set(0,'defaultAxesFontSize',14);

%  GP, RN, SN, PU, CN, DN, 
% [181 187 210 179 177 189]
% [232 224 226 216 214 230]
% BGind = [232 224 226 216 214 230];
BGind = [187 210 179 177 189];
QSM_L = meanStat(:,BGind,:);
QSM_L = cat(2,QSM_L,mean(QSM_L,2));
QSM_R = meanStat(:,BGind+1,:);
QSM_R = cat(2,QSM_R,mean(QSM_R,2));
QSM_LR = (QSM_L + QSM_R)/2 - repmat(meanStat(:,REFind,:),[1 length(BGind)+1 1]);

medianROI_Rtemp = zeros(length(BGind)+1, length(QSMfile_list));
medianROI_Ptemp = zeros(length(BGind)+1, length(QSMfile_list));
medianROI_Stemp = zeros(length(BGind)+1, length(QSMfile_list));
for rr = 1:length(BGind)+1
    for ii = 1:length(QSMfile_list)
        data = QSM_LR(:,rr,ii); 
        [R,P] = corr(ageList(:),data(:),'Type','Pearson');
        medianROI_Rtemp(rr,ii) = R;
        medianROI_Ptemp(rr,ii) = P;
        F = fit(ageList(:),data(:),'poly1');
        medianROI_Stemp(rr,ii) = F.p1;
    end
end

figure('position', [100 100 750 300]);
v = violinplot(medianROI_Rtemp(1:end-1,:));
v(1).ViolinColor = [0.4940    0.1840    0.5560];
for ii = 2:4
    v(ii).ViolinColor = [0.8500    0.3250    0.0980];
end
for ii = 5:8
    v(ii).ViolinColor = [0    0.4470    0.7410];
end
for ii = 9:12
    v(ii).ViolinColor = [0.4660    0.6740    0.1880];
end
hold on; 
m = plot(1:12, medianROI_Rtemp(end,:), 'o-', 'Color', [0 0 0], ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5]/2);
ax = gca;
ax.XTick = 1:length(QSM_name);
ax.XTickLabel = QSM_name(:);
ax.XTickLabelRotation = 30;
ylabel('Corr. Coef. in BG'); % ylim([0.2 0.7]);
xlim([0 13]); ylim([0.2 0.8]);
legend(m, {'Mean of BG ROIs'}, 'box','off','location','best');
export_fig([img_path '/Age-BG_R'], '-png','-transparent'); % close;

BG_age_R = medianROI_Rtemp(end,:)';

%% with Iron

set(0,'defaultAxesFontSize',14);

% [GP RN SN PU DN CN TH WM(ACR)]
Iron_mean = [21.30 19.48 18.46 13.32 10.35 9.28 4.76]; % 4.24
Iron_SD = [3.49 6.86 6.52 3.43 4.86 2.14 1.16]; % 4.24

ROIind = [181 187 210 179 189 177 208]; % 4.24
% ROIind = [232 224 226 216 230 214 222 131];

% Iron_mean = [205 nan nan 160 nan 105 50 47 29 46 36];
% Iron_SD = [39 nan nan 37 nan 27 12 11 10 11 9];
% ROIind = [181 187 210 179 189 177 208 131 111 133 135];

QSM_L = meanStat(:,ROIind,:); % 131
QSM_R = meanStat(:,ROIind+1,:); % 131
QSM_LR = (QSM_L + QSM_R)/2 - repmat(meanStat(:,REFind,:),[1 length(Iron_mean) 1]);
QSM_mean = squeeze(nanmean(QSM_LR,1));
QSM_SD = squeeze(nanstd(QSM_LR,1));

figure('position', [100 0 1200 900]);
R_wIron = zeros(length(QSM_mean),1);
for ii = 1:12
    subaxis(3,4,ii,'SpacingVert',0.11)
%     errorbarxy(Iron_mean, QSM_mean(:,ii), Iron_SD, QSM_SD(:,ii),...
%         {'ko',[0 0.4470 0.7410],[0 0.4470 0.7410]}); hold on;
    %     plot(Iron_mean, QSM_mean(:,ii),'o','Color',[0 0.4470 0.7410],...
    %         'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2);
    for xx = 1:length(Iron_mean)
        errorbarxy(Iron_mean(xx), QSM_mean(xx,ii), Iron_SD(xx), QSM_SD(xx,ii),{'ko',[0 0.4470 0.7410],[0 0.4470 0.7410]});
        hold on;
        if xx == 8
            s(xx) = scatter(Iron_mean(xx), QSM_mean(xx,ii), 80, 'ko', 'filled'); hold on;
        else
            s(xx) = scatter(Iron_mean(xx), QSM_mean(xx,ii), 80, 'o', 'filled'); hold on;
        end
    end
    hold on;
    validInd = ~isnan(Iron_mean);
    [F,gof] = fit(Iron_mean(validInd)', QSM_mean(validInd,ii),'poly1');
    h = plot(F);
    h.Color = [0.5 0.5 0.5];
    [R,P] = corr(Iron_mean(validInd)', QSM_mean(validInd,ii),'Type','Pearson');
    R_wIron(ii) = R;
    
    text_posy = 0.11;
    
    if P < 0.0001
        text(2,text_posy,sprintf('r = %.3f\np < 0.0001', R));
    elseif P < 0.001
        text(2,text_posy,sprintf('r = %.3f\np < 0.001', R));
    else
        text(2,text_posy,sprintf('r = %.3f\np = %.3f', R, P));
    end
    
    title(strrep(strrep(QSM_name(ii),'_','-'),'-meanEcho',''));
    if ii > 8
        xlabel('mg Iron /100g');
    else
        xlabel('');
    end
    if mod(ii,4) == 1
        ylabel('Susceptibility (ppm)');
    else
        ylabel('');
    end
    ylim([-0.04 0.13]);
    legend('off');
end
export_fig([img_path '/Iron'], '-png','-transparent'); % close;

%% plot Iron example

set(0,'defaultAxesFontSize',14);

ii = 1;
figure('position', [100 100 500 400]);
for xx = 1:length(Iron_mean)
    errorbarxy(Iron_mean(xx), QSM_mean(xx,ii), Iron_SD(xx), QSM_SD(xx,ii),{'ko',[0 0.4470 0.7410],[0 0.4470 0.7410]});
    hold on;
    if xx == 8
        s(xx) = scatter(Iron_mean(xx), QSM_mean(xx,ii), 80, 'ko', 'filled'); hold on;
    else
        s(xx) = scatter(Iron_mean(xx), QSM_mean(xx,ii), 80, 'o', 'filled'); hold on;
    end
end
% GP RN SN PU DN CN TH WM(ACR)
legend(s,{'GP','RN','SN','PU','DN','CN','TH'},'box','off','location','best');
title(QSM_name{ii});
xlabel('mg Iron /100g Weight');
ylabel('Susceptibility (ppm)'); ylim([-0.05 0.1]);
export_fig([img_path '/Iron-example'], '-png','-transparent'); % close;

%% with Iron-Age equation

set(0,'defaultAxesFontSize',14);

% [WM(ACR) GP CD PU]
QSM_L = medianROI(:,[131 181 177 179],:);
QSM_R = medianROI(:,[131 181 177 179]+1,:);
QSM_LR = (QSM_L + QSM_R)/2 - repmat(meanStat(:,REFind,:),[1 4 1]);
% QSM_LR = reshape(QSM_LR,[],10);

Iron_WM = 3.95*(1-exp(-0.10*ageList))+0.31;
Iron_GP = 21.41*(1-exp(-0.09*ageList))+0.37;
Iron_CD = 9.66*(1-exp(-0.05*ageList))+0.33;
Iron_PU = 14.62*(1-exp(-0.04*ageList))+0.46;

Iron_age = [Iron_WM Iron_GP Iron_CD Iron_PU];

clear s;
figure('position', [100 0 1200 900]);
R_wIronAge = zeros(length(QSM_mean),1);
for ii = 1:12
    subaxis(3,4,ii,'SpacingVert',0.11)
    
    for nn = 1:4
        s(nn) = plot(Iron_age(:,nn), QSM_LR(:,nn,ii), '+'); hold on;
    end
    xlim([0 30]); 
    QSM_temp = QSM_LR(:,:,ii);
    [F,~] = fit(Iron_age(:), QSM_temp(:),'poly1');
    f = plot(F);
    f.Color = [0.5 0.5 0.5];
    [R,P] = corr(Iron_age(:), QSM_temp(:),'Type','Pearson');
    R_wIronAge(ii) = R;
    if P < 0.0001
        text(2,0.12,sprintf('r = %.3f\np < 0.0001', R));
    elseif P < 0.001
        text(2,0.12,sprintf('r = %.3f\np < 0.001', R));
    else
        text(2,0.12,sprintf('r = %.3f\np = %.3f', R, P));
    end
    if ii == 12
        legend(s,{'WM','GP','CN','PU'},...
                'box','off','location','southeast');
    else
        legend('off');
    end
    title(QSM_name{ii});
    if ii > 8
        xlabel('mg Iron /100g');
    else
        xlabel('')
    end
    if mod(ii,4) == 1
        ylabel('Susceptibility (ppm)');
    else
        ylabel('')
    end
    ylim([-0.05 0.15]);
end
% legend('WM','GP','CD','PU');

export_fig([img_path '/Iron_age'], '-png','-transparent'); % close;

%% plot Iron example

% set(0,'defaultAxesFontSize',12);

ii = 1;
figure('position', [100 100 300 400]);
clear s
for nn = 1:4
    s(nn) = plot(Iron_age(:,nn), QSM_LR(:,nn,ii), '+'); hold on;
end
QSM_temp = QSM_LR(:,:,ii);
[F,~] = fit(Iron_age(:), QSM_temp(:),'poly1');
xlim([0 30]); ylim([-0.05 0.15]);
f = plot(F); f.Color = [0.5 0.5 0.5];
[R,P] = corr(Iron_age(:), QSM_temp(:),'Type','Pearson');
if P < 0.0001
    text(2,0.12,sprintf('r = %.3f\np < 0.0001', R));
elseif P < 0.001
    text(2,0.12,sprintf('r = %.3f\np < 0.001', R));
else
    text(2,0.12,sprintf('r = %.3f\np = %.3f', R, P));
end
legend(s,{'WM','GP','CN','PU'},...
                'box','off','location','southeast');
title(QSM_name{ii});
xlabel('mg Iron /100g Weight');
ylabel('Susceptibility (ppm)');
export_fig([img_path '/wIronAge-' QSM_name{ii}], ...
    '-png','-transparent'); % close;

%% iron correlation

set(0,'defaultAxesFontSize',14);

figure('position', [100 100 750 300]);
s1 = plot(1:12, R_wIron, 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2); hold on;
% s2 = plot(1:12, R_wIronAge, 'o', ...
%     'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.8500 0.3250 0.0980]/2); hold on;
ax = gca;
ax.XTick = 1:12;
ax.XTickLabel = QSM_name(:);
ax.XTickLabelRotation = 30;
ax.Box = 'off';
ylabel('Corr. Coef. with Iron'); ylim([0.9 1]);
xlim([0 13]);
% legend({'Mean Iron','Iron(age)'},'box','off','location','best');
export_fig([img_path '/Iron_R'], '-png','-transparent'); % close;

%% save with and without white matter

figure('position', [100 100 300 300]);
scatter(R1, R2, 80, 'o', 'filled'); hold on;
xlabel('R with WM ROI'); ylabel('R without WM ROI');

[Rmethod,Pmethod] = corr(R1, R2,'Type','Pearson');
text(0.96, 0.925, sprintf('r = %.3f\np<0.0001', Rmethod));
legend({'QSM methods'}, 'box', 'off', 'location', 'best');

export_fig([img_path '/R_iron_WM'], '-png','-transparent'); % close;

%% save

T_Corr = table(QSM_name', R_wIron, R_wIronAge, BG_age_R);

save([matout_root '/QSMcomp_all.mat'], 'T_Corr', '-append');

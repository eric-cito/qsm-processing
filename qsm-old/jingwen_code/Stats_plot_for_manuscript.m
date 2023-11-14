clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath('/home/jyao3/010_MATLAB_Utils/');
addpath('/home/jyao3/010_MATLAB_Utils/violinplot');
addpath('/home/jyao3/010_MATLAB_Utils/export_fig');

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220103.xlsx','Sheet','NoRep');
statusList = [T.status];

subjList = T.b_num; subjList(~strcmp(statusList,'HC')) = [];
examList = T.t_num; examList(~strcmp(statusList,'HC')) = [];
sexList = T.sex; sexList(~strcmp(statusList,'HC')) = [];
ageList = T.age; ageList(~strcmp(statusList,'HC')) = [];

matout_root = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/Data';

Mind = strcmp(sexList,'M');
Find = strcmp(sexList,'F');

fprintf('Age mean %.3f SD %.3f \n', mean(ageList), std(ageList));
fprintf('M %i F %i \n', sum(Mind), sum(Find));

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

%% plot ROI example

set(0,'defaultAxesFontSize',14);

% GP, RN, SN, PU, DN, CN
% [181 187 210 179 189 177];

ROIind = 181;
figure('position', [100 100 1500 500]);
% data range
data = squeeze(mean(squeeze(medianROI(:,ROIind:ROIind+1,:)),2)) - squeeze(medianROI(:,212,:));
ymin = min(data(:)); ymax = max(data(:)); yrange = ymax - ymin;
ymin = ymin - 0.01*yrange; ymax = ymax + 0.1*yrange;
for ii = 1:10

    subplot(2,5,ii)
    data = mean(squeeze(medianROI(:,ROIind:ROIind+1,ii)),2) - squeeze(medianROI(:,212,ii)); %
    % err = mean(squeeze(stdROI(:,ROIind:ROIind+1,ii)),2);
    plot(ageList(Mind), data(Mind), 'o', ...
        'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2); hold on;
    plot(ageList(Find), data(Find), 'o', ...
        'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.8500 0.3250 0.0980]/2); hold on;
    % errorbar(ageList(Mind), data(Mind), err(Mind), 'o'); hold on;
    % errorbar(ageList(Find), data(Find), err(Find), 'o');
    ylim([ymin ymax]); % [-0.05 0.2]
    xlim([20 80]);
    
    [F,~] = fit(ageList(:),data(:),'poly1');
    h = plot(F);
    h.Color = [0.5 0.5 0.5];
    [R,P] = corrcoef(ageList(:),data(:));
    if ii == 1
        if P(1,2) < 0.0001
            legend({'Male', 'Female', sprintf('r = %.3f\np < 0.0001', R(1,2))},...
                'box','off','location','best');
        elseif P(1,2) < 0.001
            legend({'Male', 'Female', sprintf('r = %.3f\np < 0.001', R(1,2))},...
                'box','off','location','best');
        else
            legend({'Male', 'Female', sprintf('r = %.3f\np = %.3f', R(1,2), P(1,2))},...
                'box','off','location','best');
        end
    else
        if P(1,2) < 0.0001
            legend(h,{sprintf('r = %.3f\np < 0.0001', R(1,2))},...
                'box','off','location','best');
        elseif P(1,2) < 0.001
            legend(h,{sprintf('r = %.3f\np < 0.001', R(1,2))},...
                'box','off','location','best');
        else
            legend(h,{sprintf('r = %.3f\np = %.3f', R(1,2), P(1,2))},...
                'box','off','location','best');
        end
    end
    
    if mod(ii,5) == 1
        ylabel('CN Susc. (ppm)'); % Putamen
    else
        ylabel('');
    end
    if ii > 5
        xlabel('Age (years)');
    else
        xlabel('');
    end
    title(strrep(strrep(QSMfile_list(ii),'_','-'),'-meanEcho',''));
end

export_fig(['Age_' ROIname{ROIind}], '-png','-transparent'); 

%% output age correlation stats

meanROI_R = zeros(length(ROIlist), length(QSMfile_list));
meanROI_P = zeros(length(ROIlist), length(QSMfile_list));
meanROI_slope = zeros(length(ROIlist), length(QSMfile_list));
medianROI_R = zeros(length(ROIlist), length(QSMfile_list));
medianROI_P = zeros(length(ROIlist), length(QSMfile_list));
medianROI_slope = zeros(length(ROIlist), length(QSMfile_list));
medianROI_intercept = zeros(length(ROIlist), length(QSMfile_list));
stdROI_R = zeros(length(ROIlist), length(QSMfile_list));
stdROI_P = zeros(length(ROIlist), length(QSMfile_list));
stdROI_slope = zeros(length(ROIlist), length(QSMfile_list));
for rr = 1:length(ROIlist)
    for ii = 1:10
        % mean
        data = squeeze(meanROI(:,rr,ii)) - squeeze(meanROI(:,212,ii)); 
        [R,P] = corr(ageList(:),data(:),'Type','Pearson');
        meanROI_R(rr,ii) = R;
        meanROI_P(rr,ii) = P;
        try
            F = fit(ageList(~isnan(data)),data(~isnan(data)),'poly1');
            meanROI_slope(rr,ii) = F.p1;
        catch
            meanROI_slope(rr,ii) = nan;
        end
        % median
        data = squeeze(medianROI(:,rr,ii)) - squeeze(medianROI(:,212,ii)); 
        [R,P] = corr(ageList(:),data(:),'Type','Pearson');
        medianROI_R(rr,ii) = R;
        medianROI_P(rr,ii) = P;
        try
            F = fit(ageList(~isnan(data)),data(~isnan(data)),'poly1');
            medianROI_slope(rr,ii) = F.p1;
            medianROI_intercept(rr,ii) = F.p0;
        catch
            medianROI_slope(rr,ii) = nan;
            medianROI_intercept(rr,ii) = nan;
        end
        % std
        data = squeeze(stdROI(:,rr,ii)); 
        [R,P] = corr(ageList(:),data(:),'Type','Pearson');
        stdROI_R(rr,ii) = R;
        stdROI_P(rr,ii) = P;
        try
            F = fit(ageList(~isnan(data)),data(~isnan(data)),'poly1');
            stdROI_slope(rr,ii) = F.p1;
        catch
            stdROI_slope(rr,ii) = nan;
        end
    end
end

% save('QSM_age_fit.mat', 'medianROI_slope', 'medianROI_intercept');

%% plot heatmap - median

QSMname = strrep(QSMfile_list,'_','-');
QSMname = strrep(QSMname,'QSM-','');
QSMname = strrep(QSMname,'-meanEcho','');

data = [medianROI_R];
data(isnan(data)) = 0;
cmap = flip(colormap_RWG([-1 1],100));

figure('position', [100 100 750 300]);
imagesc(data(1:209,:)',[-1 1]); colormap(cmap); colorbar; hold on;
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
export_fig('Age-Corr-Median-R', '-png','-transparent'); % close;

r = ones(100, 1);
g = linspace(1, 0, 100)'; 
cmap = [r g g];

data = [medianROI_P];
data(isnan(data)) = 1;

figure('position', [100 100 750 300]);
imagesc(data(1:209,:)', [0 0.05]); colormap(flip(cmap)); colorbar; hold on;
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
title('QSM-Age Correlation P-value');
export_fig('Age-Corr-Median-P', '-png','-transparent'); % close;

%% BG correlation coefficients

set(0,'defaultAxesFontSize',14);

% GP, RN, SN, PU, CN % DN, 
BGind = [181 182 187 188 210 211 179 180 177 178]; %  189 190
QSM_L = medianROI(:,BGind(1:2:end),:);
QSM_L = cat(2,QSM_L,mean(QSM_L,2));
QSM_R = medianROI(:,BGind(2:2:end),:);
QSM_R = cat(2,QSM_R,mean(QSM_R,2));
QSM_LR = (QSM_L + QSM_R)/2 - repmat(medianROI(:,212,:),[1 length(BGind)/2+1 1]);

medianROI_Rtemp = zeros(length(BGind)/2+1, length(QSMfile_list));
medianROI_Ptemp = zeros(length(BGind)/2+1, length(QSMfile_list));
medianROI_Stemp = zeros(length(BGind)/2+1, length(QSMfile_list));
for rr = 1:length(BGind)/2+1
    for ii = 1:10
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
for ii = 2:5
    v(ii).ViolinColor = [0    0.4470    0.7410];
end
for ii = 6:8
    v(ii).ViolinColor = [0.8500    0.3250    0.0980];
end
for ii = 9:10
    v(ii).ViolinColor = [0.4660    0.6740    0.1880];
end
hold on; 
m = plot(1:10, medianROI_Rtemp(end,:), 'o-', 'Color', [0 0 0], ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5]/2);
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30;
ylabel('Corr. Coef. in BG'); % ylim([0.2 0.7]);
xlim([0 11]); ylim([0.2 0.7]);
legend(m, {'Mean of BG ROIs'}, 'box','off','location','best');
export_fig('Age-BG_R', '-png','-transparent'); % close;

%% with Iron

set(0,'defaultAxesFontSize',14);

% [GP RN SN PU DN CN TH WM(ACR)]
Iron_mean = [21.30 19.48 18.46 13.32 10.35 9.28 4.76 4.24];
Iron_SD = [3.49 6.86 6.52 3.43 4.86 2.14 1.16 0.88];

QSM_L = medianROI(:,[181 187 210 179 189 177 208 131],:);
QSM_R = medianROI(:,[181 187 210 179 189 177 208 131]+1,:);
QSM_LR = (QSM_L + QSM_R)/2 - repmat(medianROI(:,212,:),[1 8 1]);
QSM_mean = squeeze(nanmean(QSM_LR,1));
QSM_SD = squeeze(nanstd(QSM_LR,1));

figure('position', [100 100 1200 500]);
R_wIron = zeros(1,length(QSM_mean));
for ii = 1:10
    subplot(2,5,ii)
    errorbarxy(Iron_mean, QSM_mean(:,ii), Iron_SD, QSM_SD(:,ii),...
        {'ko',[0 0.4470 0.7410],[0 0.4470 0.7410]}); hold on;
    plot(Iron_mean, QSM_mean(:,ii),'o','Color',[0 0.4470 0.7410],...
        'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2);
    hold on;
    [F,gof] = fit(Iron_mean', QSM_mean(:,ii),'poly1');
    h = plot(F);
    h.Color = [0.5 0.5 0.5];
    [R,P] = corr(Iron_mean', QSM_mean(:,ii),'Type','Pearson');
    R_wIron(ii) = R;
    if P < 0.0001
        text(2,0.12,sprintf('r = %.3f\np < 0.0001', R));
    elseif P < 0.001
        text(2,0.12,sprintf('r = %.3f\np < 0.001', R));
    else
        text(2,0.12,sprintf('r = %.3f\np = %.3f', R, P));
    end
    
    title(strrep(strrep(QSMfile_list(ii),'_','-'),'-meanEcho',''));
    if ii > 5
        xlabel('mg Iron /100g');
    else
        xlabel('');
    end
    if ii == 1 || ii == 6
        ylabel('Susceptibility (ppm)');
    else
        ylabel('');
    end
    ylim([-0.05 0.15]);
    legend('off');
end
export_fig('Iron', '-png','-transparent'); % close;

%% plot Iron example

set(0,'defaultAxesFontSize',14);

ii = 6;
figure('position', [100 100 300 500]);
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
legend(s,{'GP','RN','SN','PU','DN','CN','TH','WM'},'box','off','location','best');
% title(strrep(strrep(QSMfile_list(ii),'_','-'),'-meanEcho',''));
xlabel('mg Iron /100g Weight');
ylabel('Susceptibility (ppm)'); ylim([-0.05 0.15]);
export_fig(['Iron-' strrep(strrep(QSMfile_list{ii},'_','-'),'-meanEcho','')], ...
    '-png','-transparent'); % close;

%% with Iron-Age equation

set(0,'defaultAxesFontSize',14);

% [WM(ACR) GP CD PU]
QSM_L = medianROI(:,[131 181 177 179],:);
QSM_R = medianROI(:,[131 181 177 179]+1,:);
QSM_LR = (QSM_L + QSM_R)/2 - repmat(medianROI(:,212,:),[1 4 1]);
% QSM_LR = reshape(QSM_LR,[],10);

Iron_WM = 3.95*(1-exp(-0.10*ageList))+0.31;
Iron_GP = 21.41*(1-exp(-0.09*ageList))+0.37;
Iron_CD = 9.66*(1-exp(-0.05*ageList))+0.33;
Iron_PU = 14.62*(1-exp(-0.04*ageList))+0.46;

Iron_age = [Iron_WM Iron_GP Iron_CD Iron_PU];

clear s;
figure('position', [100 100 1500 500]);
R_wIronAge = zeros(1,length(QSM_mean));
for ii = 1:10
    subplot(2,5,ii)
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
    if ii == 6
        legend(s,{'WM','GP','CN','PU'},...
                'box','off','location','best');
    else
        legend('off');
    end
    title(strrep(strrep(QSMfile_list(ii),'_','-'),'-meanEcho',''));
    if ii > 5
        xlabel('mg Iron /100g');
    else
        xlabel('')
    end
    if mod(ii,5) == 1
        ylabel('Susceptibility (ppm)');
    else
        ylabel('')
    end
    ylim([-0.05 0.15]);
end
% legend('WM','GP','CD','PU');

export_fig('Iron_age', '-png','-transparent'); % close;

%% plot Iron example

set(0,'defaultAxesFontSize',12);

ii = 6;
figure('position', [100 100 300 300]);
clear s
for nn = 1:4
    s(nn) = plot(Iron_age(:,nn), QSM_LR(:,nn,ii), '+'); hold on;
end
xlim([0 30]); ylim([-0.05 0.2]);
QSM_temp = QSM_LR(:,:,ii);
[F,~] = fit(Iron_age(:), QSM_temp(:),'poly1');
f = plot(F);
[R,P] = corrcoef(Iron_age(:), QSM_temp(:));
R_wIronAge(ii) = R(1,2);
if P(1,2) < 0.0001
    legend([f s],{sprintf('r = %.3f\np < 0.0001', R(1,2)),'WM','GP','CN','PU'},...
        'box','off','location','best');
elseif P(1,2) < 0.001
    legend([f s],{sprintf('r = %.3f\np < 0.001', R(1,2)),'WM','GP','CN','PU'},...
        'box','off','location','best');
else
    legend([f s],{sprintf('r = %.3f\np = %.3f', R(1,2), P(1,2)),'WM','GP','CN','PU'},...
        'box','off','location','best');
end
title(strrep(strrep(QSMfile_list(ii),'_','-'),'-meanEcho',''));
xlabel('mg Iron /100g Weight');
ylabel('Susceptibility (ppm)'); ylim([-0.05 0.2]);
export_fig(['wIronAge-' strrep(strrep(QSMfile_list{ii},'_','-'),'-meanEcho','')], ...
    '-png','-transparent'); % close;

%% iron correlation

set(0,'defaultAxesFontSize',14);

figure('position', [100 100 750 300]);
s1 = plot(1:10, R_wIron, 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2); hold on;
s2 = plot(1:10, R_wIronAge, 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.8500 0.3250 0.0980]/2); hold on;
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30;
ax.Box = 'off';
ylabel('Corr. Coef. with Iron'); ylim([0.8 1]);
xlim([0 11]);
legend({'Mean Iron','Iron(age)'},'box','off','location','best');
export_fig('Iron_R', '-png','-transparent'); % close;

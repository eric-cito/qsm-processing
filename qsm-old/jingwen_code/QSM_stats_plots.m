clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath('/home/jyao3/010_MATLAB_Utils/');
addpath('/home/jyao3/010_MATLAB_Utils/violinplot');

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList.xlsx');
statusList = [T.status];

subjList = T.b_num; subjList(~strcmp(statusList,'HC')) = [];
examList = T.t_num; examList(~strcmp(statusList,'HC')) = [];
sexList = T.sex; sexList(~strcmp(statusList,'HC')) = [];
ageList = T.age; ageList(~strcmp(statusList,'HC')) = [];

matout_root = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/Data';

Mind = strcmp(sexList,'M');
Find = strcmp(sexList,'F');

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
for ROIind = [177 179 181 187 189 210 208 115]
    figure('position', [100 100 1500 500]);
    for ii = 1:10
        subplot(2,5,ii)
        data = mean(squeeze(medianROI(:,ROIind:ROIind+1,ii)),2) - squeeze(medianROI(:,212,ii)); %
        err = mean(squeeze(stdROI(:,ROIind:ROIind+1,ii)),2);
        errorbar(ageList(Mind), data(Mind), err(Mind), 'o'); hold on;
        errorbar(ageList(Find), data(Find), err(Find), 'o');
        [F,gof] = fit(ageList(:),data(:),'poly1');
        h = plot(F);
        h.Color = [0.5 0.5 0.5];
        [R,P] = corrcoef(ageList(:),data(:));
        legend({'Male', 'Female', sprintf('fit R %.4f P %.4f', R(1,2), P(1,2))},...
            'box','off','location','best');
        title(strrep(strrep(QSMfile_list(ii),'_','-'),'-meanEcho',''));
        ylim([-0.1 0.2]); 
        ylabel('Susceptibility (ppm)');
        xlabel('Age');
        if mod(ii,5) == 1
            ylabel(ROIname{ROIind}); 
        end
    end
    export_fig([ROIname{ROIind} '_refCSF'], '-png','-transparent'); % close;
end

%% plot ROI example

set(0,'defaultAxesFontSize',14);

ROIind = 179; % [177 179 181 187 189 210 208 115]
ii = 1;
figure('position', [100 100 300 300]);
data = mean(squeeze(medianROI(:,ROIind:ROIind+1,ii)),2) - squeeze(medianROI(:,212,ii)); %
err = mean(squeeze(stdROI(:,ROIind:ROIind+1,ii)),2);
errorbar(ageList(Mind), data(Mind), err(Mind), 'o'); hold on;
errorbar(ageList(Find), data(Find), err(Find), 'o');
ylim([-0.05 0.2]); % [-0.05 0.2]
xlim([20 80]);

[F,gof] = fit(ageList(:),data(:),'poly1');
h = plot(F);
h.Color = [0.5 0.5 0.5];
[R,P] = corrcoef(ageList(:),data(:));
legend({'Male', 'Female', sprintf('R %.2e\nP %.2e', R(1,2), P(1,2))},...
    'box','off','location','best');
ylabel('Putamen Susc. (ppm)'); % Putamen
xlabel('Age (years)');
title(strrep(strrep(QSMfile_list(ii),'_','-'),'-meanEcho',''));

export_fig([ROIname{ROIind} '-' strrep(strrep(QSMfile_list{ii},'_','-'),'-meanEcho','')], ...
    '-png','-transparent'); % close;

%% output age correlation stats

meanROI_R = zeros(length(ROIlist), length(QSMfile_list));
meanROI_P = zeros(length(ROIlist), length(QSMfile_list));
meanROI_slope = zeros(length(ROIlist), length(QSMfile_list));
medianROI_R = zeros(length(ROIlist), length(QSMfile_list));
medianROI_P = zeros(length(ROIlist), length(QSMfile_list));
medianROI_slope = zeros(length(ROIlist), length(QSMfile_list));
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
        catch
            medianROI_slope(rr,ii) = nan;
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

%% with Iron

set(0,'defaultAxesFontSize',12);

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
    errorbarxy(Iron_mean, QSM_mean(:,ii), Iron_SD, QSM_SD(:,ii),{'ko','k','k'});
    hold on;
    [F,gof] = fit(Iron_mean', QSM_mean(:,ii),'poly1');
    f = plot(F);
    [R,P] = corr(Iron_mean', QSM_mean(:,ii),'Type','Pearson');
    R_wIron(ii) = R;
    legend(f,{sprintf('R %.2e\nP %.2e', R, P)},...
        'box','off','location','best');
    title(strrep(strrep(QSMfile_list(ii),'_','-'),'-meanEcho',''));
    if ii > 5
        xlabel('mg Iron /100g Weight');
    else
        xlabel('');
    end
    if ii == 1 || ii == 6
        ylabel('Susceptibility (ppm)');
    else
        ylabel('');
    end
    ylim([-0.05 0.2]);
end
export_fig('wIron', '-png','-transparent'); % close;

%% plot Iron example

set(0,'defaultAxesFontSize',10);

ii = 6;
figure('position', [100 100 300 300]);
for xx = 1:length(Iron_mean)
    errorbarxy(Iron_mean(xx), QSM_mean(xx,ii), Iron_SD(xx), QSM_SD(xx,ii),{'ko','k','k'});
    hold on;
    s(xx) = scatter(Iron_mean(xx), QSM_mean(xx,ii), 60, 'o', 'filled'); hold on;
end
% GP RN SN PU DN CN TH WM(ACR)
[F,gof] = fit(Iron_mean', QSM_mean(:,ii),'poly1');
f = plot(F);
legend(s,{'GP','RN','SN','PU','DN','CN','TH','WM'},'box','off','location','best');
% [R,P] = corr(Iron_mean', QSM_mean(:,ii),'Type','Pearson');
% legend(f,{sprintf('R %.2e\nP %.2e', R, P)},...
%     'box','off','location','best');
title(strrep(strrep(QSMfile_list(ii),'_','-'),'-meanEcho',''));
xlabel('mg Iron /100g Weight');
ylabel('Susceptibility (ppm)'); ylim([-0.05 0.2]);
export_fig(['wIron-' strrep(strrep(QSMfile_list{ii},'_','-'),'-meanEcho','')], ...
    '-png','-transparent'); % close;

%% with Iron-Age equation

set(0,'defaultAxesFontSize',10);

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

figure('position', [100 100 1500 500]);
R_wIronAge = zeros(1,length(QSM_mean));
for ii = 1:10
    subplot(2,5,ii)
    for nn = 1:4
        s(nn) = plot(Iron_age(:,nn), QSM_LR(:,nn,ii), '+'); hold on;
    end
    xlim([0 30]); ylim([-0.05 0.2]);
    QSM_temp = QSM_LR(:,:,ii);
    [F,~] = fit(Iron_age(:), QSM_temp(:),'poly1');
    f = plot(F);
    [R,P] = corrcoef(Iron_age(:), QSM_temp(:));
    R_wIronAge(ii) = R(1,2);
    legend(f,{sprintf('R %.2e\nP %.2e', R(1,2), P(1,2))},...
        'box','off','location','best');
    title(strrep(strrep(QSMfile_list(ii),'_','-'),'-meanEcho',''));
    xlabel('mg Iron /100g Fresh Weight');
    ylabel('Susceptibility (ppm)'); ylim([-0.05 0.2]);
end
% legend('WM','GP','CD','PU');

% export_fig('wIron_age', '-png','-transparent'); % close;

%% plot Iron example

set(0,'defaultAxesFontSize',10);

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
legend([f s],{sprintf('R %.2e\nP %.2e', R(1,2), P(1,2)),'WM','GP','CD','PU'},...
    'box','off','location','best');
title(strrep(strrep(QSMfile_list(ii),'_','-'),'-meanEcho',''));
xlabel('mg Iron /100g Weight');
ylabel('Susceptibility (ppm)'); ylim([-0.05 0.2]);
export_fig(['wIronAge-' strrep(strrep(QSMfile_list{ii},'_','-'),'-meanEcho','')], ...
    '-png','-transparent'); % close;

%% plot

data = [meanROI_R medianROI_R stdROI_R];
data(isnan(data)) = 0;
cmap = flip(colormap_RWG([-1 1],100));

figure('position', [100 100 500 700]);
imagesc(data,[-1 1]); colormap(cmap); colorbar; hold on;
plot([10.5 10.5],[0 208],'k');
plot([20.5 20.5],[0 208],'k');
ax = gca;
ax.XTick = [5 15 25];
ax.XTickLabel = {'Mean','Median','SD'};
ax.YTick = [1 109 177];
ax.YTickLabel = {'Cortical GM','WM','Basal Ganglia'};
title('Correlation Coefficient');
export_fig('Corr-Pearson-R', '-png','-transparent'); % close;

data = [meanROI_slope medianROI_slope stdROI_slope]*1e3;
data(isnan(data)) = 0;
cmap = flip(colormap_RWG([-2 2],100));

figure('position', [100 100 500 700]);
imagesc(data,[-1 1]); colormap(cmap); colorbar; hold on;
plot([10.5 10.5],[0 208],'k');
plot([20.5 20.5],[0 208],'k');
ax = gca;
ax.XTick = [5 15 25];
ax.XTickLabel = {'Mean','Median','SD'};
ax.YTick = [1 109 177];
ax.YTickLabel = {'Cortical GM','WM','Basal Ganglia'};
title('Fitting Slope (1000x)');
export_fig('Corr-Slope', '-png','-transparent'); % close;

r = ones(100, 1);
g = linspace(1, 0, 100)'; 
cmap = [r g g];

data = [meanROI_P medianROI_P stdROI_P];
data(isnan(data)) = 1;

figure('position', [100 100 500 700]);
imagesc(data, [0 0.05]); colormap(flip(cmap)); colorbar; hold on;
plot([10.5 10.5],[0 208],'k');
plot([20.5 20.5],[0 208],'k');
ax = gca;
ax.XTick = [5 15 25];
ax.XTickLabel = {'Mean','Median','SD'};
ax.YTick = [1 109 177];
ax.YTickLabel = {'Cortical GM','WM','Basal Ganglia'};
title('Correlation P-value');
export_fig('Corr-Pearson-P', '-png','-transparent'); % close;

%% plot heatmap - median

QSMname = strrep(QSMfile_list,'_','-');
QSMname = strrep(QSMname,'QSM-','');
QSMname = strrep(QSMname,'-meanEcho','');

data = [medianROI_R];
data(isnan(data)) = 0;
cmap = flip(colormap_RWG([-1 1],100));

figure('position', [100 100 500 700]);
imagesc(data(:,:),[-1 1]); colormap(cmap); colorbar; hold on;
plot([0 11],[109.5 109.5],'k');
plot([0 11],[177.5 177.5],'k');
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30;
ax.YTick = floor([(1+109)/2 (109+177)/2 (177+212)/2]);
ax.YTickLabelRotation = 90;
ax.YTickLabel = {'Cortical GM','WM','BG'};
title('QSM-Age Correlation Coefficient');
export_fig('Corr-Median-R', '-png','-transparent'); % close;

data = [medianROI_slope]*1e3;
data(isnan(data)) = 0;
cmap = flip(colormap_RWG([-2 2],100));

figure('position', [100 100 500 700]);
imagesc(data(:,:),[-1 1]); colormap(cmap); colorbar; hold on;
plot([0 11],[109.5 109.5],'k');
plot([0 11],[177.5 177.5],'k');
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30;
ax.YTick = floor([(1+109)/2 (109+177)/2 (177+212)/2]);
ax.YTickLabelRotation = 90;
ax.YTickLabelRotation = 90;
ax.YTickLabel = {'Cortical GM','WM','BG'};
title('QSM-Age Fitting Slope');
export_fig('Corr-Median-Slope', '-png','-transparent'); % close;

r = ones(100, 1);
g = linspace(1, 0, 100)'; 
cmap = [r g g];

data = [medianROI_P];
data(isnan(data)) = 1;

figure('position', [100 100 500 700]);
imagesc(data(:,:), [0 0.05]); colormap(flip(cmap)); colorbar; hold on;
plot([0 11],[109.5 109.5],'k');
plot([0 11],[177.5 177.5],'k');
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30;
ax.YTick = floor([(1+109)/2 (109+177)/2 (177+212)/2]);
ax.YTickLabelRotation = 90;
ax.YTickLabel = {'Cortical GM','WM','BG'};
title('QSM-Age Correlation P-value');
export_fig('Corr-Median-P', '-png','-transparent'); % close;

%% ranking

[R,ind] = sort(medianROI_slope(:,3),'descend');
P = medianROI_P(ind,3);
ROI = ROIname(ind);

num = 1:212; num = num(ind)';

T = table(ROI, R, P, num);

disp(T(1:20,:));

figure;
plot(1:212, sort(medianROI_R(:,1), 'descend'),'k-'); hold on;
plot(1:212, sort(medianROI_R(:,[2 3 4 5]), 'descend'),'+-');
plot(1:212, sort(medianROI_R(:,[6 7 8]), 'descend'),'-');
plot(1:212, sort(medianROI_R(:,[9 10]), 'descend'),'--');
xlim([1 100]);
axis tight;
legend(QSMname(:));
ylabel('Corr. Coef.');
xlabel('ROI');
export_fig('sortR', '-png','-transparent'); % close;

%% BG correlation coefficients

set(0,'defaultAxesFontSize',14);

BGind = [181 182 187 188 210 211 179 180 189 190 177 178];
QSM_L = medianROI(:,BGind(1:2:end),:);
QSM_R = medianROI(:,BGind(2:2:end),:);
QSM_LR = (QSM_L + QSM_R)/2 - repmat(medianROI(:,212,:),[1 6 1]);

medianROI_Rtemp = zeros(length(BGind)/2, length(QSMfile_list));
medianROI_Ptemp = zeros(length(BGind)/2, length(QSMfile_list));
for rr = 1:length(BGind)/2
    for ii = 1:10
        data = QSM_LR(:,rr,ii); 
        [R,P] = corr(ageList(:),data(:),'Type','Pearson');
        medianROI_Rtemp(rr,ii) = R;
        medianROI_Ptemp(rr,ii) = P;
    end
end

figure('position', [100 100 900 300]);
violinplot(medianROI_Rtemp(:,:));
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30;
ylabel('Corr. Coef. in BG'); ylim([0.2 0.7]);
xlim([0 11]);
export_fig('BG_R', '-png','-transparent'); % close;

figure('position', [100 100 500 300]);
s1 = plot(1:10, R_wIron, 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2); hold on;
s2 = plot(1:10, R_wIronAge, 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.8500 0.3250 0.0980]/2); hold on;
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30;
ylabel('Corr. Coef. with Iron'); ylim([0.8 1]);
xlim([0 11]);
legend({'Mean Iron','Iron(age)'},'box','off','location','best');
export_fig('Iron_R', '-png','-transparent'); % close;

%% scaling iLSQR and QSMGAN

X = squeeze(meanROI(1,:,3))';
Y = squeeze(meanROI(1,:,6))';
indnan = isnan(X) | isnan(Y);

X = X(~indnan); Y = Y(~indnan);

figure;
plot(X, Y, '.'); hold on;
[F,gof] = fit(X, Y,'poly1');
plot(F);
[R,P] = corrcoef(X, Y);
legend({'data', sprintf('fit R %.4f P %.4f slope %.4f intercept %.4f', ...
    R(1,2), P(1,2), F.p1, F.p2)},...
    'box','off','location','best');
xlabel('iLSQR');
ylabel('QSMGAN');

export_fig('iLSQR-QSMGAN', '-png','-transparent'); % close;

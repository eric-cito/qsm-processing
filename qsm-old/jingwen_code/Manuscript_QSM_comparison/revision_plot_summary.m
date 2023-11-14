clc; clear;

set(0,'defaultAxesFontSize',14);
addpath('/home/jyao3/010_MATLAB_Utils/spider_plot');

%% QSM info

QSM_name = {'iLSQR' 'QSIP' 'SSTGV' 'SSTV' 'STARQSM' 'FANSI' 'HDQSM' 'MEDI' ...
    'QSMGAN' 'QSMnet+' 'xQSM' 'iQSM'};

%% load data

matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data';
load([matout_root '/QSMcomp_all.mat']);

%% correlation

figure;
subplot(231);
plot_corr(T_Corr.BG_age_R, T_Corr.R_wIron);
subplot(232);
plot_corr(T_Corr.BG_age_R, T_HD.ES_all_PM);
subplot(233);
plot_corr(T_Corr.BG_age_R, T_HD.ES_all_EM);
subplot(234);
plot_corr(T_Corr.R_wIron, T_HD.ES_all_PM);
subplot(235);
plot_corr(T_Corr.R_wIron, T_HD.ES_all_EM);
subplot(236);
plot_corr(T_HD.ES_all_PM, T_HD.ES_all_EM);

%% overall performance

T = [T_COSMOS T_Corr(:,2:end) T_HD(:,2:end)];
T = T(:,[11 15 4 7 5 6]);

Tmin = min(T{:,:},[],1);
Trange = max(T{:,:},[],1) - min(T{:,:},[],1);
Tnorm = T;
Tnorm{:,:} = (Tnorm{:,:} - repmat(Tmin,[12 1]))./repmat(Trange,[12 1]);

plotColor = zeros(12,3);
plotColor(1,:) = [0.4940    0.1840    0.5560];
plotColor(2:4,:) = repmat([0.8500    0.3250    0.0980],[3 1]);
plotColor(5:8,:) = repmat([0    0.4470    0.7410],[4 1]);
plotColor(9:12,:) = repmat([0.4660    0.6740    0.1880],[4 1]);
plotStyle = {'-','-','--',':','-','--',':','-.','-','--',':','-.'};

selInd = [1];

figure('position', [100 100 600 500]);
spider_plot(Tnorm{selInd,:},'AxesLabels',{'HC-PM','HC-EM','COSMOS-SSIM','Corr-Age','Corr-Iron','Visual'}, ...
    'AxesDisplay','none','FillOption','on','MarkerSize',10,'AxesLabelsEdge','none', ...
    'Color',plotColor(selInd,:),'LineStyle',plotStyle(selInd),'LabelFontSize',14,...
    'AxesLimits',repmat([0;1],[1 size(Tnorm,2)]));
legend(QSM_name(selInd),'location','northeast','box','off','FontSize',12);

% export_fig('Summ1', '-png','-transparent'); % close;

selInd = [2:4];

figure('position', [100 100 600 500]);
spider_plot(Tnorm{selInd,:},'AxesLabels',{'HC-PM','HC-EM','COSMOS-SSIM','Corr-Age','Corr-Iron','Visual'}, ...
    'AxesDisplay','none','FillOption','on','MarkerSize',10,'AxesLabelsEdge','none', ...
    'Color',plotColor(selInd,:),'LineStyle',plotStyle(selInd),'LabelFontSize',14,...
    'AxesLimits',repmat([0;1],[1 size(Tnorm,2)]));
legend(QSM_name(selInd),'location','northeast','box','off','FontSize',12);

% export_fig('Summ2', '-png','-transparent'); % close;

selInd = [5:8];

figure('position', [100 100 600 500]);
spider_plot(Tnorm{selInd,:},'AxesLabels',{'HC-PM','HC-EM','COSMOS-SSIM','Corr-Age','Corr-Iron','Visual'}, ...
    'AxesDisplay','none','FillOption','on','MarkerSize',10,'AxesLabelsEdge','none', ...
    'Color',plotColor(selInd,:),'LineStyle',plotStyle(selInd),'LabelFontSize',14,...
    'AxesLimits',repmat([0;1],[1 size(Tnorm,2)]));
legend(QSM_name(selInd),'location','northeast','box','off','FontSize',12);

% export_fig('Summ3', '-png','-transparent'); % close;

selInd = [9:12];

figure('position', [100 100 600 500]);
spider_plot(Tnorm{selInd,:},'AxesLabels',{'HC-PM','HC-EM','COSMOS-SSIM','Corr-Age','Corr-Iron','Visual'}, ...
    'AxesDisplay','none','FillOption','on','MarkerSize',10,'AxesLabelsEdge','none', ...
    'Color',plotColor(selInd,:),'LineStyle',plotStyle(selInd),'LabelFontSize',14,...
    'AxesLimits',repmat([0;1],[1 size(Tnorm,2)]));
legend(QSM_name(selInd),'location','northeast','box','off','FontSize',12);

% export_fig('Summ4', '-png','-transparent'); % close;

%% all methods

selInd = [1 3 5 12];
plotStyle = {'-','-','-','-'};

figure('position', [100 100 600 500]);
spider_plot(Tnorm{selInd,:},'AxesLabels',{'HC-PM','HC-EM','COSMOS-SSIM','Corr-Age','Corr-Iron','Visual'}, ...
    'AxesDisplay','none','FillOption','on','MarkerSize',10,'AxesLabelsEdge','none', ...
    'Color',plotColor(selInd,:),'LineStyle',plotStyle,'LabelFontSize',14,...
    'AxesLimits',repmat([0;1],[1 size(Tnorm,2)]));
legend(QSM_name(selInd),'location','northeast','box','off','FontSize',12);

% export_fig('Summ5', '-png','-transparent'); % close;

%% function

function [] = plot_corr(A, B)

scatter(A, B);
[R,P] = corr(A, B);
fprintf('r %.3f p %.3f\n', R, P);

end
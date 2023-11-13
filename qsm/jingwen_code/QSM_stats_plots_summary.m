clc; clear;

set(0,'defaultAxesFontSize',14);
addpath('/home/jyao3/010_MATLAB_Utils/spider_plot');

%% QSM info

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

QSMorder = [3 10 1 2 4 5 8 9 6 7];
QSMfile_list = QSMfile_list(QSMorder);
QSMname = strrep(QSMfile_list,'_','-');

%% overall performance

T = readtable('/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/summarize.xlsx');
T.Properties.RowNames = T.Var1;
T = T([2 3 5 6 7 9],2:end);
T{'Visual',:} = -T{'Visual',:};

T{'SSIM',:} = [0.932752807548321;0.942135128340641;0.924763110013671;0.923193571837909;...
    0.920339928206209;0.915461092413069;0.930331439384648;0.930215043115122;...
    0.956621542143198;0.942349362155985]';
T{'CorrAge',:} = [0.5849    0.5795    0.6201    0.5295    0.5901    0.5674    0.6012    0.5964    0.4749    0.6582];
T{'CorrIron',:} = [0.9420    0.9463    0.9048    0.9395    0.9417    0.9605    0.9614    0.9575    0.9382    0.9156];
% T{'NMSE',:} = [0.979399706931571;0.731487962737004;0.976372392791937;0.967994235529792;1.15519168676460;1.33291126607560;1.03339829484722;1.06528165316647;0.467867600455459;0.623944248064062]';
T{'HD-PMHC',:} = [2.6120    2.4054    2.1693    2.3464    2.5235    2.3464    2.6120    2.4644    1.4609    2.0512];
T{'HD-EMHC',:} = [4.6932    4.3464    4.5848    4.5848    4.5848    4.1730    4.4331    4.2163    3.1324    3.9562];

Tmin = min(T{:,:},[],2);
Trange = max(T{:,:},[],2) - min(T{:,:},[],2);
Tnorm = T;
Tnorm{:,:} = (Tnorm{:,:} - repmat(Tmin,[1 10]))./repmat(Trange,[1 10]);

plotColor = zeros(10,3);
plotColor(1,:) = [0.4940    0.1840    0.5560];
plotColor(2:5,:) = repmat([0    0.4470    0.7410],[4 1]);
plotColor(6:8,:) = repmat([0.8500    0.3250    0.0980],[3 1]);
plotColor(9:10,:) = repmat([0.4660    0.6740    0.1880],[2 1]);
plotStyle = {'-','-','--',':','-.','-','--',':','-','--'};

selInd = [1];

figure('position', [100 100 600 500]);
spider_plot(Tnorm{:,selInd}','AxesLabels',{'HC-PM','HC-EM','COSMOS-SSIM','Corr-Age','Corr-Iron','Visual'}, ...
    'AxesDisplay','none','FillOption','on','MarkerSize',10,'AxesLabelsEdge','none', ...
    'Color',plotColor(selInd,:),'LineStyle',plotStyle(selInd),'LabelFontSize',14,...
    'AxesLimits',repmat([0;1],[1 size(Tnorm,1)]));
legend(QSMname(selInd),'location','northeast','box','off','FontSize',12);

export_fig('Summ1', '-png','-transparent'); % close;

selInd = [2:5];

figure('position', [100 100 600 500]);
spider_plot(Tnorm{:,selInd}','AxesLabels',{'HC-PM','HC-EM','COSMOS-SSIM','Corr-Age','Corr-Iron','Visual'}, ...
    'AxesDisplay','none','FillOption','on','MarkerSize',10,'AxesLabelsEdge','none', ...
    'Color',plotColor(selInd,:),'LineStyle',plotStyle(selInd),'LabelFontSize',14,...
    'AxesLimits',repmat([0;1],[1 size(Tnorm,1)]));
legend(QSMname(selInd),'location','northeast','box','off','FontSize',12);

export_fig('Summ2', '-png','-transparent'); % close;

selInd = [6:8];

figure('position', [100 100 600 500]);
spider_plot(Tnorm{:,selInd}','AxesLabels',{'HC-PM','HC-EM','COSMOS-SSIM','Corr-Age','Corr-Iron','Visual'}, ...
    'AxesDisplay','none','FillOption','on','MarkerSize',10,'AxesLabelsEdge','none', ...
    'Color',plotColor(selInd,:),'LineStyle',plotStyle(selInd),'LabelFontSize',14,...
    'AxesLimits',repmat([0;1],[1 size(Tnorm,1)]));
legend(QSMname(selInd),'location','northeast','box','off','FontSize',12);

export_fig('Summ3', '-png','-transparent'); % close;

selInd = [9:10];

figure('position', [100 100 600 500]);
spider_plot(Tnorm{:,selInd}','AxesLabels',{'HC-PM','HC-EM','COSMOS-SSIM','Corr-Age','Corr-Iron','Visual'}, ...
    'AxesDisplay','none','FillOption','on','MarkerSize',10,'AxesLabelsEdge','none', ...
    'Color',plotColor(selInd,:),'LineStyle',plotStyle(selInd),'LabelFontSize',14,...
    'AxesLimits',repmat([0;1],[1 size(Tnorm,1)]));
legend(QSMname(selInd),'location','northeast','box','off','FontSize',12);

export_fig('Summ4', '-png','-transparent'); % close;

selInd = [1 5 7 10];
plotStyle = {'-','-','-','-'};

figure('position', [100 100 600 500]);
spider_plot(Tnorm{:,selInd}','AxesLabels',{'HC-PM','HC-EM','COSMOS-SSIM','Corr-Age','Corr-Iron','Visual'}, ...
    'AxesDisplay','none','FillOption','on','MarkerSize',10,'AxesLabelsEdge','none', ...
    'Color',plotColor(selInd,:),'LineStyle',plotStyle,'LabelFontSize',14,...
    'AxesLimits',repmat([0;1],[1 size(Tnorm,1)]));
legend(QSMname([1 5 7 10]),'location','northeast','box','off','FontSize',12);

export_fig('Summ5', '-png','-transparent'); % close;

%% correlation

T = readtable('/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/summarize.xlsx');
T.Properties.RowNames = T.Var1;
T = T(:,2:end);
T{'Visual',:} = -T{'Visual',:};

T{'SSIM',:} = [0.932752807548321;0.942135128340641;0.924763110013671;0.923193571837909;0.920339928206209;0.915461092413069;0.930331439384648;0.930215043115122;0.956621542143198;0.942349362155985]';

[r,p] = corr(T{'HD-EMHC',:}',T{'CorrIron',:}','Type','Spearman')
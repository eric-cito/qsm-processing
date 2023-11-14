clear; clc;
warning('off');

%% add path

addpath('/home/jyao3/010_MATLAB_Utils/NIfTI');
addpath('/home/jyao3/010_MATLAB_Utils/violinplot');
addpath('/home/jyao3/010_MATLAB_Utils/export_fig');

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

img_path = '/working/lupolab/jingwen/001_QSM/temp';

%% read in list

ListFile = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/VisualList.mat';
load(ListFile);

order = T.Order;
T = sortrows(T,'Order');

%% QSM method list

QSMfile_list = {'QSM_iLSQR_meanEcho' ...
    'QSM_STARQSM_meanEcho' ...
    'QSM_FANSI_nonlinearTV_meanEcho' ...
    'QSM_HDQSM_meanEcho' ...
    'QSM_MEDI_meanEcho' ...
    'QSM_QSIP_meanEcho' ...
    'QSM_SSTGV_meanEcho' ...
    'QSM_SSTV_meanEcho' ...
    'QSM_QSMGAN_meanEcho' ...
    'QSM_QSMnet_meanEcho'};

QSMlist = zeros(size(T.Order));
for ii = 1:length(QSMfile_list)
    QSMlist(strcmp(T.QSMfile, QSMfile_list{ii})) = ii;
end

%% read in result

assess1 = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/QSM_VisualScore_JVM.xlsx';
Tresult = readtable(assess1,'Sheet','ScoringSheet','Range','B:D');

Sstreaking1 = Tresult.Var1;
Sunnatural1 = Tresult.Var2;
Snoise1 = Tresult.Var3;

assess2 = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/QSM_VisualScore_Johanna.xlsx';
Tresult = readtable(assess2,'Sheet','ScoringSheet','Range','B:D');

Sstreaking2 = Tresult.Var1;
Sunnatural2 = Tresult.Var2;
Snoise2 = Tresult.Var3;

%% histogram

figure('position', [100 100 900 400]);
subplot(131);
histogram(Sstreaking1);
xlabel('Score'); ylabel('Streaking');
subplot(132);
histogram(Sunnatural1);
xlabel('Score'); ylabel('Unnaturalness');
subplot(133);
histogram(Snoise1);
xlabel('Score'); ylabel('Noise');

export_fig([img_root 'ScoreHist_rater1'], '-png','-transparent');

figure('position', [100 100 900 400]);
subplot(131);
histogram(Sstreaking2);
xlabel('Score'); ylabel('Streaking');
subplot(132);
histogram(Sunnatural2);
xlabel('Score'); ylabel('Unnaturalness');
subplot(133);
histogram(Snoise2);
xlabel('Score'); ylabel('Noise');

export_fig([img_root 'ScoreHist_rater2'], '-png','-transparent');

%% average the two raters

Sstreaking = 0.5*(Sstreaking1 + Sstreaking2);
Sunnatural = 0.5*(Sunnatural1 + Sunnatural2);
Snoise = 0.5*(Snoise1 + Snoise2);

%% calculate mean

meanScores = zeros(10,3);
stdScores = zeros(10,3);
for ii = 1:length(QSMfile_list)
    meanScores(ii,1) = mean(Sstreaking(QSMlist == ii));
    meanScores(ii,2) = mean(Sunnatural(QSMlist == ii));
    meanScores(ii,3) = mean(Snoise(QSMlist == ii));
    stdScores(ii,1) = std(Sstreaking(QSMlist == ii));
    stdScores(ii,2) = std(Sunnatural(QSMlist == ii));
    stdScores(ii,3) = std(Snoise(QSMlist == ii));
end

meanAll = mean(meanScores,2);
maxAll = max(meanScores,[],2);

%% plot

set(0,'defaultAxesFontSize',14);

QSMname = strrep(QSMfile_list,'_','_');
QSMname = strrep(QSMname,'_meanEcho','');
QSMname = strrep(QSMname,'_nonlinearTV','');
QSMname = strrep(QSMname,'QSM_','');
QSMname = strrep(QSMname,'QSMnet','QSMnet+');

figure('position', [100 100 1500 900]);
subplot(321);
plotScore(Sstreaking, QSMlist, QSMname);
ylabel('Streaking'); 

subplot(323);
plotScore(Sunnatural, QSMlist, QSMname);
ylabel('Unnaturalness'); 

subplot(325);
plotScore(Snoise, QSMlist, QSMname);
ylabel('Noise'); 

subplot(322);
plot(1:10, meanScores(:,1), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2); hold on;
plot(1:10, meanScores(:,2), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.8500 0.3250 0.0980]/2);
plot(1:10, meanScores(:,3), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.9290 0.6940 0.1250]/2);
plot(1:10, meanAll, 'o-', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.4940 0.1840 0.5560]/2);
% hold on; plot([0 11], 1.96*[1 1], 'k--');
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30; xlim([0 11]); ylim([-0.5 3.5]);
ax.Box = 'off';
legend({'Streaking','Unnaturalness','Noise','Mean of Scores'},'box','off','location','best','FontSize',12);
ylabel('Mean Score');

export_fig([img_root 'ViusalResult'], '-png','-transparent');

%% plot functions

function [] = plotScore(Score, QSMlist, QSMname)

v = violinplot(Score, QSMlist, 'showMean', true); hold on;
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
plot([0 11], [1 1]*mean(Score), 'k:');
ax = gca;
ax.XTick = 1:10;
ax.XTickLabel = QSMname(:);
ax.XTickLabelRotation = 30;
xlim([0 11]); ylim([-0.5 3.5]);

end
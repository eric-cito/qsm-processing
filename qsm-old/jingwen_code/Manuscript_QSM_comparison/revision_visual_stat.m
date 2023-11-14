clear; clc;
warning('off');

%% add path

addpath('/home/jyao3/010_MATLAB_Utils/NIfTI');
addpath('/home/jyao3/010_MATLAB_Utils/violinplot');
addpath('/home/jyao3/010_MATLAB_Utils/export_fig');

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

img_path = '/working/lupolab/jingwen/001_QSM/temp';

%% QSM method list

QSMfile_list = {'QSM_iLSQR_meanEcho' ...
    'QSM_QSIP_meanEcho' ...
    'QSM_SSTGV_meanEcho' ...
    'QSM_SSTV_meanEcho' ...
    'QSM_STARQSM_meanEcho' ...
    'QSM_FANSI_nonlinearTV_meanEcho' ...
    'QSM_HDQSM_meanEcho' ...
    'QSM_MEDI_meanEcho' ...
    'QSM_QSMGAN_meanEcho' ...
    'QSM_QSMnet_meanEcho' ...
    'QSM_xQSM2_meanEcho' ...
    'QSM_iQSM2_meanEcho'};

%% read in result

assess1 = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/QSM_VisualScore_result.xlsx';
Tresult = readtable(assess1,'Sheet','ScoringSheet','Range','A1:F121');

Sstreaking1 = Tresult.Streaking;
Sunnatural1 = Tresult.Unnaturalness;
Snoise1 = Tresult.Noise;
QSMlist = Tresult.Method;

Tresult = readtable(assess1,'Sheet','ScoringSheet','Range','N1:P121');

Sstreaking2 = Tresult.Streaking;
Sunnatural2 = Tresult.Unnaturalness;
Snoise2 = Tresult.Noise;

%% heatmap

figure;
subplot(131);
tbl = table(Sstreaking1, Sstreaking2);
heatmap(tbl,'Sstreaking1', 'Sstreaking2');

subplot(132);
tbl = table(Sunnatural1, Sunnatural2);
heatmap(tbl,'Sunnatural1', 'Sunnatural2');

subplot(133);
tbl = table(Snoise1, Snoise2);
heatmap(tbl,'Snoise1', 'Snoise2');

%% average the two raters

Sstreaking = mean([Sstreaking1 Sstreaking2],2,'omitnan');
Sunnatural = mean([Sunnatural1 Sunnatural2],2,'omitnan');
Snoise = mean([Snoise1 Snoise2],2,'omitnan');

%% calculate mean

meanScores = zeros(12,3);
stdScores = zeros(12,3);
for ii = 1:length(QSMfile_list)
    meanScores(ii,1) = mean(Sstreaking(strcmp(QSMlist,QSMfile_list{ii})));
    meanScores(ii,2) = mean(Sunnatural(strcmp(QSMlist,QSMfile_list{ii})));
    meanScores(ii,3) = mean(Snoise(strcmp(QSMlist,QSMfile_list{ii})));
    stdScores(ii,1) = std(Sstreaking(strcmp(QSMlist,QSMfile_list{ii})));
    stdScores(ii,2) = std(Sunnatural(strcmp(QSMlist,QSMfile_list{ii})));
    stdScores(ii,3) = std(Snoise(strcmp(QSMlist,QSMfile_list{ii})));
end

meanAll = mean(meanScores,2);
maxAll = max(meanScores,[],2);

%% plot

set(0,'defaultAxesFontSize',14);

QSM_name = {'iLSQR' 'QSIP' 'SSTGV' 'SSTV' 'STARQSM' 'FANSI' 'HDQSM' 'MEDI' ...
    'QSMGAN' 'QSMnet+' 'xQSM' 'iQSM'};

QSMlist_cat = categorical(QSMlist, QSMfile_list, QSM_name);

figure('position', [100 100 1500 900]);
subplot(321);
plotScore(Sstreaking, QSMlist_cat, QSM_name);
ylabel('Streaking'); 

subplot(323);
plotScore(Sunnatural, QSMlist_cat, QSM_name);
ylabel('Unnaturalness'); 

subplot(325);
plotScore(Snoise, QSMlist_cat, QSM_name);
ylabel('Noise'); 

subplot(322);
plot(1:length(QSM_name), meanScores(:,1), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0 0.4470 0.7410]/2); hold on;
plot(1:length(QSM_name), meanScores(:,2), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.8500 0.3250 0.0980]/2);
plot(1:length(QSM_name), meanScores(:,3), 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.9290 0.6940 0.1250]/2);
plot(1:length(QSM_name), meanAll, 'o-', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5] + [0.4940 0.1840 0.5560]/2);
ax = gca;
ax.XTick = 1:length(QSM_name);
ax.XTickLabel = QSM_name(:);
% ax.XTickLabelRotation = 30; 
xlim([0 length(QSM_name)+1]); ylim([-0.5 3.5]);
ax.Box = 'off';
legend({'Streaking','Unnaturalness','Noise','Mean of Scores'},'box','off','location','best','FontSize',12);
ylabel('Mean Score');

export_fig([img_path '/ViusalResult'], '-png','-transparent');

%% plot functions

function [] = plotScore(Score, QSMlist, QSMname)

v = violinplot(Score, QSMlist, 'showMean', true); hold on;
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
plot([0 length(QSMname)+1], [1 1]*mean(Score,'omitnan'), 'k:');
ax = gca;
ax.XTick = 1:length(QSMname);
ax.XTickLabel = QSMname(:);
% ax.XTickLabelRotation = 30;
xlim([0 length(QSMname)+1]); ylim([-0.5 3.5]);

end
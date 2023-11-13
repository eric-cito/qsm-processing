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
T_orig = sortrows(T,'Order');

ListFile = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/VisualList_revision_120.mat';
load(ListFile);

order_rev = T.Order;
T_rev = sortrows(T,'Order');

ListFile = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/VisualList_revision.mat';
load(ListFile);

order_rev = T.Order;
T_rev40 = sortrows(T,'Order');

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
    'QSM_QSMnet_meanEcho' ...
    'QSM_xQSM2_meanEcho' ...
    'QSM_iQSM2_meanEcho'};

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

assess2 = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/QSM_VisualScore_revision_JL.xlsx';
Tresult = readtable(assess2,'Sheet','ScoringSheet','Range','B:D');

Sstreaking2_rev = Tresult.Var1;
Sunnatural2_rev = Tresult.Var2;
Snoise2_rev = Tresult.Var3;

%% histogram

figure('position', [100 100 900 400]);
subplot(131);
histogram(Sstreaking2); hold on;
histogram(Sstreaking2_rev);
xlabel('Score'); ylabel('Streaking');
subplot(132);
histogram(Sunnatural2); hold on;
histogram(Sstreaking2_rev);
xlabel('Score'); ylabel('Unnaturalness');
subplot(133);
histogram(Snoise2); hold on;
histogram(Snoise2_rev);
xlabel('Score'); ylabel('Noise');

export_fig([img_path 'ScoreHist_rater2'], '-png','-transparent');

%% overlap cases

T_orig.subj(strcmp(T_orig.subj, 'temp_110620_keep_GET_FMRI_RAW_FILES')) = {'temp_110620_keep'};

overlapInd_orig = nan(size(order));
overlapInd_rev = nan(size(order_rev));
k = 0;
for ii = 1:height(T_rev)
    ind = find(strcmp(T_orig.subj,T_rev.subj{ii}) & strcmp(T_orig.QSMfile,T_rev.QSMfile{ii}));
    if ~isempty(ind)
        k = k+1;
        overlapInd_rev(ii) = k;
        overlapInd_orig(ind) = k;
    end
end

[~,ind_orig] = sort(overlapInd_orig);
[~,ind_rev] = sort(overlapInd_rev);

T_orig_overlap = T_orig(ind_orig,:);
T_rev_overlap = T_rev(ind_rev,:);

%% compare two

score1 = Sunnatural2(T_orig_overlap.Order);
score2 = Sunnatural2_rev(T_rev_overlap.Order(1:100));

histogram(score1); hold on;
histogram(score2);

%% average the two raters

% Sstreaking = 0.5*(Sstreaking1 + Sstreaking2);
% Sunnatural = 0.5*(Sunnatural1 + Sunnatural2);
% Snoise = 0.5*(Snoise1 + Snoise2);

Sstreaking1 = [Sstreaking1(T_orig_overlap.Order); nan(20,1)];
Sunnatural1 = [Sunnatural1(T_orig_overlap.Order); nan(20,1)];
Snoise1 = [Snoise1(T_orig_overlap.Order); nan(20,1)];

Sstreaking = mean([Sstreaking1 Sstreaking2_rev(T_rev_overlap.Order)],2,'omitnan');
Sunnatural = mean([Sunnatural1 Sunnatural2_rev(T_rev_overlap.Order)],2,'omitnan');
Snoise = mean([Snoise1 Snoise2_rev(T_rev_overlap.Order)],2,'omitnan');

%% calculate mean

meanScores = zeros(10,3);
stdScores = zeros(10,3);
for ii = 1:length(QSMfile_list)
    meanScores(ii,1) = mean(Sstreaking(strcmp(T_rev_overlap.QSMfile,QSMfile_list{ii})));
    meanScores(ii,2) = mean(Sunnatural(strcmp(T_rev_overlap.QSMfile,QSMfile_list{ii})));
    meanScores(ii,3) = mean(Snoise(strcmp(T_rev_overlap.QSMfile,QSMfile_list{ii})));
    stdScores(ii,1) = std(Sstreaking(strcmp(T_rev_overlap.QSMfile,QSMfile_list{ii})));
    stdScores(ii,2) = std(Sunnatural(strcmp(T_rev_overlap.QSMfile,QSMfile_list{ii})));
    stdScores(ii,3) = std(Snoise(strcmp(T_rev_overlap.QSMfile,QSMfile_list{ii})));
end

meanAll = mean(meanScores,2);
maxAll = max(meanScores,[],2);

%% plot

set(0,'defaultAxesFontSize',14);

QSM_name = {'iLSQR' 'QSIP' 'SSTGV' 'SSTV' 'STARQSM' 'FANSI' 'HDQSM' 'MEDI' ...
    'QSMGAN' 'QSMnet+' 'xQSM' 'iQSM'};

QSMlist = categorical(T_rev_overlap.QSMfile, QSMfile_list([1 6 7 8 2 3 4 5 9 10 11 12]), QSM_name(1:12));

figure('position', [100 100 1500 900]);
subplot(321);
plotScore(Sstreaking, QSMlist, QSM_name);
ylabel('Streaking'); 

subplot(323);
plotScore(Sunnatural, QSMlist, QSM_name);
ylabel('Unnaturalness'); 

subplot(325);
plotScore(Snoise, QSMlist, QSM_name);
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

export_fig([img_path 'ViusalResult'], '-png','-transparent');

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
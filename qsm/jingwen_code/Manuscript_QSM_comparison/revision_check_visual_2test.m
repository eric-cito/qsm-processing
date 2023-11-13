clear; clc;
warning('off');

%% add path

addpath('/home/jyao3/010_MATLAB_Utils/NIfTI');
addpath('/home/jyao3/010_MATLAB_Utils/violinplot');
addpath('/home/jyao3/010_MATLAB_Utils/export_fig');

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

img_root = '/home/jyao3/030_QSM/img_temp/';

%% read in list

ListFile = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/VisualList.mat';
load(ListFile);

order = T.Order;
Tall = sortrows(T,'Order');

%% find the ones in the second test

ListFile = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/VisualList_revision.mat';
load(ListFile);

Tsel = T;
Tsel(contains(Tsel.QSMfile,'iQSM'),:) = [];
Tsel(contains(Tsel.QSMfile,'xQSM'),:) = [];

indSel = zeros(20,1);
for ii = 1:20
    indSel(ii) = find(contains(Tall.subj,Tsel.subj{ii}) & strcmp(Tall.QSMfile,Tsel.QSMfile{ii}));
end

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

QSMlist = zeros(size(Tall.Order));
for ii = 1:length(QSMfile_list)
    QSMlist(strcmp(Tall.QSMfile, QSMfile_list{ii})) = ii;
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
histogram(Sstreaking1); hold on;
histogram(Sstreaking1(indSel));
xlabel('Score'); ylabel('Streaking');
subplot(132);
histogram(Sunnatural1); hold on;
histogram(Sunnatural1(indSel));
xlabel('Score'); ylabel('Unnaturalness');
subplot(133);
histogram(Snoise1); hold on;
histogram(Snoise1(indSel));
xlabel('Score'); ylabel('Noise');

figure('position', [100 100 900 400]);
subplot(131);
histogram(Sstreaking2); hold on;
histogram(Sstreaking2(indSel));
xlabel('Score'); ylabel('Streaking');
subplot(132);
histogram(Sunnatural2); hold on;
histogram(Sunnatural2(indSel));
xlabel('Score'); ylabel('Unnaturalness');
subplot(133);
histogram(Snoise2); hold on;
histogram(Snoise2(indSel));
xlabel('Score'); ylabel('Noise');

%% heat map

tbl = table(Sstreaking1, Sstreaking2, Sunnatural1, Sunnatural2, Snoise1, Snoise2);

heatmap(tbl,'Snoise1','Snoise2');


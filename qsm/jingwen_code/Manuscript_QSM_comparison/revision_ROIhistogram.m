clear; clc;
warning('off');

%% Add path

addpath('/home/jyao3/030_QSM/01_Code/Manuscript_BGQSM_HD');
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

%% read data

examPath = '/data/7T_hunt/b3025/t12518/swan_qsm/HDBET_allQSM/FSseg';
QSMfile = [examPath '/QSM_iLSQR_meanEcho_reg.nii.gz'];
ROIfile = [examPath '/Seg_ANTS_manual.nii.gz'];

nii = load_nii(QSMfile);
QSMmap = double(nii.img);

nii = load_nii(ROIfile);
ROImask = double(nii.img);

%% ROI histogram

ROIdata_CN = QSMmap(ROImask >=1 & ROImask <=2);
ROIdata_PU = QSMmap(ROImask >=3 & ROImask <=4);
ROIdata_GP = QSMmap(ROImask >=5 & ROImask <=8);

figure('position', [100 0 400 400]);
histogram(ROIdata_CN, [-0.1:0.01:0.3], 'EdgeColor', 'none'); hold on;
histogram(ROIdata_PU, [-0.1:0.01:0.3], 'EdgeColor', 'none');
histogram(ROIdata_GP, [-0.1:0.01:0.3], 'EdgeColor', 'none');
legend({'CN','PU','GP'});

xlabel('Susceptibility (ppm)');
ylabel('Voxel count');

img_path = '/working/lupolab/jingwen/001_QSM/temp/';
export_fig([img_path 'ROIhist'], '-png','-transparent'); % close;

% figure;
% histogram(ROIdata_PU); hold on;
% plot(mean(ROIdata_PU)*[1 1], [0 500], '--');
% plot(median(ROIdata_PU)*[1 1], [0 500], '-');

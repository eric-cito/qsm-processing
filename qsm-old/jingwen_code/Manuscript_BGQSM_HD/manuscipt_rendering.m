clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/010_MATLAB_Utils/'));

dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/SPM_results/QSM';
img_root = '/working/lupolab/jingwen/001_QSM/temp';

%% load images

nii = load_nii(['/working/lupolab/jingwen/001_QSM/03_QSM_HD/Atlas/Seg_Subcortical.nii.gz']);
ROImask = double(nii.img);

nii = load_nii([dataPath '/qsm1fa2over3.nii.gz']);
Statmask = double(nii.img);

nii = load_nii(['/working/lupolab/jingwen/001_QSM/03_QSM_HD/Atlas/MNI_T1_brain.nii.gz']);
Brain_img = double(nii.img);
Brain_mask = double(nii.img > 0);

%% render images

StatmaskPU = Statmask.*(ROImask == 3 | ROImask == 4);

figure('position', [100 100 500 500]);

plotMask = smooth3(ROImask == 3 | ROImask == 4,'gaussian',[5 5 5],2);
hisoBG = patch(isosurface(plotMask,0.1),...
    'FaceColor',[0.9 0.9 0.9], ...
    'FaceAlpha',0.2, ...
    'EdgeColor','none');
isonormals(plotMask,hisoBG); hold on;

plotCNMask = smooth3(StatmaskPU == 1 | StatmaskPU == 3,'gaussian',[3 3 3],2);
hisoCN = patch(isosurface(plotCNMask,0.1),...
    'FaceColor', [0.8500 0.3250 0.0980], ...
    'FaceAlpha',1, ...
    'EdgeColor','none');
isonormals(plotCNMask,hisoCN);

plotCNMask = smooth3(StatmaskPU == 2 | StatmaskPU == 3,'gaussian',[3 3 3],2);
hisoCN = patch(isosurface(plotCNMask,0.1),...
    'FaceColor', [0.3010 0.7450 0.9330], ...
    'FaceAlpha',1, ...
    'EdgeColor','none');
isonormals(plotCNMask,hisoCN);

view(30,0);
axis tight equal off

lightangle(45,45);
lighting gouraud

pause(1); export_fig([img_root '/StatsPU'], '-png','-transparent');

%% render images

StatmaskPU = Statmask.*(ROImask == 1 | ROImask == 2);

figure('position', [100 100 500 500]);

plotMask = smooth3(ROImask == 1 | ROImask == 2,'gaussian',[3 3 3],2);
hisoBG = patch(isosurface(plotMask,0.1),...
    'FaceColor',[0.9 0.9 0.9], ...
    'FaceAlpha',0.2, ...
    'EdgeColor','none');
isonormals(plotMask,hisoBG); hold on;

plotCNMask = smooth3(StatmaskPU == 1 | StatmaskPU == 3,'gaussian',[3 3 3],2);
hisoCN = patch(isosurface(plotCNMask,0.1),...
    'FaceColor', [0.8500 0.3250 0.0980], ...
    'FaceAlpha',1, ...
    'EdgeColor','none');
isonormals(plotCNMask,hisoCN);

plotCNMask = smooth3(StatmaskPU == 2 | StatmaskPU == 3,'gaussian',[3 3 3],2);
hisoCN = patch(isosurface(plotCNMask,0.1),...
    'FaceColor', [0.3010 0.7450 0.9330], ...
    'FaceAlpha',1, ...
    'EdgeColor','none');
isonormals(plotCNMask,hisoCN);

view(30,0);
axis tight equal off

lightangle(45,45);
lighting gouraud

pause(1); export_fig([img_root '/StatsCN'], '-png','-transparent');
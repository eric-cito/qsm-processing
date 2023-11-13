clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/'));

segPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/SPM_results/QSM/qsm1fa2over3.nii.gz';
atlasPath = '/netopt/rhel7/fsl/data/atlases/Striatum/striatum-con-label-thr25-7sub-1mm.nii.gz';

qsmTPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/SPM_results/QSM/HDHCt_QSM.nii';
faTPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/SPM_results/QSM/HDHCt_FA.nii';

%% load data

nii = load_nii(segPath);
qsmMask = nii.img == 1 | nii.img == 3;
faMask = nii.img == 2 | nii.img == 3;

nii = load_nii(atlasPath);
segMask = nii.img;
limbMask = segMask == 1;
execMask = segMask == 2;
sensMask = segMask == 3;

nii = load_nii(qsmTPath);
qsmT = nii.img;
nii = load_nii(faTPath);
faT = nii.img;

%% show images

figure;
subplot(1,2,1);
imagesc(qsmMask(:,:,77) + 2*faMask(:,:,77)); axis tight off equal
subplot(1,2,2);
imagesc(segMask(:,:,77)); axis tight off equal

%% calculate percentage

Nstriatum = sum(segMask(:) > 0);

NQSM = [sum(segMask(:).*qsmMask(:) == 1) ...
    sum(segMask(:).*qsmMask(:) == 2) ...
    sum(segMask(:).*qsmMask(:) == 3) ...
    sum(segMask(:).*qsmMask(:) == 4) ...
    sum(segMask(:).*qsmMask(:) == 5) ...
    sum(segMask(:).*qsmMask(:) == 6) ...
    sum(segMask(:).*qsmMask(:) == 7)]';

NFA = [sum(segMask(:).*faMask(:) == 1) ...
    sum(segMask(:).*faMask(:) == 2) ...
    sum(segMask(:).*faMask(:) == 3) ...
    sum(segMask(:).*faMask(:) == 4) ...
    sum(segMask(:).*faMask(:) == 5) ...
    sum(segMask(:).*faMask(:) == 6) ...
    sum(segMask(:).*faMask(:) == 7)]';

Nall = [sum(segMask(:) == 1) ...
    sum(segMask(:) == 2) ...
    sum(segMask(:) == 3) ...
    sum(segMask(:) == 4) ...
    sum(segMask(:) == 5) ...
    sum(segMask(:) == 6) ...
    sum(segMask(:) == 7)]';

Nnone = Nall - NQSM - NFA;

%% plot barplot

figure('Position', [100 100 500 500]);
X = {'Limbic','Executive','Rostral Motor','Caudal Motor','Parietal','Occipital','Temporal'};
bar([NQSM NFA]./Nall*100); ylim([0 40]);
xticks(1:7);
xticklabels(X);
ylabel('Striatal Subregion Fraction (%)'); 
legend({'QSM t > 3','FA t > 3'}, 'box', 'off', 'location', 'best');

img_root = '/working/lupolab/jingwen/001_QSM/temp';
export_fig([img_root '/subregion'], '-png','-transparent');

%% plot histogram

% figure('position',[100 100 600 800]);
% subplot(3,1,1);
% histogram(qsmT(limbMask & qsmT ~= 0),[-6:0.5:6]); hold on
% histogram(faT(limbMask & faT ~= 0),[-6:0.5:6]);
% subplot(3,1,2);
% histogram(qsmT(execMask & qsmT ~= 0),[-6:0.5:6]); hold on;
% histogram(faT(execMask & faT ~= 0),[-6:0.5:6]);
% subplot(3,1,3);
% histogram(qsmT(sensMask & qsmT ~= 0),[-6:0.5:6]); hold on;
% histogram(faT(sensMask & faT ~= 0),[-6:0.5:6]);

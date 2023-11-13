clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/'));

dataPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/SubjImgDataBG_ANTS_iLSQR_1004.mat';
img_root = '/home/jyao3/030_QSM/img_temp';

%% Load data

load(dataPath);

indHC = strcmp({HD_BGanalysis.group},'HC');
indPM = strcmp({HD_BGanalysis.group},'PM');
indEM = strcmp({HD_BGanalysis.group},'EM');
indMan = strcmp({HD_BGanalysis.group},'Manifest');

statusList = nan(1,length(HD_BGanalysis));

age = [HD_BGanalysis.age];
sex = {HD_BGanalysis.sex};
CAG = [HD_BGanalysis.CAG];
volList = zeros(3,length(HD_BGanalysis));

statusList(indHC) = 1;
statusList(indPM) = 2;
statusList(indEM) = 3;
statusList(indMan) = 3;
statusList(CAG < 36) = nan;

%% demographics

% CAPS = [HD_BGanalysis.CAPS];
CAPS = age.*(CAG-35.5);
CAPS(isnan(statusList)) = nan;

CAG = [HD_BGanalysis.CAG];
TMS = [HD_BGanalysis.TMS];
DCL = [HD_BGanalysis.DCL];
TFC = [HD_BGanalysis.TFC];

ageOnset = zeros(1,length(HD_BGanalysis));
for ii = 1:length(HD_BGanalysis)
    if isnan(CAG(ii))
        ageOnset(ii) = nan;
    else
        ageOnset(ii) = medianYearOnset(CAG(ii));
    end
end

AOO = ageOnset;
YTO = AOO - age;

indPM = statusList == 2;

fprintf('=== PM === \n');
printDemo(indPM, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

indEM = statusList == 3;

fprintf('=== EM === \n');
printDemo(indEM, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

indPMfar = indPM & YTO > 15;
indPMnear = indPM & YTO <= 15;

indEM1 = (indEM | indMan) & DCL < 4;
indEM2 = (indEM | indMan) & DCL >= 4;

fprintf('=== PM far === \n');
printDemo(indPMfar, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

fprintf('=== PM near === \n');
printDemo(indPMnear, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

fprintf('=== EM DCL < 4 === \n');
printDemo(indEM1, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

fprintf('=== EM DCL >= 4 === \n');
printDemo(indEM2, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC);

%% load QSM data

nii = load_nii(['/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/allQSMnorm.nii.gz']);
allQSM = double(nii.img);

nii = load_nii(['/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/Seg_atlas.nii.gz']);
ROImask = double(nii.img);

nii = load_nii(['/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/MNI_T1.nii']);
Brain_mask = double(nii.img > 0);

%% plot

figure;
plotBGMask = smooth3(Brain_mask > 0);
hisoBG = patch(isosurface(plotBGMask,0.1),...
    'FaceColor',[0.75 1 1], ...
    'FaceAlpha',0.5, ...
    'EdgeColor','none');
isonormals(plotBGMask,hisoBG); hold on;

plotCNMask = smooth3(ROImask == 1);
hisoCN = patch(isosurface(plotCNMask,0.1),...
    'FaceColor',[1 0 1], ...
    'EdgeColor','none');
isonormals(plotCNMask,hisoCN);

view(30,15);
axis tight equal off

lightangle(45,45);
lighting gouraud

%% voxelwise analysis

QSMlist = reshape(allQSM,[],size(allQSM,4));

indValid = find(Brain_mask > 0);

tbl = table(age',sex',statusList',zeros(size(age))',...
    'VariableNames',{'Age','Sex','Group','QSM'});
tbl.Sex = categorical(tbl.Sex);
tbl.Group = categorical(tbl.Group);

tList = zeros(size(QSMlist,1),4);
pList = zeros(size(QSMlist,1),4);
f = waitbar(0,'1','Name','Processing voxel ...');
for ii = 1:length(indValid)
    tbl.QSM = QSMlist(indValid(ii),:)';
    mdl = fitlm(tbl);
    tList(indValid(ii),:) = mdl.Coefficients{{'Age','Sex_M','Group_2','Group_3'},'tStat'};
    pList(indValid(ii),:) = mdl.Coefficients{{'Age','Sex_M','Group_2','Group_3'},'pValue'};
    waitbar(ii/length(indValid),f,sprintf('%.3f%% done.',ii/length(indValid)*100))
end

%% save to matrix

tImage = reshape(tList,[size(Brain_mask) 4]);
pImage = reshape(pList,[size(Brain_mask) 4]);

nii = load_nii(['/working/lupolab/jingwen/001_QSM/03_QSM_HD/allQSMnorm.nii.gz']);
nii.hdr.dime.dim(5) = 4;
nii.img = tImage;
save_nii(nii, ['/working/lupolab/jingwen/001_QSM/03_QSM_HD/tStat.nii.gz']);
nii.img = pImage;
save_nii(nii, ['/working/lupolab/jingwen/001_QSM/03_QSM_HD/pVal.nii.gz']);

%% functions

function [] = printDemo(indEM, age, sex, CAG, CAPS, AOO, YTO, TMS, DCL, TFC)

fprintf('Age mean %.3f SD %.3f \n', mean(age(indEM)), std(age(indEM)));
fprintf('M %i F %i \n', sum(strcmp(sex(indEM),'M')), sum(strcmp(sex(indEM),'F')));
fprintf('CAG mean %.3f SD %.3f \n', mean(CAG(indEM)), std(CAG(indEM)));
fprintf('CAG range %.3f - %.3f \n', min(CAG(indEM)), max(CAG(indEM)));
fprintf('CAPS median %.3f range %.3f - %.3f \n', median(CAPS(indEM)), min(CAPS(indEM)), max(CAPS(indEM)));
fprintf('AOO median %.3f range %.3f - %.3f \n', median(AOO(indEM)), min(AOO(indEM)), max(AOO(indEM)));
fprintf('YTO median %.3f range %.3f - %.3f \n', median(YTO(indEM)), min(YTO(indEM)), max(YTO(indEM)));
fprintf('TMS median %.3f range %.3f - %.3f \n', median(TMS(indEM)), min(TMS(indEM)), max(TMS(indEM)));
fprintf('DCL median %.3f range %.3f - %.3f \n', median(DCL(indEM)), min(DCL(indEM)), max(DCL(indEM)));
fprintf('TFC median %.3f range %.3f - %.3f \n', median(TFC(indEM)), min(TFC(indEM)), max(TFC(indEM)));

end

function [ageOnset] = medianYearOnset(CAG)

syms f(x)
f(x) = (1+exp(pi/sqrt(3)*(-21.54-exp(9.56-0.146*CAG)+x)./(sqrt(35.55+exp(17.72-0.327*CAG))))).^-1 - 0.5;

tmp = vpasolve(f);
ageOnset = double(tmp);

end
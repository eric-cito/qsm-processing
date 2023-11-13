clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/031_HD_NDM'));

FSresult_root = '/working/lupolab/jingwen/002_HD_NDM/AllT1N4';
matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/';

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220209.xlsx','Sheet','NoRep');
T(T.CAG < 36,:) = [];
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status_reclass];
ageList = [T.age];
sexList = double(strcmp(T.sex,'F'));

indHC = strcmp(statusList,'HC');
indPM = strcmp(statusList,'PM');
indEM = strcmp(statusList,'EM') | strcmp(statusList,'Manifest');

groupList = nan(size(ageList));
groupList(indHC) = 1;
groupList(indPM) = 2;
groupList(indEM) = 3;

R = [ageList sexList indHC indPM indEM];
names = {'Age','Sex','HC', 'PM', 'EM'};

save([matout_root '/Data/Covariates.mat'], 'R', 'names');

%% read in ROI index -  QSM atlas

ROItxt = fileread('/working/lupolab/jingwen/090_Brain_ROIs_QSM_Atlas/ROIfromWord.txt');
ROItxt = strsplit(ROItxt, '\n');
ROItxt(cellfun('isempty', ROItxt)) = [];
ROItxt(strcmp(ROItxt, ' ')) = [];

ROIlist = cell(length(ROItxt),2);

for ii = 1:length(ROItxt)
    temp = strsplit(ROItxt{ii},' ');
    temp(cellfun('isempty', temp)) = [];
    if length(temp) > 1
        ROIlist{ii,1} = [temp{1:end-1}];
        ROIlist{ii,2} = str2double(temp{end});
    end
end

% keep the BG ones
ROIlist = ROIlist(177:190,:);
disp('QSM atlas:');
disp(ROIlist);

%% read in ROI index - FS

ROIexcel = '/working/lupolab/jingwen/002_HD_NDM/FSRegionList.xlsx';
T = readtable(ROIexcel);
ROIname = T.Var1;
ROIindex = T.Var2;

ROIlistFS = cell(length(ROIname),2);

for ii = 1:length(ROIname)
    ROIlistFS{ii,1} = ROIname{ii};
    ROIlistFS{ii,2} = ii;
end

% keep the BG ones
ROIlistFS = ROIlistFS([71 80 72 81 73 82],:);
disp('FreeSurfer:');
disp(ROIlistFS);

ROIlistFSL = ROIlistFS;
ROIlistFSL(:,2) = num2cell([11 50 12 51 13 52]');
disp('FSL:');
disp(ROIlistFSL);

%% Loop through subjects and clean up ROI files

for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    fprintf('%s is %s \n', subjPath, statusList{ii});
    
    exam_id = [subjList{ii} '_' examList{ii}];
    
    QSMfile_root = [subjPath '/swan_qsm/HDBET_allQSM/'];
    ROIfile_atlas = [QSMfile_root '/FSseg/QSM_atlas_ROI_noVent.nii.gz'];
    ROIfile_FS = [QSMfile_root '/FSseg/FSseg.nii.gz'];
    ROIfile_FSL = [QSMfile_root '/FSseg/firstseg.nii.gz'];
    
    ROIfile_atlasN = [QSMfile_root '/FSseg/Seg_atlas.nii.gz'];
    ROIfile_FSN = [QSMfile_root '/FSseg/Seg_FS.nii.gz'];
    ROIfile_FSLN = [QSMfile_root '/FSseg/Seg_FSL.nii.gz'];
    
    % QSM atlas
    nii = load_nii(ROIfile_atlas);
    QSM_ROI = zeros(size(nii.img));
    QSM_ROI(nii.img == 191) = 1;
    QSM_ROI(nii.img == 192) = 2;
    QSM_ROI(nii.img == 193) = 3;
    QSM_ROI(nii.img == 194) = 4;
    QSM_ROI(nii.img == 195) = 5;
    QSM_ROI(nii.img == 196) = 6;
    QSM_ROI(nii.img == 197 | nii.img == 199) = 7;
    QSM_ROI(nii.img == 198 | nii.img == 200) = 8;
    QSM_ROI(nii.img == 201) = 9;
    QSM_ROI(nii.img == 202) = 10;
    QSM_ROI(nii.img == 203) = 11;
    QSM_ROI(nii.img == 204) = 12;
    nii.img = QSM_ROI;
    save_nii(nii, ROIfile_atlasN);
    
    % FreeSurfer atlas
    nii = load_nii(ROIfile_FS);
    QSM_ROI = zeros(size(nii.img));
    QSM_ROI(nii.img == 71) = 1;
    QSM_ROI(nii.img == 80) = 2;
    QSM_ROI(nii.img == 72) = 3;
    QSM_ROI(nii.img == 81) = 4;
    QSM_ROI(nii.img == 73) = 5;
    QSM_ROI(nii.img == 82) = 6;
    nii.img = QSM_ROI;
    save_nii(nii, ROIfile_FSN);
    
    % FSL atlas
    nii = load_nii(ROIfile_FSL);
    QSM_ROI = zeros(size(nii.img));
    QSM_ROI(nii.img == 11) = 1;
    QSM_ROI(nii.img == 50) = 2;
    QSM_ROI(nii.img == 12) = 3;
    QSM_ROI(nii.img == 51) = 4;
    QSM_ROI(nii.img == 13) = 5;
    QSM_ROI(nii.img == 52) = 6;
    nii.img = QSM_ROI;
    save_nii(nii, ROIfile_FSLN);
    
end

%% Loop through subjects and concatenate registered QSM files

cmd_str = '3dTcat ';
for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    fprintf('%s is %s \n', subjPath, statusList{ii});
    
    QSMfile_root = [subjPath '/swan_qsm/HDBET_allQSM/MNIreg/'];
    QSMfile_path = [QSMfile_root 'QSM_iLSQR_meanEcho_MNInorm.nii.gz'];
    
    if exist(QSMfile_path,'file') > 0
        cmd_str = [cmd_str ' ' QSMfile_path];
    end
    
end
cmd_str = [cmd_str ' -prefix ' matout_root 'Data/allQSMnorm.nii.gz'];
system(cmd_str);

%% Loop through subjects and save registered & normalized QSM files
% ANTS registration individual files and normalize to lateral ventricles

VentMask_path = [matout_root '/Data/mask_lat_vent.nii.gz'];

cat_str = '3dTcat ';
for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    fprintf('%s is %s \n', subjPath, statusList{ii});
    
    QSMfile_root = [subjPath '/swan_qsm/HDBET_allQSM/'];
    QSMfile_reg_path = [QSMfile_root 'ANTSreg/QSM_iLSQR_MNI.nii.gz'];
    newfile_path = [matout_root '/ANTSreg_QSM_subj/' subjList{ii} '_' examList{ii} '_QSM_iLSQR_MNI.nii'];
    
    if exist(newfile_path, 'file') == 0
        cmd = sprintf('3dmaskave -q -mask %s %s', VentMask_path, QSMfile_reg_path);
        [~,out] = system(cmd);
        normVal = strsplit(out);
        normVal = normVal{end-1};
        fprintf(' - Ventrical QSM = %s\n', normVal);
        
        cmd = sprintf('3dcalc -a %s -expr a-%s -prefix %s', QSMfile_reg_path, normVal, newfile_path);
        [~,out_3dcalc] = system(cmd);
    end
    
    if exist(newfile_path,'file') > 0
        cat_str = [cat_str ' ' newfile_path];
    end
    
end
cat_str = [cat_str ' -prefix ' matout_root '/Data/allQSM_ANTS.nii'];
system(cat_str);

%% Loop through subjects and concatenate registered QSM files - ANTS registration

cmd_str = '3dTcat ';
for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    fprintf('%s is %s \n', subjPath, statusList{ii});
    
    QSMfile_root = [subjPath '/swan_qsm/HDBET_allQSM/'];
    QSMfile_path = [QSMfile_root 'FSseg/QSM_iLSQR_meanEcho_reg.nii.gz'];
    QSMfile_reg_path = [QSMfile_root 'ANTSreg/QSM_iLSQR_MNI.nii.gz'];
    MNIfile_path = [QSMfile_root 'ANTSreg/MNI_T1.nii.gz'];
    MNItrans1_path = [QSMfile_root 'ANTSreg/T1_MNI_1Warp.nii.gz'];
    MNItrans2_path = [QSMfile_root 'ANTSreg/T1_MNI_0GenericAffine.mat'];
    
    if exist(QSMfile_reg_path,'file') == 0
        cmd = sprintf('antsApplyTransforms -d 3 -i %s -r %s -t %s -t %s -o %s', ...
            QSMfile_path, MNIfile_path, MNItrans1_path, MNItrans2_path, QSMfile_reg_path);
        system(cmd);
    end
    
    if exist(QSMfile_reg_path,'file') > 0
        cmd_str = [cmd_str ' ' QSMfile_reg_path];
    end
    
end
cmd_str = [cmd_str ' -prefix ' matout_root '/Data/allQSM_ANTS.nii.gz'];
system(cmd_str);

%% MNI space

nii = load_nii([matout_root '/Data/Whole_brain_segmentation_DGM.nii']);
QSM_ROI = zeros(size(nii.img));
QSM_ROI(nii.img == 191) = 1;
QSM_ROI(nii.img == 192) = 2;
QSM_ROI(nii.img == 193) = 3;
QSM_ROI(nii.img == 194) = 4;
QSM_ROI(nii.img == 195) = 5;
QSM_ROI(nii.img == 196) = 6;
QSM_ROI(nii.img == 197 | nii.img == 199) = 7;
QSM_ROI(nii.img == 198 | nii.img == 200) = 8;
QSM_ROI(nii.img == 201) = 9;
QSM_ROI(nii.img == 202) = 10;
QSM_ROI(nii.img == 203) = 11;
QSM_ROI(nii.img == 204) = 12;
nii.img = QSM_ROI;
save_nii(nii, [matout_root '/Data/Seg_atlas.nii.gz']);

%% MNI space - new atlas segmentations from PD & cerebellum atlases
% http://nist.mni.mcgill.ca/multi-contrast-pd25-atlas/
% https://github.com/DiedrichsenLab/cerebellar_atlases/tree/master/Diedrichsen_2009

nii = load_nii(['/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/' ...
    'PD_subcortical.nii.gz']);
QSM_ROI = zeros(size(nii.img));
QSM_ROI(nii.img == 7) = 1; % L caudate
QSM_ROI(nii.img == 8) = 2; % R caudate
QSM_ROI(nii.img == 9) = 3; % L putamen
QSM_ROI(nii.img == 10) = 4; % R putamen
QSM_ROI(nii.img == 11) = 5; % L GPe
QSM_ROI(nii.img == 12) = 6; % R GPe
QSM_ROI(nii.img == 13) = 7; % L GPi
QSM_ROI(nii.img == 14) = 8; % R GPi
QSM_ROI(nii.img == 15) = 9; % L thalamus
QSM_ROI(nii.img == 16) = 10; % R thalamus
QSM_ROI(nii.img == 1) = 11; % L red nucleus
QSM_ROI(nii.img == 2) = 12; % R red nucleus
QSM_ROI(nii.img == 3) = 13; % L substantia nigra
QSM_ROI(nii.img == 4) = 14; % R substantia nigra
QSM_ROI(nii.img == 5) = 15; % L subthalamic nucleus
QSM_ROI(nii.img == 6) = 16; % R subthalamic nucleus

nii = load_nii(['/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/' ...
    'Cerebellum_atlas_MNI.nii.gz']);
QSM_ROI(nii.img == 29) = 17; % L dentate
QSM_ROI(nii.img == 30) = 18; % R dentate

nii.img = QSM_ROI;
save_nii(nii, ['/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data/' ...
    'Seg_Subcortical.nii.gz']);
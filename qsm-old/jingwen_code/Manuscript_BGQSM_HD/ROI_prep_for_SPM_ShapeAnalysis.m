clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/031_HD_NDM'));

FSresult_root = '/working/lupolab/jingwen/002_HD_NDM/AllT1N4';
matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/';

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220809.xlsx','Sheet','NoRep');
T(T.CAG < 36,:) = [];
T(strcmp(T.status_reclass,'MISSING'),:) = [];
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status_reclass];
ageList = [T.age];
sexList = double(strcmp(T.sex,'F'));

indHC = strcmp(statusList,'HC');
indPM = strcmp(statusList,'PM');
indEM = strcmp(statusList,'EM') | strcmp(statusList,'Manifest');

%% save covariates

R = [ageList sexList indHC indPM indEM];
names = {'Age','Sex','HC', 'PM', 'EM'};

save([matout_root '/SPM_results/Covariates.mat'], 'R', 'names');

R(38,:) = [];
save([matout_root '/SPM_results/Covariates_DTI.mat'], 'R', 'names');

R = [ageList sexList indHC ~indHC];
names = {'Age','Sex','HC', 'HD'};

save([matout_root '/SPM_results/Covariates_2group.mat'], 'R', 'names');

R(38,:) = [];
save([matout_root '/SPM_results/Covariates_DTI_2group.mat'], 'R', 'names');

%% Loop through subjects and copy ROI files

mkdir([matout_root '/ANTSreg_ROI_subj/']);
mkdir([matout_root '/ANTSreg_ROI_subj/HC']);
mkdir([matout_root '/ANTSreg_ROI_subj/PM']);
mkdir([matout_root '/ANTSreg_ROI_subj/EM']);
mkdir([matout_root '/ANTSreg_ROI_subj/Manifest']);
for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    fprintf('%s is %s \n', subjPath, statusList{ii});
    
    QSMfile_root = [subjPath '/swan_qsm/HDBET_allQSM/'];
    ROI_path = [QSMfile_root 'FSseg/Seg_ANTS.nii.gz'];
    SA_path = [matout_root '/ANTSreg_ROI_subj/' statusList{ii} '/' subjList{ii} '_' examList{ii} '_ROIs.nii.gz'];
    
    copyfile(ROI_path, SA_path);
    
end

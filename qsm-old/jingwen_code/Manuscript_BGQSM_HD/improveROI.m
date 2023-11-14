clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/031_HD_NDM'));

matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/';

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220809.xlsx','Sheet','NoRep');
T(strcmp(T.status_reclass,'MISSING'),:) = [];
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status_reclass];

%% Loop through subjects

for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii} '/swan_qsm/HDBET_allQSM/FSseg/'];
    
    if exist([subjPath 'Seg_ANTS_manual.nii.gz'],'file') == 0
        disp(['Redo ROI for ' subjPath]);
        cd(subjPath);
        system('bash /home/jyao3/030_QSM/01_Code/Manuscript_BGQSM_HD/DN_ROI.sh');
    end
    
end
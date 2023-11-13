clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220209.xlsx','Sheet','NoRep');
subjList = [T.b_num];
examList = [T.t_num];

%% tissue segmentation and fix QSMatlas ROIs

FS_root = '/working/lupolab/jingwen/002_HD_NDM/AllT1N4';

for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    Segfile_root = [subjPath '/swan_qsm/HDBET_allQSM/FSseg/'];
    
    exam_id = [subjList{ii} '_' examList{ii}];
    fprintf('Exam %s \n', exam_id);
    
    T1file = [FS_root '/' exam_id '_T1_N4.nii'];
    maskFile = [Segfile_root 'brain_mask_reg.nii.gz'];
    brainFile = [Segfile_root 'T1_brain.nii.gz'];
    fastFile = [Segfile_root 'T1_brain_seg_0.nii.gz'];
    
    if exist(brainFile,'file') ~= 2
        
        fprintf(' - Brain extraction \n');
        cmd = sprintf('3dcalc -a %s -b %s -expr a*b -prefix %s', ...
            T1file, maskFile, [Segfile_root '/T1_brain.nii.gz']);
        system(cmd);
    end
    
    if exist(fastFile,'file') ~= 2
        cmd = sprintf('fast -t 3 -g -o %s %s', [Segfile_root '/T1_brain'], ...
            [Segfile_root '/T1_brain.nii.gz']);
        system(cmd);
    end
    
    atlasFile = [Segfile_root 'QSM_atlas_ROI_reg.nii.gz'];
    atlasFixFile = [Segfile_root 'QSM_atlas_ROI_noVent.nii.gz'];
    
    if exist(atlasFixFile,'file') ~= 2
        cmd = sprintf('3dcalc -a %s -b %s -expr ''a*(1-b)'' -prefix %s', ...
            atlasFile, fastFile, atlasFixFile);
        system(cmd);
    end
        
end

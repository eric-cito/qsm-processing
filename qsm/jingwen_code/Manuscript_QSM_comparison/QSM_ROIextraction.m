clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','NoRep');
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status];
ageList = [T.age];

matout_root = '/working/lupolab/jingwen/001_QSM/04_QSM_Comparison_revision/';
mkdir(matout_root);

%% read in ROI index

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

%% Loop through subjects

for ii = 1:length(subjList)
    
    status = statusList{ii};
    fprintf('%s is %s \n', subjList{ii}, status);
    
    tic
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    
    % extract ROI values and save in .mat files
    exam_id = [subjList{ii} '_' examList{ii}];
    QSMfile_root = [subjPath '/swan_qsm/HDBET_allQSM/FSseg'];
    QSMfile_list = {'QSM_iLSQR_meanEcho_reg' ...
        'QSM_STARQSM_meanEcho_reg' ...
        'QSM_FANSI_nonlinearTV_meanEcho_reg' ...
        'QSM_HDQSM_meanEcho_reg' ...
        'QSM_MEDI_meanEcho_reg' ...
        'QSM_QSIP_meanEcho_reg' ...
        'QSM_SSTGV_meanEcho_reg' ...
        'QSM_SSTV_meanEcho_reg' ...
        'QSM_QSMGAN_meanEcho_reg' ...
        'QSM_QSMnet_meanEcho_reg'};
    ROIfile = [QSMfile_root '/QSM_atlas_ROI_noVent.nii.gz'];
    CSFfile = [QSMfile_root '/latven_mask_reg.nii.gz'];
    Maskfile = [QSMfile_root '/brain_mask_reg.nii.gz'];
    
    if exist([matout_root '/' exam_id '_erode_noVent.mat'],'file') == 0
        [QSMdata, QSMstats] = ...
            stats_QSM(QSMfile_root, QSMfile_list, ROIfile, CSFfile, ...
            Maskfile, ROIlist, 1);
        QSMstats(1).SubjName = subjList{ii};
        QSMstats(1).ExamName = examList{ii};
        save([matout_root '/' exam_id '_erode_noVent.mat'], 'QSMdata', 'QSMstats');
    end
    
    toc
end

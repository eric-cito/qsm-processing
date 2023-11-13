clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/031_HD_NDM'));

matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data';

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','NoRep');
T(T.CAG < 36,:) = [];
T(strcmp(T.status_reclass,'MISSING'),:) = [];
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status_reclass];
sexList = [T.sex];
ageList = [T.age];

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

%% read in ROI index - ANTS

ROIlist_ANTS = cell(18,1);
ROIlist_ANTS{1} = 'Caudate_L'; ROIlist_ANTS{2} = 'Caudate_R';
ROIlist_ANTS{3} = 'Putamen_L'; ROIlist_ANTS{4} = 'Putamen_R';
ROIlist_ANTS{5} = 'GPe_L'; ROIlist_ANTS{6} = 'GPe_R';
ROIlist_ANTS{7} = 'GPi_L'; ROIlist_ANTS{8} = 'GPi_R';
ROIlist_ANTS{9} = 'Thalamus_L'; ROIlist_ANTS{10} = 'Thalamus_R';
ROIlist_ANTS{11} = 'RN_L'; ROIlist_ANTS{12} = 'RN_R';
ROIlist_ANTS{13} = 'SN_L'; ROIlist_ANTS{14} = 'SN_R';
ROIlist_ANTS{15} = 'STN_L'; ROIlist_ANTS{16} = 'STN_R';
ROIlist_ANTS{17} = 'DN_L'; ROIlist_ANTS{18} = 'DN_R';

%% QSM file list

QSMfile_list = { ...
    'QSM_iLSQR' ...
    'QSM_STARQSM' ...
    'QSM_FANSI_nonlinearTV' ...
    'QSM_HDQSM' ...
    'QSM_MEDI' ...
    'QSM_QSIP' ...
    'QSM_SSTGV' ...
    'QSM_SSTV' ...
    'QSM_QSMGAN' ...
    'QSM_QSMnet' ...
    'QSM_xQSM2' ...
    'QSM_iQSM2'};

%% Loop through subjects

% load([matout_root '/SubjImgData_QSMcomp.mat'],'DataStruct');

% clear DataStruct;
% DataStruct(length(subjList)) = struct;

flag_erode = 0;

for ii = 1:length(subjList)
    
    fprintf('%s is %s \n', subjList{ii}, statusList{ii});
    
    DataStruct(ii).subjID = subjList{ii};
    DataStruct(ii).examID = examList{ii};
    DataStruct(ii).group = statusList{ii};
    DataStruct(ii).age = ageList(ii);
    DataStruct(ii).sex = sexList{ii};
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii} ...
        '/swan_qsm/HDBET_allQSM/'];
    
    CSF_path = [subjPath 'FSseg/latven_mask_reg.nii.gz'];
    nii = load_nii(CSF_path);
    VENmask = double(nii.img) > 0;
    
    brain_path = [subjPath 'FSseg/brain_mask_reg.nii.gz'];
    nii = load_nii(brain_path);
    brain_mask = double(nii.img) > 0;
    
    ROI_path = [subjPath 'FSseg/QSM_atlas_ROI_noVent.nii.gz'];
    nii = load_nii(ROI_path);
    QSM_ROI = double(nii.img);
    
    ROI_path = [subjPath 'FSseg/Seg_ANTS_manual.nii.gz'];
    nii = load_nii(ROI_path);
    ANTS_ROI = double(nii.img);
    
    for qq = 1:length(QSMfile_list)
        
        fprintf(' QSMfile %s \n', QSMfile_list{qq});
        
        QSM_path = [subjPath 'FSseg/' QSMfile_list{qq} '_meanEcho_reg.nii.gz'];
        nii = load_nii(QSM_path);
        QSMmap = double(nii.img);
        
        [QSMstats] = ...
            revision_QSMstat(QSMmap, QSM_ROI, VENmask, ...
            brain_mask, ROIlist, flag_erode);
        
        [Q] = BGstats(ROIlist_ANTS, ANTS_ROI, QSMmap, flag_erode);
        
        DataStruct(ii).(QSMfile_list{qq}) = [QSMstats; Q];
    end
    
end

%% save data

save([matout_root '/SubjImgData_QSMcomp_noerode.mat'],'DataStruct');

%% check 0 and nan and inf values

ROIvalue = zeros(length(subjList),1);
for ii = 1:length(subjList)
    
    fprintf('%s is %s \n', subjList{ii}, statusList{ii});
    
    ROIvalue(ii) = DataStruct(ii).QSM_iLSQR.ROIstd(212);
    
end

%% functions

function [Q] = BGstats(ROIlist, QSM_ROI, QSMmap, flag_erode, offset)

if nargin < 5
    offset = 0;
end

% median and mad QSM for BG ROIs
Q = table;
for rr = 1:length(ROIlist)
    ROImask = QSM_ROI == rr;
    [ROImedian, ROImad, ROImean, ROIstd] = ROIstats(ROImask, QSMmap, offset, flag_erode);
    ROIname = ROIlist(rr,1);
    T_roi = table(ROIname, ROImean, ROIstd, ROImedian, ROImad);
    Q = [Q; T_roi];
end

ROIname = {'GP_L'};
ROImask = (QSM_ROI == 5) | (QSM_ROI == 7);
[ROImedian, ROImad, ROImean, ROIstd] = ROIstats(ROImask, QSMmap, offset, flag_erode);
T_roi = table(ROIname,ROImean,ROIstd,ROImedian,ROImad);
Q = [Q; T_roi];

ROIname = {'GP_R'};
ROImask = (QSM_ROI == 6) | (QSM_ROI == 8);
[ROImedian, ROImad, ROImean, ROIstd] = ROIstats(ROImask, QSMmap, offset, flag_erode);
T_roi = table(ROIname,ROImean,ROIstd,ROImedian,ROImad);
Q = [Q; T_roi];

end

function [ROImedian, ROImad, ROImean, ROIstd] = ROIstats(ROImask, QSMmap, offset, flagErode)

se = strel('sphere',1);

if flagErode
    ROImask = imerode(ROImask,se);
end
ROIdata = QSMmap(ROImask);

ROIdata(ROIdata == 0) = [];
ROIdata(isnan(ROIdata)) = [];
ROIdata(isinf(ROIdata)) = [];

ROImean = mean(ROIdata-offset);
ROIstd = std(ROIdata-offset);
ROImedian = median(ROIdata-offset);
ROImad = mad(ROIdata-offset,1);

end

clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/031_HD_NDM'));

matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data';

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20230213.xlsx','Sheet','NoRep');
T(T.CAG < 36,:) = [];
T(strcmp(T.status_reclass2,'MISSING'),:) = [];
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status_reclass2];
sexList = [T.sex];

TMSList = [T.TMS];
DCLList = [T.DCL];
TFCList = [T.TFC];
ageList = [T.age];
cagList = [T.CAG];
DARTList = [T.DART_ToC];
FlankerList = [T.Flanker_Score];
MatchList = [T.Match_CorrectTotal];
SetShiftList = [T.SetShift_Shift_Score];
DARTTimingList = [T.DART_Total_AvgReaction];

% CAPS = 100*ageList.*(cagList-30)/627;
% AOO = CAPS*45/100;
% YTO = AOO - ageList;

%% read in ROI index -  QSM atlas

ROIlist = cell(18,1);
ROIlist{1} = 'Caudate_L'; ROIlist{2} = 'Caudate_R';
ROIlist{3} = 'Putamen_L'; ROIlist{4} = 'Putamen_R';
ROIlist{5} = 'GPe_L'; ROIlist{6} = 'GPe_R';
ROIlist{7} = 'GPi_L'; ROIlist{8} = 'GPi_R';
ROIlist{9} = 'Thalamus_L'; ROIlist{10} = 'Thalamus_R';
ROIlist{11} = 'RN_L'; ROIlist{12} = 'RN_R';
ROIlist{13} = 'SN_L'; ROIlist{14} = 'SN_R';
ROIlist{15} = 'STN_L'; ROIlist{16} = 'STN_R';
ROIlist{17} = 'DN_L'; ROIlist{18} = 'DN_R';

%% Loop through subjects

clear HD_BGanalysis;
HD_BGanalysis(length(subjList)) = struct;

flag_erode = 1;
FSresult_root = '/working/lupolab/jingwen/002_HD_NDM/AllT1N4';
figure('Position', [100 100 800 800]);

for ii = 1:length(subjList)
    
    fprintf('%s is %s \n', subjList{ii}, statusList{ii});
    
    HD_BGanalysis(ii).subjID = subjList{ii};
    HD_BGanalysis(ii).examID = examList{ii};
    HD_BGanalysis(ii).group = statusList{ii};
    HD_BGanalysis(ii).age = ageList(ii);
    HD_BGanalysis(ii).sex = sexList{ii};
    HD_BGanalysis(ii).CAG = cagList(ii);
    HD_BGanalysis(ii).TMS = TMSList(ii);
    HD_BGanalysis(ii).DCL = DCLList(ii);
    HD_BGanalysis(ii).TFC = TFCList(ii);
    HD_BGanalysis(ii).DART = DARTList(ii);
    HD_BGanalysis(ii).Flanker = FlankerList(ii);
    HD_BGanalysis(ii).Match = MatchList(ii);
    HD_BGanalysis(ii).SetShift = SetShiftList(ii);
    HD_BGanalysis(ii).DARTtiming = DARTTimingList(ii);
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii} ...
        '/swan_qsm/HDBET_allQSM/'];
    
    QSM_path = [subjPath 'FSseg/QSM_iQSM2_meanEcho_reg.nii.gz'];
    nii = load_nii(QSM_path);
    QSMsubj = double(nii.img);
    CSF_path = [subjPath 'FSseg/latven_mask_reg.nii.gz'];
    nii = load_nii(CSF_path);
    ROImask = double(nii.img) > 0;
    
    % median ventricle QSM
    [CSFmedian, ~] = ROIstats(ROImask, QSMsubj, 0, 0);
    fprintf('CSF median: %.4f\n', CSFmedian);
    
    ROI_path = [subjPath 'FSseg/Seg_ANTS_manual.nii.gz'];
    nii = load_nii(ROI_path);
    ROImask = double(nii.img);
    
    % QC plot
%     subplot(1,2,1); imagesc(squeeze(QSMsubj(:,120,:)), [-0.15 0.15]); axis off equal tight;
%     subplot(1,2,2); imagesc(squeeze(ROImask(:,120,:))); axis off equal tight;
%     drawnow;
    
    % median and mad QSM for BG ROIs
    [Q] = BGstats(ROIlist, ROImask, QSMsubj, CSFmedian, flag_erode);
    
    % FA
    DTIfile = [subjPath 'DTI/FA_regT1.nii.gz'];
    if exist(DTIfile,'file') > 0
        nii = load_nii(DTIfile); DTIsubj = double(nii.img); DTIsubj(DTIsubj == 0) = nan;
        [DTI_FA] = BGstats(ROIlist, ROImask, DTIsubj, 0, flag_erode);
    else
        DTI_FA = nan(size(Q));
    end
    
    % MD
    DTIfile = [subjPath 'DTI/MD_regT1.nii.gz'];
    if exist(DTIfile,'file') > 0
        nii = load_nii(DTIfile); DTIsubj = double(nii.img); DTIsubj(DTIsubj == 0) = nan;
        [DTI_MD] = BGstats(ROIlist, ROImask, DTIsubj, 0, flag_erode);
    else
        DTI_MD = nan(size(Q));
    end
    
    % RD
    DTIfile = [subjPath 'DTI/RD_regT1.nii.gz'];
    if exist(DTIfile,'file') > 0
        nii = load_nii(DTIfile); DTIsubj = double(nii.img); DTIsubj(DTIsubj == 0) = nan;
        [DTI_RD] = BGstats(ROIlist, ROImask, DTIsubj, 0, flag_erode);
    else
        DTI_RD = nan(size(Q));
    end
    
    % intracranial volume
    exam_id = [subjList{ii} '_' examList{ii}];
    if isfolder([FSresult_root '/FSresults/FSresult_' exam_id '_T1_N4/'])
        segPath = [FSresult_root '/FSresults/FSresult_' exam_id '_T1_N4/'];
        
    else
        segPath = [FSresult_root '/FSresult_' exam_id '_T1_N4/'];
    end
    [segname, segindex, segstats, eTIV] = load_segstats('aseg.stats',segPath);
    
    % volume
    voxvol = nii.hdr.dime.pixdim(2)*nii.hdr.dime.pixdim(3)*nii.hdr.dime.pixdim(4); % mm3
    V = zeros(1,length(ROIlist));
    for rr = 1:length(ROIlist)
        V(rr) = sum(ROImask == rr,'all')*voxvol;
    end
    V = V/eTIV*1500; % normalized by total intracranial volume
    
    SubjImgData = table(ROIlist, V', Q(1,:)', Q(2,:)', DTI_FA(1,:)', DTI_MD(1,:)', DTI_RD(1,:)', ...
        'VariableName', ...
        {'ROI' 'Volume' 'QSMmedian' 'QSMmad' 'FA' 'MD' 'RD'});
    
    HD_BGanalysis(ii).imData = SubjImgData;
    
end

%% save data

save([matout_root '/SubjImgDataBG_ANTS_iQSM_0213_erode.mat'],'HD_BGanalysis');

%% functions

function [Q] = BGstats(ROIlist, QSM_ROI, QSMmap, offset, flagErode)

if nargin < 5
    flagErode = 0;
end

if nargin < 4
    offset = 0;
end

% median and mad QSM for BG ROIs
Q = zeros(2,length(ROIlist));
for rr = 1:length(ROIlist)
    ROImask = QSM_ROI == rr;
    [ROImedian, ROImad] = ROIstats(ROImask, QSMmap, offset, flagErode);
    Q(:,rr) = [ROImedian; ROImad];
end

end

function [ROImedian, ROImad] = ROIstats(ROImask, QSMmap, offset, flagErode)

se = strel('sphere',1);

if flagErode
    ROImask = imerode(ROImask,se);
end
ROIdata = QSMmap(ROImask);

ROIdata(ROIdata == 0) = [];
ROIdata(isnan(ROIdata)) = [];
ROIdata(isinf(ROIdata)) = [];

ROImedian = median(ROIdata-offset);
ROImad = mad(ROIdata-offset,1);

end

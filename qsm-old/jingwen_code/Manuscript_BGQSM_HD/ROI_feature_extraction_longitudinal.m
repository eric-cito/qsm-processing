clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/031_HD_NDM'));

matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data';

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20230213.xlsx','Sheet','Rep');

%% extract demo info

subjList = [T.b_num];
examList = [T.t_num];
DOBunique = unique(T.dob);

clear HDsubj_long;
HDsubj_long(length(DOBunique)) = struct;
for ii = 1:length(DOBunique)
    indSubj = T.dob == DOBunique(ii);
    
    HDsubj_long(ii).bnum(1) = subjList(indSubj & T.scan_num == 1);
    HDsubj_long(ii).tnum(1) = examList(indSubj & T.scan_num == 1);
    HDsubj_long(ii).scandate(1) = T.scan_date(indSubj & T.scan_num == 1);
    HDsubj_long(ii).scanday(1) = 0;
    
    HDsubj_long(ii).age(1) = T.age(indSubj & T.scan_num == 1);
    
    if sum(indSubj & T.scan_num == 2) > 0
        HDsubj_long(ii).bnum(2) = subjList(indSubj & T.scan_num == 2);
        HDsubj_long(ii).tnum(2) = examList(indSubj & T.scan_num == 2);
        HDsubj_long(ii).scandate(2) = T.scan_date(indSubj & T.scan_num == 2);
        HDsubj_long(ii).scanday(2) = days(HDsubj_long(ii).scandate(2) - HDsubj_long(ii).scandate(1));
        HDsubj_long(ii).age(2) = T.age(indSubj & T.scan_num == 2);
    else
        HDsubj_long(ii).bnum(2) = {''};
        HDsubj_long(ii).tnum(2) = {''};
        HDsubj_long(ii).scandate(2) = {''};
        HDsubj_long(ii).scanday(2) = nan;
        HDsubj_long(ii).age(2) = nan;
    end
    
    if sum(indSubj & T.scan_num == 3) > 0
        HDsubj_long(ii).bnum(3) = subjList(indSubj & T.scan_num == 3);
        HDsubj_long(ii).tnum(3) = examList(indSubj & T.scan_num == 3);
        HDsubj_long(ii).scandate(3) = T.scan_date(indSubj & T.scan_num == 3);
        HDsubj_long(ii).scanday(3) = days(HDsubj_long(ii).scandate(3) - HDsubj_long(ii).scandate(1));
        HDsubj_long(ii).age(3) = T.age(indSubj & T.scan_num == 3);
    else
        HDsubj_long(ii).bnum(3) = {''};
        HDsubj_long(ii).tnum(3) = {''};
        HDsubj_long(ii).scandate(3) = {''};
        HDsubj_long(ii).scanday(3) = nan;
        HDsubj_long(ii).age(3) = nan;
    end
    
    HDsubj_long(ii).status = T.status_reclass(indSubj & T.scan_num == 1);
    HDsubj_long(ii).sex = T.sex(indSubj & T.scan_num == 1);
    HDsubj_long(ii).cag = T.cag_2(indSubj & T.scan_num == 1);
end

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

%% Loop through subjects and extract ROI data

flagErode = 1;

FSresult_root = '/working/lupolab/jingwen/002_HD_NDM/AllT1N4';

for ii = 1:length(DOBunique)
    
    indSubj = T.dob == DOBunique(ii);
    
    tp1_path = T.QSM_path{indSubj & T.scan_num == 1};
    fprintf('TP1: %s\n', tp1_path);
    
    % intracranial volume
    exam_id = [subjList{indSubj & T.scan_num == 1} '_' examList{indSubj & T.scan_num == 1}];
    if isfolder([FSresult_root '/FSresults/FSresult_' exam_id '_T1_N4/'])
        segPath = [FSresult_root '/FSresults/FSresult_' exam_id '_T1_N4/'];
    else
        segPath = [FSresult_root '/FSresult_' exam_id '_T1_N4/'];
    end
    [segname, segindex, segstats, eTIV] = load_segstats('aseg.stats',segPath);
    
    [Q1, V1] = extractQSM(tp1_path, ROIlist, flagErode);
    V1 = V1/eTIV*1500; % normalized by total intracranial volume
    
    if sum(indSubj & T.scan_num == 2) > 0
        tp2_path = T.QSM_path{indSubj & T.scan_num == 2};
        fprintf('TP2: %s\n', tp2_path);
        [Q2, V2] = extractQSM(tp2_path, ROIlist, flagErode);
    else
        Q2 = nan(size(Q1));
        V2 = nan(size(V1));
    end
    V2 = V2/eTIV*1500; % normalized by total intracranial volume
    
    if sum(indSubj & T.scan_num == 3) > 0
        tp3_path = T.QSM_path{indSubj & T.scan_num == 3};
        fprintf('TP3: %s\n', tp3_path);
        [Q3, V3] = extractQSM(tp3_path, ROIlist, flagErode);
    else
        Q3 = nan(size(Q1));
        V3 = nan(size(V1));
    end
    V3 = V3/eTIV*1500; % normalized by total intracranial volume
    
    SubjImgData = table([ROIlist; {'ST_L'}; {'ST_R'}], V1', Q1(1,:)', Q1(2,:)', ...
        V2', Q2(1,:)', Q2(2,:)', V3', Q3(1,:)', Q3(2,:)', ...
        'VariableName', ...
        {'ROI' 'VolumeTP1' 'QSMmedianTP1' 'QSMmadTP1' ...
        'VolumeTP2' 'QSMmedianTP2' 'QSMmadTP2' ...
        'VolumeTP3' 'QSMmedianTP3' 'QSMmadTP3'});
    
    HDsubj_long(ii).imData = SubjImgData;
    
end

%% save data

save([matout_root '/SubjImgDataBG_Longitudinal_1228_iQSM_erode.mat'],'HDsubj_long');

%% functions

function [Q, V] = extractQSM(tp1_path, ROIlist, flagErode)

QSM_path = [tp1_path 'FSseg/QSM_iQSM2_meanEcho_reg.nii.gz'];
nii = load_nii(QSM_path);
QSMsubj = double(nii.img);
CSF_path = [tp1_path 'FSseg/latven_mask_reg.nii.gz'];
nii = load_nii(CSF_path);
CSFmask = double(nii.img) > 0;

% median ventricle QSM
[CSFmedian, ~] = ROIstats(CSFmask, QSMsubj, 0, 0);
fprintf('CSF median: %.4f\n', CSFmedian);

ROI_path = [tp1_path 'FSseg/Seg_ANTS_manual.nii.gz'];
nii = load_nii(ROI_path);
ROImask = double(nii.img);

% QC plot
subplot(1,2,1); imagesc(squeeze(QSMsubj(:,120,:)), [-0.15 0.15]); axis off equal tight;
subplot(1,2,2); imagesc(squeeze(ROImask(:,120,:)) + squeeze(CSFmask(:,120,:))); axis off equal tight;
drawnow;

% median and mad QSM for BG ROIs
[Q] = BGstats(ROIlist, ROImask, QSMsubj, CSFmedian, flagErode);

% volume
voxvol = nii.hdr.dime.pixdim(2)*nii.hdr.dime.pixdim(3)*nii.hdr.dime.pixdim(4); % mm3
V = zeros(1,length(ROIlist));
for rr = 1:length(ROIlist)
    V(rr) = sum(ROImask == rr,'all')*voxvol;
end
V(rr+1) = sum(ROImask == 1 | ROImask == 3,'all')*voxvol;
V(rr+2) = sum(ROImask == 2 | ROImask == 4,'all')*voxvol;

end

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

ROImask = QSM_ROI == 1 | QSM_ROI == 3;
[ROImedian, ROImad] = ROIstats(ROImask, QSMmap, offset, flagErode);
Q(:,length(ROIlist)+1) = [ROImedian; ROImad];

ROImask = QSM_ROI == 2 | QSM_ROI == 4;
[ROImedian, ROImad] = ROIstats(ROImask, QSMmap, offset, flagErode);
Q(:,length(ROIlist)+2) = [ROImedian; ROImad];

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

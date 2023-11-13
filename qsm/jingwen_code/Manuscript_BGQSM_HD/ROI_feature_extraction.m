clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/031_HD_NDM'));

FSresult_root = '/working/lupolab/jingwen/002_HD_NDM/AllT1N4';
matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD';

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220809.xlsx','Sheet','NoRep');
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status_reclass];
sexList = [T.sex];

TMSList = [T.TMS];
DCLList = [T.DCL];
TFCList = [T.TFC];
ageList = [T.age];
cagList = [T.CAG];
CAPS = 100*ageList.*(cagList-30)/627;
AOO = CAPS*45/100;
YTO = AOO - ageList;

% QSMfile_list = {'QSM_FANSI_nonlinearTV_meanEcho' ...
%     'QSM_HDQSM_meanEcho' ...
%     'QSM_iLSQR_meanEcho' ...
%     'QSM_MEDI_meanEcho' ...
%     'QSM_QSIP_meanEcho' ...
%     'QSM_QSMGAN_meanEcho' ...
%     'QSM_QSMnet_meanEcho' ...
%     'QSM_SSTGV_meanEcho' ...
%     'QSM_SSTV_meanEcho' ...
%     'QSM_STARQSM_meanEcho'};

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
disp(ROIlistFS);

%% Loop through subjects

clear HD_BGanalysis V;
HD_BGanalysis(length(subjList)) = struct;

DTIflag = 1;

for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    fprintf('%s is %s \n', subjPath, statusList{ii});
    
    HD_BGanalysis(ii).subjID = subjList{ii};
    HD_BGanalysis(ii).examID = examList{ii};
    HD_BGanalysis(ii).group = statusList{ii};
    HD_BGanalysis(ii).age = ageList(ii);
    HD_BGanalysis(ii).sex = sexList{ii};
    HD_BGanalysis(ii).CAG = cagList(ii);
    HD_BGanalysis(ii).CAPS = CAPS(ii);
    HD_BGanalysis(ii).AOO = AOO(ii);
    HD_BGanalysis(ii).YTO = YTO(ii);
    HD_BGanalysis(ii).TMS = TMSList(ii);
    HD_BGanalysis(ii).DCL = DCLList(ii);
    HD_BGanalysis(ii).TFC = TFCList(ii);
    
    exam_id = [subjList{ii} '_' examList{ii}];
    
    QSMfile_root = [subjPath '/swan_qsm/HDBET_allQSM/'];
    ROIfile_atlas = [QSMfile_root '/FSseg/QSM_atlas_ROI_noVent.nii.gz'];
    ROIfile_FS = [QSMfile_root '/FSseg/FSseg.nii.gz'];
    CSFfile = [QSMfile_root '/FSseg/latven_mask_reg.nii.gz'];
    QSMfile = [QSMfile_root '/FSseg/QSM_iLSQR_meanEcho_reg.nii.gz'];
    
    % load the maps
    nii = load_nii(CSFfile);
    QSM_LatVen = double(nii.img);
    nii = load_nii(QSMfile);
    QSMmap = double(nii.img);
    
    % load the ROI files
    nii = load_nii(ROIfile_atlas);
    QSM_ROI_atlas = double(nii.img);
    QSM_ROI_atlas(QSM_LatVen > 0) = 0;
    nii = load_nii(ROIfile_FS);
    QSM_ROI_FS = double(nii.img);
    QSM_ROI_FS(QSM_LatVen > 0) = 0;
    
    % intracranial volume
    segPath = [FSresult_root '/FSresults/FSresult_' exam_id '_T1_N4/'];
    [segname, segindex, segstats, eTIV] = load_segstats('aseg.stats',segPath);
    
    % QSM Atlas BG volume
    nii = load_nii(ROIfile_atlas);
    voxvol = nii.hdr.dime.pixdim(2)*nii.hdr.dime.pixdim(3)*nii.hdr.dime.pixdim(4); % mm3
    V = zeros(1,length(ROIlist));
    for rr = 1:length(ROIlist)
        V(rr) = sum(QSM_ROI_atlas == ROIlist{rr,2},'all')*voxvol;
    end
    
    ROIname = [ROIlist(:,1); {'SubstantiaNigra_L'}];
    ROImask = (QSM_ROI_atlas == 197) | (QSM_ROI_atlas == 199);
    V(end+1) = sum(ROImask,'all');
    ROIname = [ROIname; {'SubstantiaNigra_R'}];
    ROImask = (QSM_ROI_atlas == 198) | (QSM_ROI_atlas == 200);
    V(end+1) = sum(ROImask,'all');
    
    % FreeSurfer BG volume
    nii = load_nii(ROIfile_FS);
    voxvol = nii.hdr.dime.pixdim(2)*nii.hdr.dime.pixdim(3)*nii.hdr.dime.pixdim(4); % mm3
    tmp = length(V);
    for rr = 1:length(ROIlistFS)
        V(tmp+rr) = sum(QSM_ROI_FS == ROIlistFS{rr,2},'all')*voxvol;
    end
    
    V = V/eTIV*1500; % normalized by total intracranial volume
    
    % median ventricle QSM
    ROImask = QSM_LatVen > 0;
    [CSFmedian, ~] = ROIstats(ROImask, QSMmap, 0, 0);
    
    % median and mad QSM for BG ROIs
    [Q] = BGstats(ROIlist, QSM_ROI_atlas, ROIlistFS, QSM_ROI_FS, QSMmap, CSFmedian);
    sizeQ = size(Q);
    
    if DTIflag
        % DTI values for BG ROIs
        DTIfile = [QSMfile_root '/DTI/FA_regT1.nii.gz'];
        [DTI_FA] = DTIstats(ROIlist, QSM_ROI_atlas, ROIlistFS, QSM_ROI_FS, DTIfile, sizeQ);
        
        DTIfile = [QSMfile_root '/DTI/MD_regT1.nii.gz'];
        [DTI_MD] = DTIstats(ROIlist, QSM_ROI_atlas, ROIlistFS, QSM_ROI_FS, DTIfile, sizeQ);
        
        DTIfile = [QSMfile_root '/DTI/RD_regT1.nii.gz'];
        [DTI_RD] = DTIstats(ROIlist, QSM_ROI_atlas, ROIlistFS, QSM_ROI_FS, DTIfile, sizeQ);
        
        DTIfile = [QSMfile_root '/DTI/ODI_regT1.nii.gz'];
        [DTI_ODI] = DTIstats(ROIlist, QSM_ROI_atlas, ROIlistFS, QSM_ROI_FS, DTIfile, sizeQ);
        
        ROInameList = [ROIname; ROIlistFS(:,1)];
        SubjImgData = table(ROInameList, V', Q(1,:)', Q(2,:)', DTI_FA(1,:)', DTI_FA(2,:)', ...
            DTI_MD(1,:)', DTI_MD(2,:)', DTI_RD(1,:)', DTI_RD(2,:)', DTI_ODI(1,:)', DTI_ODI(2,:)', ...
            'VariableName', ...
            {'ROI' 'NormVolume' 'QSMmedian' 'QSMmad' 'FAmedian' 'FAmad' ...
            'MDmedian' 'MDmad' 'RDmedian' 'RDmad' 'ODImedian' 'ODImad'});
    else
        ROInameList = [ROIname; ROIlistFS(:,1)];
        SubjImgData = table(ROInameList, V', Q(1,:)', Q(2,:)', ...
            'VariableName', ...
            {'ROI' 'NormVolume' 'QSMmedian' 'QSMmad'});
    end
    
    HD_BGanalysis(ii).imData = SubjImgData;
    
end

%% save data

save([matout_root '/SubjImgDataBG_iLSQR.mat'],'HD_BGanalysis');

%% local QSM stats

load([matout_root '/SubjImgDataBG_SSTGV.mat']);

% age regression in HC
qsmMat = zeros(length(HD_BGanalysis),size(HD_BGanalysis(1).imData,1));
for ii = 1:length(HD_BGanalysis)
    qsmMat(ii,:) = HD_BGanalysis(ii).imData.QSMmedian;
end

age = [HD_BGanalysis.age];
indHC = strcmp({HD_BGanalysis.group},'HC');

qsmThresh = zeros(size(qsmMat));
for ii = 1:size(qsmMat,2)
    mdl = fitlm(age(indHC)',qsmMat(indHC,ii));
    [ypred,yci] = predict(mdl,age');
    qsmThresh(:,ii) = yci(:,2);
end

for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    fprintf('%s is %s \n', subjPath, statusList{ii});
    
    exam_id = [subjList{ii} '_' examList{ii}];
    
    QSMfile_root = [subjPath '/swan_qsm/HDBET_allQSM/'];
    ROIfile_atlas = [QSMfile_root '/FSseg/QSM_atlas_ROI_noVent.nii.gz'];
    ROIfile_FS = [QSMfile_root '/FSseg/FSseg.nii.gz'];
    CSFfile = [QSMfile_root '/FSseg/latven_mask_reg.nii.gz'];
    QSMfile = [QSMfile_root '/FSseg/QSM_iLSQR_meanEcho_reg.nii.gz'];
    
    % load the maps
    nii = load_nii(CSFfile);
    QSM_LatVen = double(nii.img);
    nii = load_nii(QSMfile);
    QSMmap = double(nii.img);
    
    % load the ROI files
    nii = load_nii(ROIfile_atlas);
    QSM_ROI_atlas = double(nii.img);
    QSM_ROI_atlas(QSM_LatVen > 0) = 0;
    nii = load_nii(ROIfile_FS);
    QSM_ROI_FS = double(nii.img);
    QSM_ROI_FS(QSM_LatVen > 0) = 0;
    
    % median ventricle QSM
    ROImask = QSM_LatVen > 0;
    [CSFmedian, ~] = ROIstats(ROImask, QSMmap, 0, 1);
    
    % median and mad QSM for BG ROIs
    [Q] = BGstats_II(ROIlist, QSM_ROI_atlas, ROIlistFS, QSM_ROI_FS, QSMmap, CSFmedian, qsmThresh(ii,:));
    sizeQ = size(Q);
    
    SubjImgDataII = table(Q(1,:)', Q(2,:)', ...
        'VariableName', ...
        {'QSMIImedian' 'QSMIImad'});
    
    HD_BGanalysis(ii).imData = [HD_BGanalysis(ii).imData SubjImgDataII];
    
end

%% save data

save([matout_root '/SubjImgDataBG_iLSQR.mat'],'HD_BGanalysis');

%% functions

function [Q] = DTIstats(ROIlist, QSM_ROI, ROIlistFS, QSM_ROI_FS, DTIfile, sizeQ)

if exist(DTIfile,'file') == 2
    nii = load_nii(DTIfile); DTImap = double(nii.img);
    [Q] = BGstats(ROIlist, QSM_ROI, ROIlistFS, QSM_ROI_FS, DTImap, 0);
else
    Q = nan(sizeQ);
end

end

function [Q] = BGstats_II(ROIlist, QSM_ROI, ROIlistFS, QSM_ROI_FS, QSMmap, CSFmedian, qsmThresh)

% median and mad QSM for BG ROIs - high QSM regions
Q = zeros(2,length(ROIlist));
for rr = 1:length(ROIlist)
    ROImask = QSM_ROI == ROIlist{rr,2};
    ROImask = ROImask & (QSMmap > qsmThresh(rr)+CSFmedian);
    [ROImedian, ROImad] = ROIstats(ROImask, QSMmap, CSFmedian, 0);
    Q(:,rr) = [ROImedian; ROImad];
end

ROImask = (QSM_ROI == 197) | (QSM_ROI == 199);
ROImask = ROImask & (QSMmap > qsmThresh(length(ROIlist)+1)+CSFmedian);
[ROImedian, ROImad] = ROIstats(ROImask, QSMmap, CSFmedian, 0);
Q(:,end+1) = [ROImedian; ROImad];

ROImask = (QSM_ROI == 198) | (QSM_ROI == 200);
ROImask = ROImask & (QSMmap > qsmThresh(length(ROIlist)+2)+CSFmedian);
[ROImedian, ROImad] = ROIstats(ROImask, QSMmap, CSFmedian, 0);
Q(:,end+1) = [ROImedian; ROImad];

for rr = 1:length(ROIlistFS)
    ROImask = QSM_ROI_FS == ROIlistFS{rr,2};
    ROImask = ROImask & (QSMmap > qsmThresh(length(ROIlist)+2+rr)+CSFmedian);
    [ROImedian, ROImad] = ROIstats(ROImask, QSMmap, CSFmedian, 0);
    Q(:,length(ROIlist) + 2 + rr) = [ROImedian; ROImad];
end

end

function [Q] = BGstats(ROIlist, QSM_ROI, ROIlistFS, QSM_ROI_FS, QSMmap, CSFmedian)

% median and mad QSM for BG ROIs
Q = zeros(2,length(ROIlist));
for rr = 1:length(ROIlist)
    ROImask = QSM_ROI == ROIlist{rr,2};
    [ROImedian, ROImad] = ROIstats(ROImask, QSMmap, CSFmedian, 1);
    Q(:,rr) = [ROImedian; ROImad];
end

ROImask = (QSM_ROI == 197) | (QSM_ROI == 199);
[ROImedian, ROImad] = ROIstats(ROImask, QSMmap, CSFmedian, 1);
Q(:,end+1) = [ROImedian; ROImad];

ROImask = (QSM_ROI == 198) | (QSM_ROI == 200);
[ROImedian, ROImad] = ROIstats(ROImask, QSMmap, CSFmedian, 1);
Q(:,end+1) = [ROImedian; ROImad];

for rr = 1:length(ROIlistFS)
    ROImask = QSM_ROI_FS == ROIlistFS{rr,2};
    [ROImedian, ROImad] = ROIstats(ROImask, QSMmap, CSFmedian, 1);
    Q(:,length(ROIlist) + 2 + rr) = [ROImedian; ROImad];
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

ROImedian = mean(ROIdata-offset);
ROImad = std(ROIdata-offset,1);

end

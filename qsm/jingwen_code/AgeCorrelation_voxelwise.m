clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath('/home/jyao3/010_MATLAB_Utils/');
addpath('/home/jyao3/010_MATLAB_Utils/NIfTI');
addpath('/home/jyao3/010_MATLAB_Utils/export_fig');

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220103.xlsx','Sheet','All');
statusList = [T.status];

subjList = T.b_num; 
examList = T.t_num; 
sexList = T.sex; 
ageList = T.age; 

% subjList(~strcmp(statusList,'HC')) = [];
% examList(~strcmp(statusList,'HC')) = [];
% sexList(~strcmp(statusList,'HC')) = [];
% ageList(~strcmp(statusList,'HC')) = [];

Mind = strcmp(sexList,'M');
Find = strcmp(sexList,'F');

fprintf('Age mean %.3f SD %.3f \n', mean(ageList), std(ageList));
fprintf('M %i F %i \n', sum(Mind), sum(Find));

%% Specify QSM names

QSMfile_list = {'QSM_FANSI_nonlinearTV_meanEcho' ...
    'QSM_HDQSM_meanEcho' ...
    'QSM_iLSQR_meanEcho' ...
    'QSM_MEDI_meanEcho' ...
    'QSM_QSIP_meanEcho' ...
    'QSM_QSMGAN_meanEcho' ...
    'QSM_QSMnet_meanEcho' ...
    'QSM_SSTGV_meanEcho' ...
    'QSM_SSTV_meanEcho' ...
    'QSM_STARQSM_meanEcho'};
QSMatlas = '/working/lupolab/QSM_Atlas_MNI_toolbox/QSM_atlas_blackbackground_MNI_1mm_v1.nii';

savePath = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/NIFTI_HC';

%% loop through all subjects

for qq = 1:length(QSMfile_list)
    
    disp(['Processing ' QSMfile_list{qq}]);
    
%     if exist([savePath '/' QSMfile_list{qq} '_MNIall.nii.gz'],'file') == 2; continue; end
    
    for ii = 1:length(subjList)
        examPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
        QSMfile_root = [examPath '/swan_qsm/HDBET_allQSM/'];
        
        disp([' - Subject ' subjList{ii}]);
        
        % register to MNI space
        QSMfile_MNI = [QSMfile_root '/MNIreg/' QSMfile_list{qq} '_MNI.nii.gz'];
        if exist(QSMfile_MNI,'file') ~= 2
            QSMfile = [QSMfile_root '/' QSMfile_list{qq} '.nii.gz'];
            QSMfile_orient = [QSMfile_root '/MNIreg' QSMfile_list{qq} '_orient.nii.gz'];
            
            cmd = sprintf('3dresample -input %s -prefix %s -orient LPI',...
                QSMfile, QSMfile_orient);
            system(cmd);
        
            cmdstr = ['applywarp --ref=' QSMatlas ...
                ' --in=' QSMfile_orient ...
                ' --warp=' QSMfile_root '/MNIreg/QSM_regNonlin_warp.nii.gz' ...
                ' --out=' QSMfile_MNI ];
            system(cmdstr);
        end
        
        % normalize to lateral ventricle median
        QSMfile_MNInorm = [QSMfile_root '/MNIreg/' QSMfile_list{qq} '_MNInorm.nii.gz'];
        if exist(QSMfile_MNInorm,'file') ~= 2
            CSFfile = [QSMfile_root '/QSMseg_latven.nii.gz'];
            [RefValue] = calcLatVenMedian(CSFfile, QSMfile);
        
            cmdstr = ['3dcalc -a ' QSMfile_MNI ...
                ' -expr a-' num2str(RefValue) ...
                ' -prefix ' QSMfile_MNInorm ];
            system(cmdstr);
        end
        
        % load QSM image
%         nii = load_nii(QSMfile_MNInorm);
%         QSMmap = double(nii.img);
%         QSMmap_all(:,:,:,ii) = QSMmap;
        
    end
    
    % save as nifti
%     disp([' - Saving QSMmaps ' QSMfile_list{qq}]);
%     niftiwrite(QSMmap_all, [savePath '/' QSMfile_list{qq} '_MNIall.nii'], 'Compressed',true);
%     clear QSMmap_all
end

%% load all QSM and voxel-wise correlation

for qq = 1:length(QSMfile_list)
    
    if 0 % exist([savePath '/' QSMfile_list{qq} '_R.nii.gz'],'file') == 2
        continue;
    end
    
    % load nifti
    nii = load_nii([savePath '/' QSMfile_list{qq} '_MNIall.nii.gz']);
    QSMmap_all = double(nii.img);
    
    % load brain mask
    nii = load_nii([savePath '/QSM_atlas_mask.nii.gz']);
    mask = double(nii.img);
    
    % vectorize things
    Nsubj = length(ageList);
    sizeMat = size(mask);
    QSMmap_vec = reshape(QSMmap_all,[],Nsubj);
    
    % index of valid voxels
    ind = find(mask > 0);
    QSMmap_vec = QSMmap_vec(ind,:);
    
    % loop through all voxels for fitting
    Rvec = zeros(1,size(QSMmap_vec,1));
    Pvec = zeros(1,size(QSMmap_vec,1));
    Svec = zeros(1,size(QSMmap_vec,1));
    Ivec = zeros(1,size(QSMmap_vec,1));
    parfor vv = 1:size(QSMmap_vec,1)
        
        if mod(vv,10000) == 0
            fprintf('processing voxel %i of %i: %.1f p.u. \n',...
                vv,size(QSMmap_vec,1),vv/size(QSMmap_vec,1)*100);
        end
        
        data = QSMmap_vec(vv,:)';
        validInd = ~isnan(data) & (data ~= 0);
        [R,P] = corr(ageList(validInd),data(validInd),'Type','Pearson');
        Rvec(vv) = R;
        Pvec(vv) = P;
        F = fit(ageList(validInd),data(validInd),'poly1');
        Svec(vv) = F.p1;
        Ivec(vv) = F.p2;
    end
    
    % back to 3D image
    Rimg = zeros(sizeMat); Rimg(ind) = Rvec;
    Pimg = zeros(sizeMat); Pimg(ind) = Pvec;
    Simg = zeros(sizeMat); Simg(ind) = Svec;
    Iimg = zeros(sizeMat); Iimg(ind) = Ivec;
    
    % save as nifti
    nii = load_nii([savePath '/QSM_dummy.nii.gz']);
    nii.img = Rimg; save_nii(nii,[savePath '/' QSMfile_list{qq} '_R.nii.gz']);
    nii.img = Pimg; save_nii(nii,[savePath '/' QSMfile_list{qq} '_P.nii.gz']);
    nii.img = Simg; save_nii(nii,[savePath '/' QSMfile_list{qq} '_Slope.nii.gz']);
    nii.img = Iimg; save_nii(nii,[savePath '/' QSMfile_list{qq} '_Intercept.nii.gz']);
 
end

%% functions

function [RefValue] = calcLatVenMedian(CSFfile, QSMfile)

nii = load_nii(QSMfile);
QSMmap = double(nii.img);

nii = load_nii(CSFfile);
QSM_LatVen = double(nii.img);

ROIdata = nonzeros(QSMmap(QSM_LatVen > 0));
ROIdata(isnan(ROIdata)) = [];

RefValue = median(ROIdata);

end
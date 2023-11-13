clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/031_HD_NDM'));

matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data';

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','NoRep');
statusList = nan(1,length(T.status_reclass));
statusList(strcmp(T.status_reclass,'HC')) = 1;
statusList(strcmp(T.status_reclass,'PM')) = 2;
statusList(strcmp(T.status_reclass,'EM')) = 3;

subjList = T.b_num; subjList(isnan(statusList)) = [];
examList = T.t_num; examList(isnan(statusList)) = [];
statusList(isnan(statusList)) = [];

indHC = statusList == 1;
indPM = statusList == 2; 
indEM = statusList == 3;

%% Specify QSM names

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
QSMatlas = '/working/lupolab/QSM_Atlas_MNI_toolbox/QSM_atlas_blackbackground_MNI_1mm_v1.nii';

matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data';

%% loop through all subjects

for qq = [11 12] % 1:length(QSMfile_list)
    
    disp(['Processing ' QSMfile_list{qq}]);
    
    for ii = 1:length(subjList)
        examPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
        QSMfile_root = [examPath '/swan_qsm/HDBET_allQSM/'];
        
        disp([' - Subject ' subjList{ii}]);
        
        % clean up
        % cmd = sprintf('mv %s/MNIregQSM* %s/MNIreg/', QSMfile_root, QSMfile_root);
        % system(cmd);
        
        % register to MNI space
        QSMfile = [QSMfile_root '/' QSMfile_list{qq} '_meanEcho.nii.gz'];
        QSMfile_MNI = [QSMfile_root '/MNIreg/' QSMfile_list{qq} '_meanEcho_MNI.nii.gz'];
        if exist(QSMfile_MNI,'file') ~= 2
            QSMfile_orient = [QSMfile_root '/MNIreg/MNIreg' QSMfile_list{qq} '_meanEcho_orient.nii.gz'];
            
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
        QSMfile_MNInorm = [QSMfile_root '/MNIreg/' QSMfile_list{qq} '_meanEcho_MNInorm.nii.gz'];
        % %%delete(QSMfile_MNInorm);
        if exist(QSMfile_MNInorm,'file') ~= 2
            CSFfile = [QSMfile_root '/QSMseg_latven.nii.gz'];
            [RefValue] = calcLatVenMedian([QSMfile_root '/FSseg/latven_mask_reg.nii.gz'], ...
                [QSMfile_root '/FSseg/' QSMfile_list{qq} '_meanEcho_reg.nii.gz']);
        
            cmdstr = ['3dcalc -a ' QSMfile_MNI ...
                ' -expr a-' num2str(RefValue) ...
                ' -prefix ' QSMfile_MNInorm ];
            system(cmdstr);
        end
        
        % load QSM image
        nii = load_nii(QSMfile_MNInorm);
        QSMmap = double(nii.img);
        QSMmap_all(:,:,:,ii) = QSMmap;
        
    end
    
    QSMmap_HC = mean(QSMmap_all(:,:,:,indHC),4);
    QSMmap_PM = mean(QSMmap_all(:,:,:,indPM),4);
    QSMmap_EM = mean(QSMmap_all(:,:,:,indEM),4);
    
    % save as nifti
    disp([' - Saving QSMmaps ' QSMfile_list{qq}]);
    nii.img = QSMmap_HC;
    save_nii(nii,[matout_root '/' QSMfile_list{qq} '_MNI_HCmean.nii.gz']);
    nii.img = QSMmap_PM;
    save_nii(nii,[matout_root '/' QSMfile_list{qq} '_MNI_PMmean.nii.gz']);
    nii.img = QSMmap_EM;
    save_nii(nii,[matout_root '/' QSMfile_list{qq} '_MNI_EMmean.nii.gz']);
    clear QSMmap_all
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


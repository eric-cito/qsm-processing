clear; clc;
warning('off');

%% Add path

addpath('/home/jyao3/030_QSM/01_Code/Manuscript_BGQSM_HD');
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

%% read in list
T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','NoRep');
subjList = [T.b_num];
examList = [T.t_num];

%% Loop through subjects

for ii = [37 45 50 62:length(subjList)]
    
    examPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    
    if exist([examPath '/swan_qsm/HDBET_allQSM/Phase.nii.gz'], 'file') == 2
        
        % QSM register to QSM atlas
        try
            exam_id = [subjList{ii} '_' examList{ii}];
            regfile_root = [examPath '/swan_qsm/HDBET_allQSM/MNIreg'];
            mkdir(regfile_root);
            
            if exist([regfile_root '/QSM_atlas_ROI.nii.gz'],'file') ~= 2
                fprintf(' - Registrating exam %s \n', exam_id);
                reg_QSM([examPath '/swan_qsm/HDBET_allQSM'], regfile_root);
            end
        catch
            disp('!!! Something went wrong!');
        end
        
%         delete([output_root '/QSM_atlas_ROI.nii.gz']);
%         delete([output_root '/QSM_atlas_ROI_reg.nii.gz']);
%         delete([output_root '/QSM_atlas_ROI_noVent.nii.gz']);
        
        % ROI register to T1
        output_root = [examPath '/swan_qsm/HDBET_allQSM/FSseg/'];
        exam_id = [subjList{ii} '_' examList{ii}];
        T1file = [output_root '/' exam_id '_T1_N4.nii'];
        if exist([output_root '/QSM_atlas_ROI.nii.gz'], 'file') == 0
            disp(['Register atlas ROI ' subjList{ii} ' ' examList{ii}]);
            Maskfile = [regfile_root '/QSM_atlas_ROI.nii.gz'];
            if exist([output_root '/QSM_atlas_ROI.nii.gz'],'file') == 0
                cmdstr = ['3dresample' ...
                    ' -input ' Maskfile ...
                    ' -master ' output_root '/../QSM_iLSQR_meanEcho.nii.gz' ...
                    ' -prefix ' output_root '/QSM_atlas_ROI.nii.gz'];
                system(cmdstr);
            end
        end
        if exist([output_root '/QSM_atlas_ROI_reg.nii.gz'], 'file') == 0
            Maskfile = [output_root '/QSM_atlas_ROI.nii.gz'];
            if exist([output_root '/QSM_atlas_ROI_reg.nii.gz'],'file') == 0
                cmdstr = ['/netopt/rhel7/fsl/bin/flirt' ...
                    ' -in ' Maskfile ...
                    ' -ref ' T1file ...
                    ' -out ' output_root '/QSM_atlas_ROI_reg.nii.gz' ...
                    ' -applyxfm -init ' output_root '/QSM2T1F.mat -interp nearestneighbour'];
                system(cmdstr);
            end
        end
        
        % fix atlas ROIs by removing CSF
        segfile_root = [output_root '/FAST/'];
        atlasFile = [output_root 'QSM_atlas_ROI_reg.nii.gz'];
        atlasFixFile = [output_root 'QSM_atlas_ROI_noVent.nii.gz'];
        fastFile = [segfile_root 'T1_brain_seg_0.nii.gz'];
        if exist(atlasFixFile,'file') ~= 2
            disp(['Fix atlas ROI ' subjList{ii} ' ' examList{ii}]);
            cmd = sprintf('3dcalc -a %s -b %s -expr ''a*(1-b)'' -prefix %s', ...
                atlasFile, fastFile, atlasFixFile);
            system(cmd);
        end
        
    end
    
end

%% subject list from COSMOS

mainPath = '/working/lupolab/eason/DL_QSM/data/7T_cosmos';
subjFolders = dir(mainPath);
subjFolders(~[subjFolders.isdir]) = [];
subjFolders(~cellfun(@isempty,strfind({subjFolders.name},'.'))) = [];
subjFolders(cellfun(@isempty,strfind({subjFolders.name},'volunteer'))) = [];

%% Loop through subjects

for ii = 1:length(subjFolders)
    
    subjPath = [mainPath '/' subjFolders(ii).name '/scan1'];
    exam_root = [subjPath '/swan_qsm/COSMOSmask_allQSM'];
    
    if exist([exam_root '/QSM_COSMOS.nii.gz'], 'file') == 2
        
        try
            exam_id = [subjList{ii} '_' examList{ii}];
            regfile_root = [exam_root '/MNIreg'];
            mkdir(regfile_root);
            
            if exist([regfile_root '/QSM_atlas_ROI.nii.gz'],'file') ~= 2
                fprintf(' - Registrating exam %s \n', exam_id);
                reg_QSM(exam_root, regfile_root);
            end
        catch
            disp('!!! Something went wrong!');
        end
        
    end
    
end

%% for failed ones

% copy T1 and MNI-T1 to MNIreg folder
% cp /working/lupolab/QSM_Atlas_MNI_toolbox/MNI152_T1_1mm_brain.nii Atlas_T1.nii.gz
% resample T1 to QSM_orient
% 3dresample -input T1.nii.gz -master QSM_orient.nii.gz -prefix T1_qsm.nii.gz
% resample brain_mask_HD to LPI as QSM_orient
% 3dresample -input ../brain_mask_HD.nii.gz -prefix brain_mask.nii.gz -orient LPI
% apply brain mask to T1 -> T1_brain
% 3dcalc -a brain_mask.nii.gz -b T1_qsm.nii.gz -expr a*b -prefix T1_brain.nii.gz

% register T1
% flirt -in T1_brain.nii.gz -ref Atlas_T1.nii.gz -out T1_LinReg.nii.gz
% -omat T1_LinReg.mat -cost corratio -dof 12 -interp trilinear
% fnirt -v --ref=Atlas_T1.nii.gz --in=T1_brain.nii.gz --aff=T1_LinReg.mat
% --intmod=global_linear --iout=T1_regNonlin.nii.gz --fout=T1_regNonlinDF.nii.gz
% --cout=T1_regNonlin_warp.nii.gz --jacrange=-1,100

% apply to QSM atlas
% invwarp -v --ref=T1_brain.nii.gz --warp=T1_regNonlin_warp.nii.gz
% --out=T1_regNonlin_warp_inv.nii.gz
% applywarp --ref=T1_brain.nii.gz --in=QSM_atlas.nii.gz
% --warp=T1_regNonlin_warp_inv.nii.gz --out=QSM_atlas_ROI.nii.gz
% --interp=nn
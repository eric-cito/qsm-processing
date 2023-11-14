clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/031_HD_NDM'));

FSresult_root = '/working/lupolab/jingwen/002_HD_NDM/AllT1N4';
matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/';

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20221129.xlsx','Sheet','NoRep');
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status_reclass];

%% Loop through subjects and recon DTI files

addpath(genpath('/home/mmorrison1/Documents/recon_mb/'));
addpath('/netopt/share/lib/local/brain/matlab/');

for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    
    exam_id = [subjList{ii} '_' examList{ii}];
    % fprintf('Exam %s \n', exam_id);
    
    [DTIpath] = findDTIfolder(subjPath);
    
    if isempty(DTIpath)
        fprintf('Exam %s DTI not done \n', exam_id);
        % in MATLAB raw folder
        % change pfiles name to EXXXX_multiband_dti_noddi & _.ref.dat & _.vrgf.dat
        % copy one E.DCM file to raw folder
        cd(['/data/7T_hunt/' subjList{ii} '/orig_raw/dti_raw']);
        
        dtiFiles = dir(pwd);
        dtiFiles(~contains({dtiFiles.name}, '_multiband_dti_noddi_')) = [];
        if ~isempty(dtiFiles)
            dtiBase = strsplit(dtiFiles(1).name,'_multiband_dti_noddi_');
            dtiBase = dtiBase{1};
            fprintf(' - DTI base : %s \n', dtiBase);
        end
        
        % recon
        if isfolder(['/data/7T_hunt/' subjList{ii} '/orig_raw/dti_raw/dicom']) == 0
            fprintf(' - Recon DTI \n');
            recon_multiband_all_parfor(dtiBase, 'E.DCM')
        end
        
        % dti processing
        if exist(['/data/7T_hunt/' subjList{ii} '/orig_raw/dti_raw/processed_dti/' dtiBase '_TOPUP_FA.nii.gz'],'file') == 0
            fprintf(' - Process DTI \n');
            cmd = sprintf('/working/lupolab/jingwen/process_MB_DTI_rad_debug dicom %s processed_dti', dtiBase);
            system(cmd);
            cmd = sprintf('mv processed_dti/data_topup_ec_* processed_dti/grid_io_eddy/');
            system(cmd);
            cmd = sprintf('mv processed_dti/*.out processed_dti/grid_io_eddy/');
            system(cmd);
        end
    end
    
end

%% Loop through subjects and copy DTI files

dsi_root = '/working/lupolab/jingwen/004_HD_DSI/data';

for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    
    [DTIpath] = findDTIfolder(subjPath);
    fprintf(' - DTI path: %s \n', DTIpath);
    
    dtiFiles = dir(DTIpath);
    dtiFiles(~contains({dtiFiles.name}, 'TOPUP_FA.nii.gz')) = [];
    if ~isempty(dtiFiles)
        dtiBase = strsplit(dtiFiles.name,'_FA.nii.gz');
        dtiBase = dtiBase{1};
        fprintf(' - DTI base name : %s \n', dtiBase);
    end
    
    processingPath = [dsi_root '/' subjList{ii} '/' examList{ii}];
    if exist([processingPath '/' dtiBase '_V1.nii.gz'],'file') == 0
        mkdir(processingPath);
        cmd_str = sprintf('cp %s/*_TOPUP_* %s', ...
            DTIpath, processingPath);
        system(cmd_str);
    end
    
end

%% Loop through subjects and perform DTIKIT

for ii = 1:length(subjList)
    
    subjPath = [dsi_root '/' subjList{ii} '/' examList{ii}];
    dtiFiles = dir(subjPath);
    dtiFiles(~contains({dtiFiles.name}, 'TOPUP_FA.nii.gz')) = [];
    if ~isempty(dtiFiles)
        dtiBase = strsplit(dtiFiles.name,'_FA.nii.gz');
        dtiBase = dtiBase{1};
        fprintf(' - DTI base : %s \n', dtiBase);
    end
    
    if exist([subjPath '/' dtiBase '_dtitk.nii.gz'],'file') == 0 && ~isempty(dtiFiles)
        cmd_str = sprintf(['export PATH=/netopt/rhel7/versions/dtitk/current/bin:${PATH}; ' ...
            'export PATH=/netopt/rhel7/versions/dtitk/current/utilities:${PATH}; ' ...
            '/netopt/rhel7/versions/dtitk/current/scripts/fsl_to_dtitk %s'], ...
            [subjPath '/' dtiBase]);
        system(cmd_str);
    end
    
    % copy HC to a new folder
    if strcmp(statusList{ii}, 'HC') && exist([subjPath '/' dtiBase '_dtitk.nii.gz'],'file') ...
            && exist([dsi_root '/../HC_DTI_template/' dtiBase '_dtitk.nii.gz'],'file') == 0
        copyfile([subjPath '/' dtiBase '_dtitk.nii.gz'], [dsi_root '/../HC_DTI_template/']);
    end
    
end

%% make DTI template from HC data

% create subject list txt

% cd([dsi_root '/../HC_DTI_template/']);
% cmd_str = sprintf(['export DTITK_ROOT=/netopt/rhel7/versions/dtitk/current; ' ...
%     'export PATH=/netopt/rhel7/versions/dtitk/current/bin:${PATH}; ' ...
%     'export PATH=/netopt/rhel7/versions/dtitk/current/utilities:${PATH}; ' ...
%     'export PATH=/netopt/rhel7/versions/dtitk/current/scripts:${PATH}; ' ...
%     'dti_template_bootstrap %s %s'], ...
%     '/home/jyao3/013_Templates/ixi_aging_template_v3.0/template/ixi_aging_template.nii.gz', ...
%     'subj.txt');
% system(cmd_str);
% 
% cmd_str = sprintf(['export DTITK_ROOT=/netopt/rhel7/versions/dtitk/current; ' ...
%     'export PATH=/netopt/rhel7/versions/dtitk/current/bin:${PATH}; ' ...
%     'export PATH=/netopt/rhel7/versions/dtitk/current/utilities:${PATH}; ' ...
%     'export PATH=/netopt/rhel7/versions/dtitk/current/scripts:${PATH}; ' ...
%     'dti_affine_population mean_initial.nii.gz subj.txt EDS 3; ' ...
%     'TVtool -in mean_affine3.nii.gz -tr; ' ...
%     'BinaryThresholdImageFilter mean_affine3_tr.nii.gz mask.nii.gz 0.01 100 1 0; ']);
% system(cmd_str);
% 
% cmd_str = sprintf(['export DTITK_ROOT=/netopt/rhel7/versions/dtitk/current; ' ...
%     'export PATH=/netopt/rhel7/versions/dtitk/current/bin:${PATH}; ' ...
%     'export PATH=/netopt/rhel7/versions/dtitk/current/utilities:${PATH}; ' ...
%     'export PATH=/netopt/rhel7/versions/dtitk/current/scripts:${PATH}; ' ...
%     'dti_diffeomorphic_population mean_affine3.nii.gz subj_aff.txt mask.nii.gz 0.002; ' ...
%     'dti_warp_to_template_group subj.txt mean_diffeomorphic_initial6.nii.gz 1.0 1.0 1.0']);
% system(cmd_str);

%% warp individual subjects to tensor

templatePath = [dsi_root '/../HC_DTI_template/mean_diffeomorphic_initial6.nii.gz'];

for ii = 1:length(subjList)
    
    exam_id = [subjList{ii} '_' examList{ii}];
    fprintf('Exam %s is %s \n', exam_id, statusList{ii});
    
    subjPath = [dsi_root '/' subjList{ii} '/' examList{ii}];
    dtiFiles = dir(subjPath);
    dtiFiles(~contains({dtiFiles.name}, 'TOPUP_FA.nii.gz')) = [];
    if ~isempty(dtiFiles)
        dtiBase = strsplit(dtiFiles.name,'_TOPUP_FA.nii.gz');
        dtiBase = dtiBase{1};
        fprintf(' - DTI base : %s \n', dtiBase);
    end
    
    % rigid registration to template
    if exist([subjPath '/' dtiBase '_TOPUP_dtitk_aff_diffeo.nii.gz'],'file') == 0
        if strcmp(statusList{ii}, 'HC') ...
                && exist([dsi_root '/../HC_DTI_template/' dtiBase '_TOPUP_dtitk_aff_diffeo.nii.gz'],'file')
            copyfile([dsi_root '/../HC_DTI_template/' dtiBase '_TOPUP_dtitk_aff_diffeo.nii.gz'], ...
                [subjPath '/' dtiBase '_TOPUP_dtitk_aff_diffeo.nii.gz']);
        else
            cmd_str = sprintf(['export DTITK_ROOT=/netopt/rhel7/versions/dtitk/current; ' ...
                'export PATH=/netopt/rhel7/versions/dtitk/current/bin:${PATH}; ' ...
                'export PATH=/netopt/rhel7/versions/dtitk/current/utilities:${PATH}; ' ...
                'export PATH=/netopt/rhel7/versions/dtitk/current/scripts:${PATH}; ' ...
                'dti_rigid_reg %s %s_TOPUP_dtitk.nii.gz EDS 4 4 4 0.01; '...
                'TVtool -in %s -out %s/final_temp_tr.nii.gz -tr; ' ...
                'BinaryThresholdImageFilter %s/final_temp_tr.nii.gz %s/mask.nii.gz 0.01 100 1 0; ' ...
                'dti_diffeomorphic_reg %s %s_TOPUP_dtitk_aff.nii.gz %s/mask.nii.gz 1 6 0.002'], ...
                templatePath, [subjPath '/' dtiBase], ...
                templatePath, subjPath, ...
                subjPath, subjPath, ...
                templatePath, [subjPath '/' dtiBase], subjPath);
            system(cmd_str);
        end
    end
    
end

%% Copy DTI to QSM folder /swan_qsm/HDBET_allQSM/DTI/

for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    
    [DTIpath] = findDTIfolder(subjPath);
    fprintf(' - DTI path: %s \n', DTIpath);
    
    dtiFiles = dir(DTIpath);
    dtiFiles(~contains({dtiFiles.name}, 'TOPUP_FA.nii.gz')) = [];
    if ~isempty(dtiFiles)
        dtiBase = strsplit(dtiFiles.name,'_TOPUP_FA.nii.gz');
        dtiBase = dtiBase{1};
        fprintf(' - DTI base : %s \n', dtiBase);
    end
    
    newPath = [subjPath '/swan_qsm/HDBET_allQSM/DTI/'];
    if exist([newPath '/RD.nii.gz'],'file') == 0
        mkdir(newPath);
        copyfile([DTIpath '/' dtiBase '_TOPUP_FA.nii.gz'], [newPath 'FA.nii.gz']);
        copyfile([DTIpath '/' dtiBase '_TOPUP_MD.nii.gz'], [newPath 'MD.nii.gz']);
        % copyfile([DTIpath '/' dtiBase '_TOPUP_RD.nii.gz'], [newPath 'RD.nii.gz']);
    end
    
end

%% Loop through subjects and register DTI to T1

for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    
    exam_id = [subjList{ii} '_' examList{ii}];
    T1file = [subjPath '/swan_qsm/HDBET_allQSM/FSseg/' exam_id '_T1_N4.nii'];
    DTIfile_root = [subjPath '/swan_qsm/HDBET_allQSM/DTI/'];
    
    FA = [DTIfile_root '/FA.nii.gz'];
    
    if exist([DTIfile_root '/RD_regT1.nii.gz'],'file') ~= 2 && exist(FA,'file') == 2
        fprintf(' - Registrating exam %s \n', exam_id);
        
        % register FA to T1
        copyfile(T1file, [DTIfile_root '/T1.nii.gz']);
        b0file = [DTIfile_root '/FA.nii.gz'];
        b0reg = [DTIfile_root '/FA_regT1.nii.gz'];
        b0regmat = [DTIfile_root '/b0_regT1.mat'];
        cmd = sprintf('flirt -in %s -ref %s -out %s -omat %s -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90',...
            b0file, T1file, b0reg, b0regmat);
        system(cmd);
        cmdstr = ['bash /home/jyao3/030_QSM/01_Code/convertMat.sh'...
            ' ' DTIfile_root '/b0_regT1.mat'...
            ' > ' DTIfile_root '/b0_regT1F.mat'];
        system(cmdstr);
        b0regmat = [DTIfile_root '/b0_regT1F.mat'];
        
        % register MD to T1
        MOVreg = [DTIfile_root '/MD_regT1.nii.gz'];
        MOV = [DTIfile_root '/MD.nii.gz'];
        cmd = sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s',...
            MOV, T1file, MOVreg, b0regmat);
        system(cmd);
        % register RD to T1
%         MOVreg = [DTIfile_root '/RD_regT1.nii.gz'];
%         MOV = [DTIfile_root '/RD.nii.gz'];
%         cmd = sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s',...
%             MOV, T1file, MOVreg, b0regmat);
%         system(cmd);
        
    end
    
end

%% helper function

function [DTIpath] = findDTIfolder(subjPath)

% find b0 and TOPUP_FA files
b0file = [subjPath '/NODDI/processed_dti/b0.nii.gz'];
FAorig = dir([subjPath '/NODDI/processed_dti/']);
FAorig(~contains({FAorig.name},'TOPUP_FA.nii.gz')) = [];
if exist(b0file, 'file') ~= 2
    b0file = [subjPath '/NODDI/b0.nii.gz'];
    FAorig = dir([subjPath '/NODDI/']);
    FAorig(~contains({FAorig.name},'TOPUP_FA.nii.gz')) = [];
    if exist(b0file, 'file') ~= 2
        b0file = [subjPath '/NODDI_melanie/b0.nii.gz'];
        FAorig = dir([subjPath '/NODDI_melanie/']);
        FAorig(~contains({FAorig.name},'TOPUP_FA.nii.gz')) = [];
        if exist(b0file, 'file') ~= 2
            b0file = [subjPath '/dti/processed_dti/b0.nii.gz'];
            FAorig = dir([subjPath '/dti/processed_dti/']);
            FAorig(~contains({FAorig.name},'TOPUP_FA.nii.gz')) = [];
            if exist(b0file, 'file') ~= 2
                b0file = [subjPath '/../NODDI/processed_dti/b0.nii.gz'];
                FAorig = dir([subjPath '/../NODDI/processed_dti/']);
                FAorig(~contains({FAorig.name},'TOPUP_FA.nii.gz')) = [];
                if exist(b0file, 'file') ~= 2
                    b0file = [subjPath '/../orig_raw/dti_raw/processed_dti/b0.nii.gz'];
                    FAorig = dir([subjPath '/../orig_raw/dti_raw/processed_dti']);
                    FAorig(~contains({FAorig.name},'TOPUP_FA.nii.gz')) = [];
                end
            end
        end
    end
end

if isempty(FAorig)
    DTIpath = '';
else
    DTIpath = FAorig.folder;
end

end

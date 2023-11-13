clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

%% read in list
T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220103.xlsx');
subjList = [T.b_num];
statusList = [T.status];

%% subject list from 7T_hunt

subjFolders = dir('/data/7T_hunt');
subjFolders(~[subjFolders.isdir]) = [];
subjFolders(~cellfun(@isempty,strfind({subjFolders.name},'.'))) = [];

%% Loop through subjects

latex = {};

for ii = 1:length(subjFolders)
    
    ind = strcmp(subjList, subjFolders(ii).name);
    ind = find(ind);
    
    if isempty(ind)
        continue
    end
    
    status = statusList{ind(1)};
    
    fprintf('%s is %s \n', subjFolders(ii).name, status);
    
%     if ~strcmp(status,'HC')
%         continue
%     end
    
    subjPath = ['/data/7T_hunt/' subjFolders(ii).name];
    examFolders = dir(subjPath);
    examFolders(~[examFolders.isdir]) = [];
    examFolders(~cellfun(@isempty,strfind({examFolders.name},'.'))) = [];
    
    for ee = 1:length(examFolders)
        examPath = [subjPath '/' examFolders(ee).name];
        if exist([examPath '/swan_qsm/HDBET_allQSM/Phase.nii.gz'], 'file') == 2
            
            try
                
                exam_id = [subjFolders(ii).name '_' examFolders(ee).name];
                regfile_root = [examPath '/swan_qsm/HDBET_allQSM/MNIreg'];
                mkdir(regfile_root);
                
                fprintf(' - Registrating exam %s \n', exam_id);
                
                if exist([regfile_root '/QSM_atlas_ROI.nii.gz'],'file') ~= 2
                    reg_QSM([examPath '/swan_qsm/HDBET_allQSM'], regfile_root);
                end
                
            catch
                disp('!!! Something went wrong!');
            end
            
            % save images for QC
            QCfolder = [regfile_root '/QC_images'];
            mkdir(QCfolder);
            
            QSMname = [regfile_root '/QSM_shift.nii.gz'];
            ROIname = [regfile_root '/QSM_atlas_ROI.nii.gz'];
            
            ATLASname = [regfile_root '/QSM_atlas.nii'];
            REGname = [regfile_root '/QSM_regNonlin.nii.gz'];
            
            if exist([QCfolder '/QSMatlas_subj.png'],'file') ~= 2 ...
                    && exist(REGname,'file') 
                targetSize = [300 300 148];
                reg_QSM_plotQC(QSMname, ROIname, targetSize, 0);
                export_fig([QCfolder '/QSMsubj_ROI'], '-png'); close;
                reg_QSM_plotQC(ATLASname, REGname, targetSize, 0);
                export_fig([QCfolder '/QSMatlas_subj'], '-png'); close;
            end
            latex_cell = ScreenShotLatexCell(exam_id, QCfolder);
            latex(end+1:end+length(latex_cell)) = latex_cell;
        end
        
    end
    
end

%% create QC report

path2reportfolder = '/home/jyao3/030_QSM/01_Code/QCreport_template';

% save as txt
fid = fopen([path2reportfolder '/ScreenShots.txt'], 'wt');
fprintf(fid, '%s\n',char(latex)');
fclose(fid);

% create report
cd(path2reportfolder);
cmd = 'buf_size=2000000 pdflatex QCreport.tex;';
system(cmd);

% move to working dir
save_root = '/working/lupolab/jingwen/001_QSM/';
copyfile('QCreport.pdf',[save_root 'QCreport_QSMreg.pdf']);

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
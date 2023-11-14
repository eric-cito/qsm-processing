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

%% process T1

for ii = 1:length(subjList)
    
    examPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    
    if exist([examPath '/swan_qsm/HDBET_allQSM/FSseg/T1_brain.nii.gz'], 'file') ~= 0
        disp([subjList{ii} ' ' examList{ii} ' already processed']);
    else
        disp(['Processing ' subjList{ii} ' ' examList{ii}]);
        
        % create T1 nifti
        if exist([examPath '/swan_qsm/HDBET_allQSM/FSseg/' subjList{ii} '_' examList{ii} '_T1_N4.nii'], 'file') == 0
            disp(' -- Processing T1N4');
            [T1file] = findT1N4(examPath);
            % system(['rm -r ' examPath '/FSseg/']);
            mkdir([examPath '/swan_qsm/HDBET_allQSM/FSseg/']);
            if ~isempty(T1file)
                copyfile(T1file, [examPath '/swan_qsm/HDBET_allQSM/FSseg/' subjList{ii} '_' examList{ii} '_T1_N4.nii']);
            else
                % create N4 T1 if not exist
                disp('  -- Processing T1');
                findT1(examPath);
            end
        end
        
        % brain extraction
        T1_path = [examPath '/swan_qsm/HDBET_allQSM/FSseg/' subjList{ii} '_' examList{ii} '_T1_N4.nii'];
        T1brain_path = [examPath '/swan_qsm/HDBET_allQSM/FSseg/T1_brain.nii.gz'];
        mask_path = [examPath '/swan_qsm/HDBET_allQSM/FSseg/brain_mask_reg.nii.gz'];
        if exist(T1brain_path,'file') == 0
            disp(' -- Brain extraction');
            cmd_str = sprintf('3dcalc -a %s -b %s -expr a*b -prefix %s', ...
                T1_path, mask_path, T1brain_path);
            system(cmd_str);
        end
    end
end

%% FAST segmentation & lateral ventricle

for ii = 1:length(subjList)
    
    examPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii} '/swan_qsm/HDBET_allQSM/FSseg/'];
    segfile_root = [examPath '/FAST/'];
    
    % delete([examPath '/latven_mask_reg.nii.gz']);
    if exist([examPath '/latven_mask_reg.nii.gz'], 'file') ~= 0
        disp([subjList{ii} ' ' examList{ii} ' already processed']);
    else
        disp(['Processing ' subjList{ii} ' ' examList{ii}]);
        
        if exist([examPath '/T1_brain_seg_0.nii.gz'],'file')
            mkdir(segfile_root);
            system(sprintf('mv %s/T1_brain_* %s', examPath, segfile_root));
        end
        
        if exist([segfile_root '/T1_brain_seg_0.nii.gz'],'file') ~= 2
            % tissue segmentation
            mkdir(segfile_root);
            cmd = sprintf('fast -t 3 -g -o %s -v %s', ...
                [segfile_root '/T1_brain'], ...
                [examPath '/T1_brain.nii.gz']);
            system(cmd);
        end
        
        if exist([segfile_root '/T1_brain_latven.nii.gz'],'file') == 0
            % erosion of brain mask and crop the lower half
            cmd = sprintf('3dmask_tool -input %s -prefix %s -dilate_result -10', ...
                [examPath '/brain_mask_reg.nii.gz'], ...
                [segfile_root '/brain_mask_erosion.nii.gz']);
            system(cmd);
            cmd = sprintf('3dcalc -a %s -expr ''step(k-70)*step(110-k)'' -prefix %s', ...
                [segfile_root '/brain_mask_erosion.nii.gz'], ...
                [segfile_root '/zcrop.nii.gz']);
            system(cmd);
            % erosion of CSF mask and crop
            cmd = sprintf('3dmask_tool -input %s -prefix %s -dilate_result -1', ...
                [segfile_root '/T1_brain_seg_0.nii.gz'], ...
                [segfile_root '/T1_brain_csf_erosion.nii.gz']);
            system(cmd);
            cmd = sprintf('3dcalc -a %s -b %s -c %s -expr a*b*c -prefix %s', ...
                [segfile_root '/T1_brain_csf_erosion.nii.gz'], ...
                [segfile_root '/brain_mask_erosion.nii.gz'], ...
                [segfile_root '/zcrop.nii.gz'], ...
                [segfile_root '/T1_brain_csf_mask.nii.gz']);
            system(cmd);
            % calculate cluster and threshold cluster size for final mask
            cmd = sprintf('cluster -i %s -t 1 --osize=%s', ...
                [segfile_root '/T1_brain_csf_mask.nii.gz'], ...
                [segfile_root '/T1_brain_csf_cluster.nii.gz']);
            system(cmd);
            cmd = sprintf('3dcalc -a %s -expr ''step(a-900)'' -prefix %s', ...
                [segfile_root '/T1_brain_csf_cluster.nii.gz'], ...
                [segfile_root '/T1_brain_latven.nii.gz']);
            system(cmd);
        end
        % copy to final folder
        copyfile([segfile_root '/T1_brain_latven.nii.gz'], [examPath '/latven_mask_reg.nii.gz']);
        
    end
    
    % plot for QC
    if exist([segfile_root '/latven.png'],'file') ~= 2 ...
            && exist([segfile_root '/T1_brain_latven.nii.gz'],'file')
        nii = load_nii([segfile_root '/T1_brain_latven.nii.gz']);
        LatVen = double(nii.img);
        figure;
        isosurface(LatVen,0.5);
        export_fig([segfile_root '/latven'], '-png'); % close;
    end
end

%% FreeSurfer segmentation

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','NoRep');
subjList = [T.b_num];
examList = [T.t_num];

FSresult_root = '/working/lupolab/jingwen/002_HD_NDM/AllT1N4';

cmd_str = ['export FREESURFER_HOME=/netopt/freesurfer-5.3.0; ' ...
    'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
    'export SUBJECTS_DIR=/working/lupolab/jingwen/002_HD_NDM/AllT1N4; ' ...
    'cd $SUBJECTS_DIR; '];
for ii = 1:length(subjList)
    
    examPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    exam_id = [subjList{ii} '_' examList{ii}];
    T1file = [examPath '/swan_qsm/HDBET_allQSM/FSseg/' exam_id '_T1_N4.nii'];
    
    % copy to FS path
    if exist([FSresult_root '/' exam_id '_T1_N4.nii'],'file') == 0
        copyfile(T1file, [FSresult_root '/' exam_id '_T1_N4.nii']);
    end
    
    % FS segmentation
    if exist([FSresult_root '/FSresults/FSresult_' exam_id '_T1_N4/stats/aseg.stats'],'file') == 0
        fprintf('- FreeSurfer processing: %s \n', exam_id);
        cmd_str = sprintf('%s recon-all -i %s -s %s -all &', ...
            cmd_str, [exam_id '_T1_N4.nii'], ...
            ['/FSresults/FSresult_' exam_id '_T1_N4']);
        %         cmd_str = sprintf('%s recon-all -s %s -all -no-isrunning &', ...
        %             cmd_str, ['FSresult_' exam_id '_T1_N4']);
    end
    
end
disp(cmd_str);

%% helper functions

function [T1file] = findT1N4(examPath)

% find N4 T1 file
T1file = dir([examPath '/images/']);
T1file(~contains({T1file.name},'t1v_N4.nii') & ~contains({T1file.name},'t1_N4.nii')) = [];
T1file(contains({T1file.name},'ct') | contains({T1file.name},'ewc') | contains({T1file.name},'iy') ...
    | contains({T1file.name},'icent') | contains({T1file.name},'wc')) = [];
if isempty(T1file)
    T1file = dir([examPath '/']);
    T1file(~contains({T1file.name},'t1v_N4.nii') & ~contains({T1file.name},'t1_N4.nii')) = [];
    T1file(cellfun(@isempty,(regexp({T1file.name},'t[0123456789]'))) & ~contains({T1file.name},'T1w')) = [];
    if isempty(T1file)
        T1file = dir([examPath '/../images/']);
        T1file(~contains({T1file.name},'t1v_N4.nii') & ~contains({T1file.name},'t1_N4.nii')) = [];
        T1file(cellfun(@isempty,(regexp({T1file.name},'t[0123456789]'))) & ~contains({T1file.name},'T1w')) = [];
        if isempty(T1file)
            T1file = dir([examPath '/../']);
            T1file(~contains({T1file.name},'t1v_N4.nii') & ~contains({T1file.name},'t1_N4.nii')) = [];
            T1file(cellfun(@isempty,(regexp({T1file.name},'t[0123456789]'))) & ~contains({T1file.name},'T1w')) = [];
            if isempty(T1file)
                fprintf(' - T1_N4 file Not found! \n');
            else
                T1file = [examPath '/../' T1file(1).name];
            end
        else
            T1file = [examPath '/../images/' T1file(1).name];
        end
    else
        T1file = [examPath '/' T1file(1).name];
    end
else
    T1file = [examPath '/images/' T1file(1).name];
end

end

function [t1Folder] = findT1(subjPath)

% find T1 dicom
examFolder = dir([subjPath '/']);
examFolder(cellfun(@isempty,(regexp({examFolder.name},'E[0123456789]')))) = [];
examPath = [subjPath '/' examFolder.name '/'];
if isempty(examFolder)
    examFolder = dir([subjPath '/../']);
    examFolder(cellfun(@isempty,(regexp({examFolder.name},'E[0123456789]')))) = [];
    examPath = [subjPath '/../' examFolder.name '/'];
end

if isempty(examFolder)
    fprintf(' - No T1 DICOM found! \n');
else
    t1Folder = dir(examPath);
    t1Folder(contains({t1Folder.name},'.')) = [];
    for jj = 1:length(t1Folder)
        [~,out] = system(sprintf('dicom_info %s',[examPath '/' t1Folder(jj).name]));
        if contains(out,'SPGR')
            t1Folder = [examPath '/' t1Folder(jj).name];
            break
        end
    end
    fprintf(' - T1 folder %s \n', t1Folder);
    
    % convert to nifti
    if exist([subjPath '/images/T1w.nii.gz'],'file') == 0
        mkdir([subjPath '/images/']);
        cmd = sprintf('dcm2niix -z y -o %s -f T1w %s', [subjPath '/images/'], t1Folder);
        system(cmd);
    end
    
    % N4 correction
    N4corr([subjPath '/images/T1w'],3);
end

end
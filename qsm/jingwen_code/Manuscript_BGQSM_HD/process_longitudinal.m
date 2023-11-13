clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/031_HD_NDM'));

matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/';

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20230213.xlsx','Sheet','Rep');
subjList = [T.b_num];
examList = [T.t_num];

%% resgiter to 1st timepoint

DOBunique = unique(T.dob);

for ii = 1:length(DOBunique)
    indSubj = T.dob == DOBunique(ii);
    
    tp1_path = T.QSM_path{indSubj & T.scan_num == 1};
    fprintf('TP1: %s\n', tp1_path);
    
    tp1_T1_path = [tp1_path 'FSseg/T1_brain.nii.gz'];
    tp1_ROI_path = [tp1_path 'FSseg/Seg_ANTS_manual.nii.gz'];
    
    if sum(indSubj & T.scan_num == 2) > 0
        tp2_path = T.QSM_path{indSubj & T.scan_num == 2};
        exam_id = [subjList{indSubj & T.scan_num == 2} '_' examList{indSubj & T.scan_num == 2}];
        fprintf('TP2: %s, ID: %s\n', tp2_path, exam_id);
        registerTP(tp2_path, tp1_T1_path, tp1_ROI_path, exam_id);
    end
    
    if sum(indSubj & T.scan_num == 3) > 0
        tp3_path = T.QSM_path{indSubj & T.scan_num == 3};
        exam_id = [subjList{indSubj & T.scan_num == 3} '_' examList{indSubj & T.scan_num == 3}];
        fprintf('TP3: %s, ID: %s\n', tp3_path, exam_id);
        registerTP(tp3_path, tp1_T1_path, tp1_ROI_path, exam_id);
    end
end

%% inverse register the ROIs and latvent masks

DOBunique = unique(T.dob);

for ii = 1:length(DOBunique)
    indSubj = T.dob == DOBunique(ii);
    
    tp1_path = T.QSM_path{indSubj & T.scan_num == 1};
    fprintf('TP1: %s\n', tp1_path);
    
    tp1_T1_path = [tp1_path 'FSseg/T1_brain.nii.gz'];
    tp1_ROI_path = [tp1_path 'FSseg/Seg_ANTS_manual.nii.gz'];
    
    if sum(indSubj & T.scan_num == 2) > 0
        tp2_path = T.QSM_path{indSubj & T.scan_num == 2};
        exam_id = [subjList{indSubj & T.scan_num == 2} '_' examList{indSubj & T.scan_num == 2}];
        fprintf('TP2: %s, ID: %s\n', tp2_path, exam_id);
        inversewarpMasks(tp1_path, tp2_path);
    end
    
    if sum(indSubj & T.scan_num == 3) > 0
        tp3_path = T.QSM_path{indSubj & T.scan_num == 3};
        exam_id = [subjList{indSubj & T.scan_num == 3} '_' examList{indSubj & T.scan_num == 3}];
        fprintf('TP3: %s, ID: %s\n', tp3_path, exam_id);
        inversewarpMasks(tp1_path, tp3_path);
    end
end

%% function

function [] = inversewarpMasks(tp1_path, subjPath)

T1native_path = [subjPath 'FSseg/T1_brain.nii.gz'];
ROInative_path = [subjPath 'FSseg/Seg_ANTS.nii.gz'];
ROInativeClean_path = [subjPath 'FSseg/Seg_ANTS_manual.nii.gz'];
CSFnative_path = [subjPath 'FSseg/FAST/T1_brain_seg_0.nii.gz'];
Latventnative_path = [subjPath 'FSseg/latven_mask_reg.nii.gz'];

delete(Latventnative_path);

ROItp1_path = [tp1_path 'FSseg/Seg_ANTS_manual.nii.gz'];
trans1_path = [subjPath 'FSseg/TP1/T1_affine2TP1_1InverseWarp.nii.gz'];
trans2_path = [subjPath 'FSseg/TP1/T1_affine2TP1_0GenericAffine.mat'];
CSFtp1_path = [tp1_path 'FSseg/latven_mask_reg.nii.gz'];

if exist(ROInativeClean_path,'file') == 0
    cmd = sprintf(['antsApplyTransforms -d 3 ' ...
        '-i %s -r %s -t [%s, 1] -t %s -o %s -n NearestNeighbor'], ...
        ROItp1_path, T1native_path, trans2_path, trans1_path, ROInative_path);
    system(cmd);
    cmd = sprintf('3dcalc -a %s -b %s -expr ''a*(1-b)'' -prefix %s', ...
        ROInative_path, CSFnative_path, ROInativeClean_path);
    system(cmd);
end

if exist(Latventnative_path,'file') == 0
    cmd = sprintf(['antsApplyTransforms -d 3 ' ...
        '-i %s -r %s -t [%s, 1] -t %s -o %s -n NearestNeighbor'], ...
        CSFtp1_path, T1native_path, trans2_path, trans1_path, Latventnative_path);
    system(cmd);
end

end

function [] = registerTP(tp_path, tp1_T1_path, tp1_ROI_path, exam_id)

tp_ROI_path = [tp_path 'FSseg/Seg_ANTS_manual.nii.gz'];
if exist(tp_ROI_path,'file') == 0
    fprintf('   ## No segmentation yet, processing \n');
    
    % copy TP1 files
    mkdir([tp_path 'FSseg/TP1/']);
    copyfile(tp1_T1_path, [tp_path 'FSseg/TP1/T1_TP1.nii.gz']);
    copyfile(tp1_ROI_path, [tp_path 'FSseg/TP1/ROI_TP1.nii.gz']);
    
    % prepare data
    T1_path = [tp_path 'FSseg/' exam_id '_T1_N4.nii'];
    T1brain_path = [tp_path 'FSseg/T1_brain.nii.gz'];
    mask_path = [tp_path 'FSseg/brain_mask_reg.nii.gz'];
    if exist(T1brain_path,'file') == 0
        cmd_str = sprintf('3dcalc -a %s -b %s -expr a*b -prefix %s', ...
            T1_path, mask_path, T1brain_path);
        system(cmd_str);
    end
    
    % ANTS registration
    if exist([tp_path 'FSseg/TP1/T1_affine2TP1_warped.nii.gz'],'file') == 0
        cmd_str = ['rootPath="' tp_path 'FSseg/";'...
            'regPath="${rootPath}TP1/";'...
            't1brainPath="${rootPath}T1_brain.nii.gz";'...
            'outputPath="${rootPath}TP1/T1_affine2TP1_";'...
            'atlasPath="${rootPath}TP1/T1_TP1.nii.gz";'...
            'cp $atlasPath $regPath;'...
            'cp $t1brainPath $regPath;'...
            'antsRegistration --dimensionality 3 --float 0 '...
            '--output [${outputPath},${outputPath}warped.nii.gz] '...
            '--interpolation Linear '...
            '--winsorize-image-intensities [0.005,0.995] '...
            '--use-histogram-matching 0 '...
            '--initial-moving-transform [$atlasPath,$t1brainPath,1] '...
            '--transform Rigid[0.1] '...
            '--metric MI[$atlasPath,$t1brainPath,1,32,Regular,0.25] '...
            '--convergence [1000x500x250x100,1e-6,10] '...
            '--shrink-factors 8x4x2x1 '...
            '--smoothing-sigmas 3x2x1x0vox '...
            '--transform Affine[0.1] '...
            '--metric MI[$atlasPath,$t1brainPath,1,32,Regular,0.25] '...
            '--convergence [1000x500x250x100,1e-6,10] '...
            '--shrink-factors 8x4x2x1 '...
            '--smoothing-sigmas 3x2x1x0vox '...
            '--transform SyN[0.1,3,0] '...
            '--metric CC[$atlasPath,$t1brainPath,1,4] '...
            '--convergence [100x70x50x20,1e-6,10] '...
            '--shrink-factors 8x4x2x1 '...
            '--smoothing-sigmas 3x2x1x0vox '...
            '--verbose'];
        system(cmd_str);
    else
        disp(['   Registration existed for ' tp_path]);
    end
    
end

end
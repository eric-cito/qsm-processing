clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/031_HD_NDM'));

FSresult_root = '/working/lupolab/jingwen/002_HD_NDM/AllT1N4';
matout_root = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/';

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','NoRep');
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status_reclass];

%% TRK processing

DTIpath_list = cell(length(subjList),1);

for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    
    [DTIpath] = findDTIfolder(subjPath);
    fprintf(' - DTI path: %s \n', DTIpath);
    DTIpath_list{ii} = DTIpath;
    
    cd(DTIpath);
    mkdir('TRK');
    cd('TRK');
    
    % transpose and construct bvec
    if exist('dti_multiB.bvec','file') == 0
        system('1dtranspose ../data.bvec dti_T.bvec');
        system('1dtranspose ../data.bval dti_T.bval');
        system('1dcat dti_T.bvec dti_T.bval > dti_multiB.bvec');
    end
    % recon tensor
    if exist('dti_tensor.nii','file') == 0
        system('/netopt/dtk/dti_recon ../data_topup_ec.nii.gz dti -gm dti_multiB.bvec 1 -b 2000 -b0 7 -ot nii');
    end
    % brain masking
    if exist('dti_dwi_brain_mask.nii.gz','file') == 0
        system('bet2 dti_dwi.nii dti_dwi_brain -f 0.5 -g 0 -m');
    end
    % track
    if exist('track_tmp.trk','file') == 0
        system('/netopt/dtk/dti_tracker dti track_tmp.trk -at 35 -m dti_dwi_brain_mask.nii.gz -it nii');
    end
    % interpolate
    if exist('dti.trk','file') == 0
        nii = load_nii('dti_b0.nii');
        stepsize = min(nii.hdr.dime.pixdim(2:4));
        system(sprintf('/netopt/dtk/spline_filter track_tmp.trk %f dti.trk', stepsize));
    end
    % register
    if exist('dti_reg.trk','file') == 0
        copyfile([subjPath '/swan_qsm/HDBET_allQSM/FSseg/T1_brain.nii.gz'],'T1_brain.nii.gz');
        T1file = 'T1_brain.nii.gz';
        b0file = 'dti_b0.nii';
        regfile = 'reg.mat';
        system(sprintf('flirt -v -in %s -ref %s -out %s -omat %s', ...
            b0file, T1file, 'b0_reg.nii.gz', regfile));
        system(sprintf('/netopt/dtk/track_transform dti.trk %s -src %s -ref %s -reg %s', ...
            'dti_reg.trk', b0file, T1file, regfile));
    end
end

%% register brain atlases

for ii = 2 % :length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    
    [DTIpath] = findDTIfolder(subjPath);
    fprintf(' - DTI path: %s \n', DTIpath);
    
    cd([DTIpath '/TRK']);
    mkdir('ROIs');
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii} ...
        '/swan_qsm/HDBET_allQSM/'];
    
    % WM
    T1native_path = 'T1_brain.nii.gz';
    ROInative_path = 'JHU_WM_ROIs.nii.gz';
    MNImask_path = '/working/lupolab/jingwen/004_HD_DSI/Atlas/JHU/JHU-ICBM-labels-1mm.nii.gz';
    MNItrans1_path = [subjPath 'ANTSreg/T1_MNI_1InverseWarp.nii.gz'];
    MNItrans2_path = [subjPath 'ANTSreg/T1_MNI_0GenericAffine.mat'];
    
    if exist(ROInative_path,'file') == 0
        cmd = sprintf(['antsApplyTransforms -d 3 ' ...
            '-i %s -r %s -t [%s, 1] -t %s -o %s -n NearestNeighbor'], ...
            MNImask_path, T1native_path, MNItrans2_path, MNItrans1_path, ROInative_path);
        system(cmd);
    end
    
    ROIlist = {'CST_R','CST_L','CP_R','CP_L','PLIC_R','PLIC_L','PCR_R','PCR_L','ACR_R','ACR_L'};
    ROInumber = [7 8 15 16 19 20 27 28 23 24];
    for rr = 1:length(ROIlist)
        cmd = sprintf('3dcalc -a JHU_WM_ROIs.nii.gz -expr ''equals(a,%i)'' -prefix ROIs/%s.nii.gz', ...
            ROInumber(rr), ROIlist{rr});
        system(cmd);
        cmd = sprintf('3dmask_tool -input ROIs/%s.nii.gz -prefix ROIs/%s_d5.nii.gz -dilate_input 5', ...
            ROIlist{rr}, ROIlist{rr});
        system(cmd);
    end
    
    % GM Cort
    ROInative_path = 'HavardOxford_GM_ROIs.nii.gz';
    MNImask_path = '/working/lupolab/jingwen/004_HD_DSI/Atlas/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz';
    if exist(ROInative_path,'file') == 0
        cmd = sprintf(['antsApplyTransforms -d 3 ' ...
            '-i %s -r %s -t [%s, 1] -t %s -o %s -n NearestNeighbor'], ...
            MNImask_path, T1native_path, MNItrans2_path, MNItrans1_path, ROInative_path);
        system(cmd);
    end
    
    for rr = 1:length(ROIlist)
        cmd = sprintf('3dcalc -a %s -expr ''equals(a,7)+equals(a,17)'' -prefix ROIs/PrPoC.nii.gz', ...
            ROInative_path);
        system(cmd);
        cmd = sprintf('3dmask_tool -input ROIs/PrPoC.nii.gz -prefix ROIs/PrPoC_d5.nii.gz -dilate_input 5');
        system(cmd);
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

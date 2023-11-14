function [] = reg_QSM(input_data_path, output_data_path)

%% add path

warning('off');
addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));

%% set up code/data paths

% Add FSL path
% FSL_PATH = '/netopt/rhel7/fsl/';
% addpath(genpath(FSL_PATH));
% setenv('FSLDIR', '/netopt/rhel7/fsl/');
% fsldir = getenv('FSLDIR');
% fsldirmpath = sprintf('%s/etc/matlab',fsldir);
% path(path, fsldirmpath);
% setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
% clear fsldir fsldirmpath;

% detect iLSQR QSM file
if exist([input_data_path '/QSM_iLSQR_meanEcho.nii.gz'],'file') ~= 2
    error('Cannot find iLSQR QSM file in the input path.');
else
    subjQSMfile = [input_data_path '/QSM_iLSQR_meanEcho.nii.gz'];
end

% detect mask file
if exist([input_data_path '/brain_mask_HD.nii.gz'],'file') ~= 2
    error('Cannot find HD mask file in the input path.');
else
    subjMaskfile = [input_data_path '/brain_mask_HD.nii.gz'];
end

%% Load individual QSM image (NIFTI) and remove background

nii = load_untouch_nii(subjQSMfile);
Img_ind_org = nii.img;
nii = load_untouch_nii(subjMaskfile);
mask = nii.img;

shift = min(Img_ind_org(:));
    
if exist([output_data_path '/QSM_shift.nii.gz'],'file') ~= 2
    Img_ind = Img_ind_org - shift.*ones(size(Img_ind_org));
    Img_ind = Img_ind.*mask + shift.*ones(size(Img_ind_org));
    nii.img = Img_ind;
    save_untouch_nii(nii, [output_data_path '/QSM_shift.nii.gz']);
end

if exist([output_data_path '/QSM_orient.nii.gz'],'file') ~= 2
        cmd = sprintf('3dresample -input %s -prefix %s -orient LPI',...
            [output_data_path '/QSM_shift.nii.gz'],...
            [output_data_path '/QSM_orient.nii.gz']);
        system(cmd);
%     cmd = sprintf('fslreorient2std %s %s',...
%         [output_data_path '/QSM_shift.nii.gz'],...
%         [output_data_path '/QSM_orient.nii.gz']);
%     system(cmd);
end

atlasQSM = '/working/lupolab/QSM_Atlas_MNI_toolbox/QSM_atlas_blackbackground_MNI_1mm_v1.nii';
atlasMask = '/working/lupolab/QSM_Atlas_MNI_toolbox/QSM_atlas_mask.nii.gz';
atlasROI = '/working/lupolab/QSM_Atlas_MNI_toolbox/Whole_brain_segmentation_DGM.nii';

copyfile(atlasQSM, [output_data_path '/QSM_atlas.nii']);

%% Run FSL linear registration

disp(' ## Linear registration');

if exist([output_data_path '/QSM_regLin.nii.gz'],'file') ~= 2      
    dof = 12;
    cmdstr = ['/netopt/rhel7/fsl/bin/flirt -v' ...
        ' -in ' output_data_path '/QSM_orient.nii.gz' ...
        ' -ref ' atlasQSM ...
        ' -out ' output_data_path '/QSM_regLin.nii.gz' ...
        ' -omat ' output_data_path '/QSM_regLin.mat' ...
        ' -cost corratio -searchrz -90 90 -dof ' num2str(dof) ' -interp trilinear'];
    system(cmdstr);
end

%% Run nonlinear registration

disp(' ## Nonlinear registration');

if exist([output_data_path '/QSM_regNonlin.nii.gz'],'file') ~= 2
    cmdstr = ['bash /home/jyao3/030_QSM/01_Code/convertMat.sh'...
        ' ' output_data_path '/QSM_regLin.mat'...
        ' > ' output_data_path '/QSM_regLinF.mat'];
    system(cmdstr);
    cmdstr = ['/netopt/rhel7/fsl/bin/fnirt -v ' ...
        ' --ref=' atlasQSM ...
        ' --in=' output_data_path '/QSM_orient.nii.gz' ...
        ... % ' --inmask=' input_data_path '/brain_mask_HD.nii.gz' ...
        ... % ' --refmask=' atlasMask ...
        ' --aff=' output_data_path '/QSM_regLinF.mat' ...
        ' --intmod=global_linear'  ...
        ' --iout=' output_data_path '/QSM_regNonlin.nii.gz' ...
        ' --fout=' output_data_path '/QSM_regNonlinDF.nii.gz' ...
        ' --cout=' output_data_path '/QSM_regNonlin_warp.nii.gz' ...
        ' --jacrange=-1,100'];
    system(cmdstr);
end

%% Inverse Transformation

disp(' ## Inverse Transformation ');

if exist([output_data_path '/QSM_regNonlin_warp_inv.nii.gz'],'file') ~= 2
    cmdstr=['invwarp -v --ref=' [output_data_path '/QSM_orient.nii.gz']...
        ' --warp=' output_data_path '/QSM_regNonlin_warp.nii.gz' ...
        ' --out=' output_data_path '/QSM_regNonlin_warp_inv.nii.gz'];
    system(cmdstr);
end

if exist([output_data_path '/QSM_atlas_ROI.nii.gz'],'file') ~= 2
    cmdstr = ['applywarp --ref=' [output_data_path '/QSM_orient.nii.gz']...
        ' --in=' atlasROI ...
        ' --warp=' output_data_path '/QSM_regNonlin_warp_inv.nii.gz' ...
        ' --out=' output_data_path '/QSM_atlas_ROI.nii.gz' ...
        ' --interp=nn'];
    system(cmdstr);
end

end

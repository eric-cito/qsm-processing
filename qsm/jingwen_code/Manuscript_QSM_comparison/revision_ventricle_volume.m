clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));

img_root = '/working/lupolab/jingwen/001_QSM/temp';

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','NoRep');
subjList = [T.b_num];
examList = [T.t_num];
statusList = T.status_reclass;

%% Loop through subjects

Vcsf = zeros(length(subjList),1);
for ii = 1:length(subjList)
    
    fprintf('%s \n', subjList{ii});
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii} ...
        '/swan_qsm/HDBET_allQSM/FSseg/'];
    
    CSF_path = [subjPath 'latven_mask_reg.nii.gz'];
    nii = load_nii(CSF_path);
    VENmask = double(nii.img) > 0;
    
    voxSize = nii.hdr.dime.pixdim(2:4);
    Vcsf(ii) = sum(VENmask,'all')*voxSize(1)*voxSize(2)*voxSize(3);
    
end

%% Ventricle ROI size

set(0,'defaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1);

figure('position', [100 100 400 400]);
histogram(Vcsf(strcmp(statusList,'HC'))/1000,[0:2:50]); hold on;
histogram(Vcsf(~strcmp(statusList,'HC'))/1000,[0:2:50]);
xlabel('Ventricle ROI volume (ml)');
ylabel('Count');
legend({'Healthy control','HD'});
export_fig([img_root '/Ventricle_volume'], '-png','-transparent'); % close;

median(Vcsf(strcmp(statusList,'HC'))/1000)
median(Vcsf(~strcmp(statusList,'HC'))/1000)

clear; clc;
warning('off');

%% Add path

addpath('/home/jyao3/030_QSM/01_Code/Manuscript_BGQSM_HD');
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

img_path = '/working/lupolab/jingwen/001_QSM/temp';

%% read data

examPath = '/data/7T_hunt/b4473/t12317/swan_qsm/HDBET_allQSM';

ImgFile = [examPath '/Magni.nii.gz'];
nii = load_nii(ImgFile);
Magni = flip(rot90(double(nii.img(:,:,:,1))),2);

ImgFile = [examPath '/Phase.nii.gz'];
nii = load_nii(ImgFile);
Phase = flip(rot90(double(nii.img(:,:,:,2))),2);

%% plot image

xrange = [38:273]; yrange = [56:247];
xrange2 = [100:170]; yrange2 = [110:200];
slice = 75;

plot_save_img(Magni,xrange,yrange,slice,[0 7000]);
export_fig([img_path '/Magni_1'], '-png','-transparent'); close;

plot_save_img(Magni,xrange2,yrange2,slice,[0 7000]);
export_fig([img_path '/Magni_2'], '-png','-transparent'); close;

plot_save_img(Phase,xrange,yrange,slice,[-pi pi]);
export_fig([img_path '/Phase_1'], '-png','-transparent'); close;

plot_save_img(Phase,xrange2,yrange2,slice,[-pi pi]);
export_fig([img_path '/Phase_2'], '-png','-transparent'); close;

%% plot QSM

QSMfile_list = {
    'QSM_iLSQR_meanEcho' ...
    'QSM_QSIP_meanEcho' ...
    'QSM_SSTGV_meanEcho' ...
    'QSM_SSTV_meanEcho' ...
    'QSM_STARQSM_meanEcho' ...
    'QSM_FANSI_nonlinearTV_meanEcho' ...
    'QSM_HDQSM_meanEcho' ...
    'QSM_MEDI_meanEcho' ...
    'QSM_QSMGAN_meanEcho' ...
    'QSM_QSMnet_meanEcho' ...
    'QSM_xQSM2_meanEcho' ...
    'QSM_iQSM2_meanEcho'};

ImgFile = [examPath '/QSMseg_latven.nii.gz'];
nii = load_nii(ImgFile);
LatVent = flip(rot90(double(nii.img)),2);

QSMimage1 = [];
QSMimage2 = [];
for ii = 1:length(QSMfile_list)
    
    ImgFile = [examPath '/' QSMfile_list{ii} '.nii.gz'];
    nii = load_nii(ImgFile);
    QSMmap = flip(rot90(double(nii.img)),2);
    
    if ii == 9
        QSMmap = QSMmap/0.5684;
    end
    
    qsmCSF = median(QSMmap(LatVent > 0));
    disp(qsmCSF);
    
    mask = QSMmap ~= 0;
    QSMmap = (QSMmap - qsmCSF).*mask;
    QSMimage1(:,:,ii) = QSMmap(xrange,yrange,slice);
    QSMimage2(:,:,ii) = QSMmap(xrange2,yrange2,slice);
    
%     plot_save_img(QSMmap,xrange,yrange,slice,[-0.15 0.15]);
%     export_fig([img_path '/' QSMfile_list{ii} '_1'], '-png','-transparent'); close;
%     
%     plot_save_img(QSMmap,xrange2,yrange2,slice,[-0.15 0.15]);
%     export_fig([img_path '/' QSMfile_list{ii} '_2'], '-png','-transparent'); close;
    
end

%% plot montage

figure('position', [100 0 800 800]);
montage(QSMimage1, 'DisplayRange', [-0.15 0.15], 'Size', [2 6]); colormap gray;
axis equal tight off; pause(1);
export_fig([img_path '/QSMexam1'], '-png','-transparent'); close;

figure('position', [100 0 800 800]);
montage(QSMimage2, 'DisplayRange', [-0.15 0.15], 'Size', [2 6]); colormap gray;
axis equal tight off; pause(1);
export_fig([img_path '/QSMexam2'], '-png','-transparent'); close;

%% read data - MNI reg

examPath = '/data/7T_hunt/b3025/t12518/swan_qsm/HDBET_allQSM/MNIreg';

ImgFile = [examPath '/QSM_atlas.nii'];
nii = load_nii(ImgFile);
Atlas = flip(rot90(double(nii.img)),2);

ImgFile = [examPath '/QSM_iLSQR_meanEcho_MNInorm.nii.gz'];
nii = load_nii(ImgFile);
MNIqsm = flip(rot90(double(nii.img)),2);

ImgFile = [examPath '/QSM_atlas_ROI.nii.gz'];
nii = load_nii(ImgFile);
AtlasROI = flip(rot90(double(nii.img)),2);

%% plot

figure('position', [100 0 800 800]);
montage(Atlas(:,:,50:10:120), 'DisplayRange', [-0.15 0.15], 'Size', [2 4]); colormap gray;
axis equal tight off; pause(1);
export_fig([img_path '/QSM_MNIatlas'], '-png','-transparent');

figure('position', [100 0 800 800]);
montage(MNIqsm(:,:,50:10:120), 'DisplayRange', [-0.15 0.15], 'Size', [2 4]); colormap gray;
axis equal tight off; pause(1);
export_fig([img_path '/QSM_MNIreg'], '-png','-transparent');

%% plot QSM_MNI averaged

QSMfile_list = {
    'QSM_iLSQR' ...
    'QSM_QSIP' ...
    'QSM_SSTGV' ...
    'QSM_SSTV' ...
    'QSM_STARQSM' ...
    'QSM_FANSI_nonlinearTV' ...
    'QSM_HDQSM' ...
    'QSM_MEDI' ...
    'QSM_QSMGAN' ...
    'QSM_QSMnet' ...
    'QSM_xQSM2' ...
    'QSM_iQSM2'};

xrange = [15:203]; yrange = [15:167];
xrange2 = [65:118]; yrange2 = [50:130];
slice = 75;

examPath = '/working/lupolab/jingwen/001_QSM/03_QSM_HD/Data';

ImgFile = [examPath '/../Atlas/MNI_T1_brain_mask.nii.gz'];
nii = load_nii(ImgFile);
mask = flip(rot90(double(nii.img)),2);

QSMimage1 = [];
QSMimage2 = [];
for ii = 1:length(QSMfile_list)
    
    ImgFile = [examPath '/' QSMfile_list{ii} '_MNI_HCmean.nii.gz'];
    nii = load_nii(ImgFile);
    QSMmap = flip(rot90(double(nii.img)),2).*mask;
    
    if ii == 9
        QSMmap = QSMmap/0.5684;
    end

    QSMimage1(:,:,ii) = QSMmap(xrange,yrange,slice);
    QSMimage2(:,:,ii) = QSMmap(xrange2,yrange2,slice);
    
end

figure('position', [100 0 800 800]);
montage(QSMimage1, 'DisplayRange', [-0.15 0.15], 'Size', [2 6]); colormap gray;
axis equal tight off; pause(1);
export_fig([img_path '/QSM_MNI_HC1'], '-png','-transparent'); % close;

figure('position', [100 0 800 800]);
montage(QSMimage2, 'DisplayRange', [-0.15 0.15], 'Size', [2 6]); colormap gray;
axis equal tight off; pause(1);
export_fig([img_path '/QSM_MNI_HC2'], '-png','-transparent'); % close;

%% helper function

function [] = plot_save_img(Img,xrange,yrange,slice,irange)

figure('position', [100 0 800 800]);
imagesc(Img(xrange,yrange,slice), irange); colormap gray;
axis equal tight off
pause(1);

end

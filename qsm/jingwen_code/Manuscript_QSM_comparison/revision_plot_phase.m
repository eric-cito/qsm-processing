clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/031_HD_NDM'));

img_path = '/working/lupolab/jingwen/001_QSM/temp';

%% read in subject list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','NoRep');
T(T.CAG < 36,:) = [];
T(strcmp(T.status_reclass,'MISSING'),:) = [];
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status_reclass];

%% Loop through subjects

rng(777);
SubjInd = [1 21 3 4 5]; % randi(length(subjList), 1, 5); 

subjPath = ['/data/7T_hunt/' subjList{SubjInd(5)} '/' examList{SubjInd(5)} ...
    '/swan_qsm/HDBET_allQSM/'];
nii = load_nii([subjPath '/brain_mask_HD.nii.gz']);
brain_mask = flip(rot90(double(nii.img)),2);
[crop] = cropImg(brain_mask);

sliceAx = 60:10:100;
sliceSag = 130:10:170;

TotalField_plot = [];
LocalField_plot = [];
TotalField_plotSag = [];
LocalField_plotSag = [];
for ii = 1:5
    
    subjPath = ['/data/7T_hunt/' subjList{SubjInd(ii)} '/' examList{SubjInd(ii)} ...
        '/swan_qsm/HDBET_allQSM/'];
    
    nii = load_nii([subjPath '/totalField_meanEcho.nii.gz']);
    TotalField = flip(rot90(double(nii.img)),2);
    
    nii = load_nii([subjPath '/localField_meanEcho.nii.gz']);
    LocalField = flip(rot90(double(nii.img)),2);
    
    TotalField_plot(:,:,(ii-1)*5+1:ii*5) = -TotalField(crop.X,crop.Y,sliceAx);
    LocalField_plot(:,:,(ii-1)*5+1:ii*5) = -LocalField(crop.X,crop.Y,sliceAx);
    
    TotalField_plotSag(:,:,(ii-1)*5+1:ii*5) = -flip(permute(TotalField(crop.X,sliceSag,crop.Z), [3 1 2]),1);
    LocalField_plotSag(:,:,(ii-1)*5+1:ii*5) = -flip(permute(LocalField(crop.X,sliceSag,crop.Z), [3 1 2]),1);
    
end

figure;
montage(TotalField_plot/(2*pi*0.0165*42.58*7)*1e3, 'DisplayRange', [-200 200], 'Size', [5 5]); colormap gray;
export_fig([img_path '/TotalField_plot'], '-png','-transparent'); close;

figure;
montage(LocalField_plot/(0.0165*42.58*7)*1e3, 'DisplayRange', [-40 40], 'Size', [5 5]); colormap gray;
export_fig([img_path '/LocalField_plot'], '-png','-transparent'); close;

figure;
montage(TotalField_plotSag/(2*pi*0.0165*42.58*7)*1e3, 'DisplayRange', [-200 200], 'Size', [5 5]); colormap gray;
export_fig([img_path '/TotalField_plotSag'], '-png','-transparent'); close;

figure;
montage(LocalField_plotSag/(0.0165*42.58*7)*1e3, 'DisplayRange', [-40 40], 'Size', [5 5]); colormap gray;
export_fig([img_path '/LocalField_plotSag'], '-png','-transparent'); close;

%% visual inspection

for ii = 40:49 % 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii} ...
        '/swan_qsm/HDBET_allQSM/'];
    
    nii = load_nii([subjPath '/localField_meanEcho.nii.gz']);
    LocalField = flip(rot90(double(nii.img)),2);
    LocalField_plot = -LocalField(crop.X,crop.Y,sliceAx);
    LocalField_plotSag = -flip(permute(LocalField(crop.X,sliceSag,crop.Z), [3 1 2]),1);
    
    subplot(211);
    montage(LocalField_plot/(0.0165*42.58*7)*1e3, 'DisplayRange', [-40 40], 'Size', [1 5]); colormap gray;
    subplot(212);
    montage(LocalField_plotSag/(0.0165*42.58*7)*1e3, 'DisplayRange', [-40 40], 'Size', [1 5]); colormap gray;
    title(subjList{ii});
    pause; 
    
end

%% function

function [crop] = cropImg(brain_mask)

[x, y, z] = ind2sub(size(brain_mask), find(brain_mask > 0));
gap = 5; % set a magin of 5 pixels
x1 = max(1, min(x)-gap); x2 = min(size(brain_mask,1), max(x)+gap);
y1 = max(1, min(y)-gap); y2 = min(size(brain_mask,2), max(y)+gap);
z1 = max(1, min(z)-gap); z2 = min(size(brain_mask,3), max(z)+gap);
if mod(x2 - x1, 2) == 0
    if x2 == size(brain_mask,1)
        x1 = x1 - 1;
    else
        x2 = x2 + 1;
    end
end
if mod(y2 - y1, 2) == 0
    if y2 == size(brain_mask,2)
        y1 = y1 - 1;
    else
        y2 = y2 + 1;
    end
end
if mod(z2 - z1, 2) == 0 && z2 == size(brain_mask,3)
    z1 = z1 - 1;
elseif mod(z2 - z1, 2) == 0
    z2 = z2 + 1;
end

crop.X = x1:x2;
crop.Y = y1:y2;
crop.Z = z1:z2;

end
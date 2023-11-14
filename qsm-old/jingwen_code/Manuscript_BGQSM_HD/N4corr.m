function N4corr(niiRoot,Nclasses)

%clear all; close all; clc;

addpath(genpath('/home/jyao3/031_HD_NDM/02_Library/N4_bias_correction'));

%% Obtain head mask
% I recommend to copy this snippet of code and use with all your images since it
% greatly accelerates bias field computation and improves results


if exist([niiRoot '.nii.gz'],'file') == 2
    niiFileName = [niiRoot '.nii.gz'];
elseif exist([niiRoot '.nii'],'file') == 2
    niiFileName = [niiRoot '.nii.gz'];
end
FileStruct=load_untouch_nii(niiFileName);
image=FileStruct.img;

%%  Triangle thresholding
level = triangle_th(hist(mat2gray(image(image>0)),256)',256);
mask = mat2gray(image) > level;

% Median filter to reduce noise
mask = medfilt3(mask);

% Fill holes
mask = imfill3(mask);

% Identify connected components
CC = bwconncomp(mask);

% Measure volume of each object
numPixels = cellfun(@numel,CC.PixelIdxList);

% Identify largest object (should be the head)
[biggest,idx] = max(numPixels);

% Extract head mask
newMask = zeros(size(mask));
newMask(CC.PixelIdxList{idx}) = 1;
 
mask = logical(newMask);

% imagine(mask);
%% Run iterative BCFCM+N4
params.N4BiasCorEXE = '/home/jyao3/031_HD_NDM/02_Library/N4_bias_correction/N4BiasFieldCorrection';
params.mask = mask;

if Nclasses == 3
[biasCorrectedImg,biasField,U] = BCFCM_N4(image,params); 
% imagine(image,biasCorrectedImg);
% imagine(U(:,:,:,3),U(:,:,:,2),U(:,:,:,1));

elseif Nclasses == 4
[biasCorrectedImg,biasField,U] = BCFCM_N4_2(image,params); 
% imagine(image,biasCorrectedImg);
% imagine(U(:,:,:,4),U(:,:,:,3),U(:,:,:,2),U(:,:,:,1));
end 

%% Write out biacorrected image

image2 = double(image);
FileStruct.img = mat2gray(biasCorrectedImg).*max(image2(:));
Filename = [niiRoot '_t1v_N4'];

save_untouch_nii(FileStruct, Filename);

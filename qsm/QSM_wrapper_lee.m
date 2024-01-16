% QSM_processing 
% Runs QSM processing via SEPIA. This code simply prepared data for this
% Written by Jingwen
% Modified by Melanie
% Reworked by Lee

clear all

%% Settings
input_data_path = '/Users/lee/data/pda440';
output_data_path = '/Users/lee/data/pda440/processed/';
correctFilter = false;

loc_dcm2niix = '/opt/homebrew/bin/dcm2niix';
loc_fsl = '/Users/lee/binaries/fsl/share/fsl/bin/';
%% add path
dir_this = fileparts(mfilename('fullpath'));

warning('off');
addpath(genpath2(dir_this, '.git'));
addpath(loc_fsl)

%% set up data paths
ID1 = extractBefore(input_data_path,'/qsm');
ptid = extractAfter(ID1,'clin/');
ptid = extractBefore(ptid,'_no.consent.yet-addpost');

cd(input_data_path);


%% Step 1: sort dicoms
myEchos = OrganiseDicoms(input_data_path);


%% Step 2 Correct offset for Re/Im images and overwrite
CorrectOffsetForReImAndCreateNiftis(input_data_path, myEchos, correctFilter, loc_dcm2niix);

%% Step 3 Create Phase from Real + Im
CreatePhase(input_data_path, myEchos)

%% Combine mag & phase echoes
fileLocator = FileLocator(input_data_path);
ConcatImages(fileLocator, 'mag',  myEchos, fileLocator.GetMagnitude_AllEchos());
ConcatImages(fileLocator, 'phase',  myEchos, fileLocator.GetPhase_AllEchos());


%% Prepare header file for processing

% error 'CODE DOES NOT WORK(!) TE files are in qsm_dicoms/echo*/real/*.json Manual override below'
% 
% te1 = 0.2;
% te2 = 0.22628;
% te3 = 0.25256;
% te4 = 0.27884;

%% Find Echo times (0018,0081)

echoes = zeros(1,length(myEchos)); %create empty vector to enter echos

for iEcho = myEchos
    fname = fileLocator.GetBasicTypeWithSuffix(iEcho, 'real', 'json'); 
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    val = jsondecode(str);
    echoes(:,iEcho)= val.EchoTime;
end

%te1 = echoes(1,1);
%te2 = echoes(1,2);
%te3 = echoes(1,3);
%te4 = echoes(1,4);
%te5 = echoes(1,5);

save('echoes.mat', 'echoes')

% te1 = 5.44e-3; % s 
% te2 = 12.0e-3; % s 
% te3 = 18.6e-3; % s 
% te4 = 25.2e-3;% 


% set up dimensions
% spatial resolution of the data, in mm
File = load_untouch_nii(fileLocator.GetBasicType(myEchos(1), 'real'));
header.voxelSize = File.hdr.dime.pixdim([2,3,4]); % nifti - pixdim 2,3,4 have the dimensions. 
% image matrix size
%header.matrixSize = [224,224,60]; %for Phillips clinical QSM
header.matrixSize = File.hdr.dime.dim([2,3,4]);
% set up parameters
header.b0 = 3;                  % magnetic field strength, in Tesla
header.b0dir = [0;0;1];         % main magnetic field direction, [x,y,z]
header.CF = header.b0 * 42.58 * 1e6;       % imaging frequency, in Hz (B0*gyromagnetic_ratio)
header.te = echoes;  % echo time for each GRE image, in seconds
header.delta_TE = echoes(1,2) - echoes(1,1);      % echo spacing, in seconds

% save to mat file
save([input_data_path '/sepia_header.mat'],'header');

%% set up algorithm parameters

%%% Please make sure the file names are the same as used below %%%%%%%%%%%%
input.magnitudeFile     = [char(fileLocator.GetMagnitude_AllEchos())];
input.phaseFile         = [char(fileLocator.GetPhase_AllEchos())];
%%% Check name - end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.headerFile        = [input_data_path '/sepia_header.mat'];

opts.writeLog           = 1;
opts.isGPU              = 0;
%opts.BETmethod          = 'HD-BET';
opts.BETmethod          = 'BET';
opts.iLSQR              = 1;
opts.QSMGAN             = 0;
opts.All                = 0;

%% run

QSMfile = dir(output_data_path);
QSMfile([QSMfile.isdir]) = [];
QSMfile(cellfun(@isempty,strfind({QSMfile.name},'QSM_'))) = [];

%if (opts.All && length(QSMfile) < 10) || (~opts.All && length(QSMfile) < 2)
    [Tcomp] = QSM_processing_CPU(input, output_data_path, opts);
%else
    %fprintf(' >> Already have all QSM maps, skip QSM processing \n');
%end
%
%
%
%
%
%
%
%
%
%
%
%% (optional) convert to dicom
error ABORT
FileStruct = load_untouch_nii('QSM_iLSQR_meanEcho.nii.gz');
img = -double(FileStruct.img); %place negative sign in front if contrast is bad
img = img + 0.15;
img(img < 0) = 0; img(img > 0.3) = 0.3;
img = rot90(int16(img/0.3*4096));
%img = int16(4096) - img;
%img = -img;
imagesc(img(:,:,1)); colormap gray;


%% load 1 dicom, change info, write dicom; 4-D SEQUENCE 
cd ..
cd('dicoms')
dcmfile = '1.3.46.670589.11.71065.5.0.7220.2022011815310387114.dcm';
metadata_orig = dicominfo(dcmfile);

uid = dicomuid;
metadata.SliceThickness = FileStruct.hdr.dime.pixdim(4);
metadata.PixelSpacing = FileStruct.hdr.dime.pixdim(2:3);
metadata.ImageOrientationPatient = metadata_orig.ImageOrientationPatient;
metadata.PatientIdentityRemoved = 'NO';
z_base = FileStruct.hdr.hist.qoffset_z;
metadata.WindowWidth = 3276; %1732
metadata.WindowCenter = 1638; %1982
metadata.ProtocolName = 'QSM_Processed3D';
metadata.SeriesDescription = 'QSM_Processed3D';
metadata.SeriesInstanceUID = uid;
metadata.SpacingBetweenSlices = metadata_orig.SpacingBetweenSlices;

mkdir(['QSM_DICOM']);
cd('QSM_DICOM')

%set the dicom dict
%dicomdict('set', '/data/morrison/data/dystonia/Starr_targeting/dicom-dict_copy.txt');
dicomdict('set','/netopt/rhel7/versions/MATLAB/R2021a/toolbox/images/iptformats/dicom-dict.txt')

for i=1:size(img,3)
    metadata.InstanceNumber = i;
    metadata.Private_2001_100a = i; %slicenumber
    metadata.SliceLocation = z_base + (i-1)*metadata.SliceThickness;
    Filename = sprintf('IM%04i.dcm',i);
    %dicomwrite(img(:,:,i),Filename, metadata, 'CreateMode', 'copy'); 
    dicomwrite(img(:,:,i),Filename, metadata, 'WritePrivate', true, 'CreateMode', 'Copy', 'UseMetadataBitDepths', true) 
    %dicomwrite(img(:,:,i),Filename, metadata, 'WritePrivate', true, 'CreateMode', 'Copy') 
    %dicomwrite(img(:,:,i), Filename, 'VR', 'explicit', 'ObjectType', 'MR Image Storage', ...
    %'WritePrivate', true, metadata);
end 

%CHECK IN MRSC

%% load 1 dicom, change info, write dicom; USING 3-D PD SEQUENCE
cd ..
cd('dicoms_PD')
dcmfile = '1.3.46.670589.11.71065.5.0.7220.2022011813443679500.dcm';
metadata = dicominfo(dcmfile);

uid = dicomuid;
metadata.SliceThickness = FileStruct.hdr.dime.pixdim(4);
metadata.PixelSpacing = FileStruct.hdr.dime.pixdim(2:3);
metadata.Columns = size(img,1);
metadata.Rows = size(img,2);
metadata.ImageOrientationPatient = metadata_orig.ImageOrientationPatient;
metadata.PatientIdentityRemoved = 'NO';
z_base = FileStruct.hdr.hist.qoffset_z;
metadata.WindowWidth = 4096; % 1732;
metadata.WindowCenter = 2048; % 1982;
metadata.ProtocolName = 'QSM_Processed';
metadata.SeriesDescription = 'QSM_Processed';
metadata.SeriesInstanceUID = uid;
metadata.EchoTrainLength = 1;
metadata.SpacingBetweenSlices = metadata_orig.SpacingBetweenSlices;

cd('/data/morrison/data/dystonia/Targeting/PHI_VS_0122/output_QSM');
mkdir('QSM_DICOM_testJY');
cd('QSM_DICOM_testJY')

%set the dicom dict
%dicomdict('set', '/data/morrison/data/dystonia/Starr_targeting/dicom-dict_copy.txt');
dicomdict('set','/netopt/rhel7/versions/MATLAB/R2021a/toolbox/images/iptformats/dicom-dict.txt')

for i=1:size(img,3)
    metadata.InstanceNumber = i;
    metadata.Private_2001_100a = i; %slicenumber
    metadata.SliceLocation = z_base + (i-1)*metadata.SliceThickness;
    Filename = sprintf('IM%04i.dcm',i);
    %dicomwrite(img(:,:,i),Filename, metadata, 'CreateMode', 'copy'); 
    dicomwrite(img(:,:,i),Filename, metadata, 'WritePrivate', true, 'CreateMode', 'Copy', 'UseMetadataBitDepths', true) 
    %dicomwrite(img(:,:,i),Filename, metadata, 'WritePrivate', true, 'CreateMode', 'Copy') 
    %dicomwrite(img(:,:,i), Filename, 'VR', 'explicit', 'ObjectType', 'MR Image Storage', ...
    %'WritePrivate', true, metadata);
end 

%% Jingwen test - write DICOM
cd ..
cd('dicoms')
dcmfile = '1.3.46.670589.11.71065.5.0.7220.2022011815310387114.dcm';
metadata = dicominfo(dcmfile);

uid = dicomuid;
metadata.PatientIdentityRemoved = 'NO';
z_base = FileStruct.hdr.hist.qoffset_z;
metadata.ProtocolName = 'QSM_Processed3D';
metadata.SeriesDescription = 'QSM_Processed3D';
metadata.SeriesInstanceUID = uid;
metadata.EchoTrainLength = 1;

cd('/data/morrison/data/dystonia/Targeting/PHI_VS_0122/output_QSM');
mkdir('QSM_DICOM_testJY');
cd('QSM_DICOM_testJY')

%set the dicom dict
%dicomdict('set', '/data/morrison/data/dystonia/Starr_targeting/dicom-dict_copy.txt');
dicomdict('set','/netopt/rhel7/versions/MATLAB/R2021a/toolbox/images/iptformats/dicom-dict.txt')

for i=1:size(img,3)
    metadata.InstanceNumber = i;
    metadata.Private_2001_100a = i; %slicenumber
    metadata.SliceLocation = z_base + (i-1)*metadata.SliceThickness;
    Filename = sprintf('IM%04i.dcm',i);
    dicomwrite(img(:,:,i),Filename, metadata, 'CreateMode', 'copy'); 
    % dicomwrite(img(:,:,i),Filename, metadata, 'WritePrivate', true, 'CreateMode', 'Copy', 'UseMetadataBitDepths', true) 
    % dicomwrite(img(:,:,i),Filename, metadata, 'WritePrivate', true, 'CreateMode', 'Copy') 
    % dicomwrite(img(:,:,i), Filename, 'VR', 'explicit', 'ObjectType', 'MR Image Storage', ...
    % 'WritePrivate', true, metadata);
end 

%% send to pacs

system('send_to_pacs --in_dir /data/morrison/data/dystonia/Targeting/PHI_VS_0122/output_QSM/QSM_DICOM_PDdcm --no_reid') 
%check in /data/dicom_mb/export/PACS/


% % %move phase data to separate folder
% listing = dir('/data/morrison/data/dystonia/Starr_targeting/PHI_J/qsm1')
% for j = 3:(size(listing,1)-1)
%     meta = dicominfo(listing(j,1).name);
%     if meta.EchoNumbers == 1
%         syscmd = ["mv " + listing(j,1).name + " echo1"];
%         system(syscmd);
%     end 
% end 
%%
% CHANGE - save path
% DCMout_root = '';
% 
% % CHANGE - path to qsm nifti
% QSMpath = '.nii.gz';
% nii = load_untouch_nii(QSMpath);
% img = double(nii.img);
%             
% % rescale to show -0.15 to 0.15
% img = img + 0.15;
% img(img < 0) = 0; img(img > 0.3) = 0.3;
% img = rot90(int16(img/0.3*4096));
% % NOTE - check the orientation of outcome DICOM
%             
% mkdir([DCMout_root '/QSM_DICOM']);
%             
% 
% meta.SeriesInstanceUID = uid;
% meta.SliceThickness = nii.hdr.dime.pixdim(4);
% meta.PixelSpacing = nii.hdr.dime.pixdim(2:3);
% meta.ImageOrientationPatient = [1 0 0 0 1 0];
% 
% % CHANGE - specify name
% meta.PatientName.FamilyName = '';
% meta.PatientName.GivenName = '';
%             
% z_base = nii.hdr.hist.qoffset_z;
%             
% for ss = 1:size(img,3)
%     meta.InstanceNumber = ss;
%     meta.SliceLocation = z_base + (ss-1)*meta.SliceThickness;
%     dcm_filename = sprintf('%s/QSM_DICOM/IM%04i.dcm',DCMout_root,ss);
%     dicomwrite(img(:,:,ss),dcm_filename,meta,...
%         'CreateMode', 'copy');
% end
% 

function ConcatImages(fileLocator, imgType, echoes, saveTo)
    %ConcatImagesAndZip concats images for each echo into a single file
    concated = [];
    for iEcho = echoes
        File = load_untouch_nii(fileLocator.GetBasicType(iEcho, imgType));
        concated = cat(4, concated, File.img);
    end

    File.img = concated;
    File.hdr.dime.dim(5) = length(echoes);
    save_untouch_nii(File, saveTo);
end


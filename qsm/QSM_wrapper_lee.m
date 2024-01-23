% QSM_processing 
% Runs QSM processing via SEPIA. This code simply prepared data for this
% Written by Jingwen
% Modified by Melanie
% Reworked by Lee

% To use, call like so
% matlab QSM_wrapper_lee --input-directory /path/to/my/dir/
% where that /dir/ has a subfolder containing dicoms, like /path/to/my/dir/dicoms/


clear all


%% INPUT PARSING

% Ensure the correct number of arguments is provided
if numel(varargin) ~= 2 || ~strcmpi(varargin{1}, '--input-directory')
    error('Usage: QSM_wrapper_lee --input-directory <directory_path>');
end

% Parse command line arguments
inputDirectoryIndex = 2;

if inputDirectoryIndex > numel(varargin)
    error('Error: Missing or invalid input directory argument');
end

input_data_path = varargin{inputDirectoryIndex};

%% Settings
%input_data_path = '/data/morrison/wip/lee/PDa447/';% '/Users/lee/data/pda440';
output_data_path = [input_data_path, 'processed/']; %'/Users/lee/data/pda440/processed/';
correctFilter = false;
philipsTrueGEFalse = false;
expectRealImaginary = false;

loc_dcm2niix = 'dcm2niix';% '/opt/homebrew/bin/dcm2niix';
loc_fsl = '/netopt/rhel7/fsl/bin/';%/Users/lee/binaries/fsl/share/fsl/bin/';
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

if philipsTrueGEFalse
    error NOT IMPLEMENTED SEE FUNCTION BELOW
    error ALSO CHECK PHASE IMAGES - THEY SEEM ALREADY UNWRAPPED OR SOMETHING??
    CreatePhaseMag_Philips

else
    if expectRealImaginary
        %% Step 2 Correct offset for Re/Im images and overwrite
        CorrectOffsetForReImAndCreateNiftis(input_data_path, myEchos, correctFilter, loc_dcm2niix);
        
        %% Step 3 Create Phase from Real + Im
        CreatePhase(input_data_path, myEchos);
    else
        ConvertPhaseMag(input_data_path, myEchos, loc_dcm2niix);
    end

    %% Combine mag & phase echoTimes
    fileLocator = FileLocator(input_data_path);
    ConcatImages(fileLocator, 'mag',  myEchos, fileLocator.GetMagnitude_AllEchos());
    ConcatImages(fileLocator, 'phase',  myEchos, fileLocator.GetPhase_AllEchos());

    
end



%% Run Sepia
RunSepia(fileLocator, myEchos, input_data_path, output_data_path)

%% (optional) convert to dicom
%SendToPACS()


function ConcatImages(fileLocator, imgType, echoTimes, saveTo)
    %ConcatImagesAndZip concats images for each echo into a single file
    concated = [];
    for iEcho = echoTimes
        File = load_untouch_nii(fileLocator.GetBasicType(iEcho, imgType));
        concated = cat(4, concated, File.img);
    end

    File.img = concated;
    File.hdr.dime.dim(5) = length(echoTimes);
    save_untouch_nii(File, saveTo);
end




function SendToPACS()
    error("NOT IMPLEMENTED")
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
end
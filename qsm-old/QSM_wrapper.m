% QSM_processing 
% Written by Jingwen
% Modified by Melanie

%% add path

warning('off');
addpath(genpath('/data/morrison/scripts/qsm/jingwen_code'));
addpath('/data/morrison/scripts/matlab_snippets/nifti_tools');
addpath('/data/morrison/scripts/qsm');

%% set up data paths

input_data_path = '/data/morrison/data/parkinsons/retro_clin/PDa420_no.consent.yet-addpost/qsm';
ID1 = extractBefore(input_data_path,'/qsm');
ptid = extractAfter(ID1,'clin/');
ptid = extractBefore(ptid,'_no.consent.yet-addpost');

output_data_path = [input_data_path '/processed_QSM'];
cd(input_data_path);


%% Remove phase chop in Im/Re images

%Step 1: sort dicoms

cd('qsm_dicoms');
system('mkdir echo1 echo2 echo3 echo4 echo5');
myDir = uigetdir; %gets directory
myFiles = dir(fullfile(myDir,'*.dcm')); %gets all wav files in struct
cd(myDir)
for k = 1:length(myFiles)
    baseFileName = myFiles(k).name;
    %fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', baseFileName);
    struct = dicominfo(baseFileName);
    if struct.EchoNumbers==1
        system(['mv ' baseFileName ' ../echo1']);
    elseif struct.EchoNumbers==2
        system(['mv ' baseFileName ' ../echo2']);
    elseif struct.EchoNumbers==3
        system(['mv ' baseFileName ' ../echo3']);
    elseif struct.EchoNumbers==4
        system(['mv ' baseFileName ' ../echo4']);
    elseif struct.EchoNumbers==5
        system(['mv ' baseFileName ' ../echo5']);
    end
end

cd ..
system(['rm -r ' myDir]); %removed 1.?? empty folder
myEchos = dir(fullfile(pwd,'e*'));

for i= 1:length(myEchos)
cd(myEchos(i).name)
system('mkdir real imaginary mag');
myFiles = dir(fullfile(pwd,'*.dcm')); %gets all wav files in struct
    for k = 1:length(myFiles)
        baseFileName = myFiles(k).name;
        %fullFileName = fullfile(myFolder, baseFileName);
        fprintf(1, 'Now reading %s\n', baseFileName);
        struct = dicominfo(baseFileName);

        % GE stores this information as a signed short of the 
        % Private Image Type(0043,102F) tag. 
        % The values 0, 1, 2, 3 correspond to magnitude, phase, real, and imaginary (respectively).
        if struct.Private_0043_102f==0 %mag
            system(['mv ' baseFileName ' mag']);
        elseif struct.Private_0043_102f==2 %real
            system(['mv ' baseFileName ' real']);
        elseif struct.Private_0043_102f==3 %imaginary
            system(['mv ' baseFileName ' imaginary']);
        end
    end
cd ..
end


%Step 2 Correct offset for Re/Im images and overwrite

for i= 1:length(myEchos)
    cd(myEchos(i).name)

    cd('real')
    myFiles = dir(fullfile(pwd,'*.dcm')); %gets all wav files in struct
    for k=1:length(myFiles)
    real1=dicomread(myFiles(k).name);real1_offset=real1(1,1);
    real1_offset
    real1=real1-real1_offset;
    metadata = dicominfo(myFiles(k).name);
    dicomwrite(real1,myFiles(k).name, metadata, 'WritePrivate', true, 'CreateMode', 'Copy', 'UseMetadataBitDepths', true) 
    end
    system('dcm2niix *')
    system(['mv real*.nii ' ptid '_qsm_e' num2str(i) '_real_fixed.nii']);
    system(['mv real*.json ' ptid '_qsm_e' num2str(i) '_real_fixed.json']);

    cd .. %back to echo level

    cd('imaginary')
    myFiles = dir(fullfile(pwd,'*.dcm')); %gets all wav files in struct
    for k=1:length(myFiles)
    imag1=dicomread(myFiles(k).name);imag1_offset=imag1(1,1);
    imag1_offset
    imag1=imag1-imag1_offset;
    metadata = dicominfo(myFiles(k).name); %09/19/23: will have to store metadata
    dicomwrite(imag1,myFiles(k).name, metadata, 'WritePrivate', true, 'CreateMode', 'Copy', 'UseMetadataBitDepths', true) 
    end
    system('dcm2niix *')
    system(['mv imag*.nii ' ptid '_qsm_e' num2str(i) '_imaginary_fixed.nii']);
    system(['mv imag*.json ' ptid '_qsm_e' num2str(i) '_imaginary_fixed.json']);

    cd .. %back to echo level
    cd .. %back to all echoes level

end


%Step 3 Create Phase from Real + Im

for i= 1:length(myEchos)
    cd(myEchos(i).name)
    cd('real')
    File = load_untouch_nii([ptid '_qsm_e' num2str(i) '_real_fixed.nii']);
    Re = File.img;
    cd ..
    cd('imaginary')
    File = load_untouch_nii([ptid '_qsm_e' num2str(i) '_imaginary_fixed.nii']);
    Im = File.img;

    cd ../../../ %should be at the qsm level

    complex1=complex(double(Re),double(Im));
    for m=1:2:size(complex1,3)
    complex1(:,:,m) = -1*complex1(:,:,m); 
    end
    phase=angle(complex1);
    File.img = double(phase);
    File.hdr.dime.cal_max= pi;
    File.hdr.dime.cal_min=-pi;
    File.hdr.dime.datatype=16;
    File.hdr.dime.bitpix=16;
    save_untouch_nii(File,[ptid '_qsm_e' num2str(i) '_phase_fixed.nii'])
    cd('qsm_dicoms')
end




%% Combine mag & phase echoes

File = load_untouch_nii([ptid '_qsm_e1.nii']);
Mag1 = File.img;
File = load_untouch_nii([ptid '_qsm_e2.nii']);
Mag2 = File.img;
File = load_untouch_nii([ptid '_qsm_e3.nii']);
Mag3 = File.img;
File = load_untouch_nii([ptid '_qsm_e4.nii']);
Mag4 = File.img;
File = load_untouch_nii([ptid '_qsm_e5.nii']);
Mag5 = File.img;

Magni = cat(4, Mag1, Mag2, Mag3, Mag4, Mag5);
File.img = Magni;
File.hdr.dime.dim(5) = 5;
save_untouch_nii(File, 'Magni');

File = load_untouch_nii([ptid '_qsm_e1_phase_fixed.nii']);
Ph1 = File.img;
File = load_untouch_nii([ptid '_qsm_e2_phase_fixed.nii']);
Ph2 = File.img;
File = load_untouch_nii([ptid '_qsm_e3_phase_fixed.nii']);
Ph3 = File.img;
File = load_untouch_nii([ptid '_qsm_e4_phase_fixed.nii']);
Ph4 = File.img;
File = load_untouch_nii([ptid '_qsm_e5_phase_fixed.nii']);
Ph5 = File.img;

Phase = cat(4, Ph1, Ph2, Ph3, Ph4, Ph5);
File.img = Phase;
File.hdr.dime.dim(5) = 5;
save_untouch_nii(File, 'Phase');


system('gzip Magni.nii');
system('gzip Phase.nii');

%% Prepare header file for processing

% set up echo times (0018,0081) units=e-3 seconds

echoes = zeros(1,5); %create empty vector to enter echos

for echo = 1:5
fname = [ptid '_qsm_e' num2str(echo) '.json']; 
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);
echoes(:,echo)= val.EchoTime;
end

te1 = echoes(1,1);
te2 = echoes(1,2);
te3 = echoes(1,3);
te4 = echoes(1,3);
te5 = echoes(1,3);

save('echoes.mat', 'echoes')

% te1 = 5.44e-3; % s 
% te2 = 12.0e-3; % s 
% te3 = 18.6e-3; % s 
% te4 = 25.2e-3;% 


% set up dimensions
% spatial resolution of the data, in mm
header.voxelSize = [1,1,1];  
% image matrix size
%header.matrixSize = [224,224,60]; %for Phillips clinical QSM
header.matrixSize = [File.hdr.dime.dim(2) File.hdr.dime.dim(3) File.hdr.dime.dim(4)];
% set up parameters
header.b0 = 3;                  % magnetic field strength, in Tesla
header.b0dir = [0;0;1];         % main magnetic field direction, [x,y,z]
header.CF = header.b0*42.58*1e6;       % imaging frequency, in Hz (B0*gyromagnetic_ratio)
header.te = [te1,te2,te3];  % echo time for each GRE image, in second
header.delta_TE = te2-te1;      % echo spacing, in second

% save to mat file
save([input_data_path '/sepia_header.mat'],'header');

%% set up algorithm parameters

%%% Please make sure the file names are the same as used below %%%%%%%%%%%%
input.magnitudeFile     = [input_data_path '/Magni.nii'];
input.phaseFile         = [input_data_path '/Phase.nii'];
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

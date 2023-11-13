% QSM_processing 
% Written by Jingwen
% Modified by Melanie

%% add path

warning('off');
addpath(genpath('/home/lreid/qsm/qsm/jingwen_code'));
addpath('/data/morrison/scripts/matlab_snippets/nifti_tools');
addpath('/home/lreid/qsm');

%% set up data paths

dir_top = '/data/morrison/wip/lee/nov6/PDa434_no.consent.yet-addpost/'; % uigedir % must end with /s
dir_dicoms_in = strcat(dir_top, 'qsm_dicoms/brain_ax_3d_swi_v1/');
ID1 = extractBefore(dir_dicoms_in,'qsm');
ptid = extractAfter(ID1,'clin/');
ptid = extractBefore(ptid,'_no.consent.yet-addpost');

dir_out_top = [dir_top '/processed_QSM/'];
dir_organised_dicoms = strcat(dir_out_top, 'raw/');
dir_offset_fixed = strcat(dir_out_top, 'offset_fixed/');
dir_phase_mag = strcat(dir_out_top, 'phase_mag/');

loc_phase = strcat(dir_phase_mag, 'Phase.nii.gz');
loc_magnitude = strcat(dir_phase_mag, 'Magni.nii.gz');

%% Execution options
noEchoesExpected = 10;

mkdir(dir_out_top)
cd(dir_dicoms_in);

SortDicoms(dir_dicoms_in, dir_organised_dicoms, noEchoesExpected)

locs_phase = [];
locs_mag = [];
for iEcho=1:noEchoesExpected
    locs_phase = [locs_phase, CreatePhase(dir_organised_dicoms, dir_out_top, iEcho)];
    locs_mag = [locs_mag, DicomToNiiFromType('mag', dir_organised_dicoms, dir_out_top, iEcho)];
end

ConcatNiftis(locs_mag, loc_phase);
ConcatNiftis(locs_phase, loc_magnitude);

% %% Combine mag & phase echoes

% File = load_untouch_nii([ptid '_qsm_e1.nii']);
% Mag1 = File.img;
% File = load_untouch_nii([ptid '_qsm_e2.nii']);
% Mag2 = File.img;
% File = load_untouch_nii([ptid '_qsm_e3.nii']);
% Mag3 = File.img;
% File = load_untouch_nii([ptid '_qsm_e4.nii']);
% Mag4 = File.img;
% File = load_untouch_nii([ptid '_qsm_e5.nii']);
% Mag5 = File.img;

% Magni = cat(4, Mag1, Mag2, Mag3, Mag4, Mag5);
% File.img = Magni;
% File.hdr.dime.dim(5) = 5;
% save_untouch_nii(File, 'Magni');

% File = load_untouch_nii([ptid '_qsm_e1_phase_fixed.nii']);
% Ph1 = File.img;
% File = load_untouch_nii([ptid '_qsm_e2_phase_fixed.nii']);
% Ph2 = File.img;
% File = load_untouch_nii([ptid '_qsm_e3_phase_fixed.nii']);
% Ph3 = File.img;
% File = load_untouch_nii([ptid '_qsm_e4_phase_fixed.nii']);
% Ph4 = File.img;
% File = load_untouch_nii([ptid '_qsm_e5_phase_fixed.nii']);
% Ph5 = File.img;

% Phase = cat(4, Ph1, Ph2, Ph3, Ph4, Ph5);
% File.img = Phase;
% File.hdr.dime.dim(5) = 5;
% save_untouch_nii(File, 'Phase');


% system('gzip Magni.nii');
% system('gzip Phase.nii');

return
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
save([dir_dicoms_in '/sepia_header.mat'],'header');

%% set up algorithm parameters

%%% Please make sure the file names are the same as used below %%%%%%%%%%%%
input.magnitudeFile     = [dir_dicoms_in '/Magni.nii'];
input.phaseFile         = [dir_dicoms_in '/Phase.nii'];
%%% Check name - end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.headerFile        = [dir_dicoms_in '/sepia_header.mat'];

opts.writeLog           = 1;
opts.isGPU              = 0;
%opts.BETmethod          = 'HD-BET';
opts.BETmethod          = 'BET';
opts.iLSQR              = 1;
opts.QSMGAN             = 0;
opts.All                = 0;

%% run

QSMfile = dir(dir_out_top);
QSMfile([QSMfile.isdir]) = [];
QSMfile(cellfun(@isempty,strfind({QSMfile.name},'QSM_'))) = [];

%if (opts.All && length(QSMfile) < 10) || (~opts.All && length(QSMfile) < 2)
    [Tcomp] = QSM_processing_CPU(input, dir_out_top, opts);
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




function SortDicoms(input_data_dir, dir_out_top, noEchoesExpected)
    % Organises dicoms by echo number
    for i = 1:noEchoesExpected
        mkdir([GetEchoDir(dir_out_top,i, 'real')])
        mkdir([GetEchoDir(dir_out_top,i, 'imaginary')])
        mkdir([GetEchoDir(dir_out_top,i, 'mag')])
    end

    myFiles = dir(fullfile(input_data_dir,'*.dcm')); %gets all dicom files in dir

    if isempty(myFiles)
        error(strcat("No files found in ", input_data_dir))
    end

    parfor k = 1:length(myFiles)
        CopyDicomIfNotFound(string(myFiles(k).name), dir_out_top);
    end
end

function echoNum = CopyDicomIfNotFound(baseFileName, dir_out_top)
        fprintf(1, 'Now reading %s\n', baseFileName);
        dicomHeader = dicominfo(baseFileName);
        echoNum = dicomHeader.EchoNumbers;
        copyTo = GetEchoDir(dir_out_top, echoNum, GetGREImageType(dicomHeader));

        if ~isfile(copyTo)
            copyfile(baseFileName, copyTo, 'f');  
        end

        function type = GetGREImageType(dicomHeader)
            % GE stores this information as a signed short of the 
            % Private Image Type(0043,102F) tag. 
            % The values 0, 1, 2, 3 correspond to magnitude, phase, real, and imaginary (respectively).
            if dicomHeader.Private_0043_102f==0 %mag
                type = 'mag';
            elseif dicomHeader.Private_0043_102f==2 %real
                type = 'real';
            elseif dicomHeader.Private_0043_102f==3 %imaginary
                type = 'imaginary';
            else
                error(['Unexpected code for tag (0043,102F): ', string(dicomHeader.Private_0043_102f)])
            end
        end

    end



function directory = GetEchoDir(dir_out_top, echoNo, type)
    % Returns the directory holding dicoms for a specific echo
    % type: provide real/imaginary/mag as a string
    directory = strcat(dir_out_top,  'echo', string(echoNo), '/', type, '/');
end


function OffsetIntensitiesToTopLeftPixel(locFrom, locTo)
    % Threadsafe so long as no I/O touching the files at hand
    real1 = dicomread(locFrom);
    real1_offset = real1(1,1);
    real1 = real1 - real1_offset;
    metadata = dicominfo(locFrom);

    dicomwrite(real1, locTo, metadata, 'WritePrivate', true, 'CreateMode', 'Copy', 'UseMetadataBitDepths', true) 
end

function FixOffset(patientID, dir_in_top, dir_out_top, echoNumber, type)


    dirIn = GetEchoDir(dir_in_top, echoNumber, type);
    dirOut = GetEchoDir(dir_out_top, echoNumber, type);

    myFiles = dir(fullfile(strcat(dirIn, '*.dcm')));

    parfor k=1:length(myFiles)
        fn = myFiles(k).name
        locTo = strcat(dirOut, fn);
        if ~isfile(locTo)
            locFrom = strcat(dirIn, fn);
            OffsetIntensitiesToTopLeftPixel(locFrom, locTo)
        end
    end

end

function loc_nii = DicomToNiiFromType(type, dir_in_top, dir_out_top, echoNumber)
    echoNumber = string(echoNumber);
    dir_dicoms = GetEchoDir(dir_in_top, echoNumber, type);
    loc_nii = GetLocRawNii(dir_out_top, echoNumber, type);
    DicomToNifti(dir_dicoms, loc_nii, type);    
end

function DicomToNifti(dir_dicoms, loc_nii, type)
    cwd = pwd();
    cd(dir_dicoms)
    system('dcm2niix *')

    for loc = dir(dir_dicoms)
        if endsWith(loc, ".nii") && contains(loc, strcat(type, "*.nii"))
            movefile(loc, loc_nii)
            break
        end
    end


    %display(strcat("mv ", type,"*.nii ", loc_nii))
    %system(strcat("mv ", type,"*.nii ", loc_nii));
    %system(['mv ' type '*.json ' patientID '_qsm_e' num2str(echoNumber) '_' type '_fixed.json']);
    system('rm *.json');

    if ~isfile(loc_nii)
        error(strcat("Could not find converted file for type ", type, " in ", dir_dicoms))
    end    

    cd(cwd);
end

function loc = GetLocRawNii(dir_out_top, echoNumber, type)
    loc = strcat(dir_out_top, "echo", echoNumber, "_", type, ".nii");
end


function loc_save_phase = CreatePhase(dir_in_top, dir_out_top, echoNumber)
    %cd(myEchos(i).name)
    %cd('real')
    echoNumber = string(echoNumber);
    
    loc_real = DicomToNiiFromType('real', dir_in_top, dir_out_top, echoNumber);
    loc_imagin = DicomToNiiFromType('imaginary', dir_in_top, dir_out_top, echoNumber);


    File = ReadNifti(loc_real);
    Re = File.img;
    File = ReadNifti(loc_imagin);
    Im = File.img;

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
    loc_save_phase = GetLocRawNii(dir_out_top, echoNumber, 'phase');
    save_untouch_nii(File,convertStringsToChars(loc_save_phase))    
end


function images = ReadNiftis(locations)
    images = [];
    for loc = locations
        images = [images, ReadNifti(loc)];
    end
end


function nifti = ReadNifti(loc)
    % hack to get around nifti loading in dependency
    nifti = load_untouch_nii(convertStringsToChars(loc));
end

function File = ConcatNiftis(filenames, loc_saveTo)
    % Concatenates 3D niftis along the 4th dimension
    % Images assumed to have the same header information

    files = ReadNiftis(filenames);

    imgData = [];
    for file = files
        imgData = cat(4,imgData, file.img);
    end

    File = files(1);
    File.img = imgData;
    File.hdr.dime.dim(5) = 5;
    save_untouch_nii(File, loc_saveTo);

end

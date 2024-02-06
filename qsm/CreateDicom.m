function CreateDicom(dir_templateDicom, loc_nifti,dir_out)
    %CREATEDICOM Creates a series of dicoms

    FileStruct = load_untouch_nii(loc_nifti);
    img = double(FileStruct.img);
    imgDims = FileStruct.hdr.dime.pixdim(2:4); % first is dud in nii

    % Convert/Scale to uint16
    img = ScaleIntensities(img, -0.4, 0.4);

    LoadDicomDict()

    % Get a dicom UID
    uid = dicomuid();
    metadata_orig = GetFirstDicom(dir_templateDicom);

    mkdir(dir_out);
    numSlices = size(img,3);
    for i=1:numSlices
        CreateDicom(i, numSlices);
    end

    function img = ScaleIntensities(img, min, max)
        % ScaleIntensities Alters image intensities to 0 - 65535 (UInt16
        % max)

        % Crop intensity range
        img(img < min) = min; 
        img(img > max) = max;

        img = (img - min)/ (max-min);

        %img = img * 65535;
    end

    function CreateDicom(iSlice, numSlices)
        metadata.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
        metadata.ImageType = 'DERIVED\SECONDARY\OTHER';
        metadata.StudyDate = metadata_orig.StudyDate;
        metadata.SeriesDate = metadata_orig.SeriesDate;
        metadata.AcquisitionDate = metadata_orig.AcquisitionDate;
        metadata.ContentDate = metadata_orig.ContentDate;
        metadata.StudyTime = metadata_orig.StudyTime;
        metadata.SeriesTime = metadata_orig.SeriesTime;
        metadata.AcquisitionTime = metadata_orig.AcquisitionTime;
        metadata.ContentTime = metadata_orig.ContentTime;
        metadata.Modality = 'MR';


        metadata.SliceThickness = imgDims(3);
        metadata.PixelSpacing = imgDims(1:2);
        metadata.ImageOrientationPatient = metadata_orig.ImageOrientationPatient;
        metadata.ImagePositionPatient = metadata_orig.ImagePositionPatient;
        metadata.FrameOfReferenceUID = metadata_orig.FrameOfReferenceUID;
        metadata.PatientIdentityRemoved = 'NO';
        %metadata.WindowWidth = 3276; %1732
        %metadata.WindowCenter = 0; %1982
        metadata.StudyInstanceUID = metadata_orig.StudyInstanceUID;
        metadata.ProtocolName = 'QSM_Processed3D';
        metadata.SeriesDescription = 'QSM_Processed3D';
        metadata.SeriesInstanceUID = uid;
        metadata.SpacingBetweenSlices = metadata_orig.SpacingBetweenSlices;

        metadata.InstanceNumber = iSlice;
        metadata.NumberOfSlices = numSlices;
        metadata.SliceLocation = iSlice * imgDims(3);
        metadata.SeriesNumber = metadata_orig.SeriesNumber * 100 + 1;

        filename = [dir_out '/' sprintf('IM%04i.dcm',i)];

        % VOXELS
        % Be careful not to LR flip here
        % Niftis are stored with a different data order to dicoms so 
        % the raw data must be flipped. There's also a rot90 here though
        % I'm not sure on its cause as that was inherited
        voxels = fliplr(rot90(img(:,:,iSlice))); % Check for flipping if altering this

        dicomwrite(voxels, filename, metadata, 'WritePrivate', true, 'CreateMode', 'Copy', 'UseMetadataBitDepths', true)
    end

end

function template = GetFirstDicom(dir)
    allDicoms = ls(dir);
    loc_templateDicom = allDicoms(1);
    template = dicominfo(loc_templateDicom);
end


function LoadDicomDict()
    dir = fileparts(mfilename('fullpath'));
    dicomdict('set', [dir '/dicom-dict.txt']);
end

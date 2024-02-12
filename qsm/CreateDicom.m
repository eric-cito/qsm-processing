function CreateDicom(dir_templateDicom, loc_nifti,dir_out)
    %CREATEDICOM Creates a series of dicoms

    FileStruct = load_untouch_nii(loc_nifti);
    img = double(FileStruct.img);
    imgDims = FileStruct.hdr.dime.pixdim(2:4); % first is dud in nii
    pixelDims = FileStruct.hdr.dime.imgdim(2:4);

    % Convert/Scale to uint16
    img = ScaleIntensities(img, -0.4, 0.4);

    LoadDicomDict()

    % Create shared vars
    uid = dicomuid();
    today = datetime('now','format', 'yyyyMMdd');
    nowTime = datetime('now','format', 'HHmmss.00');
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

        img = (img - min) / (max-min);

        img = uint16(img * 65535);
    end

    function CreateDicom(iSlice, numSlices)
        % -- Patient 
        metadata.PatientIdentityRemoved = 'NO';
        % -- Clinical Trial Subject not included

        CopyFieldIfFound(metadata_orig, metadata, "StudyID") % Not req   
        % -- General Image
        metadata.ImageType = 'DERIVED\SECONDARY\OTHER'; % Not Req
        CopyFieldIfFound(metadata_orig, metadata, "AcquisitionDate"); % Not Req
        metadata.ContentDate = today;
        metadata.ContentTime = nowTime;

        % -- Image Plane
        metadata.SliceLocation = metadata_orig.SliceLocation + (iSlice-1) * pixelDims(3); % not req
        metadata.PixelSpacing = pixelDims(1:2);


        % -- Image Pixel
        metadata.SamplesPerPixel = 1; % 3 if RGB
        metadata.PhotometricInterpretation = 'MONOCHROME2';
        % Planar configuration is req if RGB
        metadata.Rows = imgDims(1);
        metadata.Cols = imgDims(2);
        metadata.PixelAspectRatio = pixelDims(1:2);
        metadata.BitsAllocated = 16;
        metadata.BitsStored = 16;
        metadata.HighBit = 15;
        metadata.PixelRepresentation = 0; % unsigned

        % -- MR Image
        CopyFieldIfFound(metadata_orig, metadata, "RepetitionTime");
        CopyFieldIfFound(metadata_orig, metadata, "InversionTime");

        % -- Overlay plane 
        % -- This has required values but this seems an error? If there's
        % no overlay I dont' know what should be here

        % -- VOI LUT
        metadata.WindowCenter = 0;
        metadata.WindowWidth = 65535;

        % -- SOP Common
        % Specific characted set is presumed set by the binary writer
        metadata.InstanceCreationDate = today; % Not req
        metadata.InstanceCreationTime = nowTime; % Not req
        metadata.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
        metadata.SOPInstanceUID = UP TO HERE ON SOP COMMON UNDER MR here: https://dicom.innolitics.com/ciods/mr-image/sop-common/00080018




        required = [

        % -- Patient
        "PatientName";
        "PatientID";
        "PatientBirthDate";
        "PatientSex";
        % For children, responsible person role is req if ResponsiblePerson
        % is set. I've assumed it is not

        
        % -- Clinical Trial Subject not included
        

        % -- General Study
        "StudyInstanceUID";
        "SeriesDate";
        "SeriesTime";
        "AccessionNumber";
        "ReferringPhysicianName";

        % -- Patient Study has no req tags
        % -- Clinical Trial Study not included
        % -- Clinical Trial Series not included
        % -- Frame of Reference
        "FrameOfReferenceUID";
        "PositionReferenceIndicator";
        
        % -- General Equipment
        "Manufacturer";

        % -- General Image
        "PatientOrientation"; % likely empty

        % -- Image Plane
        "ImagePositionPatient";
        "ImageOrientationPatient";

        % -- Image Pixel req all set above

        % -- Container Identifier
        "ContainerIdentifier";
        "IssuerOfTheContainerIdentifierSequence";% Req, empty if unknown
        "ContainerTypeCodeSequence"; % Req, empty if unknown
        "SpecimenDescriptionSequence";

        % -- MR Image
        "ScanningSequence";
        "SequenceVariant"; 
        "ScanOptions"; % Req, empty if unknown
        "MRAcquisitionType"; % Req, empty if unknown
        "EchoTime";
        "EchoTrainLength";
        

        


         

        % metadata. = metadata_orig.
        % metadata. = metadata_orig.
        % metadata. = metadata_orig.



        metadata.StudyDate = metadata_orig.StudyDate;
        metadata.StudyTime = metadata_orig.StudyTime;
        
        metadata.AcquisitionTime = metadata_orig.AcquisitionTime;
        metadata.ContentDate = metadata_orig.ContentDate;
        metadata.ContentTime = metadata_orig.ContentTime;

        

        
        
        metadata.PatientSex = metadata_orig.PatientSex;
        
        metadata.StudyNumber = metadata_orig.StudyNumber;
        metadata.PatientOrientation = metadata_orig.PatientOrientation;

        
        metadata.Modality = 'MR';

        metadata.SliceThickness = pixelDims(3);
        metadata.SpacingBetweenSlices = pixelDims(3);
        
        metadata.ImageOrientationPatient = metadata_orig.ImageOrientationPatient;
        metadata.ImagePositionPatient = metadata_orig.ImagePositionPatient;
        metadata.FrameOfReferenceUID = metadata_orig.FrameOfReferenceUID;
        
        
        %metadata.WindowWidth = 3276; %1732
        %metadata.WindowCenter = 0; %1982
        metadata.StudyInstanceUID = metadata_orig.StudyInstanceUID;
        metadata.ProtocolName = 'QSM_Processed3D';
        metadata.SeriesDescription = 'QSM_Processed3D';
        metadata.SeriesInstanceUID = uid;
        metadata.SpacingBetweenSlices = metadata_orig.SpacingBetweenSlices;

        metadata.InstanceNumber = iSlice;
        metadata.NumberOfSlices = numSlices;
        
        metadata.SeriesNumber = metadata_orig.SeriesNumber * 1;
        metadata.InStackPositionNumber = iSlice;
        metadata.StackID = metadata_orig.StackID;

        filename = [dir_out '/' sprintf('IM%04i.dcm',i)];

        % VOXELS
        % Be careful not to LR flip here
        % Niftis are stored with a different data order to dicoms so 
        % the raw data must be flipped. There's also a rot90 here though
        % I'm not sure on its cause as that was inherited
        voxels = fliplr(rot90(img(:,:,iSlice))); % Check for flipping if altering this

        dicomwrite(voxels, filename, metadata, 'WritePrivate', true, 'CreateMode', 'Copy', 'UseMetadataBitDepths', true)
    end

    function CopyFieldOrDefault(src, dest, field, valIfNotFound)
        if isfield(dest, field)
            val = src.getfield(field);
        else
            val = valIfNotFound;
        end

        dest.setfield(field, val)
    end

    function CopyFieldIfFound(src, dest, field)
        if isfield(dest, field)
            dest.setfield(field, src.getfield(field))
        end
    end

end

function template = GetFirstDicom(dicomDir)
    allDicoms = [dir(fullfile(dicomDir,'**/*.DCM')); dir(fullfile(dicomDir,'**/*.dcm'))];
    loc_templateDicom = [allDicoms(1).folder '/' allDicoms(1).name];
    template = dicominfo(loc_templateDicom);
end


function LoadDicomDict()
    dir = fileparts(mfilename('fullpath'));
    dicomdict('set', [dir '/dicom-dict.txt']);
end

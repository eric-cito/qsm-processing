function CreateDicom(loc_templateDicom, loc_nifti,dir_out)
    %CREATEDICOM Creates a series of dicoms
    FileStruct = load_untouch_nii(loc_nifti);
    img = double(FileStruct.img);
    imgDims = FileStruct.hdr.dime.pixdim(2:4); % first is dud in nii

    % Convert/Scale to uint16
    img = ScaleIntensities(img, -0.3, 0.3);

    LoadDicomDict()

    % Get a dicom UID
    uid = dicomuid();
    metadata_orig = dicominfo(loc_templateDicom);

    mkdir(dir_out);
    numSlices = size(img,3);
    for i=1:numSlices
        CreateDicom(i, numSlices);
    end

    function img = ScaleIntensities(img, min, max)
        % ScaleIntensities Alters image intensities to 0 - 65535 (UInt16
        % max)

        % Crop intensity range
        img(img < min) = -min; 
        img(img > max) = max;

        img = (img - min)/ (max-min);

        img = img * 65535;
    end

    function CreateDicom(iSlice, numSlices)
        metadata.SliceThickness = imgDims(3);
        metadata.PixelSpacing = imgDims(1:2);
        metadata.ImageOrientationPatient = metadata_orig.ImageOrientationPatient;
        metadata.PatientIdentityRemoved = 'NO';
        %metadata.WindowWidth = 3276; %1732
        %metadata.WindowCenter = 0; %1982
        metadata.ProtocolName = 'QSM_Processed3D';
        metadata.SeriesDescription = 'QSM_Processed3D';
        metadata.SeriesInstanceUID = uid;
        metadata.SpacingBetweenSlices = metadata_orig.SpacingBetweenSlices;

        metadata.InstanceNumber = iSlice;
        metadata.NumberOfSlices = numSlices;
        metadata.SliceLocation = iSlice * imgDims(3);

        filename = sprintf('IM%04i.dcm',i);
        dicomwrite(img(:,:,i), filename, metadata, 'WritePrivate', true, 'CreateMode', 'Copy', 'UseMetadataBitDepths', true)
    end

end


function LoadDicomDict()
    dir = fileparts(mfilename('fullpath')).Directory;
    dicomdict('set', [dir 'dicom-dist.txt']);
end

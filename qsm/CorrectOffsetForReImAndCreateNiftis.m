function CorrectOffsetForReImAndCreateNiftis(topDir, myEchos, correctFilter, loc_dcm2niix)
    % WARNING: This is dangerous upside-down code that will overwrite dicoms only
    % to then convert them to nifti, rather than converting niftis and
    % correcting them in place. It should be refactored into two methods -
    % one which converts to nifti, and another that will correct them

    fileLocator = FileLocator(topDir);

    for iEcho= 1:length(myEchos)

        dir_dest = fileLocator.GetEchoDir(iEcho);
        if exist(dir_dest,"dir")
            % pre clean to avoid duplicates etc
            system(strcat("rm -r ", dir_dest));
        end

        mkdir(dir_dest)
        
        ConvertAndDefilter(iEcho, 'real')
        ConvertAndDefilter(iEcho, 'imaginary')
    end

    function ConvertAndDefilter(iEcho, imgType)
        dicomDir = fileLocator.GetDicomDir(iEcho, imgType);
        if correctFilter
            CorrectGEFilter(dicomDir);
        end

        system(strcat(loc_dcm2niix, " -o ", dir_dest, " ", ' -f "TEMP" ', " ", dicomDir));
        system(strcat("mv ", dir_dest, "TEMP*.nii ", fileLocator.GetBasicType(iEcho, imgType)));
        system(strcat("mv ", dir_dest, "TEMP*.json ", fileLocator.GetBasicTypeWithSuffix(iEcho, imgType, "json")));
    end
    
    function CorrectGEFilter(directory)
        % Attempts to remove the GE filter applied to dicoms, overwriting those
        % dicoms
        myFiles = dir(fullfile(directory,'*.dcm')); %gets all wav files in struct
        parfor k=1:length(myFiles)
            loc = [myFiles(k).folder '/' myFiles(k).name];
            real1 = dicomread(loc);
            real1_offset=real1(1,1);
            
            real1 = real1-real1_offset;
            metadata = dicominfo(loc);
            dicomwrite(real1, loc, metadata, 'WritePrivate', true, 'CreateMode', 'Copy', 'UseMetadataBitDepths', true) 
        end
    end

end

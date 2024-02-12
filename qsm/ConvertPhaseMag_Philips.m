function ConvertPhaseMag_Philips(topDir,myEchos, loc_dcm2niix)
    %ConvertPhaseMag Converts dicoms for phase and magnitude into nifti
    fileLocator = FileLocator(topDir);

    for iEcho= 1:length(myEchos)
        CreatePhaseForEcho(myEchos(iEcho))
    end

    function CreatePhaseForEcho(echo)
        Convert(echo, 'phase');
        Convert(echo, 'mag');
    end

    function Convert(iEcho, imgType)
        dicomDir = fileLocator.GetDicomDir(iEcho, imgType);
        dir_dest = fileLocator.GetEchoDir(iEcho);
        if ~exist("dir_dest","dir")
            mkdir(dir_dest)
        end

        system(strcat(loc_dcm2niix, " -o ", dir_dest, " ", ' -f "TEMP" ', " ", dicomDir));
        system(strcat("mv ", dir_dest, "TEMP*.nii ", fileLocator.GetBasicType(iEcho, imgType)));
        system(strcat("mv ", dir_dest, "TEMP*.json ", fileLocator.GetBasicTypeWithSuffix(iEcho, imgType, "json")));
    end
end
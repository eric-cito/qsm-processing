function ConvertPhaseMag_Philips(fileLocator,myEchos, loc_dcm2niix)
    %ConvertPhaseMag Converts dicoms for phase and magnitude into nifti

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

% % function myEchos = CreatePhaseMag_Philips(fileLocator, myEchos, loc_dcm2niix)
% % 
% % 
% %     for iEcho= 1:length(myEchos)
% %         CreatePhaseForEcho(myEchos(iEcho))
% %     end
% % 
% % 
% %         dirTmp = [tempname,'/']; mkdir(dirTmp); c =
% %         onCleanup(@()rmdir(dirTmp,'s')); % delete temp dir on function
% %         exit
% % 
% % 
% % 
% %         File.img = double(phase); save_untouch_nii(File,
% %         fileLocator.GetPhase(echo))
% % 
% %         File.img = double(mag); save_untouch_nii(File,
% %         fileLocator.GetMagnitude(echo))
% % 
% %         system(["dcm2niix", dir_dicoms, dirTmp]) system(["mv", dirTmp
% %         "*_ph.nii*", fileLocator.GetPhase_AllEchos()]) system(["mv",
% %         dirTmp "*.nii*", fileLocator.GetPhase_AllEchos()])
% % 
% % 
% % end
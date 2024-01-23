function myEchos = CreatePhaseMag_Philips(fileLocator, dir_dicoms)

    error NOT IMPLEMENTED AS ONLY ONE ECHO EXISTS IN CURRENT DATA


        dirTmp = [tempname,'/'];
        mkdir(dirTmp);
        c = onCleanup(@()rmdir(dirTmp,'s')); % delete temp dir on function exit

        

        File.img = double(phase);
        save_untouch_nii(File, fileLocator.GetPhase(echo))
    
        File.img = double(mag);
        save_untouch_nii(File, fileLocator.GetMagnitude(echo))
        
        system(["dcm2niix", dir_dicoms, dirTmp])
        system(["mv", dirTmp "*_ph.nii*", fileLocator.GetPhase_AllEchos()])
        system(["mv", dirTmp "*.nii*", fileLocator.GetPhase_AllEchos()])
        

end
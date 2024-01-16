function CreatePhase(topDir,myEchos)
    fileLocator = FileLocator(topDir);

    for iEcho= 1:length(myEchos)
        CreatePhaseForEcho(myEchos(iEcho))
    end

    function CreatePhaseForEcho(echo)

        File = load_untouch_nii(fileLocator.GetReal(echo));
        Re = File.img;
        File = load_untouch_nii(fileLocator.GetImaginary(echo));
        Im = File.img;

        complex1=complex(double(Re),double(Im));
        % In GE, every second slice is effectively inverted due to phase wrap
        for m=1:2:size(complex1,3)
            complex1(:,:,m) = -1*complex1(:,:,m); 
        end

        phase = angle(complex1);
        mag = abs(complex1);
    
        File.hdr.dime.cal_max= pi;
        File.hdr.dime.cal_min=-pi;
        File.hdr.dime.datatype=16;
        File.hdr.dime.bitpix=16;
    
        File.img = double(phase);
        save_untouch_nii(File, fileLocator.GetPhase(echo))
    
        File.img = double(mag);
        save_untouch_nii(File, fileLocator.GetMagnitude(echo))
    end
end
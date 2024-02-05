function RunSepia(fileLocator, myEchos, input_data_path, output_data_path)
    %% Prepare header file for processing
    echoTimes = ReadEchoTimes(fileLocator, myEchos);
    save([input_data_path 'echoes.mat'], 'echoTimes')
    
    % set up dimensions
    % spatial resolution of the data, in mm
    File = load_untouch_nii(fileLocator.GetPhase(myEchos(1)));
    header.voxelSize = File.hdr.dime.pixdim([2,3,4]); % nifti - pixdim 2,3,4 have the dimensions. 
    % image matrix size
    header.matrixSize = File.hdr.dime.dim([2,3,4]);
    % set up parameters
    header.b0 = 3;                  % magnetic field strength, in Tesla
    header.b0dir = [0;0;1];         % main magnetic field direction, [x,y,z]
    header.CF = header.b0 * 42.58 * 1e6;       % imaging frequency, in Hz (B0*gyromagnetic_ratio)
    header.te = echoTimes;  % echo time for each GRE image, in seconds
    header.delta_TE = echoTimes(1,2) - echoTimes(1,1);      % echo spacing, in seconds
    
    % save to mat file
    save([input_data_path 'sepia_header.mat'],'header');
    
    %% set up algorithm parameters
    
    %%% Please make sure the file names are the same as used below %%%%%%%%%%%%
    input.magnitudeFile     = [char(fileLocator.GetMagnitude_AllEchos())];
    input.phaseFile         = [char(fileLocator.GetPhase_AllEchos())];
    %input.magnitudeFile     = '/data/morrison/wip/lee/pda403/mag_allEchos_denoised.nii';
    %input.phaseFile         = '/data/morrison/wip/lee/pda403/phase_allEchos_denoised.nii';
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
    if ~exist("output_data_path","dir")
        mkdir(output_data_path)
    end
    
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
end

function echoTimes = ReadEchoTimes(fileLocator, myEchos)
    echoTimes = zeros(1,length(myEchos)); %create empty vector to enter echos
    
    for iEcho = myEchos
        fname = fileLocator.GetBasicTypeWithSuffix(iEcho, 'real', 'json'); 
        if ~exist(fname, "file")
            fname = fileLocator.GetBasicTypeWithSuffix(iEcho, 'phase', 'json'); 
        end
        fid = fopen(fname); 
        raw = fread(fid,inf); 
        str = char(raw'); 
        fclose(fid); 
        val = jsondecode(str);
        echoTimes(:,iEcho)= val.EchoTime;
    end
end

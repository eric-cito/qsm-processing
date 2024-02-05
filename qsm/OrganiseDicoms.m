function allEchoes = OrganiseDicoms(topDir)
%ORGANISEDICOMS Organises dicoms by type (real/im/etc) and echo number
%   Only designed for QSM use. If topDir is empty, it will be removed

fileLocator = FileLocator(topDir);

myFiles = [dir(fullfile(topDir,'**/*.DCM')); dir(fullfile(topDir,'**/*.dcm'))];

disp("Reading Dicoms...")
[echoes, imgTypes] = ReadDicoms(myFiles);
disp("Reading Dicoms completed")

allEchoes = unique(echoes);
foundImgTypes = unique(imgTypes);

CreateDirectories(allEchoes, foundImgTypes);

disp("Moving Dicoms...")
MoveFiles()

disp("Cleaning up")
RemoveEmptyDirs(topDir)

disp("Dicom Organisation Complete")

function MoveFiles()
    %MoveFiles moves files from their original to their organised locations
    for k = 1:length(myFiles)
        f = myFiles(k);
        src = strcat(f.folder, "/", f.name); % important this is a true string or strcmp does not work properly
        destDir = fileLocator.GetDicomDir(echoes(k), imgTypes(k));
        dest = string(strcat(destDir, f.name));
        if ~strcmp(src, dest)
            movefile(src, dest);
        end
    end
end

function [echoes, imgTypes] = ReadDicoms(myFiles)
    % READDICOMS reads all dicoms and returns their echo number and imgType
    numFiles = length(myFiles);
    
    echoes = ones(1,numFiles);
    imgTypes = strings(1,numFiles);
    
    parfor index = 1:numFiles
        f = myFiles(index);
        [echoNum, imageType] = ReadDicom(strcat(f.folder, '/', f.name));
        echoes(index) = echoNum;
        imgTypes(index) = imageType;
    end
end


function dirs = CreateDirectories(allEchoes, foundImgTypes)
    %CreateDirectories creates required destination directories 
    dirs = [];

    for echo = allEchoes
        for imgType = foundImgTypes
            d = fileLocator.GetDicomDir(echo, imgType);
            system(strcat("mkdir -p ", d));

            dirs = [dirs d];
        end
    end    
end

end

function [echoNum, imageType] = ReadDicom(path)
    %ReadDicom Reads echo number and image type (mag/phase/real/imaginary)
    %from a dicom file on disk

    dcmInfo = dicominfo(path);
    
    echoNum = dcmInfo.EchoNumbers;

    if isfield(dcmInfo, 'Private_0043_102f')
        % GE stores this information as a signed short of the 
        % Private Image Type(0043,102F) tag. 
        % The values 0, 1, 2, 3 correspond to magnitude, phase, real, and imaginary (respectively).
        if dcmInfo.Private_0043_102f==0 %mag
            imageType = 'mag';
        elseif dcmInfo.Private_0043_102f==1
            imageType = 'phase';
        elseif dcmInfo.Private_0043_102f==2 %real
            imageType = 'real';
        elseif dcmInfo.Private_0043_102f==3 %imaginary
            imageType = 'imaginary';
        end

    else
        % Philips
        if contains(dcmInfo.ImageType, "PHASE MAP")
            imageType = 'phase';
        else
            imageType = 'mag';
        end
    end
end
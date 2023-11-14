%separates qsm dicoms into phase and magnitude
%written by melanie 05/02/23

system('mkdir echo1 echo2 echo3');
myDir = uigetdir; %gets directory
myFiles = dir(fullfile(myDir,'*.dcm')); %gets all wav files in struct
for k = 1:length(myFiles)
    baseFileName = myFiles(k).name;
    %fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', baseFileName);
    struct = dicominfo(baseFileName);
    if struct.EchoNumbers==1
        system(['mv ' baseFileName ' echo1']);
    elseif struct.EchoNumbers==2
        system(['mv ' baseFileName ' echo2']);
    elseif struct.EchoNumbers==3
        system(['mv ' baseFileName ' echo3']);
    end
end


system('mkdir mag phase');
myDir = uigetdir; %gets directory
myFiles = dir(fullfile(myDir,'*.dcm')); %gets all wav files in struct
for k = 1:length(myFiles)
    baseFileName = myFiles(k).name;
    %fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', baseFileName);
    struct = dicominfo(baseFileName);
    if struct.Private_0043_102f==0 %mag
        system(['mv ' baseFileName ' mag']);
    elseif struct.Private_0043_102f==1
        system(['mv ' baseFileName ' phase']);
    end
end

addpath('/data/morrison/scripts/matlab_snippets/nifti_tools/');
FileStruct = load_untouch_nii('phase_DBS_GAD_PRE-OP_20221114124910_3_e3_ph.nii');
Phase = FileStruct.img;
[x, y, z] = size(Phase);
new_phase=Phase;
new_phase(:,:,2:2:z)=Phase(:,:,2:2:z).*(-1);
FileStruct.img=new_phase;
save_untouch_nii(FileStruct, 'new_phase');


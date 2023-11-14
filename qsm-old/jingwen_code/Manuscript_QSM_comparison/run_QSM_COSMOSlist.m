% clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));
addpath(genpath('/working/lupolab/eason/DL_QSM/code_matlab/cosmos'));
addpath('/netopt/share/lib/local/brain/matlab/');
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));

%% subject list from COSMOS

mainPath = '/working/lupolab/eason/DL_QSM/data/7T_cosmos';
subjFolders = dir(mainPath);
subjFolders(~[subjFolders.isdir]) = [];
subjFolders(~cellfun(@isempty,strfind({subjFolders.name},'.'))) = [];
subjFolders(cellfun(@isempty,strfind({subjFolders.name},'volunteer'))) = [];

%% Loop through subjects

for ii = 1:length(subjFolders)
    
    fprintf('Processing %s \n', subjFolders(ii).name);
    
    subjPath = [mainPath '/' subjFolders(ii).name '/scan1'];
    
    try
        exam_root = [subjPath '/swan_qsm/COSMOSmask_allQSM'];
        mkdir(exam_root);
        
        input_data_path = [subjPath '/swan_qsm'];
        output_data_path = exam_root;
        pfile_path = [subjPath '/scan1_JY'];
        
        if exist([output_data_path '/QSM_COSMOS.nii.gz'],'file') == 0
            if exist([subjPath '/swan_qsm/HDBET_allQSM/Magni.nii.gz'],'file') ~= 0
                copyfile([subjPath '/swan_qsm/HDBET_allQSM/Magni.nii.gz'],[output_data_path '/Magni.nii.gz']);
                copyfile([subjPath '/swan_qsm/HDBET_allQSM/Phase.nii.gz'],[output_data_path '/Phase.nii.gz']);
                copyfile([subjPath '/swan_qsm/HDBET_allQSM/brain_mask.nii.gz'],[output_data_path '/brain_mask.nii.gz']);
                copyfile([subjPath '/swan_qsm/HDBET_allQSM/QSM_COSMOS.nii.gz'],[output_data_path '/QSM_COSMOS.nii.gz']);
                cd(input_data_path);
            else
                [all_TE, ~, ~] = QSM_recon('', 'scan1', input_data_path, pfile_path);
                cd(input_data_path);
                set_SEPIA_header;
                header.te = all_TE*1e-3;
                header.delta_TE = all_TE(2)*1e-3 - all_TE(1)*1e-3;
                save('sepia_header.mat','header');
            end
        end
        
        %% set up parameters
        
        input.magnitudeFile     = [input_data_path '/Magni.nii.gz'];
        input.phaseFile         = [input_data_path '/Phase.nii.gz'];
        input.headerFile        = [input_data_path '/sepia_header.mat'];
        
        opts.writeLog           = 1;
        opts.isGPU              = 1;
        opts.BETmethod          = 'HD-BET';
        opts.iLSQR              = 1;
        opts.QSMGAN             = 1;
        opts.All                = 1;
        
        %% copy COSMOS brain masks to the folder
        
        if exist([output_data_path '/brain_mask_HD.nii.gz'],'file') == 0
            for scanInd = 1:3
                cmd = sprintf('i2nii -o %s -z y %s', ...
                    output_data_path, ...
                    [mainPath '/' subjFolders(ii).name '/scan' num2str(scanInd) '/scan' num2str(scanInd) '_binarymask2.idf']);
                system(cmd);
            end
            % get transformation matrices
            load([mainPath '/' subjFolders(ii).name '/scan1_spm_tm.mat']);
            % register the brain masks and take the intersection
            brain_mask1 = read_idf_image([mainPath '/' subjFolders(ii).name '/scan1/scan1_binarymask2']);
            brain_mask2 = read_idf_image([mainPath '/' subjFolders(ii).name '/scan2/scan2_binarymask2']);
            brain_mask3 = read_idf_image([mainPath '/' subjFolders(ii).name '/scan3/scan3_binarymask2']);
            brain_mask1 = permute(brain_mask1.img, [2,1,3]);
            brain_mask2 = imwarp_spm(brain_mask1, permute(brain_mask2.img, [2,1,3]), tm21);
            brain_mask3 = imwarp_spm(brain_mask1, permute(brain_mask3.img, [2,1,3]), tm31);
            brain_mask = brain_mask1 .* brain_mask2 .* brain_mask3;
            brain_mask(isnan(brain_mask)) = 0; % remove nan caused by interp3
            brain_mask = brain_mask > 0.5; % binarize the transformed mask
            
            % save as nii
            nii = load_untouch_nii([output_data_path '/brain_mask.nii.gz']);
            nii.img = permute(brain_mask,[2,1,3]);
            save_untouch_nii(nii,[output_data_path '/brain_mask_COSMOS.nii.gz']);
            
            % change to float
            cmd = sprintf('3dcalc -a %s -datum float -expr a -prefix %s', ...
                [output_data_path '/brain_mask_COSMOS.nii.gz'], ...
                [output_data_path '/brain_mask_HD.nii.gz']);
            system(cmd);
            
        end
        
        %% run
        
        QSMfile = dir(output_data_path);
        QSMfile([QSMfile.isdir]) = [];
        QSMfile(cellfun(@isempty,strfind({QSMfile.name},'QSM_'))) = [];
        
        system(['rm -r ' output_data_path '/iQSM2']);
        system(['rm ' output_data_path '/QSM_iQSM2_meanEcho.nii.gz']);
        if 1 % (opts.All && length(QSMfile) < 16) || (~opts.All && length(QSMfile) < 2)
            [Tcomp] = QSM_processing(input, output_data_path, opts);
        else
            fprintf(' >> Already have all QSM maps, skip QSM processing \n');
        end
    catch
        disp('!!! Something went wrong!');
    end
    
    %% plot images
    QSM_file = [output_data_path '/QSM_iQSM2_meanEcho.nii.gz'];
    nii = load_nii(QSM_file);
    
    figure;
    montage(nii.img(:,:,10:20:end),'DisplayRange',[-0.15 0.15]);
    drawnow;
end
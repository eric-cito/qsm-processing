clear; clc;
warning('off');

%% Add path

addpath('/home/jyao3/030_QSM/01_Code/Manuscript_BGQSM_HD');
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

%% read in list
T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','Rep');
subjList = [T.b_num];
examList = [T.t_num];

%% process QSM

for ii = 1:length(subjList) % [15 33]
    
    examPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    
%     system(['rm ' examPath '/swan_qsm/HDBET_allQSM/QSM_iQSM2_meanEcho.nii.gz']);
%     system(['rm -r ' examPath '/swan_qsm/HDBET_allQSM/iQSM2']);
    if exist([examPath '/swan_qsm/HDBET_allQSM/QSM_iQSM2_meanEcho.nii.gz'], 'file') ~= 0
        disp([subjList{ii} ' ' examList{ii} ' already processed'])
    else
        disp(['Processing ' subjList{ii} ' ' examList{ii}])
        
        try
            exam_id = [subjList{ii} '_' examList{ii}];
            exam_root = [examPath '/swan_qsm/HDBET_allQSM'];
            mkdir(exam_root);
            
            fprintf(' - Processing exam %s \n', exam_id);
            
            run_QSM('', [examPath '/swan_qsm'], exam_root, 1);
        catch
            disp('!!! Something went wrong!');
        end
        
        % organize folder - put all idf in folder
        mkdir([examPath '/swan_qsm/idf_echoes']);
        cmd = sprintf('mv %s/*echo*idf %s', [examPath '/swan_qsm'], [examPath '/swan_qsm/idf_echoes']);
        system(cmd);
        cmd = sprintf('mv %s/*echo*int2 %s', [examPath '/swan_qsm'], [examPath '/swan_qsm/idf_echoes']);
        system(cmd);
        cmd = sprintf('mv %s/*echo*real %s', [examPath '/swan_qsm'], [examPath '/swan_qsm/idf_echoes']);
        system(cmd);
        
    end
end

%% register QSM to T1

QSMfile_list = {'QSM_iLSQR_meanEcho','QSM_FANSI_nonlinearTV_meanEcho',...
    'QSM_QSMGAN_meanEcho','QSM_SSTV_meanEcho','QSM_HDQSM_meanEcho',...
    'QSM_SSTGV_meanEcho','QSM_STARQSM_meanEcho','QSM_MEDI_meanEcho',...
    'QSM_QSIP_meanEcho','QSM_QSMnet_meanEcho','QSM_iQSM2_meanEcho','QSM_xQSM2_meanEcho'};

for ii = 1:length(subjList)
    
    examPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    output_root = [examPath '/swan_qsm/HDBET_allQSM/FSseg'];
    exam_id = [subjList{ii} '_' examList{ii}];
    T1file = [output_root '/' exam_id '_T1_N4.nii'];
    
%     delete([output_root '/../QSM_iQSM2_meanEcho.nii.gz']);
%     delete([output_root '/QSM_iQSM2_meanEcho_reg.nii.gz']);
    
    if exist([output_root '/QSM_iQSM2_meanEcho_reg.nii.gz'], 'file') ~= 0
    else
        disp(['Processing T1reg ' subjList{ii} ' ' examList{ii}]);
        
        % register Magni to T1 and apply to QSM maps
        if exist([output_root '/Magni_reg.nii.gz'],'file') == 0
            dof = 6;
            cmdstr = ['/netopt/rhel7/fsl/bin/flirt' ...
                ' -in ' output_root '/../Magni.nii.gz' ...
                ' -ref ' T1file ...
                ' -out ' output_root '/Magni_reg.nii.gz' ...
                ' -omat ' output_root '/QSM2T1.mat' ...
                ' -cost corratio -searchrz -90 90 -dof ' num2str(dof) ' -interp trilinear'];
            system(cmdstr);
        end
        if exist([output_root '/QSM2T1F.mat'],'file') == 0
            cmdstr = ['bash /home/jyao3/030_QSM/01_Code/convertMat.sh'...
                ' ' output_root '/QSM2T1.mat'...
                ' > ' output_root '/QSM2T1F.mat'];
            system(cmdstr);
        end
        
        if exist([output_root '/../QSM_iQSM2_meanEcho.nii.gz'],'file') == 0
            cmd = sprintf('3dresample -master %s -prefix %s -input %s', ...
                [output_root '/../QSM_iLSQR_meanEcho.nii.gz'], ...
                [output_root '/../QSM_iQSM2_meanEcho.nii.gz'], ...
                [output_root '/../iQSM2/iQSM_echo_fitted.nii']);
            system(cmd);
        end
        
        for qq = 1:length(QSMfile_list)
            if exist([output_root '/' QSMfile_list{qq} '_reg.nii.gz'],'file') == 0
                cmdstr = ['/netopt/rhel7/fsl/bin/flirt' ...
                    ' -in ' output_root '/../' QSMfile_list{qq} '.nii.gz' ...
                    ' -ref ' T1file ...
                    ' -out ' output_root '/' QSMfile_list{qq} '_reg.nii.gz' ...
                    ' -applyxfm -init ' output_root '/QSM2T1F.mat -interp trilinear'];
                system(cmdstr);
            end
        end
        Maskfile = [output_root '/../brain_mask_HD.nii.gz'];
        if exist([output_root '/brain_mask_reg.nii.gz'],'file') == 0
            cmdstr = ['/netopt/rhel7/fsl/bin/flirt' ...
                ' -in ' output_root '/../brain_mask_HD.nii.gz' ...
                ' -ref ' T1file ...
                ' -out ' output_root '/brain_mask_reg.nii.gz' ...
                ' -applyxfm -init ' output_root '/QSM2T1F.mat -interp nearestneighbour'];
            system(cmdstr);
        end
    end
    
end

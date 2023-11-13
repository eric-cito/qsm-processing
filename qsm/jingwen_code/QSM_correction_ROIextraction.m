clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220909.xlsx','Sheet','NoRep');
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status];
ageList = [T.age];

matout_root = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/Data';

QSMfile_list = {'QSM_FANSI_nonlinearTV_meanEcho' ...
    'QSM_HDQSM_meanEcho' ...
    'QSM_iLSQR_meanEcho' ...
    'QSM_MEDI_meanEcho' ...
    'QSM_QSIP_meanEcho' ...
    'QSM_QSMGAN_meanEcho' ...
    'QSM_QSMnet_meanEcho' ...
    'QSM_SSTGV_meanEcho' ...
    'QSM_SSTV_meanEcho' ...
    'QSM_STARQSM_meanEcho' ...
    };

%% subject list from 7T_hunt

subjFolders = dir('/data/7T_hunt');
subjFolders(~[subjFolders.isdir]) = [];
subjFolders(~cellfun(@isempty,strfind({subjFolders.name},'.'))) = [];

%% read in ROI index

ROItxt = fileread('/working/lupolab/jingwen/090_Brain_ROIs_QSM_Atlas/ROIfromWord.txt');
ROItxt = strsplit(ROItxt, '\n');
ROItxt(cellfun('isempty', ROItxt)) = [];
ROItxt(strcmp(ROItxt, ' ')) = [];

ROIlist = cell(length(ROItxt),2);

for ii = 1:length(ROItxt)
    temp = strsplit(ROItxt{ii},' ');
    temp(cellfun('isempty', temp)) = [];
    if length(temp) > 1
        ROIlist{ii,1} = [temp{1:end-1}];
        ROIlist{ii,2} = str2double(temp{end});
    end
end

%% Loop through subjects

latex = {};

for ii = 1:length(subjFolders)
    
    ind = strcmp(subjList, subjFolders(ii).name);
    ind = find(ind);
    
    if isempty(ind)
        continue
    end
    
    status = statusList{ind(1)};
    
    fprintf('%s is %s \n', subjFolders(ii).name, status);
    
    %     if ~strcmp(status,'HC')
    %         continue
    %     end
    
    subjPath = ['/data/7T_hunt/' subjFolders(ii).name];
    examFolders = dir(subjPath);
    examFolders(~[examFolders.isdir]) = [];
    examFolders(~cellfun(@isempty,strfind({examFolders.name},'.'))) = [];
    
    for ee = 1:length(examFolders)
        examPath = [subjPath '/' examFolders(ee).name];
        if exist([examPath '/swan_qsm/HDBET_allQSM/Phase.nii.gz'], 'file') == 2
            
            % segmenting lateral ventricles as reference tissue
            if exist([examPath '/swan_qsm/HDBET_allQSM/QSMseg_latven.nii.gz'],'file') ~= 2
                try
                    exam_id = [subjFolders(ii).name '_' examFolders(ee).name];
                    segfile_root = [examPath '/swan_qsm/HDBET_allQSM/FASTseg'];
                    mkdir(segfile_root);
                    
                    fprintf(' - Segmenting exam %s \n', exam_id);
                    
                    segLatVen_QSM(examPath, segfile_root);
                    
                    % plot for QC
                    if exist([segfile_root '/QSMseg_latven.png'],'file') ~= 2
                        nii = load_nii([segfile_root '/QSMseg_latven.nii.gz']);
                        LatVen = double(nii.img);
                        isosurface(LatVen,0.5);
                        export_fig([segfile_root '/QSMseg_latven'], '-png'); close;
                    end
                    latex_cell = ScreenShotLatexCell(exam_id, segfile_root);
                    latex(end+1:end+length(latex_cell)) = latex_cell;
                    
                    % copy to one dir up
                    if exist([segfile_root '/../QSMseg_latven.nii.gz'],'file') ~= 2
                        copyfile([segfile_root '/QSMseg_latven.nii.gz'], ...
                            [segfile_root '/../QSMseg_latven.nii.gz']);
                    end
                    
                catch
                    disp('!!! Something went wrong during segmentation!');
                end
            end
            
            % extract ROI values and save in .mat files
            % try
            exam_id = [subjFolders(ii).name '_' examFolders(ee).name];
            if exist([matout_root '/' exam_id '_erode.mat'],'file') ==2
                load([matout_root '/' exam_id '_erode.mat'], 'QSMstats');
            end
            if 1 % exist([matout_root '/' exam_id '_erode.mat'],'file') ~= 2 ...
                %  || height(QSMstats(1).QSMtable) < 213
                QSMfile_root = [examPath '/swan_qsm/HDBET_allQSM/'];
                QSMfile_list = {'QSM_FANSI_nonlinearTV_meanEcho' ...
                    'QSM_HDQSM_meanEcho' ...
                    'QSM_iLSQR_meanEcho' ...
                    'QSM_MEDI_meanEcho' ...
                    'QSM_QSIP_meanEcho' ...
                    'QSM_QSMGAN_meanEcho' ...
                    'QSM_QSMnet_meanEcho' ...
                    'QSM_SSTGV_meanEcho' ...
                    'QSM_SSTV_meanEcho' ...
                    'QSM_STARQSM_meanEcho'};
                ROIfile = [QSMfile_root '/MNIreg/QSM_atlas_ROI.nii.gz'];
                CSFfile = [QSMfile_root '/QSMseg_latven.nii.gz'];
                Maskfile = [QSMfile_root '/brain_mask_HD.nii.gz'];
                
                [QSMdata, QSMstats] = ...
                    stats_QSM(QSMfile_root, QSMfile_list, ROIfile, CSFfile, ...
                    Maskfile, ROIlist, 1);
                QSMstats(1).SubjName = subjFolders(ii).name;
                QSMstats(1).ExamName = examFolders(ee).name;
                save([matout_root '/' exam_id '_erode.mat'], 'QSMdata', 'QSMstats');
            end
            
            %             % catch
            %                 disp('!!! Something went wrong during value extraction!');
            %             end
            
        end
        
    end
    
end

%% Load QSM slopes

for nn = 1:length(QSMfile_list)
    QSMslopeFile = ['/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/NIFTI_HC/' ...
        QSMfile_list{nn} '_Slope.nii.gz'];
    nii = load_nii(QSMslopeFile);
    QSMslope(:,:,:,nn) = double(nii.img);
end

%% Loop through subjects - age corrected QSM

medianAge = 40; %  median(ageList);

for ii = 1:length(subjList)
    
    age = ageList(ii);
    examPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    exam_id = [subjList{ii} '_' examList{ii}];
    
    fprintf('Processing %s \n', exam_id);
    
    if exist([matout_root '/' exam_id '_voxCorr.mat'],'file') ~= 2
        QSMfile_root = [examPath '/swan_qsm/HDBET_allQSM'];
        
        ROIfile = ['/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/NIFTI_HC/' ...
            'Whole_brain_segmentation_DGM.nii'];
        Maskfile = ['/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/NIFTI_HC/' ...
            'MNI_GMWM.nii.gz'];
        
        [QSMdata, QSMstats] = ...
            stats_QSM_voxCorr(QSMfile_root, QSMfile_list, ROIfile, ...
            Maskfile, ROIlist, 1, QSMslope, age, medianAge);
        QSMstats(1).SubjName = subjList{ii};
        QSMstats(1).ExamName = examList{ii};
        save([matout_root '/' exam_id '_voxCorr.mat'], 'QSMdata', 'QSMstats');
    end
    
end

%% create QC report

path2reportfolder = '/home/jyao3/030_QSM/01_Code/QCreport_template';

% save as txt
fid = fopen([path2reportfolder '/ScreenShots.txt'], 'wt');
fprintf(fid, '%s\n',char(latex)');
fclose(fid);

% create report
cd(path2reportfolder);
% cmd = 'buf_size=2000000 pdflatex QCreport.tex;';
% system(cmd);

% move to working dir
save_root = '/working/lupolab/jingwen/001_QSM/';
copyfile('QCreport.pdf',[save_root 'QCreport_QSMseg.pdf']);

cd(save_root);
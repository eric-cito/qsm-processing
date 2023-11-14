clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

%% read in list
T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList.xlsx');
subjList = [T.b_num];
statusList = [T.status];

%% subject list from 7T_hunt

subjFolders = dir('/data/7T_hunt');
subjFolders(~[subjFolders.isdir]) = [];
subjFolders(~cellfun(@isempty,strfind({subjFolders.name},'.'))) = [];

%% Loop through subjects

latex = {};

for ii = 73:length(subjFolders) % [73:length(subjFolders)] % 7 22 63
    
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
        if exist([examPath '/swan_qsm'], 'dir') == 7
            
            try
                
                exam_id = [subjFolders(ii).name '_' examFolders(ee).name];
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
            
            % save images for QC
            QCfolder = [exam_root '/QC_images'];
            mkdir(QCfolder);
            
            QSMfile = dir(exam_root);
            QSMfile([QSMfile.isdir]) = [];
            QSMfile(cellfun(@isempty,strfind({QSMfile.name},'QSM_'))) = [];
            
            if isempty(QSMfile)
                disp('!!! No QSM to plot!');
            else
                for nn = 1:length(QSMfile)
                    
                    QSMname = [exam_root '/' QSMfile(nn).name];
                    C = strsplit(QSMfile(nn).name,'.');
                    savename = C{1};
                    
                    if exist([QCfolder '/' savename '.png'],'file') ~= 2
                        targetSize = [300 300 148];
                        
                        nii = load_untouch_nii(QSMname);
                        QSM = nii.img;
                        currSize = size(QSM);
                        padSize = (targetSize - currSize)/2;
                        
                        QSM_Ax = padarray(flip(rot90(QSM(:,:,74),1)),padSize([2 1]));
                        QSM_Cor = padarray(flip(rot90(squeeze(QSM(150,:,:)),1)),padSize([3 2]));
                        QSM_Sag = padarray(flip(rot90(squeeze(QSM(:,150,:)),1)),padSize([3 1]));
                        
                        fig = figure('position',[100 100 500 500]);
                        imshow([QSM_Ax; QSM_Cor; QSM_Sag],[-0.15 0.15]); colormap gray;
                        export_fig([QCfolder '/' savename], '-png'); close;
                    end
                end
                
                latex_cell = ScreenShotLatexCell(exam_id, QCfolder);
                latex(end+1:end+length(latex_cell)) = latex_cell;
            end
            
        end
        
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
cmd = 'buf_size=2000000 pdflatex QCreport.tex;';
system(cmd);

% move to working dir
save_root = '/working/lupolab/jingwen/001_QSM/';
copyfile('QCreport.pdf',[save_root 'QCreport_QSM.pdf']);
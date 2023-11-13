clear; clc;
warning('off');

%% Add path

addpath(genpath('/home/jyao3/030_QSM/01_Code'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/NIfTI'));
addpath(genpath('/home/jyao3/010_MATLAB_Utils/export_fig'));

%% read in list

T = readtable('/working/lupolab/jingwen/001_QSM/QSM_SubjectList_20220209.xlsx','Sheet','NoRep');
subjList = [T.b_num];
examList = [T.t_num];
statusList = [T.status];
ageList = [T.age];

matout_root = '/working/lupolab/jingwen/001_QSM/02_QSM_HC_age/Data';

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

%% Loop through subjects and copy DTI files to /swan_qsm/HDBET_allQSM/DTI/

for ii = 32 % 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    
    exam_id = [subjList{ii} '_' examList{ii}];
    fprintf('Exam %s \n', exam_id);
    
    QSMfile = [subjPath '/swan_qsm/HDBET_allQSM/Magni_brain.nii.gz'];
    
    % find b0 and TOPUP_FA files
    b0file = [subjPath '/NODDI/processed_dti/b0.nii.gz'];
    FAorig = dir([subjPath '/NODDI/processed_dti/']);
    FAorig(~contains({FAorig.name},'TOPUP_FA.nii.gz')) = [];
    if exist(b0file, 'file') ~= 2
        b0file = [subjPath '/NODDI/b0.nii.gz'];
        FAorig = dir([subjPath '/NODDI/']);
        FAorig(~contains({FAorig.name},'TOPUP_FA.nii.gz')) = [];
        if exist(b0file, 'file') ~= 2
            b0file = [subjPath '/NODDI_melanie/b0.nii.gz'];
            FAorig = dir([subjPath '/NODDI_melanie/']);
            FAorig(~contains({FAorig.name},'TOPUP_FA.nii.gz')) = [];
            if exist(b0file, 'file') ~= 2
                b0file = [subjPath '/dti/processed_dti/b0.nii.gz'];
                FAorig = dir([subjPath '/dti/processed_dti/']);
                FAorig(~contains({FAorig.name},'TOPUP_FA.nii.gz')) = [];
                if exist(b0file, 'file') ~= 2
                    b0file = [subjPath '/../NODDI/processed_dti/b0.nii.gz'];
                    FAorig = dir([subjPath '/../NODDI/processed_dti/']);
                    FAorig(~contains({FAorig.name},'TOPUP_FA.nii.gz')) = [];
                    if exist(b0file, 'file') ~= 2
                        b0file = [subjPath '/../orig_raw/dti_raw/processed_dti/b0.nii.gz'];
                        FAorig = dir([subjPath '/../orig_raw/dti_raw/processed_dti']);
                        FAorig(~contains({FAorig.name},'TOPUP_FA.nii.gz')) = [];
                    end
                end
            end
        end
    end
    
    if isempty(FAorig) || length(FAorig) > 1
        disp('Check FA file in NODDI directory!');
        continue;
    end
    
    if exist(b0file, 'file') == 2
        
        % copy files from DTI folder
        QSMfile_root = [subjPath '/swan_qsm/HDBET_allQSM/'];
        DTIfile_root = [subjPath '/swan_qsm/HDBET_allQSM/DTI/'];
        mkdir(DTIfile_root);
        
        FAfile = [DTIfile_root '/FA.nii.gz'];
        b0reg = [DTIfile_root '/b0_reg.nii.gz'];
        b0regmat = [DTIfile_root '/b0_reg.mat'];
        FAreg = [DTIfile_root '/FA_reg.nii.gz'];
        FAregmat = [DTIfile_root '/FA_reg.mat'];
        
        %         cmd = sprintf('rm %s/FA_reg.nii.gz',QSMfile_root);
        %         system(cmd);
        %         cmd = sprintf('rm %s/DTI/*reg*',QSMfile_root);
        %         system(cmd);
        
        % copy files
        FApath = [FAorig.folder '/' FAorig.name];
        if exist([DTIfile_root '/FA.nii.gz'],'file') ~= 2
            copyfile(FApath, [DTIfile_root '/FA.nii.gz']);
        end
        Fpath = strrep(FApath,'TOPUP_FA.nii.gz','TOPUP_FA_alignT1.nii.gz');
        if exist([DTIfile_root '/FA_regT1.nii.gz'],'file') ~= 2 && exist(Fpath,'file') == 2
            copyfile(Fpath, [DTIfile_root '/FA_regT1.nii.gz']);
        end
        Fpath = strrep(FApath,'TOPUP_FA.nii.gz','TOPUP_MD_alignT1.nii.gz');
        if exist([DTIfile_root '/MD_regT1.nii.gz'],'file') ~= 2 && exist(Fpath,'file') == 2
            copyfile(Fpath, [DTIfile_root '/MD_regT1.nii.gz']);
        elseif exist([DTIfile_root '/MD.nii.gz'],'file') ~= 2
            Fpath = strrep(FApath,'TOPUP_FA.nii.gz','TOPUP_MD.nii.gz');
            copyfile(Fpath, [DTIfile_root '/MD.nii.gz']);
        end
        Fpath = strrep(FApath,'TOPUP_FA.nii.gz','TOPUP_RD_alignT1.nii.gz');
        if exist([DTIfile_root '/RD_regT1.nii.gz'],'file') ~= 2 && exist(Fpath,'file') == 2
            copyfile(Fpath, [DTIfile_root '/RD_regT1.nii.gz']);
        elseif exist([DTIfile_root '/RD.nii.gz'],'file') ~= 2
            Fpath = strrep(FApath,'TOPUP_FA.nii.gz','TOPUP_RD.nii.gz');
            if exist(Fpath,'file') == 2
                copyfile(Fpath, [DTIfile_root '/RD.nii.gz']);
            else
                cmd = sprintf('3dcalc -a %s -b %s -expr ''mean(a,b)'' -prefix %s', ...
                    strrep(FApath,'TOPUP_FA.nii.gz','TOPUP_L2.nii.gz'), ...
                    strrep(FApath,'TOPUP_FA.nii.gz','TOPUP_L3.nii.gz'), Fpath);
                system(cmd);
                copyfile(Fpath, [DTIfile_root '/RD.nii.gz']);
            end
        end
        Fpath = strrep(FApath,'TOPUP_FA.nii.gz','NODDI_odi_alignT1.nii.gz');
        if exist([DTIfile_root '/ODI_regT1.nii.gz'],'file') ~= 2 && exist(Fpath,'file') == 2
            copyfile(Fpath, [DTIfile_root '/ODI_regT1.nii.gz']);
        elseif exist([DTIfile_root '/ODI.nii.gz'],'file') ~= 2
            Fpath = strrep(FApath,'TOPUP_FA.nii.gz','NODDI_odi.nii.gz');
            if exist(Fpath,'file') == 2
                copyfile(Fpath, [DTIfile_root '/ODI.nii.gz']);
            end
        end
        
        if exist(FAreg,'file') ~= 2
            fprintf(' - Registrating exam %s \n', exam_id);
            
            copyfile(b0file, [DTIfile_root '/b0.nii.gz']);
            copyfile([FAorig.folder '/' FAorig.name], [DTIfile_root '/FA.nii.gz']);
            
            % register to QSM
            cmd = sprintf('flirt -in %s -ref %s -out %s -omat %s -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90',...
                b0file, QSMfile, b0reg, b0regmat);
            system(cmd);
            cmdstr = ['bash /home/jyao3/030_QSM/01_Code/convertMat.sh'...
                ' ' DTIfile_root '/b0_reg.mat'...
                ' > ' DTIfile_root '/b0_regF.mat'];
            system(cmdstr);
            b0regmat = [DTIfile_root '/b0_regF.mat'];
            cmd = sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s',...
                FAfile, QSMfile, FAreg, b0regmat);
            system(cmd);
            
            % copy to one dir up
            copyfile(FAreg, [DTIfile_root '../FA_reg.nii.gz']);
            
        end
        
    end
    
end

%% Loop through subjects and register DTI to T1

FS_root = '/working/lupolab/jingwen/002_HD_NDM/AllT1N4';

for ii = 1:length(subjList)
    
    subjPath = ['/data/7T_hunt/' subjList{ii} '/' examList{ii}];
    
    exam_id = [subjList{ii} '_' examList{ii}];
    T1file = [FS_root '/' exam_id '_T1_N4.nii'];
    DTIfile_root = [subjPath '/swan_qsm/HDBET_allQSM/DTI/'];
    
    FAreg = [DTIfile_root '/RD_regT1.nii.gz'];
    FA = [DTIfile_root '/FA.nii.gz'];
    
%     if exist([DTIfile_root '/b0_regT1F.mat'],'file') == 2
%         system(['rm ' DTIfile_root '/*_regT1*']);
%     end
    
    if exist(FAreg,'file') ~= 2 && exist(FA,'file') == 2
        fprintf('Exam %s \n', exam_id);
        fprintf(' - Registrating exam %s \n', exam_id);
        
        % register FA to T1
        copyfile(T1file, [DTIfile_root '/T1.nii.gz']);
        b0file = [DTIfile_root '/FA.nii.gz'];
        b0reg = [DTIfile_root '/FA_regT1.nii.gz'];
        b0regmat = [DTIfile_root '/b0_regT1.mat'];
        cmd = sprintf('flirt -in %s -ref %s -out %s -omat %s -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90',...
            b0file, T1file, b0reg, b0regmat);
        system(cmd);
        cmdstr = ['bash /home/jyao3/030_QSM/01_Code/convertMat.sh'...
            ' ' DTIfile_root '/b0_regT1.mat'...
            ' > ' DTIfile_root '/b0_regT1F.mat'];
        system(cmdstr);
        b0regmat = [DTIfile_root '/b0_regT1F.mat'];
        
        % register MD to T1
        MOVreg = [DTIfile_root '/MD_regT1.nii.gz'];
        MOV = [DTIfile_root '/MD.nii.gz'];
        cmd = sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s',...
            MOV, T1file, MOVreg, b0regmat);
        system(cmd);
        % register RD to T1
        MOVreg = [DTIfile_root '/RD_regT1.nii.gz'];
        MOV = [DTIfile_root '/RD.nii.gz'];
        cmd = sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s',...
            MOV, T1file, MOVreg, b0regmat);
        system(cmd);
        % register ODI to T1
        MOVreg = [DTIfile_root '/ODI_regT1.nii.gz'];
        MOV = [DTIfile_root '/ODI.nii.gz'];
        cmd = sprintf('flirt -in %s -ref %s -out %s -applyxfm -init %s',...
            MOV, T1file, MOVreg, b0regmat);
        system(cmd);
        
    end
    
end

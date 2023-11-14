function [] = run_QSM(currentPath, input_data_path, output_data_path, run_All)

% Usage
% run from terminal
% cd to the exam folder containing Exxxx under bxxxx/txxxx/
% command: 
% matlab -nodisplay -nodesktop -r "cd('/working/lupolab/jingwen/001_QSM/01_Code/');run_QSM('$(pwd)');exit;"

%% add path

warning('off');
addpath(genpath('/home/jyao3/030_QSM/01_Code'));

%% set up data paths

if nargin < 4
    run_All = 0;
end

if nargin == 0
    input_data_path = '/working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308';
    output_data_path = '/working/lupolab/jingwen/001_QSM/Test_data/b4468_t12308/output_HDBET';
elseif nargin == 1
    path_root = '/data/7T_hunt';
    % currentPath = pwd;
    fprintf('## Processing path : %s \n', currentPath);
    
    C = strsplit(currentPath,'/');
    if length(C) < 2; error('Wrong directory!'); end;
    tnumber = C{end};
    bnumber = C{end-1};
    if ~strcmp(bnumber(1),'b') && ~strcmp(bnumber(1:4),'temp')
        fprintf('bnumber : %s \n', bnumber);
        error('Wrong directory!');
    end
    if ~strcmp(tnumber(1),'t') && ~strcmp(tnumber(1:3),'for')
        fprintf('tnumber : %s \n', tnumber);
        error('Wrong directory!');
    end
    
    input_data_path = [path_root '/' bnumber '/' tnumber '/swan_qsm'];
    output_data_path = [input_data_path '/output_HDBET'];
else
    C = strsplit(input_data_path,'/');
    if length(C) < 3; error('Wrong directory!'); end;
    tnumber = C{end-1};
    bnumber = C{end-2};
    if ~strcmp(bnumber(1),'b') && ~strcmp(bnumber(1:4),'temp')
        fprintf('bnumber : %s \n', bnumber);
        error('Wrong directory!');
    end
    if ~strcmp(tnumber(1),'t') && ~strcmp(tnumber(1:3),'for')
        fprintf('tnumber : %s \n', tnumber);
        error('Wrong directory!');
    end
end

%% reconstruct QSM

if exist([input_data_path '/Magni.nii.gz'], 'file') ~= 2
    [all_TE, ~, ~] = QSM_recon(bnumber, tnumber, input_data_path);
    
    cd(input_data_path)
    set_SEPIA_header;
    
    if sum(isnan(all_TE)) == 0
        header.te = all_TE*1e-3; % ms -> s
        header.delta_TE = all_TE(2)*1e-3 - all_TE(1)*1e-3;
    end
    save('sepia_header.mat','header');
    
elseif nargin >= 3 && exist([input_data_path '/Magni.nii.gz'], 'file') ~= 2
    C = strsplit(input_data_path,'/');
    if length(C) < 3; error('Wrong directory!'); end;
    tnumber = C{end-1};
    bnumber = C{end-2};
    if ~strcmp(bnumber(1),'b') && ~strcmp(bnumber(1:4),'temp')
        fprintf('bnumber : %s \n', bnumber);
        error('Wrong directory!');
    end
    if ~strcmp(tnumber(1),'t') && ~strcmp(tnumber(1:3),'for')
        fprintf('tnumber : %s \n', tnumber);
        error('Wrong directory!');
    end
    
    [all_TE, ~, ~] = QSM_recon(bnumber, tnumber, input_data_path);
    
    cd(input_data_path)
    set_SEPIA_header;
    
    header.te = all_TE*1e-3;
    header.delta_TE = all_TE(2)*1e-3 - all_TE(1)*1e-3;
    save('sepia_header.mat','header');
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
opts.All                = run_All;

%% run

QSMfile = dir(output_data_path);
QSMfile([QSMfile.isdir]) = [];
QSMfile(cellfun(@isempty,strfind({QSMfile.name},'QSM_'))) = [];

if (opts.All) || (~opts.All && length(QSMfile) ~= 2)
    [Tcomp] = QSM_processing(input, output_data_path, opts);
else
    fprintf(' >> Already have all QSM maps, skip QSM processing \n');
end

end

% QSM_processing example

%% add path

warning('off');
addpath(genpath('/home/jyao3/030_QSM/01_Code'));

%% set up data paths

%%% Customize this %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_data_path = '/path/to/subject/folder/containing/niftis';
output_data_path = [input_data_path '/output_QSM'];
%%% Customize this - end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare header file for processing

%%% Customize this %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up echo times
te1 = 6e-3; % s
te2 = 9.5e-3; % s
te3 = 13e-3; % s
te4 = 16.5e-3;% s

% set up dimensions
% spatial resolution of the data, in mm
header.voxelSize = [0.8,0.8,1];  
% image matrix size
header.matrixSize = [300,300,148];

% set up parameters
header.b0 = 7;                  % magnetic field strength, in Tesla
%%% Customize this - end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

header.b0dir = [0;0;1];         % main magnetic field direction, [x,y,z]
header.CF = header.b0*42.58*1e6;       % imaging frequency, in Hz (B0*gyromagnetic_ratio)
header.te = [te1,te2,te3,te4];  % echo time for each GRE image, in second
header.delta_TE = te2-te1;      % echo spacing, in second

% save to mat file
save([input_data_path '/sepia_header.mat'],'header');

%% set up algorithm parameters

%%% Please make sure the file names are the same as used below %%%%%%%%%%%%
input.magnitudeFile     = [input_data_path '/Magni.nii.gz'];
input.phaseFile         = [input_data_path '/Phase.nii.gz'];
%%% Check name - end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input.headerFile        = [input_data_path '/sepia_header.mat'];

opts.writeLog           = 1;
opts.isGPU              = 1;
opts.BETmethod          = 'HD-BET';
opts.iLSQR              = 1;
opts.QSMGAN             = 1;
opts.All                = 0;

%% run

[Tcomp] = QSM_processing(input, output_data_path, opts);

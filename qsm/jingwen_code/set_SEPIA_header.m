
%% set up echo times
te1 = 6e-3; % s
te2 = 9.5e-3; % s
te3 = 13e-3; % s
te4 = 16.5e-3;% s

%% set up dimensions
% spatial resolution of the data, in mm
header.voxelSize = [0.8,0.8,1];  
% image matrix size
header.matrixSize = [300,300,148];

%% set up parameters
header.b0 = 7;                  % magnetic field strength, in Tesla
header.b0dir = [0;0;1];        % main magnetic field direction, [x,y,z]
header.CF = header.b0*42.58*1e6;       % imaging frequency, in Hz (B0*gyromagnetic_ratio)
header.te = [te1,te2,te3,te4];  % echo time for each GRE image, in second
header.delta_TE = te2-te1;      % echo spacing, in second

%% save to mat file
save('sepia_header.mat','header');
%% [chi] = Wrapper_QSM_HDQSM(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
%
% Input
% --------------
% localField    : local field map (tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% chi           : magnetic susceptibility map, in ppm
%
% Description: This is a wrapper function to access QSIP-QSM for SEPIA
%
% Jingwen Yao
% jingwen.yao@ucsf.ucsf.edu
% Date created: 25 June 2021
% Date modified: 25 June 2021
%
%
function [chi] = Wrapper_QSM_QSIP(totalField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)

sepia_universal_variables;

% add path
addpath(fullfile(SEPIA_HOME,'addons','qsm','QSIP'));
pn_pre = fullfile(SEPIA_HOME,'addons','qsm','QSIP','data');
pn_pre = [pn_pre '/'];

%% Set up parameters

method                  = algorParam.qsm.method;

% get algorithm parameters
algorParam              = check_and_set_algorithm_default(algorParam);

params.atlas_thr        = algorParam.qsm.QSIP.atlas_thr;  % 0.96;        % Atlas threshold
params.atlas_flag       = algorParam.qsm.QSIP.atlas_flag; % 1;           % set to 1 to use atlas, 0 to not use it
params.chi0             = 0.4;                              % unitless, in ppm (do not change)
params.delta            = -9.5;                             % unitless, in ppm (do not change)
params.num_iter         = algorParam.qsm.QSIP.num_iter; % = 300;          % number of CG iterations
params.bo               = headerAndExtraData.b0;            % B0 field in Tesla
params.te               = headerAndExtraData.te;            % echo time in seconds
params.voxelsize        = voxelSize;                        % voxel size

params.LPWeight         = algorParam.qsm.QSIP.LPWeight;
params.OBJWeight        = algorParam.qsm.QSIP.OBJWeight;
params.EXTWeight        = algorParam.qsm.QSIP.EXTWeight;

phase_uw = totalField;  % Rad /(2*pi*params.te*headerAndExtraData.b0*gyro);
magni = headerAndExtraData.magn;
dims = size(mask);

%% Display algorithm parameters

disp('The following parameters are being used...');

disp(['atlas_thr                    = ' num2str(params.atlas_thr)]);
disp(['atlas_flag                   = ' num2str(params.atlas_flag)]);
disp(['chi0                         = ' num2str(params.chi0)]);
disp(['delta                        = ' num2str(params.delta)]);
disp(['num_iter                     = ' num2str(params.num_iter)]);
disp(['LPWeight                     = ' num2str(params.LPWeight)]);
disp(['OBJWeight                    = ' num2str(params.OBJWeight)]);
disp(['EXTWeight                    = ' num2str(params.EXTWeight)]);

%% Prepare data

[atlas, noise_mask, est0] = run_prepare_data(pn_pre, dims);
atlas = atlas(end:-1:1, end:-1:1, end:-1:1, :);
noise_mask = noise_mask(end:-1:1, end:-1:1, end:-1:1);
est0 = est0(end:-1:1, end:-1:1, end:-1:1);

%% main

% BG field removal
[~, nfm_dipole, ~, ~, ~, ~, ~, ~] = ...
    fn_fastDipole_atlas_new(atlas, mask, magni, phase_uw, params);

% QSIP
[~, ~, ~, chi_est_ppm, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
    KOP_intsources_nfm_corr3_lambda(nfm_dipole, mask, magni, phase_uw, noise_mask, params);

chi = chi_est_ppm.*mask;

end

%% set default parameter for unspecific input

function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.qsm.QSIP.atlas_thr      = algorParam.qsm.QSIP.atlas_thr;    catch;  algorParam2.qsm.QSIP.atlas_thr 	= 0.96; end
try algorParam2.qsm.QSIP.atlas_flag     = algorParam.qsm.QSIP.atlas_flag;  	catch;  algorParam2.qsm.atlas_flag     	= 1;    end
try algorParam2.qsm.QSIP.num_iter       = algorParam.qsm.QSIP.num_iter;     catch;  algorParam2.qsm.QSIP.num_iter   = 300;  end
try algorParam2.qsm.QSIP.LPWeight       = algorParam.qsm.QSIP.LPWeight;     catch;  algorParam2.qsm.QSIP.LPWeight 	= 1e-8; end
try algorParam2.qsm.QSIP.OBJWeight      = algorParam.qsm.QSIP.OBJWeight;  	catch;  algorParam2.qsm.QSIP.OBJWeight  = 1e-6; end
try algorParam2.qsm.QSIP.EXTWeight      = algorParam.qsm.QSIP.EXTWeight;    catch;  algorParam2.qsm.QSIP.EXTWeight  = 1e10; end

end
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
% Description: This is a wrapper function to access HD-QSM for SEPIA
%
% Jingwen Yao
% jingwen.yao@ucsf.ucsf.edu
% Date created: 25 June 2021
% Date modified: 25 June 2021
%
%
function [chi] = Wrapper_QSM_SSTGV(totalField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)

sepia_universal_variables;

% add path
addpath(fullfile(SEPIA_HOME,'addons','qsm','SS_TGV_TV'));

%% Calculate kernel

N = matrixSize;
B0_dir = headerAndExtraData.b0dir;
b0 = headerAndExtraData.b0;
D = create_dipole_kernel(B0_dir, voxelSize, N, 1);

%% Generate SMV kernels and masks

min_radius = 1;
max_radius = 5;
step_size_radius = 1;

out = create_SMVkernel(totalField, mask, min_radius, max_radius, step_size_radius, N, voxelSize);

SMV_kernels = out.SMV_kernel;
SMV_inv_kernel = out.SMV_inv_kernel;
SMV_masks = out.SMV_mask;
SMV_phase = out.SMV_phase;
mask_Sharp = out.mask_eval;

%% Set up parameters

method                  = algorParam.qsm.SS.method;

% get algorithm parameters
algorParam              = check_and_set_algorithm_default(algorParam);

params.phase_unwrap     = totalField;                   % Unwrapped total phase !!! Rad

switch method
    case 'SS-TGV'
        params.alpha0           = algorParam.qsm.SSTGV.alpha0;  % Regularization param for ||Gx-v||_1
        params.alpha1           = params.alpha0/2;              % Regularization param for ||Ev||_1
        params.mu0              = algorParam.qsm.SSTGV.mu0;     % Aug Lagrangian param for z0 = Ev
        params.mu1              = params.mu0;                   % Aug Lagrangian param for z1 = Gx
        params.mu2              = params.mu0;                   % Aug Lagrangian param for z2 = HDFx
        params.maxOuterIter     = algorParam.qsm.SSTGV.maxOuterIter; % Max number of iter
        params.tol_soln         = algorParam.qsm.SSTGV.tol_soln;% Stopping criterion: RMSE change in solution
    case 'SS-TV'
        params.alpha            = algorParam.qsm.SSTV.alpha;    % Regularization param for ||Gx||_1
        params.mu1              = algorParam.qsm.SSTV.mu1;      % Aug Lagrangian param for z1 = Gx
        params.mu2              = params.mu1;                   % Aug Lagrangian param for z2 = HDFx
        params.maxOuterIter     = algorParam.qsm.SSTV.maxOuterIter; % Max number of iter
        params.tol_soln         = algorParam.qsm.SSTV.tol_soln; % Stopping criterion: RMSE change in solution
end

params.N                = matrixSize;                   % Number of voxels
params.M                = SMV_masks;                    % Mask for each reliable region
params.H                = SMV_kernels;                  % SMV kernels
params.D                = D;                            % Dipole kernel

%% Display algorithm parameters

disp('The following parameters are being used...');
switch method
    case 'SS-TGV'
        disp(['alpha0                   = ' num2str(params.alpha0)]);
        disp(['alpha1                   = ' num2str(params.alpha1)]);
        disp(['mu0,mu1,mu2              = ' num2str(params.mu0)]);
    case 'SS-TV'
        disp(['alpha                    = ' num2str(params.alpha)]);
        disp(['mu1,mu2                  = ' num2str(params.mu1)]);
end

disp(['maxOuterIter             = ' num2str(params.maxOuterIter)]);
disp(['tol_soln                 = ' num2str(params.tol_soln)]);

%% main

switch method
    case 'SS-TGV'
        out_ss = SS_TGV_QSM(params);
    case 'SS-TV'
        out_ss = SS_TV_QSM(params);
end

chi = mask_Sharp .* out_ss.x /(2*pi*b0*gyro*headerAndExtraData.te); % Rad -> ppm

end

%% set default parameter for unspecific input

function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.qsm.SSTGV.alpha0        = algorParam.qsm.SSTGV.alpha0;      catch;  algorParam2.qsm.SSTGV.alpha0        = 6e-5; end
try algorParam2.qsm.SSTGV.mu0           = algorParam.qsm.SSTGV.mu0;         catch;  algorParam2.qsm.SSTGV.mu0           = 3e-2; end
try algorParam2.qsm.SSTGV.maxOuterIter  = algorParam.qsm.SSTGV.maxOuterIter; catch; algorParam2.qsm.SSTGV.maxOuterIter  = 100; end
try algorParam2.qsm.SSTGV.tol_soln      = algorParam.qsm.SSTGV.tol_soln;    catch;  algorParam2.qsm.SSTGV.tol_soln      = 1; end

try algorParam2.qsm.SSTV.alpha          = algorParam.qsm.SSTV.alpha;        catch;  algorParam2.qsm.SSTV.alpha          = 7e-5; end
try algorParam2.qsm.SSTV.mu1            = algorParam.qsm.SSTV.mu1;          catch;  algorParam2.qsm.SSTV.mu1            = 3e-2; end
try algorParam2.qsm.SSTV.maxOuterIter   = algorParam.qsm.SSTV.maxOuterIter; catch;  algorParam2.qsm.SSTV.maxOuterIter   = 100; end
try algorParam2.qsm.SSTV.tol_soln       = algorParam.qsm.SSTV.tol_soln;     catch;  algorParam2.qsm.SSTV.tol_soln       = 1; end

end
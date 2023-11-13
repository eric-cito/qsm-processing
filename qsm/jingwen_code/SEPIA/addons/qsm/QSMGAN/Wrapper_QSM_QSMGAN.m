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
function [chi] = Wrapper_QSM_HDQSM(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)

sepia_universal_variables;

b0 = headerAndExtraData.b0;

% get algorithm parameters
algorParam      = check_and_set_algorithm_default(algorParam);

params.input    = localField; % Hz
params.weight   = headerAndExtraData.weights;   % magnitude data
params.mask     = mask;
% params.regweight        = algorParam.qsm.regweight;       % Regularization Spatially Variable Weight. Not used if not specified
params.tol_update       = algorParam.qsm.tol_update;        % Convergence Limit of the Second Stage (Default = 1.0)
params.maxOuterIterL1   = algorParam.qsm.maxOuterIterL1;    % Iterations of Fist Stage (Default = 20)
params.mu1L1            = algorParam.qsm.mu1L1;             % Gradient Consistency Weight of Fist Stage (Default = sqrt(params.mu1L2)
params.alphaL1          = algorParam.qsm.alphaL1;           % Regularization Weight of Fist Stage (Default = sqrt(params.alpha11L2)
params.maxOuterIterL2   = algorParam.qsm.maxOuterIterL2;    % Iterations of Second Stage (Default = 80)
params.mu1L2            = algorParam.qsm.mu1L2;             % Gradient Consistency Weight of Second Stage 
params.alphaL2          = algorParam.qsm.alphaL2;           % Regularization Weight of Second Stage

method      = algorParam.qsm.method;

% add path
addpath(fullfile(SEPIA_HOME,'addons','HDqsm'));

%% Calculate kernel

% mode:  0 for the continuous kernel proposed by Salomir, et al. 2003.
%        1 for the discrete kernel proposed by Milovic, et al. 2017.
%        2 for the Integrated Green function proposed by Jenkinson, et al. 2004
mode = 0;
kernel = dipole_kernel_fansi( matrixSize, voxelSize, mode );

params.kernel   = kernel;        % Dipole Kernel in the Frequency Space

%% Display algorithm parameters

disp('The following parameters are being used...');
disp(['tol_update               = ' num2str(params.tol_update)]);
% disp(['regweight                = ' num2str(params.regweight)]);
disp(['maxOuterIterL1           = ' num2str(params.maxOuterIterL1)]);
disp(['mu1L1                    = ' num2str(params.mu1L1)]);
disp(['alphaL1                  = ' num2str(params.alphaL1)]);
disp(['maxOuterIterL2           = ' num2str(params.maxOuterIterL2)]);
disp(['mu1L2                    = ' num2str(params.mu1L2)]);
disp(['alphaL2                  = ' num2str(params.alphaL2)]);

%% main

out = HDqsm(params);

chi = out.x/(b0*gyro) .* mask; % Hz -> ppm

end

%% set default parameter for unspecific input

function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

% try algorParam2.qsm.regweight   = algorParam.qsm.regweight;     catch; algorParam2.qsm.regweight    = [];  end
try algorParam2.qsm.tol_update      = algorParam.qsm.tol_update;        catch; algorParam2.qsm.tol_update       = 1.0; end
try algorParam2.qsm.maxOuterIterL2  = algorParam.qsm.maxOuterIterL2;	catch; algorParam2.qsm.maxOuterIterL2   = 280; end
try algorParam2.qsm.mu1L2           = algorParam.qsm.mu1L2;             catch; algorParam2.qsm.mu1L2            = 100 * 10^-4.785; end
try algorParam2.qsm.alphaL2         = algorParam.qsm.alphaL2;           catch; algorParam2.qsm.alphaL2          = 10^-4.785; end
try algorParam2.qsm.maxOuterIterL1  = algorParam.qsm.maxOuterIterL1;    catch; algorParam2.qsm.maxOuterIterL1   = 20; end
try algorParam2.qsm.mu1L1           = algorParam.qsm.mu1L1;             catch; algorParam2.qsm.mu1L1            = sqrt(algorParam2.qsm.mu1L2);    end
try algorParam2.qsm.alphaL1         = algorParam.qsm.alphaL1;           catch; algorParam2.qsm.alphaL1          = sqrt(algorParam2.qsm.alphaL2);  end

end
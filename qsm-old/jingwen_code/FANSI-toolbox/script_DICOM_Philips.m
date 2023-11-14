% Read and process 3T GRE data (Philips Ingenia) for QSM and R2* reconstructions.
% This example acquisition uses the following sequence parameters:
% TE1=4.0ms, δTE=3.7ms, 5 echoes, TR=37.1ms, flip angle=17°, 
% 288×288×170 matrix with 0.73×0.73×1mm3 voxels. Tacq=5:52 min. 
% The z axis points in the direction of the main field.
% For this sequence Real, imaginary, Magnitude and Phase components were saved.
%
% You may use this script to build your own full QSM pipeline. 
% Most of the functions included in this toolbox are presented as an example, with brief descriptive text.
% Please note that we provide such examples just for completition/exemplary purposes. Many of them 
% are considered in development stages, and we cannot guarantee the best possible results.
% We point to dedicated tools found in other toolboxes in the case best alternatives are available.
% This is NOT an optimized pipeline script.
%
% Last modified by Carlos Milovic in 2020.07.21

set(0,'DefaultFigureWindowStyle','docked')

%%-------------------------------------------------------------------------
%% Read Dicom data
%%-------------------------------------------------------------------------

filename = 'IM_0004';
A = dicomread(filename);
iA = dicominfo(filename);

% Read acquisition paramaters from header
sz = size(A);
N = [sz(1) sz(2) iA.Private_2001_1018]; %sz(4)/4]; % Image size, for RMPI data
ER = iA.Private_2001_1014; %10; % Number of echoes
voxel_size = [0.73 0.73 1];
B0 = 3; % Manually set field strength
gyro = 2*pi*42.58; % Manually set the gyromagnetic ratio

% Correct for Slope and Intercept
perFrameFields = fieldnames(iA.PerFrameFunctionalGroupsSequence);
for i=1:sz(4)
slope = iA.PerFrameFunctionalGroupsSequence.(perFrameFields{i}).PixelValueTransformationSequence.Item_1.RescaleSlope;
inter = iA.PerFrameFunctionalGroupsSequence.(perFrameFields{i}).PixelValueTransformationSequence.Item_1.RescaleIntercept;
B(:,:,1,i) = double(A(:,:,1,i))*slope+inter;
end
clear A slope inter

% Sepparate magnitude, phase, real and imaginary data
ImTot = squeeze(B(:,:,1,1:(sz(4)/4)));
MagTot = squeeze(B(:,:,1,(sz(4)/4+1):(2*sz(4)/4)));
PhTot = squeeze(B(:,:,1,(2*sz(4)/4+1):(3*sz(4)/4)));
ReTot = squeeze(B(:,:,1,(3*sz(4)/4+1):sz(4)));
clear B sz

% Read effective echo times
TE = zeros(1, ER); % In ms
for i=1:ER
    TE(i) = iA.PerFrameFunctionalGroupsSequence.(perFrameFields{i}).MREchoSequence.Item_1.EffectiveEchoTime;
end
clear perFrameFields

phs_scale = TE(ER)*1e-3 * gyro * B0;

% Magnitude and Phase containers. 4th dimension is echo.
MV=zeros([N ER]);
PV=zeros([N ER]);
RV=zeros([N ER]);
IV=zeros([N ER]);
for i=1:N(3)
    for j = 1:ER
    MV(:,:,i,j)=MagTot(:,:,j+(i-1)*ER);
    PV(:,:,i,j)=PhTot(:,:,j+(i-1)*ER)/1000;
    RV(:,:,i,j)=ReTot(:,:,j+(i-1)*ER);
    IV(:,:,i,j)=ImTot(:,:,j+(i-1)*ER);
    end
end
clear MagTot PhTot ReTot ImTot
CV = RV+1i*IV;

PV2 = angle(CV);
MV2 = abs(CV);
clear CV

% Data inspection
imagesc3d2(MV(:,:,:,1), N/2, 1, [90,90,0], [0,3140], 0, 'mag 1st echo');
imagesc3d2(MV2(:,:,:,1), N/2, 2, [90,90,0], [0,3140], 0, 'magComplex 1st echo');
imagesc3d2(MV(:,:,:,1)-MV2(:,:,:,1), N/2, 3, [90,90,0], [-50,50], 0, 'mag-magComplex 1st echo');

imagesc3d2(PV(:,:,:,1), N/2, 6, [90,90,0], [-3.14,3.14], 0, 'phase 1st echo');
imagesc3d2(PV2(:,:,:,1), N/2, 7, [90,90,0], [-3.14,3.14], 0, 'phaseComplex 1st echo');
imagesc3d2(angle( exp(1i*(PV(:,:,:,1)-PV2(:,:,:,1)))), N/2, 8, [90,90,0], [-0.03,0.03], 0, 'phase diff 1st echo');

save rawdata.mat MV MV2 PV PV2 N ER TE voxel_size phs_scale

% In the next sections we'll use the magnitude and phase calculated from the complex data.

%%-------------------------------------------------------------------------
%% R2* reconstruction
%%-------------------------------------------------------------------------

% Fast (log domain) reconstruction
[ M0, R2s, T2s ] = t2s( MV2, TE*1e-3, ones(size(MV2)) );
imagesc3d2(R2s, N/2, 11, [90,90,0], [0,100], 0, 'R2s [Hz]');
imagesc3d2(T2s*1000, N/2, 12, [90,90,0], [0,100], 0, 'T2s [ms]');
imagesc3d2(M0, N/2, 13, [90,90,0], [0,3140], 0, 'Mag_0 (PD)');

% Nonlinear TV regularized reconstruction
paramst2.input = MV2/max(MV2(:)); % ideally in the 0-1 range
paramst2.alpha1 = 1e-3;
paramst2.mu1 = 25*paramst2.alpha1;
paramst2.te = TE; %in ms
outt2_tv = nlT2sTV(paramst2);
imagesc3d2(outt2_tv.r2*1000, N/2, 14, [90,90,0], [0,100], 0, 'nR2s [Hz]');
imagesc3d2(1.0./(outt2_tv.r2+eps), N/2, 15, [90,90,0], [0,100], 0, 'nT2s [ms]');
imagesc3d2(outt2_tv.m0*max(MV2(:)), N/2, 16, [90,90,0], [0,3140], 0, 'nMag_0 (PD)');

% Please check also the additional functions included in this toolbox for T2*/R2* 
% mapping ('T2Reg' folder):
% nlT2sL2: Nonlinear, Tikhonov regularization
% nlT2sTVechoes: Nonlinear, with implicit TV regularization.
% nlT2sTGV: Nonlinear, with Total Generalized Variation regularization.
% nlT2sTGVechoes: Nonlinear, with implicit TGV regularization.
% Implicit regularizations are penalties applied to the projections M0.*exp(R2s*TE(t))
% instead than to the relaxometric parameters.

%%-------------------------------------------------------------------------
%% Mask creation
%%-------------------------------------------------------------------------
% Very simple mask creation, by magnitude thresholding and morphological operations
% We encourage to use dedicated tools as BET2 and others.

% First select where signal is above threshold
mask = single(outt2_tv.m0*max(MV2(:))>700);

% Use R2* refine selection
mask2 = single(abs(R2s)>3).*single(abs(R2s)<35).*mask;

% Select larger conected area
cc = bwconncomp(mask2);
numPixels = cellfun(@numel,cc.PixelIdxList); 
[biggest,idx] = max(numPixels); 
mask2 = zeros(N);
mask2(cc.PixelIdxList{idx}) = 1;
mask2 = imfill(mask2); % Fill holes

% Remove inferior 5 slices to avoid artifacts in the boundary
for i = 1:5
  mask2(:,:,i) = 0;
end

% Clean segmentation by using morphological filters
se = strel('sphere',3);
mask2=imdilate(mask2,se); % Aperture
mask2=imerode(mask2,se);
se = strel('sphere',10); % Closure
mask2=imerode(mask2,se);
mask2=imdilate(mask2,se);
se = strel('sphere',6); % Aperture
mask2=imdilate(mask2,se);
mask2=imerode(mask2,se);
se = strel('disk',2); % Erode 2 voxels to avoid data corruption in the cortical areas
mask2=imerode(mask2,se);
mask2=imerode(mask2,se); % Further Closure
mask2=imdilate(mask2,se);
cc = bwconncomp(mask2); % Just in case, select larger connected area again to remove any "island"
numPixels = cellfun(@numel,cc.PixelIdxList); 
[biggest,idx] = max(numPixels); 
mask2 = zeros(N);
mask2(cc.PixelIdxList{idx}) = 1;
mask2 = imfill(mask2);

imagesc3d2(mask,N/2, 18, [90,90,0], [0,1], 0, 'mask');
imagesc3d2(mask2,N/2, 19, [90,90,0], [0,1], 0, 'mask2');


%%-------------------------------------------------------------------------
%% Phase preprocessing: Unwrapping and Multi-echo fitting
%%-------------------------------------------------------------------------

% Phase unwrapping

% Simple solution: Best Path unwrapping - MEDI or QSMBox Toolboxes
%for i=1:ER
%[unw(:,:,:,i)] = unwrapPhase(MV2(:,:,:,i), PV2(:,:,:,i), N);
%end
%
% We recommend using SEGUE, by Anita Karsa, et al.
% http://dx.doi.org/10.1109/TMI.2018.2884093
% https://xip.uclb.com/i/software/SEGUE.html
% or ROMEO, by Barbara Dymerska, et al.
% https://github.com/korbinian90/ROMEO.jl/releases
% as more efficient alternatives for phase unwrapping.

for i=1:ER
Inputs.Mask = mask2;
Inputs.Phase = PV2(:,:,:,i);
[unw(:,:,:,i)] = SEGUE(Inputs);
end

% Avoid Laplacian-based unwrapping for multiecho processing, as they remove
% harmonic components of the phase.
% function [iFreq ] = unwrapLaplacian(iFreq_raw, matrix_size, voxel_size)
%
% The iterative approach by Schofield et al, based on the approximated laplacian solution
% often yields artifacts, but may work with moderated wrapping.
% function [ out ] = unwrap( phase, voxel_size )
%
% Gradient domain direct solvers should yield results similar to Schofield's Laplacian unwrapping.
% function [ kout uout ] = gduw( img, mask, vsize )
% function [ kout uout ] = gduw_unbound( img, vsize )
%
% Approximate better solutions to Laplacian operators may be found using the Nonlinear 
% iterative Laplacian unwrapping funtion:
% for i=1:ER
%     params = [];
%     params.input = PV2(:,:,:,i);
%     params.weight = MV2(:,:,:,i)/max(MV(:));
%     params.lambda = 1e-1;
%     params.voxel_size = voxel_size;
%     out = nluwl2(params);
%     unw(:,:,:,i) = out.psi;
% end


% Multi-echo fitting

% First estimation of field map and offset based on direct calculation
echo1 = 2; % We use two unwrapped echoes to get a first estimate of the phase offset and field map.
echo2 = 4;
% Initial estimations are made by linear fitting.
thetai = (TE(echo1)*unw(:,:,:,echo2)-TE(echo2)*unw(:,:,:,echo1))/(TE(echo1)-TE(echo2));
bi = (unw(:,:,:,echo2)-unw(:,:,:,echo1))/(TE(echo2)-TE(echo1));
imagesc3d2((thetai), N/2, 21, [90,90,0], [-1*3.14,1*3.14], 0, 'thetai');
imagesc3d2(bi*(TE(echo2)-TE(echo1)).*mask2, N/2, 22, [90,90,0], [-1*3.14,1*3.14], 0, 'bi');

% If needed, clean estimations and correct unwrapping errors

% Filter results to reduce errors
mb = medfilt3( medfilt3(bi) );
mt = medfilt3( medfilt3(thetai) );
% We are not looking for perfect estimates of b and theta at this stage. Just initial solutions as close as possible to 
% the ideal solutions, by eliminating wrapping artifacts.

% Another refinement to this estimation may be done by fitting a polynomial to the phase offset, inside the eroded mask region:
Order = 4;      % degree of 3d polyfit
mt = polyfit3D_NthOrder(mt, imerode(mask2,se), Order);
imagesc3d2(mt.*mask2, N/2, 23, [90,90,0], [-pi,pi], 0, 'polynomial fit to mt');

% With this smooth phase offset, we may recalculate the field map using all unwrapped echoes to improve CNR.
mb = zeros(N);
den = zeros(N);
for i=1:ER
   mb = mb+(unw(:,:,:,i)-mt).*MV2(:,:,:,i)/max(MV(:));
   den = den + TE(i)*MV2(:,:,:,i)/max(MV(:));
end
mb = mb./(den+eps);
imagesc3d2(mb*(TE(echo2)-TE(echo1)), N/2, 24, [90,90,0], [-1*3.14,1*3.14], 0, 'bc');

% Please note that we are ignoring conductivity related effects in the phase offset, that are reintroduced later
% in the nonlinear solver.


% Nonlinear multi-echo fit
for i=1:ER
   W(:,:,:,i) = MV2(:,:,:,i)*TE(i); % Weight to optimize CNR in the phase data
end
W = W/max(W(:));

paramsMultiEcho = [];
paramsMultiEcho.input = unw; % wrapped phase data may be used in this function
paramsMultiEcho.weight = W;
paramsMultiEcho.lambda = 3e-4; % use a small regularization weight for denoising purposes
paramsMultiEcho.mu = 1e0;
paramsMultiEcho.b0 = mb; % Use the initial field map and offset to initializate the nonlinear solver
paramsMultiEcho.theta0 = mt; % Avoid introducing wrapping artifacts.
paramsMultiEcho.TE = TE;
paramsMultiEcho.tol_update = 1e-2;
paramsMultiEcho.maxOuterIter = 500;
outc = nlme_L2gradient(paramsMultiEcho); % this function imposes smooth gradients to the phase offset.
% outc = nlme_tik(paramsMultiEcho); % use this function to impose Tikhonov constrains to the phase offset.

% Please note that these nonlinear functions give an estimation of the noise in the phase and
% complex domains. You may use these noise estimates to update the weigths

imagesc3d2(outc.theta, N/2, 25, [90,90,0], [-1*3.14,1*3.14], 0, 'thetac');
imagesc3d2(outc.b*(TE(echo2)-TE(echo1)), N/2, 26, [90,90,0], [-1*3.14,1*3.14], 0, 'bc');

imagesc3d2(outc.lnoise, N/2, 27, [90,90,0], [0,1], 0, 'lnoise');
imagesc3d2(outc.cnoise, N/2, 28, [90,90,0], [0,0.1], 0, 'cnoise');
%%-------------------------------------------------------------------------
%% Phase preprocessing: Background field removal
%%-------------------------------------------------------------------------

% You may use external functions like LBV, PDF or VSHARP. 
% Here you'll find examples on how to use the background field removal functions in this toolbox as well.

% % LBV filtering example - MEDI Toolbox
% resultLBV = LBV(outc.b,mask2,N,voxel_size);
% imagesc3d2(resultLBV, N/2, 31, [90,90,0], [-0.12,0.12], 0, 'LBV result');
% % A Polinomial fit is often required to remove large scale remnants:
% Order = 4;      % degree of 3d polyfit
% I_fitted = polyfit3D_NthOrder(resultLBV, mask2, Order); % Included in FANSI-Toolbox
% filteredLBV = resultLBV - I_fitted .* mask2;
% imagesc3d2(filteredLBV*(TE(echo2)-TE(echo1)), N/2, 32, [90,90,0], [-1.12,1.12], 0, 'filtered LBV');
%
% % PDF example - MEDI Toolbox
% [localPDF xPDFb] = PDF(outc.b, 1./(1000*MV2(:,:,:,1)+eps), mask2,N,voxel_size, [0 0 1]);
% imagesc3d2(localPDF*(TE(echo2)-TE(echo1)), N/2, 33, [90,90,-90], [-0.12,0.12], 0, 'PDF result');
% You should get similar results with the "wPDF" function included in the FANSI-Toolbox.
% See the nonlinear example below to set the parameters.


% Calculate the dipole kernel
kernel = dipole_kernel_fansi( N, voxel_size, 0 );

% Nonlinear PDF example
ppdf = [];
ppdf.maxOuterIter = 500;
ppdf.tol_update = 0.01;
ppdf.K = kernel;
ppdf.input = outc.b*(TE(echo2)-TE(echo1)); % convert to radians
ppdf.mask = mask2;
ppdf.weight = outt2_tv.m0.*mask2;
ppdf.mu1 = 0.0025;                  % gradient consistency
ppdf.alpha1 = eps;%0.0001;              % gradient L1 penalty
nout = nlPDF(ppdf);
imagesc3d2(nout.local.*mask2, N/2, 34, [90,90,0], [-1.12,1.12], 0, 'PDF result');
% If you have a mask that closelly matches all the tissues (rejects air), you may use it to preconditionate
% the solution with a first estimation of the magnetization due to those interfaces.
% ppdf.outermask = mask;
% noutp = nlPDFp(ppdf);
% If you also would like to introduce linear gradients to this model, you may use:
% noutp2 = nlPDFph(ppdf);
% Please also see the description to the "fitmodels" and "backmodel" functions.

% Multiscale Spherical Mean Value 
% 
% This function yields similar results to VSHARP, but in a more efficient computational way.
[ local ] = msmv_noextrapolation( outc.b, mask2, voxel_size,5 );
imagesc3d2(local*(TE(echo2)-TE(echo1)).*mask2, N/2, 35, [90,90,0], [-1.12,1.12], 0, 'MSMV result');

b2  = extrapolatefield_iter( outc.b, mask2,2.5e-2,250 );
[ local2 ] = msmv_noextrapolation( b2, mask2, voxel_size,5 );
imagesc3d2(b2*(TE(echo2)-TE(echo1)), N/2, 36, [90,90,0], [-1*3.14,1*3.14], 0, 'bc');
imagesc3d2(local2*(TE(echo2)-TE(echo1)).*mask2, N/2, 36, [90,90,0], [-1.12,1.12], 0, 'MSMV result');
% For a slightly better estimation (but computationaly more expensive) near the boundaries of the mask, you may try:
[ local3 ] = msmv( outc.b, mask2, voxel_size,5 );
imagesc3d2(local3*(TE(echo2)-TE(echo1)), N/2, 37, [90,90,0], [-1.12,1.12], 0, 'MSMV result');

% A SHARP deconvolution step may be required
layers = 5; % Number of layers used.
dkernel = 1-gauss_kernel( N, voxel_size, 2^(layers-1) );
dreg = 0.05; % Deconvolution regularization
msharp = ifftn(conj(dkernel).*fftn(local2)./(conj(dkernel).*(dkernel)+dreg));
imagesc3d2(msharp*(TE(echo2)-TE(echo1)), N/2, 38, [90,90,0], [-1.12,1.12], 0, 'MSHARP result');

% Alternativelly, you may incorporate the deconvolution kernel into the QSM problem:
% kernel = kernel.*dkernel; % where "kernel" is the dipole kernel.

%%-------------------------------------------------------------------------
% QSM reconstruction (FANSI Toolbox: wTV and WH-QSM(TV))
%%-------------------------------------------------------------------------

% If not done before, calculate the dipole kernel.
kernel = dipole_kernel_fansi( N, voxel_size, 0 );

% Redefine the phase-to-ppm factor
deltaTE = (TE(echo2)-TE(echo1)); % You may use also any individual echo time
phs_scale = deltaTE * gyro * B0;

% Set main parameters
params = [];
params.K = kernel;
params.input = local*deltaTE;
params.weight = outt2_tv.m0.*mask2./outc.cnoise; % or mask2./outl.cnoise;
params.weight = params.weight/max(params.weight(:); % normalization for stability


% Find optimal parameters via L-Curve analysis or Frequency Analysis

% Create frequency masks
[ m1, m2, m3 ] = create_freqmasksE( N, voxel_size, kernel );

alpha = 10.^(-(1:20)/5-1);
for ka = 1:20
    
    params.alpha1 = alpha(ka);
    out = wTV(params); 
    
    % Calculate the costs associated to each term of the optimization functional
    [ data_cost, reg_cost ] = compute_costs( out.x.*mask2, params.input.*mask2, kernel );    
    dcw(ka) = data_cost;
    rcw(ka) = reg_cost;
    
    % Calculate the Mean Amplitudes for each mask
    [e1(ka), e2(ka), e3(ka)] = compute_freqe(out.x.*Mask,m1,m2,m3);
end

% L-Curve analysis
[ Kappa ] = draw_lcurve( alpha, rcw, dcw, 11 );
[ KappaM ] = draw_lcurve_median( alpha, rcw, dcw, 12 );
% The point where the curvature is zero (from right to left) is usually more visually appealing than 
% finding the maximum curvature

% Frequency analysis
[opt23index, alpha_opt] = draw_freque(alpha,e1,e2,e3,33);
% Best results are achieved when the mask amplitudes for masks m2 and m3 are equal.

% Typically, the L-curve analysis will give a slightly over-regularized solution, whereas the frequency analysis
% yields slightly under-regularized solutions. You may manually fine tune in the range proposed by both methods.
% Alternativelly, you may use a factor between 1.5-2.0 on alpha_opt to increase regularization.

% Example optimal reconstruction
    params.alpha1 = 1.5*alpha(opt23index); % Slightly over-regularize to avoid noise amplification
    out = wTV(params); 
imagesc3d2(out.x/phs_scale, N/2, 41, [90,90,-90], [-0.1,0.1], [], 'wTV')

% You may repeat the same analysis for other reconstruction functions.
% If you don't change the input phase scale, or weight, you may use the same regularization weight.

% You may find functions in the root of the FANSI function to solve this problem with:
% - Linear or nonlinear data fidelity terms (w and nl prefix).
% - Total Variation or Total Generalized Variation (TV and TGV sufix) regularization.
% Example:
outnl = nlTGV( params ); % Nonlinear Total Generalized Variation. Use same parameters as before.


% Weak Harmonics (WH-QSM) example
%
% In the case that background field remnants may be corrupting the calculation of the susceptibility map,
% the Weak Harmonics QSM method provides with a joint estimation of those remnnats fields.
% The TV/TGV regularization parameters may be calculated using the standard FANSI methods, for simplicity.
% You may then decrease slighlty the regularization weight for WH-QSM to avoid over-smoothing.

params.mask = mask2;
params.beta = 200;
params.muh = 5; % Fine tune these parameters. Usually, beta is between 20 to 50 times larger than muh for optimal results.
                % Too large muh values may result in understimation of the remnant fields.
                % Too small muh values may result in understimation of susceptibility values.
out = WH_nlTV(params); % Weak Harmonics with nonlinear data fidelity term and TV regularization.


% Robust to Phase Inconsistencies QSM (PI-QSM)
%
% This method uses an L1-norm data fidelity term to reject extreme noise or errors in the phase, thus preventing streaking artifacts.
% Both the linear and nonlinear versions require fine-tuning of the regularization weight. Weights from FANSI or WH-QSM will not
% produce optimal results. Please run the L-curve and/or Frequency Analysis again.
% In addition, a scale factor may be needed to modify the data fidelity weight (params.weight):
lambda = 0.8;
params.weight = outt2_tv.m0.*mask2./outc.cnoise; % or mask2./outl.cnoise;
params.weight = lambda*params.weight/max(params.weight(:);
                % Smaller lambda values reject more voxels with inconsistencies. Too small values may yield inpainting effects. 
                % Too large values will give similar results to FANSI methods.
outL1nl = nlL1TV(params); % or wL1TV(params) of the linear version. TGV regularization is not yet implemented.


%% Early stopping methods
% It was shown (Polak et al) that a nonlinear functional may avoid regularization if a gradient descent algorithm is used, 
% by stopping the algorithm before it diverges. Thus, the number of iterations become the main parameter.

paramsNDI = [];
paramsNDI.K = kernel;
paramsNDI.input = local*deltaTE;
paramsNDI.weight = outt2_tv.m0.*mask2./outc.cnoise; % or mask2./outl.cnoise;
paramsNDI.weight = params.weight/max(params.weight(:); % normalization for stability

params.maxOuterIter = 300; % A gradient descent algorithm may be slow
params.show_iters = true; % Inspect the reconstructions, at each iteration. Stop with Ctrl+C if the quality degrades.
outNDI = ndi(params); 
% Modify the maxOuterIter parameter with the aparent optimal number of iterations


% A more efficient way to solve the functional is using a Conjugate Gradient descent.
params.maxOuterIter = 50; % Fewer iterations are needed
params.show_iters = true; % Inspect the reconstructions, at each iteration. Stop with Ctrl+C if the quality degrades.
outNDI = ndiCG(params); 
% Modify the maxOuterIter parameter with the aparent optimal number of iterations (typically between 12-20).

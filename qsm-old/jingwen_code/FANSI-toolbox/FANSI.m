function out = FANSI( phase, magn, spatial_res, alpha, noise, options, B0_dir )
% FANSI wildcard function. This may be used to call the included functions
% globally, for a simplified use.
%
% Parameters: 
% phase: local field map data
% magn: magnitude data
% spatial_res: Spatial resolution vector, in mm, or normalized (mean = 1).
% alpha: gradient L1 penalty, regularization weight
% noise: noise standard deviation in the complex signal
% options (structure with fields):
% options.iterations:
% options.update
% options.nonlinear: linear or nonlinear algorithm? (default = true, i.e nonlinear method)
% options.tgv: TV or TGV regularization? (default = false, i.e. TV regularization)
% options.kernel_mode: 0 for the continuous kernel proposed by Salomir, et al. 2003.
%                      1 for the discrete kernel proposed by Milovic, et al. 2017.
%                      2 for the Integrated Green function proposed by Jenkinson, et al. 2004
% options.gradient_mode: 0 to use the vector field. 
%                        1 for the L1 norm, and 
%                        2 for the L2 norm
% options.mu: ADMM Lagrange multiplier, regularization term.
%
% B0_dir: main field direction, e.g. [0 0 1]
%
% Created by Carlos Milovic, 30.03.2017
% Modified by Julio Acosta-Cabronero, 26.05.2017
% Last modified by Carlos Milovic, , 07.07.2020
isnonlinear = true;
istgv = false;
params = [];

if nargin > 5
    if isfield(options,'kernel_mode')
         kmode = options.kernel_mode;
    else
        kmode = 0;
    end
    if isfield(options,'iterations')
         params.maxOuterIter = options.iterations;
    end
    if isfield(options,'update')
         params.tol_update = options.update;
    end
    
    if isfield(options,'gradient_mode')
        gmode = options.gradient_mode;
        Gm = gradient_calc(magn,gmode); 
        Gm = max(Gm,noise); % Binary weighting not implemented here
        params.regweight = mean(Gm(:))./Gm;
    end
    if isfield(options,'nonlinear')
         isnonlinear = options.nonlinear;
    end
    if isfield(options,'tgv')
         istgv = options.tgv;
    end
else
        kmode = 0;
end

if nargin > 6
params.K = dipole_kernel_angulated( size(phase), spatial_res, B0_dir ); 
else
params.K = dipole_kernel_fansi( size(phase), spatial_res, kmode ); 
end
    

    params.input = phase;
    params.weight = magn; 
 
    params.alpha1 = alpha;            % gradient L1 penalty
    if isfield(options,'mu')
         params.mu1 = options.mu;
    end



if isnonlinear == true
    if istgv == true
        out = nlTGV(params);
    else
        out = nlTV(params);
    end
    
else
    if istgv == true
        out = wTGV(params);
    else
        out = wTV(params);
    end
        
end



end


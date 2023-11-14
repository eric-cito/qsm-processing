function res = qsip_init()
% param = init()
%
% function returns a structure with the entries that are needed for the reconstruction.
% The user MUST supply operators afterwords!
%
% See:
%	demo_SheppLoganTV demo_angio_simulation demo_Brain_2D.m
% 
% (c) Michael Lustig 2007

res.K = []; % The measurement operator (undersmapled Fourier for example)
res.XFM = []; % Sparse transform operator
res.TV = []; 	% the Total variation operator
res.data = []; % measurements to reconstruct from
res.FT0 = [];

%res.TVWeight = 0.01;	% TV penalty
%res.xfmWeight = 0.01;   % transform l1 penalty

res.Itnlim = 20;	% default number of iterations
res.gradToll = 1e-30	% step size tollerance stopping criterea (not used)

res.l1Smooth = 1e-15;	% smoothing parameter of L1 norm
res.pNorm = 1;  % type of norm to use (i.e. L1 L2 etc)

% line search parameters
res.lineSearchItnlim = 150;
res.lineSearchAlpha = 0.01; %default
%res.lineSearchAlpha = 0.001;
res.lineSearchBeta = 0.6;   %default
%res.lineSearchBeta = 0.75;

res.lineSearchT0 = 1 ; % step size to start with

%termination criteria
res.term = [];

res.Range = [];            % range of allowable susceptibility (chi1). range = upper bound - lower bound
res.Lowerbd = [];          % lower bound on susceptibility labels (chi1).

res.mag = [];
res.mag_tissue_mean = [];

res.QUADa = 0;      %coefficient for quadratic penalty on susc values outside allowed range
res.QUADb = 0;      %coefficient for quadratic penalty on susc values outside allowed range
res.QUADc = 0;      %coefficient for quadratic penalty on susc values outside allowed range

res.l1prior = 0; 

res.fmap_acq = [];    % the acquired fieldmap
res.model_err = 0;   % error between a fmap computed from the true susc map and the acquired fmap due to having an incomplete antomical model to create the true susc map (ie. missing slices outside fov)
res.LPWeight = 0;
res.KSWeight3 = 0;
res.KSWeight3_constant = [];
res.WVWeight = 0;
res.TVWeight = 0.0; 	% Weight for TV penalty
res.OUTWeight = 0;    % penalty for susceptibility outside the head (This is alpha in Eq 9 of cornell paper, not alpha^2)
res.GDWeight = 0.0;     % penalty for large gradients in susceptibility (This is beta in Eq 9 of cornell paper, not beta^2)
res.KSWeight = 0; %weight var comes from testing script.    % penalty for disagreement between fourier transform of normalized estimated susc map and fourier transform of normalized magnitude image (hypothesis is that they have similar spatial freq structure)
res.KSWeight2 = 0;   % penalty for disagreement between the estimated fmap in kspace and true fmap in kspace (ie. a fourier trans of susc estimate and mult by kernel is compared to true fmap in ks)
res.QUADWeight = 0;
res.xfmWeight = 0;

res.mask = [];
res.OUTmask = [];
res.GDmask = [];
res.scales = [];
res.chi1_max = [];
res.noise_term = [];
res.chi1_air = [];
res.chi1_tissue = [];
res.ind_gd = [];
res.chi1_mask = [];
res.c1 = [];
res.c2 = [];
res.x0 = [];
res.PFWeight = 0;
res.chi0 = [];
res.chi0_sparse = [];       %-9.05e-6; 
res.delta = [];
res.delta_sparse = [];      % 1e-6;   
res.kmask = [];
res.noise_term = [];
res.ind_gd = [];
function [x, f1, obj_vals] = fnlCg_sc_lp(params)

%-----------------------------------------------------------------------

% set line search parameters
x = params.x0;                             
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha;     
beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;
t = 1;
termCrit = params.term;              

% set other parameters
mask2 = params.mask2;
mask3 = params.mask3;
mask1 = params.mask1;
obj_vals = [];
scales = params.scales;


if params.OBJWeight

    FTXFM_fmap_acq = params.K*(mask2 .* params.fmap_acq);

else
    FTXFM_fmap_acq = 0;
end

if params.LPWeight     %calc laplacian of acquired fmap - apply 3D laplacian
    
    LpFMAP = 0 * x;
    
    fmap_acq = params.fmap_acq;
                                    
    LpFMAP(2:end-1, 2:end-1, 2:end-1) = params.lp_cc * fmap_acq(2:end-1, 2:end-1, 2:end-1) + fmap_acq(1:end-2, 2:end-1, 2:end-1) + fmap_acq(3:end, 2:end-1, 2:end-1) ...
                                        + fmap_acq(2:end-1, 1:end-2, 2:end-1) + fmap_acq(2:end-1, 3:end, 2:end-1) + fmap_acq(2:end-1, 2:end-1, 1:end-2) + fmap_acq(2:end-1, 2:end-1, 3:end);   
else
    LpFMAP = 0;
end


%% start optimization alg 
 
g0 = wGradient(x, params, mask2, LpFMAP, mask3, mask1, FTXFM_fmap_acq);

dx = -g0;

% iterations
while(1)

% backtracking line-search

	% pre-calculate values, such that it would be cheap to compute the objective many times for efficient line-search
	[FTXFMtx, FTXFMtdx, LpFTx, LpFTdx] = preobjective(x, dx, params);
    
    % note t = 0 to compute f0
    f0 = objective(FTXFMtx, FTXFMtdx, x, dx, 0, params, mask2, LpFTx, LpFTdx, LpFMAP, mask3, mask1);
        
    t = t0;
        [f1, ERRobj, RMSerr, ABSerr]  =  objective(FTXFMtx, FTXFMtdx, x, dx, t, params, mask2, LpFTx, LpFTdx, LpFMAP, mask3, mask1);
	
	lsiter = 0;             %line search iteration number

    while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)    
		lsiter = lsiter + 1;
		t = t * beta;
        [f1, ERRobj, RMSerr, ABSerr]  =  objective(FTXFMtx, FTXFMtdx, x, dx, t, params, mask2, LpFTx, LpFTdx, LpFMAP, mask3, mask1);
    end
    
	if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2
		t0 = t0 * beta;
	end 
	
	if lsiter<1
		t0 = t0 / beta;
    end
    
	%--------- uncomment for debug purposes ------------------------	
	disp(sprintf('Iter. %d   , obj: %f, L-S: %d', k, f1, lsiter));
    obj_vals = [obj_vals, f1];
    
    if ~(mod(k,10))
        iter = num2str(k/10);
        name = strcat('susc_est_ppm_',iter);
        susc_est_ppm = (params.chi0 * ones(size(mask2)) + params.delta * x) * 10^6;
        save_avw(susc_est_ppm, name,'f',scales)
    end   
	%---------------------------------------------------------------
  	
    % check if stopping criteria are met
	if (k > params.Itnlim) | (norm(dx(:)) < gradToll)  
        disp('stopping criteria: k > params.Itnlim) | (norm(dx(:)) < gradToll')
		break;
    end
    
    if f1 < termCrit
        disp('termination criteria met: total obj fn (f1) < termCrit')
        break;
    end  
        
	x = (x + t*dx);
	
    % conjugate gradient calculation
    g1 = wGradient(x, params, mask2, LpFMAP, mask3, mask1, FTXFM_fmap_acq);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  -g1 + bk* dx;
	k = k + 1;
     
end

function [FTXFMtx, FTXFMtdx, LpFTx, LpFTdx] = preobjective(x, dx, params)
% This function precalculates transforms to make line search cheap

if params.OBJWeight
    FTXFMtx = params.K*x;       
    FTXFMtdx = params.K*dx;     
else
    FTXFMtx = 0;
    FTXFMtdx = 0;
end

if params.LPWeight
    
    %if this hasn't already been calculated, then compute it. 
    if ~(params.OBJWeight)
        fmap_est = params.K*x;
        fmap_est_dx = params.K*dx;
    else
        fmap_est = FTXFMtx;
        fmap_est_dx = FTXFMtdx;
    end
    
    %calc laplacian of estimated fmap - apply 3D laplacian 
    LpFTx = 0 * x;
    LpFTdx = 0 * x;
    
    LpFTx(2:end-1, 2:end-1, 2:end-1) = params.lp_cc * fmap_est(2:end-1,2:end-1, 2:end-1) + fmap_est(1:end-2,2:end-1, 2:end-1) + fmap_est(3:end,2:end-1, 2:end-1) ...
                                        + fmap_est(2:end-1,1:end-2, 2:end-1)+ fmap_est(2:end-1,2:end-1, 1:end-2) + fmap_est(2:end-1,2:end-1, 3:end);
    LpFTdx(2:end-1, 2:end-1, 2:end-1) =  params.lp_cc * fmap_est_dx(2:end-1,2:end-1, 2:end-1) + fmap_est_dx(1:end-2,2:end-1, 2:end-1) + fmap_est_dx(3:end,2:end-1, 2:end-1) ...
                                        + fmap_est_dx(2:end-1,1:end-2, 2:end-1) + fmap_est_dx(2:end-1,3:end, 2:end-1) + fmap_est_dx(2:end-1,2:end-1, 1:end-2) + fmap_est_dx(2:end-1,2:end-1, 3:end);                                                                           
else  
    LpFTx = 0;
    LpFTdx = 0;   
    LpFMAP = 0;    
end
       

function [res, obj, RMS, ABSerr] = objective(FTXFMtx, FTXFMtdx, x, dx, t, params, mask2, LpFTx, LpFTdx, LpFMAP, mask3, mask1)
%calculated the objective function
%
%my comment: obj = F* W' *x - y
%            let yhat = F* W' *x

ABSerr = [];

p = params.pNorm;
    
if params.OBJWeight
    
    obj = FTXFMtx + t*FTXFMtdx - params.data;  

    obj = mask2 .* (FTXFMtx + t*FTXFMtdx - params.data);   %my commment: new value of yhat = current value + t* change_yhat_change_x
 
    obj = obj(:)'*obj(:);                       % ||y_hat_new - y||^2  =  ||F* W' *x - y||^2
    
else
    obj = 0;   

end
    

if params.LPWeight
    
    w =  mask1(:) .* (LpFMAP(:) - (LpFTx(:) + t*LpFTdx(:)) - params.model_err(:)); 
    LP =  (w .* conj(w) + params.l1Smooth).^(p/2);   % p = 1 for L1 norm  
else
    LP = 0;
end

if params.EXTWeight
    
    w = x + t*dx;
    
    EXT = mask3(:) .* (params.chi1_null_vol(:) - w(:));  
 
    EXT = EXT(:)'*EXT(:);        
    
else
    EXT = 0;
end

if params.OBJWeight
    RMS = sqrt(obj/sum(abs(params.data(:))>0));
    ABSerr = obj^.5;
else
    RMS = [];
end

LP = sum(LP.*params.LPWeight(:));
 
EXT = EXT * params.EXTWeight;

OBJ = obj * params.OBJWeight;

res = OBJ + (LP) + (EXT);


function grad = wGradient(x, params, mask2, LpFMAP, mask3, mask1, FTXFM_fmap_acq)

gradObj = 0;
gradLP = 0;
gradEXT = 0;

if params.OBJWeight
    gradObj = gOBJ(x, params, mask2, FTXFM_fmap_acq); 
end

if params.LPWeight
    gradLP = gLP(x,params,LpFMAP,mask1);
end

if params.EXTWeight
    gradEXT = gEXT(x, params, mask3); 
end

grad = params.OBJWeight .* gradObj + (params.LPWeight).*gradLP + (params.EXTWeight) .* gradEXT;


function gradObj = gOBJ(x, params, mask2, FTXFM_fmap_acq)
% computes the gradient of the data consistency (L2 norm squared penalty).
% note FTXFM_fmap_acq is computed at start of this code. 
    %gradObj = params.FT * (params.FT' * (mask2 .* (params.FT * (params.FT'*x)))) - (params.FT * (params.FT'*(mask2 .* params.data)));
    
    %gradObj = params.K * (mask2 .* (params.K*x)) - (params.K*(mask2 .* params.data));
    gradObj = params.K * (mask2 .* (params.K*x)) - FTXFM_fmap_acq;
    gradObj = 2*gradObj;
    
function gradLP = gLP(x,params,LpFMAP,mask1)
% This  function computes the gradient of the L1 norm of (laplacian of acq fmap - laplacian of
% estimated fmap) via several applications of the chain rule

p = params.pNorm;   %this parameter is set to 1 in init.m for the L1 norm, 2 for the L2 norm.
 
fmap_est = params.K*x;

%calc laplacian of estimated fmap - apply 3D laplacian 
LpFTx = 0 * x;

LpFTx(2:end-1, 2:end-1, 2:end-1) = params.lp_cc * fmap_est(2:end-1,2:end-1, 2:end-1) + fmap_est(1:end-2,2:end-1, 2:end-1) + fmap_est(3:end,2:end-1, 2:end-1) ...
                                        + fmap_est(2:end-1,1:end-2, 2:end-1)+ fmap_est(2:end-1,2:end-1, 1:end-2) + fmap_est(2:end-1,2:end-1, 3:end);

w =  mask1 .* (LpFMAP - LpFTx); 

w2 = w ./ ( (w .* conj(w) + params.l1Smooth).^(1-p/2) );  %this is approximately equal to w when p = 1

k = params.K*w2;

LP_k = 0 * w;
LP_k(2:end-1, 2:end-1, 2:end-1) = params.lp_cc * k(2:end-1,2:end-1, 2:end-1) + k(1:end-2,2:end-1, 2:end-1) + k(3:end,2:end-1, 2:end-1) + k(2:end-1,1:end-2, 2:end-1) + k(2:end-1,3:end, 2:end-1) + k(2:end-1,2:end-1, 1:end-2) + k(2:end-1, 2:end-1, 3:end);

gradLP = -1 * p * mask1 .* LP_k;

function gradEXT = gEXT(x, params, mask3)

gradEXT = -1*mask3 .* mask3.*(params.chi1_null_vol - x);  





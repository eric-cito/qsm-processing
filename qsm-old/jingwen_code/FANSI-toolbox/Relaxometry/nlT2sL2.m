% Nonlinear T2/T2s mapping with a Total Variation regularization.
% This uses ADMM to solve the functional.
%
% Created by Carlos Milovic in 2018.05.025
%
function out = nlT2sL2(params)
%
% input: params - structure with the following required fields:
% params.input - magnitude tensor [N x TE]
% params.alpha1 - gradient L1 penalty, regularization weight
% params.mu1 - gradient consistency weight
% params.te - echo time vector in miliseconds
%                 and the following optional fields:
% params.mu2 - fidelity consistency weight (recommended value = 1.0)
% params.maxOuterIter - maximum number of ADMM iterations (recommended = 50)
% params.tol_update - convergence limit, change rate in the solution (recommended = 1.0)
% params.delta_tol - convergence tolerance, for the Newton-Raphson solver (recommended = 1e-6)
% params.weight - data fidelity spatially variable weight (recommended = magnitude_data). Not used if not specified
% params.regweight - regularization spatially variable weight. Not used if not specified
% params.N - array size
% params.precond - preconditionate solution (for stability)
%
% output: out - structure with the following fields:
% out.m0 - magnitude at excitation
% out.r2 - 1/T2s
% out.iter - number of iterations needed
% out.time - total elapsed time (including pre-calculations)
%


tic

mu = params.mu1;
lambda = params.alpha1;

if isfield(params,'mu2')
    mu2 = params.mu2;
else
    mu2 = 1.0;
end

if isfield(params,'N')
    N = params.N;
else
    n = size(params.input);
    N = n(1:3);
end

if isfield(params,'maxOuterIter')
    num_iter = params.maxOuterIter;
else
    num_iter = 50;
end

if isfield(params,'tol_update')
    tol_update  = params.tol_update;
else
    tol_update = 1;
end


if isfield(params,'weight')
    weight = params.weight;
    
else
    weight = ones(size(params.input));
end

if isfield(params,'regweight')
    regweight = params.regweight;
    if length(size(regweight)) == 3
        regweight = repmat(regweight,[1,1,1,3]);
    end
else
    regweight = ones([N 3]);
end

if ~isfield(params,'delta_tol')
    delta_tol = 1e-6;
else
    delta_tol = params.delta_tol;
end

te = params.te;
Ne = length(te);
W = weight.*weight;
clear weight;
magn = params.input.*W;


r2 = zeros(N, 'single')+1e-6;
m0 = 1.2*magn(:,:,:,1);%zeros(N, 'single'); %no need to precond

if isfield(params,'precond')
    precond = params.precond;
else
    precond = true;
end

if precond
    z2 =  zeros(N,'single')/50; %better without precond
else
    z2 = zeros(N,'single');
end
s2 = zeros(N,'single');


%tic
for t = 1:num_iter
    % update r2 
    
    r2_prev = r2;
    r2 = max( ( mu2*(z2-s2) )./( eps + mu2 + lambda ) ,1e-5) ;
    
    r2_update = 100 * norm(r2(:)-r2_prev(:)) / norm(r2(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(r2_update)])
    
    if r2_update < tol_update
        break
    end
    
%imagesc3d2(r2, N/2, 1, [90,90,-90], [0,0.1], 0, 'r2');
    if t < num_iter
        
        % update m0
        m0 = 0;
        denom = 0;
        for echo = 1:Ne
            m0 = m0 + magn(:,:,:,echo).*exp(-te(echo)*z2);
            denom = denom + W(:,:,:,echo).*exp(-2*te(echo)*z2);
        end
        m0 = max(m0./(denom+eps),0.0);        
        
%imagesc3d2(m0, N/2, 2, [90,90,-90], [0,0.25], 0, 'm0');
%imagesc3d2(z2, N/2, 3, [90,90,-90], [0,0.1], 0, 'z2old');
        
        
        %z2 =  rhs_z2 ./ mu2 ;        
        % Newton-Raphson method for z2
        delta = inf;
        inn = 0;
        while (delta > delta_tol && inn < 50)
            inn = inn + 1;
            norm_old = norm(z2(:));
            
            up = zeros(N,'single');
            down = 0;
            for echo = 1:Ne                
                up = up - te(echo)*exp(-te(echo)*z2).*(W(:,:,:,echo).*m0.*exp(-te(echo)*z2)-magn(:,:,:,echo));
                down = down+te(echo)*te(echo)*exp(-te(echo)*z2).*( 2*W(:,:,:,echo).*m0.*exp(-te(echo)*z2)-magn(:,:,:,echo) );
            end
        %disp(sum(up(:)))
        %disp(sum(down(:)))
            
            update = (m0.*up+mu2*(z2-s2-r2))./(m0.*down+mu2+eps);
            z2 = max(z2 - update ,1e-5);
            delta = norm(update(:)) / norm_old;
        end
        disp(delta)
        
        s2 = s2 + r2 - z2;
    end
end
out.time = toc;toc

out.r2 = r2;
out.m0 = m0;
out.iter = t;


end

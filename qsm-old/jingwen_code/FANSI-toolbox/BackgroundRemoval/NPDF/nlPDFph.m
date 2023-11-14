function out = nlPDFph(params)
% Preconditioned Hybrid Projection onto Dipole Fields (wPDF), for background field removal.
% A Tikhonov term is included, for stability.
% This uses ADMM to solve the functional.
% A rough air-tissue model with linear gradients is fitted to facilitate convergence 
% and accuracy of the algorithm. An hybrid data fidelity term starts from a linear L2-norm projection,
% and is mixed with the nonlinear L2-norm projection as function of the iteration number.
%
% Parameters: params - structure with 
% Required fields:
% params.input: wrapped total field map
% params.K: dipole kernel in the frequency space
% params.mask: ROI binary mask
% Optional fields:
% params.outermask: binary mask that selects all tissues (not air, and no erosion).
% params.ph_unwrap: unwrapped total field map
% params.airmodel: precalculated magnetization field due to air-tissue interfaces.
% params.alpha1: Tikhonov regularization weight (use small values)
% params.mu1: gradient consistency weight (ADMM weight, recommended = 1.8e-3)
% params.mu2: fidelity consistency weight (ADMM weight, recommended value = 1.0)
% params.maxOuterIter: maximum number of iterations (recommended = 150)
% params.tol_update: convergence limit, update ratio of the solution (recommended = 0.1)
% params.weight: data fidelity spatially variable weight (recommended = magnitude_data). 
% params.precond: preconditionate solution by smart initialization
%
% Output: out - structure with the following fields:
% out.phi: background field map
% out.local: local field map
% out.x: external susceptibility sources
% out.iter: number of iterations needed
% out.time: total elapsed time (including pre-calculations)
%
% Modified by Carlos Milovic in 2018.10.15
% Last modified by Carlos Milovic in 2020.07.15


tic

    % Required parameters
    mask = params.mask;
    kernel = params.K;
    
    % Optional parameters
    if isfield(params,'outermask')
         outermask = params.outermask;
    else
        outermask = mask;
    end
    
    if isfield(params,'alpha1')
         lambda = params.alpha1;
    else
        lambda = 1e-5;
    end
    
    if isfield(params,'mu1')
         mu = params.mu1;
    else
        mu = 0.0018;
    end

    if isfield(params,'mu2')
         mu2 = params.mu2;
    else
        mu2 = 1.0;
    end
    
    N = size(params.input);

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
        weight = ones(N);
    end
    weight = weight.*weight;
    
    if ~isfield(params,'delta_tol')
        delta_tol = 1e-6;
    else
        delta_tol = params.delta_tol;
    end
    
    
    % Redefinition of variable for computational efficiency
    Wy = (weight.*params.input./(weight+mu2));

x = zeros(N, 'single');
z = zeros(N, 'single');
s = zeros(N, 'single');


kernel = params.K;

if isfield(params,'precond')
    precond = params.precond;
else
    precond = true;
end

background = 0;
if precond
    
    if isfield(params,'ph_unwrap')
        z2  = params.ph_unwrap;
    else
        %[z2, lapunwrap, delta] = iterlap_unwrap( params.input, mask, weight, 250 );
        z2 = unwrap(params.input.*mask, [1 1 1]);
    end
    
    phase_unwrapped = z2;
    
    imagesc3d2(mask.*z2, N/2, 106, [90,90,-90], [-30,30], 0, 'z2input-iteruw');
    imagesc3d2(mask.*z2, N/2, 107, [90,90,-90], [-6,6], 0, 'z2-lineal');
    
    
    
    if isfield(params,'airmodel')
        phi_e  = params.airmodel;
    else
        phi_e = backmodel( outer_mask, 0.0, 2 );
    end

    % LSqrt-Fit background model to unwrapped data, with linear functions
    [ background, kappa ] = fitmodels( z2, mask, weight, phi_e );
    
    %imagesc3d2(background, N/2, 100, [90,90,-90], [-6,6], 0, 'background');
    disp(['Background fitting. Kappa: ', num2str(kappa)])
    z2=phase_unwrapped-background;
    z2 = mask.*z2;
        
else
        z2 = zeros(N,'single');
end
    
s2 = zeros(N,'single');
K2 = abs(kernel).^2;

%tic
for t = 1:num_iter
    % update x : susceptibility estimate
    
    x_prev = x;
    x = mu*(1-mask).*(z-s)./(lambda+mu*(1-mask));

    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if x_update < tol_update
        break
    end
    
    if t < num_iter
        
        
        
            Fz = ( mu2*conj(kernel).*fftn(z2-s2) + mu*fftn((1-mask).*x+s) )./( mu2*K2+mu+eps );
            z = real(ifftn(Fz));
    
        
        
        dt = t/num_iter;
        
        rhs_z2 = mu2*real(ifftn(kernel.*fftn(z))+s2  );
        z2a =  rhs_z2 ./ (mu2) ;

        % Newton-Raphson method
       %if t > 100
        delta = inf;
        inn = 0;
        while (delta > delta_tol && inn < 50)
            inn = inn + 1;
            norm_old = norm(z2a(:));
            
            update = ( weight .* sin(z2a - params.input+background) + mu2*z2a - rhs_z2 ) ./ ( weight .* cos(z2a - params.input+background) + mu2 );            
        
            z2a = z2a - update;     
            delta = norm(update(:)) / norm_old;
        end        
        disp(delta)
       %else
           
           z2b =  (weight.*(phase_unwrapped-background)./(weight+mu2)) + rhs_z2./(weight + mu2);
           z2 = (1-dt)*z2b+dt*z2a;
       %end
        
            s2 = s2 + real(ifftn(kernel.*fftn(z))) - z2;
        
        
        s = s + x - z;  
    end
end
%x=x+kappa(5)*(1-outermask);
out.time = toc;toc

out.x = x;
out.phi = background+real(ifftn(kernel.*fftn((1-mask).*x)));


out.local = params.input;%(params.input - out.phi);
for i = 1:25
    out_old = out.local;
out.local = out.local + 2*pi*round( (out.phi - out.local)/(2*pi) );

if sum(abs(out_old(:)-out.local(:))) < 1
    break;
end

end
out.local = out.local-out.phi;
out.iter = t;



end

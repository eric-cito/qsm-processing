function out = nlPDFp(params)
% Preconditioned Nonlinear Projection onto Dipole Fields (wPDF), for background field removal.
% A Tikhonov term is included, for stability.
% This uses ADMM to solve the functional.
% A rough air-tissue model is fitted to facilitate convergence and accuracy of the algorithm.
%
% Parameters: params - structure with 
% Required fields:
% params.input: total field map
% params.K: dipole kernel in the frequency space
% params.mask: ROI binary mask
% Optional fields:
% params.outermask: binary mask that selects all tissues (not air, and no erosion).
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
    se = strel('ball',2,2);%strel('sphere',2);
    outermask=imdilate(outermask,se);
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


    if isfield(params,'precond')
        precond = params.precond;
    else
        precond = true;
    end
    
background = 0;
if precond
        z2 = params.input;
        
    phi_e = real(ifftn( kernel.*fftn(1-outermask) ));
    kappa = 0;
    phi0 = 0;

    a1 = sum( weight(:).*mask(:).*phi_e(:) );
    a2 = sum( weight(:).*mask(:).*phi_e(:).*phi_e(:) );
    a3 = sum( weight(:).*mask(:) );
    f1 = sum( weight(:).*mask(:).*phi_e(:).*z2(:) );
    f2 = sum( weight(:).*mask(:).*z2(:) );
    det = a1*a1-a2*a3;
    
    phi0 = (a1*f1 - a2*f2)/det;
    kappa = (-a3*f1 + a1*f2)/det;
    
    background = kappa*phi_e+phi0;
    disp(['Background fitting. Kappa: ', num2str(kappa), '   Phi0: ', num2str(phi0)])
    z2 = zeros(N,'single');    

    z2 = params.input;
    z2 =(z2-background);
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
    
        
        
        rhs_z2 = mu2*real(ifftn(kernel.*fftn(z))+s2  );
        z2 =  rhs_z2 ./ (mu2) ;

        % Newton-Raphson method
        delta = inf;
        inn = 0;
        while (delta > delta_tol && inn < 50)
            inn = inn + 1;
            norm_old = norm(z2(:));
            
            update = ( weight .* sin(z2 - params.input+background) + mu2*z2 - rhs_z2 ) ./ ( weight .* cos(z2 - params.input+background) + mu2 );            
        
            z2 = z2 - update;     
            delta = norm(update(:)) / norm_old;
        end        
        disp(delta)
        
            %z2 = Wy + mu2*real(ifftn(kernel.*fftn(z))+s2)./(weight + mu2);
        
            s2 = s2 + real(ifftn(kernel.*fftn(z))) - z2;
        
        
        s = s + x - z;  
    end
end
% Extract output values
x=x+kappa*(1-outermask);
out.time = toc;toc

out.x = x;
out.phi = real(phi0+ifftn(kernel.*fftn((1-mask).*x)));


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

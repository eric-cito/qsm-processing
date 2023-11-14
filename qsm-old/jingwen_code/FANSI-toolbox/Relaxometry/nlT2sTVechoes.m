% Nonlinear T2/T2s mapping with a Total Variation regularization.
% This uses ADMM to solve the functional.
%
% Created by Carlos Milovic in 2018.01.04
%
function out = nlT2sTVechoes(params)
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
    delta_tol = 1e-4;
else
    delta_tol = params.delta_tol;
end

te = params.te;
Ne = length(te);
%W = weight.*weight;
%clear weight;
magn = params.input;


z_dx = zeros([N Ne], 'single');
z_dy = zeros([N Ne], 'single');
z_dz = zeros([N Ne], 'single');

s_dx = zeros([N Ne], 'single');
s_dy = zeros([N Ne], 'single');
s_dz = zeros([N Ne], 'single');

r2 = zeros(N, 'single')+1e-2;
m0 = 1.2*magn(:,:,:,1);%zeros(N, 'single'); %no need to precond

if isfield(params,'precond')
    precond = params.precond;
else
    precond = true;
end

if precond
    z2 =  magn; %better without precond
else
    z2 = zeros([N Ne],'single');
end
s2 = zeros([N Ne],'single');



[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);

E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));

E1t = conj(E1);
E2t = conj(E2);
E3t = conj(E3);

EE2 = E1t .* E1 + E2t .* E2 + E3t .* E3;

%tic
for t = 1:num_iter
    % update r2 
    
        
    r2_prev = r2;
        % Newton-Raphson method for r2
        delta = inf;
        inn = 0;
        while (delta > delta_tol && inn < 20)
            inn = inn + 1;
            norm_old = norm(r2(:));
            
            up = zeros(N,'single');
            down = zeros(N,'single');
            for echo = 1:Ne                
                up = up - m0.*te(echo).*exp(-te(echo)*r2).*(m0.*exp(-te(echo)*r2)-z2(:,:,:,echo)+s2(:,:,:,echo));
                down = down+m0.*te(echo).*te(echo).*exp(-te(echo)*r2).*( 2*m0.*exp(-te(echo)*r2)-z2(:,:,:,echo)+s2(:,:,:,echo) );
            end
        %disp(sum(up(:)))
        %disp(sum(down(:)))
            
            update = up./(down+eps);
            r2 = min(max(r2 - update ,1e-3),2/te(1));
            delta = norm(update(:)) / norm_old;
        end
        disp(delta)
        
    r2_update = 100 * norm(r2(:)-r2_prev(:)) / norm(r2(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(r2_update)])
     
            up = zeros(N,'single');
            down = zeros(N,'single');
            for echo = 1:Ne                
                up = up + exp(-te(echo)*r2).*(z2(:,:,:,echo)-s2(:,:,:,echo));
                down = down+exp(-2*te(echo)*r2);
            end   
    m0 = max(up./(down+eps),0.0);
    
    
    if (r2_update < tol_update) && (t > 10)
        break
    end
    
        
    
    ll = lambda/mu;
    if t < num_iter
            for echo = 1:Ne   

    tx = E1t .* fftn(z_dx(:,:,:,echo) - s_dx(:,:,:,echo));
    ty = E2t .* fftn(z_dy(:,:,:,echo) - s_dy(:,:,:,echo));
    tz = E3t .* fftn(z_dz(:,:,:,echo) - s_dz(:,:,:,echo));
    
    z2(:,:,:,echo) = min(max(real(ifftn( ( mu * (tx + ty + tz) + mu2*fftn(m0.*exp(-te(echo)*r2)+s2(:,:,:,echo)) +fftn(magn(:,:,:,echo)))./( 1.0 + mu2 + mu * EE2 ) )),0.0),1.0);
    
    
        % update z 
        Fx = fftn(z2(:,:,:,echo));
        x_dx = real(ifftn(E1 .* Fx));
        x_dy = real(ifftn(E2 .* Fx));
        x_dz = real(ifftn(E3 .* Fx));
        
        z_dx(:,:,:,echo) = max(abs(x_dx + s_dx(:,:,:,echo)) - regweight(:,:,:,1)*ll, 0) .* sign(x_dx + s_dx(:,:,:,echo));
        z_dy(:,:,:,echo) = max(abs(x_dy + s_dy(:,:,:,echo)) - regweight(:,:,:,2)*ll, 0) .* sign(x_dy + s_dy(:,:,:,echo));
        z_dz(:,:,:,echo) = max(abs(x_dz + s_dz(:,:,:,echo)) - regweight(:,:,:,3)*ll, 0) .* sign(x_dz + s_dz(:,:,:,echo));
        
        % update s : Lagrange multiplier
        s_dx(:,:,:,echo) = s_dx(:,:,:,echo) + x_dx - z_dx(:,:,:,echo);
        s_dy(:,:,:,echo) = s_dy(:,:,:,echo) + x_dy - z_dy(:,:,:,echo);
        s_dz(:,:,:,echo) = s_dz(:,:,:,echo) + x_dz - z_dz(:,:,:,echo);
        
        
        
        s2(:,:,:,echo) = s2(:,:,:,echo) + m0.*exp(-te(echo)*r2) - z2(:,:,:,echo);
        end
    end
end
out.time = toc;toc

out.r2 = r2;
out.t2 = 1.0./(r2+eps);
out.m0 = m0;
out.iter = t;


end

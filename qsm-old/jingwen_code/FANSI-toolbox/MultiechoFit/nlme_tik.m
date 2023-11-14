function out = nlme_tik(params)
% Nonlinear multi-echo fit with Tikhonov regularization. ADMM solver.
%
% Parameters: params - structure with 
% Required fields:
% params.input: all echoes phase container (4th dimension echo index), in radians
% params.TE - echo times, in ms.
% params.lambda: Tikhonov weight
% params.mu: ADMM consistency weight.
% Optional fields:
% params.maxOuterIter: maximum number of ADMM iterations (recommended = 300)
% params.tol_update: convergence limit, change rate in the solution (recommended = 0.1)
% params.weight: magnitude data container (all echoes), rescaled to the [0,1] range.
% params.precond: preconditionate solution (for stability)
% params.b0: Initial guess of the field map
% params.theta0: Initial guess of the phase offset 
%     *** phase(TE) = TE*b0 + theta0 ***
%
% Output: out - structure with the following fields:
% out.b: field map in radians
% out.theta: phase offset in radians
% out.iter: number of iterations needed
% out.time: total elapsed time (including pre-calculations)
% out.lnoise: noise estimation in the phase domain
% out.cnoise: noise estimation in the complex image domain
%
% Modified by Carlos Milovic in 2018.02.05
% Last modified by Carlos Milovic in 2020.07.14

tic

    % Required parameters
    lambda = params.lambda;
    TE = params.TE;
    nE = length(TE);
    sz = size(params.input);
    N = [sz(1) sz(2) sz(3)];

    % Optional parameters
    if isfield(params,'mu')
         mu = params.mu;
    else
        mu = 1.0;
    end

    if isfield(params,'maxOuterIter')
        num_iter = params.maxOuterIter;
    else
        num_iter = 300;
    end
    
    if isfield(params,'tol_update')
       tol_update  = params.tol_update;
    else
       tol_update = 0.1;
    end

    if isfield(params,'weight')
        weight = params.weight;
    else
        weight = ones(size(params.input));
    end
    weight = weight.*weight;
    
    if ~isfield(params,'delta_tol')
        delta_tol = 1e-6;
    else
        delta_tol = params.delta_tol;
    end
    

    % Variable initialization
b = zeros(N, 'single');
theta = zeros(N, 'single');
z = zeros(size(params.input), 'single');
s = zeros(size(params.input), 'single');


if isfield(params,'precond')
    precond = params.precond;
else
    precond = true;
end

if precond
    if isfield(params,'b0')
        b  = params.b0;
    else
        % Estimate the field map using the difference between two unwrapped echoes
        u1 = unwrap(params.input(:,:,:,1), [1 1 1]);
        u2 = unwrap(params.input(:,:,:,2), [1 1 1]);
        b = (u2(:,:,:)-u1(:,:,:))/(TE(2)-TE(1));
    end
    
    
    if isfield(params,'theta0')
        theta  = params.theta0;
    else
        % Estimate the phase offset using the weighted difference between two unwrapped echoes
        u1 = unwrap(params.input(:,:,:,1), [1 1 1]);
        u2 = unwrap(params.input(:,:,:,2), [1 1 1]);
        theta = (TE(1)*u2(:,:,:)-TE(2)*u1(:,:,:))/(TE(1)-TE(2));
    end
    
    
    for t = 1:nE
        model = TE(t)*b+theta;
        z(:,:,:,t) = params.input(:,:,:,t);%TE(t)*b+theta;
        % Initialize this variable with a rough unwrapping of the data.
        for i = 1:5
            z(:,:,:,t) = z(:,:,:,t) + 2*pi*round( (model - z(:,:,:,t))/(2*pi) );
        end
        
    end
        
end

% Precomputation of matrix elements 
A(1) = sum( TE.*TE );
A(2) = sum( TE );
A(3) = A(2);
A(4) = nE+lambda/mu;
dA = A(1)*A(4)-A(2)*A(3);


%tic
for t = 1:num_iter
    % update b and theta
    
    b_prev = b;
    
    F1 = zeros( N );
    F2 = zeros( N );
    for tt = 1:nE
        F1 = F1 + TE(tt)*(z(:,:,:,tt)-s(:,:,:,tt));
        F2 = F2 + (z(:,:,:,tt)-s(:,:,:,tt));
    end
    % b and theta are solved jointly
    b = (A(4)*F1-A(2)*F2)/dA;
    theta = (A(1)*F2-A(3)*F1)/dA;

    b_update = 100 * norm(b(:)-b_prev(:)) / norm(b(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(b_update)])
    
    if b_update < tol_update && t > 25
        break
    end
    
    if t < num_iter
    
        rhs_z = zeros([N nE]);
        for tt = 1:nE
        rhs_z(:,:,:,tt) = mu*(s(:,:,:,tt)+TE(tt)*b+theta);
        end
        z =  rhs_z ./ (mu) ;

        % Newton-Raphson method
       %if t > 100
        delta = inf;
        inn = 0;
        while (delta > delta_tol && inn < 50)
            inn = inn + 1;
            norm_old = norm(z(:));
            
            update = ( weight .* sin(z - params.input) + mu*z - rhs_z ) ./ ( weight .* cos(z - params.input) + mu );            
        
            z = z - update;     
            delta = norm(update(:)) / norm_old;
        end        
        
        for tt = 1:nE
        s(:,:,:,tt) = s(:,:,:,tt) + TE(tt)*b+theta - z(:,:,:,tt);  
        end
        
        
    end
    
end

% Extract output values
out.lnoise = zeros(size(b));
out.cnoise = zeros(size(b));
for i = 1:nE
    y = theta+TE(i)*b;
    out.lnoise = out.lnoise + (y-params.input(:,:,:,i)).^2
    out.cnoise = out.cnoise + weight(:,:,:,i).*(abs(exp(1i*y)-exp(1i*params.input(:,:,:,i)))).^2
end
out.lnoise = sqrt(out.lnoise/nE);
out.cnoise = sqrt(out.cnoise/nE);

out.time = toc;toc
out.iter = t;

out.b = b;
out.theta = theta;

end

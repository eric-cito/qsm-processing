function out = nluwl2(params)
% Nonlinear iterative unwrapping. 
% An initial solution is found by displacing the wrapping artifacts by offsetting the input phase.
% The gradient of several shifted phases are computed and the median value is calculated.
% This is set as a regularization term in an iterative functional, where the gradient of the solution
% should be similar with this mean calculation. The data fidelity term uses a nonlinear model (complex
% image domain), to match the solution to the acquired phase.
%
% Parameters: params - structure with 
% Required fields:
% params.input: wrapped phase data
% params.lambda: regularization weight
% Optional fields:
% params.mu: consistency weight (ADMM weight)
% params.maxOuterIter: maximum number of ADMM iterations (recommended = 50)
% params.tol_update: convergence limit, change rate in the solution (recommended = 1.0)
% params.weight: data fidelity spatially variable weight (recommended = magnitude_data).
% params.precond: preconditionate solution (for stability)
%
% Output: out - structure with the following fields:
% out.psi - unwrapped phase
% out.iter - number of iterations needed
% out.time - total elapsed time (including pre-calculations)
%
% Modified by Carlos Milovic in 2017.11.23
% Last modified by Carlos Milovic in 2020.07.14


tic

lambda = params.lambda;

if isfield(params,'mu')
    mu = params.mu;
else
    mu = 1.0;
end

sz = size(params.input);
N = [sz(1) sz(2) sz(3)];

if isfield(params,'maxOuterIter')
    num_iter = params.maxOuterIter;
else
    num_iter = 150;
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

if isfield(params,'voxel_size')
    voxel_size = params.voxel_size;
else
    voxel_size = [1 1 1];
end

psi = zeros(N, 'single');
z = zeros(N, 'single');
s = zeros(N, 'single');


if isfield(params,'precond')
    precond = params.precond;
else
    precond = true;
end

if precond
    % This founds the image that generates the modified gradient by solving it
    % in a least squared sense. This is equivalent to solve a Poisson equation
    % Here a rapid solver based on the FFT is used.
    [ kuw, luw] = gduw(  params.input, ones(N), voxel_size );
    z = luw;%(kuw+luw)/2;
        
end


% Define the operators
[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);

E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));

E1t = conj(E1);
E2t = conj(E2);
E3t = conj(E3);

EE2 = E1t .* E1 + E2t .* E2 + E3t .* E3;


% Calculate the gradient of the input data
GX(:,:,:,1) = real(ifftn( E1.*fftn(params.input)));
GY(:,:,:,1) = real(ifftn( E2.*fftn(params.input)));
GZ(:,:,:,1) = real(ifftn( E3.*fftn(params.input)));

% Shift the phase data, and recalculate the gradients
temp = params.input + pi/3;
temp = angle(exp(1i*temp));
GX(:,:,:,2) = real(ifftn( E1.*fftn(temp)));
GY(:,:,:,2) = real(ifftn( E2.*fftn(temp)));
GZ(:,:,:,2) = real(ifftn( E3.*fftn(temp)));

% Shift the phase data again, and recalculate the gradients
temp = params.input + 2*pi/3;
temp = angle(exp(1i*temp));
GX(:,:,:,3) = real(ifftn( E1.*fftn(temp)));
GY(:,:,:,3) = real(ifftn( E2.*fftn(temp)));
GZ(:,:,:,3) = real(ifftn( E3.*fftn(temp)));
clear temp;

% Calculate the median value
v(:,:,:,1) = median(GX,4);
v(:,:,:,2) = median(GY,4);
v(:,:,:,3) = median(GZ,4);
% Calculate the divergence of such gradient data.
FLPhi = E1t.*fftn(v(:,:,:,1)) + E2t.*fftn(v(:,:,:,2)) + E3t.*fftn(v(:,:,:,3));
clear v GX GY GZ;
FLPhi = FLPhi./(EE2+mu/lambda); % Apply factors for computational efficiency

for t = 1:num_iter
    
    psi_prev = psi;
    Fpsi = FLPhi + (mu/lambda)*fftn(z-s)./(EE2+mu/lambda);
    psi = real(ifftn(Fpsi));
    
    p_update = 100 * norm(psi(:)-psi_prev(:)) / norm(psi(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(p_update)])
    
    if p_update < tol_update && t > 3
        break
    end
    
    if t < num_iter
    
        rhs_z = mu*(s+psi);
        z = rhs_z/mu;

        % Newton-Raphson method
        delta = inf;
        inn = 0;
        while (delta > delta_tol && inn < 10)
            inn = inn + 1;
            norm_old = norm(z(:));
            
            update = ( weight .* sin(z - params.input) + mu*z - rhs_z ) ./ ( weight .* cos(z - params.input) + mu );            
        
            z = z - update;     
            delta = norm(update(:)) / norm_old;
        end        
        
        s = s +psi - z;  
        
    end
    
end
% Extract output values
out.time = toc;toc
out.psi = psi;
out.iter = t;

end

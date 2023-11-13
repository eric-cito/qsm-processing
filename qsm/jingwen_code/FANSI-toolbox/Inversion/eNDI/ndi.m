function out = ndi(params)
% Nonlinear Dipole Inversion. Gradient Descent solver.
% Based on Polak D, et al. NMR Biomed 2020.
%
% Parameters: params - structure with 
% Required fields:
% params.input: local field map, in radians
% params.K: dipole kernel in the frequency space
% Optional fields:
% params.alpha: regularization weight (use small values for stability)
% params.maxOuterIter: maximum number of iterations (recommended = 100)
% params.weight: data fidelity spatially variable weight (recommended = magnitude_data).
% params.tau: gradient descent rate
% params.precond: preconditionate solution (for stability)
% params.show_iters: show intermediate results for each iteration.
%
% Output: out - structure with the following fields:
% out.x: calculated susceptibility map, in radians
% out.time: total elapsed time (including pre-calculations)
%
% Last modified by Carlos Milovic in 2020.07.14


tic

    % Required parameters
    kernel = params.K;

    if isfield(params,'alpha')
         alpha = params.alpha;
    else
        alpha = 1E-5;
    end
    
    if isfield(params,'tau')
         tau = params.tau;
    else
        tau = 2.0; % Accelerate it slightly. Too large values may cause fast divergence.
    end
    
    N = size(params.input);

    if isfield(params,'maxOuterIter')
        num_iter = params.maxOuterIter;
    else
        num_iter = 100;
    end
    
    if isfield(params,'weight')
        weight = params.weight;
    else
        weight = ones(N);
    end
    weight = weight.*weight;
    

    if isfield(params,'precond')
        precond = params.precond;
    else
        precond = true;
    end
    
    if precond
        x = params.weight.*params.input;
    else
        x = zeros(N, 'single');
    end



%tic
for t = 1:num_iter
    % update x : susceptibility estimate
    x_prev = x;
    phix = susc2field(kernel,x);
    x = x_prev - tau*susc2field( conj(kernel),weight.*sin(phix-params.input ) ) - tau*alpha*x;
    
    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if show_iters
        imagesc3d2(x, N/2, 31, [90,90,-90], [-1,1], [], ['NDI - Iteration: ', num2str(t), '   Update: ', num2str(x_update)] )
        drawnow;
    end

end
out.time = toc;toc
out.x = x;
out.iter = t;

end

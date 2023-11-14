function out = ndiCG(params)
% Nonlinear Dipole Inversion. Conjugate Gradient Descent solver.
% This is a faster alternative to ndi().
%
% Parameters: params - structure with 
% Required fields:
% params.input: local field map, in radians
% params.K: dipole kernel in the frequency space
% Optional fields:
% params.alpha: regularization weight (use small values for stability)
% params.maxOuterIter: maximum number of iterations (recommended = 100)
% params.weight: data fidelity spatially variable weight (recommended = magnitude_data).
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
        alpha = 1E-6;
    end
    
    N = size(params.input);

    if isfield(params,'maxOuterIter')
        num_iter = params.maxOuterIter;
    else
        num_iter = 50;
    end
    
    if isfield(params,'weight')
        weight = params.weight;
    else
        weight = ones(N);
    end
    weight = weight.*weight;
    
    if isfield(params,'show_iters')
        show_iters = params.show_iters;
    else
        show_iters = false;
    end
    
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


% Initiate with a gradient descent, with automatic step by line search.
    phix = susc2field(kernel,x);
    dx = -susc2field( conj(kernel),weight.*sin(phix-params.input ) ) - alpha*x;


    B = dx.*dx;
    A = dx.*(susc2field(conj(kernel),weight.*sin(susc2field(kernel,dx))) +alpha*dx);
    tau = sum(B(:))/(sum(A(:))+eps);

    x_prev = x;
    x = x_prev + tau*dx;
    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(0), '   Update: ', num2str(x_update)])

    s = dx;
%tic
% Start CG steps
for t = 1:num_iter
    % update x : susceptibility estimate
    x_prev = x;
    phix = susc2field(kernel,x);
    dx_prev = dx;
    dx = -susc2field( conj(kernel),weight.*sin(phix-params.input ) ) - alpha*x;
    
    betaPR = max(sum(dx(:).*(dx(:)-dx_prev(:)))/(sum(dx(:).*dx(:))+eps),0); % Automatic reset
    
    s = dx + betaPR*s;
    
    
    B = s.*dx;
    A = s.*(susc2field(conj(kernel),weight.*sin(susc2field(kernel,s))) +alpha*s);
    tau = sum(B(:))/(sum(A(:))+eps);
    
    x = x_prev + tau*s;
    
    x_update = 100 * norm(x(:)-x_prev(:)) / norm(x(:));
    disp(['Iter: ', num2str(t), '   Update: ', num2str(x_update)])
    
    if show_iters
        imagesc3d2(x, N/2, 31, [90,90,-90], [-1,1], [], ['NDICG - Iteration: ', num2str(t), '   Update: ', num2str(x_update)] )
        drawnow;
    end
end
out.time = toc;toc
out.x = x;
out.iter = t;

end

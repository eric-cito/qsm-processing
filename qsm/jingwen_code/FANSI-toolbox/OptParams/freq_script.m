% Example script of the creation and comparison of the frequency masks.

params = []; % Define the required parameters
params.input = phase_use;
params.K = kernel;


% Create the masks. Fine tune internal parameters if needed.
[m1, m2, m3 ] = create_freqmasks( spatial_res, kernel );

% Explore a range or regularization weights, and calculate the reconstructions
% with the respective mean mask amplitudes
alpha_center = 10^(-(27+10)/10);
ma = -2:0.1:2;
for k = 1:length(ma)
    
alpha(k) = (10^ma(k))*alpha_center;
params.alpha1 = alpha(k); % Update the regularization weight

outt = wTV(params); 

[e1(k), e2(k), e3(k)] = compute_freqe(outt.x.*mask_use,m1,m2,m3);

end

% Draw the mask amplitudes and zeta functions,
% and obtain the optimal regularization weight
[opt23index, alpha_opt] = draw_freque(alpha,e1,e2,e3,11);

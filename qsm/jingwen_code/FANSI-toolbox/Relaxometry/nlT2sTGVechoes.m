% Nonlinear T2/T2s mapping with a Total Generalized Variation regularization.
% This uses ADMM to solve the functional.
%
% Created by Carlos Milovic in 2019.01.23
%
function out = nlT2sTGVechoes(params)
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

mu1 = params.mu1;
    if isfield(params,'mu0')
        mu0 = params.mu0;
    else
        mu0 = 2*mu1;
    end
    
    alpha1 = params.alpha1;
    if isfield(params,'alpha0')
        alpha0 = params.alpha0;
    else
        alpha0 = 2*alpha1;
    end


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


%z_dx = zeros([N Ne], 'single');
%z_dy = zeros([N Ne], 'single');
%z_dz = zeros([N Ne], 'single');

%s_dx = zeros([N Ne], 'single');
%s_dy = zeros([N Ne], 'single');
%s_dz = zeros([N Ne], 'single');

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


    
    % Precompute gradient-related matrices
    [k1, k2, k3] = ndgrid(0:N(1)-1, 0:N(2)-1, 0:N(3)-1);
    E1 = 1 - exp(2i .* pi .* k1 / N(1));
    E2 = 1 - exp(2i .* pi .* k2 / N(2));
    E3 = 1 - exp(2i .* pi .* k3 / N(3)); 

    Et1 = conj(E1);     Et2 = conj(E2);     Et3 = conj(E3);     Kt = conj(K);
    
    E1tE1 = Et1.*E1;    E2tE2 = Et2.*E2;    E3tE3 = Et3.*E3;
    mu0_over_2_E1tE2 = mu0/2*Et1.*E2;
    mu0_over_2_E1tE3 = mu0/2*Et1.*E3;
    mu0_over_2_E2tE3 = mu0/2*Et2.*E3;

    a0 = mu2*Kt.*K;
    a0_mu1_E_sos    = a0 + mu1*(E1tE1 + E2tE2 + E3tE3);
    mu1I_mu0_E_wsos1 = mu1 + mu0*(E1tE1 + (E2tE2 + E3tE3)/2);
    mu1I_mu0_E_wsos2 = mu1 + mu0*(E1tE1/2 + E2tE2 + E3tE3/2);
    mu1I_mu0_E_wsos3 = mu1 + mu0*((E1tE1 + E2tE2)/2 + E3tE3);
    
    %% Precomputation for Cramer's Rule 
    
    a1 = a0_mu1_E_sos; a2 = mu1I_mu0_E_wsos1;  a3 = mu1I_mu0_E_wsos2;  a4 = mu1I_mu0_E_wsos3;
    a5 = -mu1*E1;       a6 = -mu1*E2;           a7 = mu0_over_2_E1tE2;  a8 = -mu1*E3;   a9 = mu0_over_2_E1tE3;  a10 = mu0_over_2_E2tE3;
    a5t = conj(a5);     a6t = conj(a6);         a7t = conj(a7);         a8t = conj(a8); a9t = conj(a9);         a10t = conj(a10);    
    
    % For x
    D11 = a2.*a3.*a4    + a7t.*a9.*a10t + a7.*a9t.*a10  - a3.*a9.*a9t   - a2.*a10.*a10t     - a4.*a7.*a7t;
    D21 = a3.*a4.*a5t   + a6t.*a9.*a10t + a7.*a8t.*a10  - a3.*a8t.*a9   - a5t.*a10.*a10t    - a4.*a6t.*a7;
    D31 = a4.*a5t.*a7t  + a6t.*a9.*a9t  + a2.*a8t.*a10  - a7t.*a8t.*a9  - a5t.*a9t.*a10     - a2.*a4.*a6t;
    D41 = a5t.*a7t.*a10t + a6t.*a7.*a9t + a2.*a3.*a8t   - a7.*a7t.*a8t  - a3.*a5t.*a9t      - a2.*a6t.*a10t;

    % For vx
    D12 = a3.*a4.*a5    + a7t.*a8.*a10t + a6.*a9t.*a10  - a3.*a8.*a9t   - a5.*a10.*a10t - a4.*a6.*a7t;
    D22 = a1.*a3.*a4    + a6t.*a8.*a10t + a6.*a8t.*a10  - a3.*a8.*a8t   - a1.*a10.*a10t - a4.*a6.*a6t;
    D32 = a1.*a4.*a7t   + a6t.*a8.*a9t  + a5.*a8t.*a10  - a7t.*a8.*a8t  - a1.*a9t.*a10  - a4.*a5.*a6t;
    D42 = a1.*a7t.*a10t + a6.*a6t.*a9t  + a3.*a5.*a8t   - a6.*a7t.*a8t  - a1.*a3.*a9t   - a5.*a6t.*a10t;

    % For vy
    D13 = a4.*a5.*a7 + a2.*a8.*a10t + a6.*a9.*a9t - a7.*a8.*a9t - a5.*a9.*a10t - a2.*a4.*a6;
    D23 = a1.*a4.*a7 + a5t.*a8.*a10t +a6.*a8t.*a9 - a7.*a8.*a8t - a1.*a9.*a10t - a4.*a5t.*a6;
    D33 = a1.*a2.*a4 + a5t.*a8.*a9t + a5.*a8t.*a9 - a2.*a8.*a8t - a1.*a9.*a9t - a4.*a5.*a5t;
    D43 = a1.*a2.*a10t + a5t.*a6.*a9t + a5.*a7.*a8t - a2.*a6.*a8t - a1.*a7.*a9t - a5.*a5t.*a10t;

    % For vz
    D14 = a5.*a7.*a10 + a2.*a3.*a8 + a6.*a7t.*a9 - a7.*a7t.*a8 - a3.*a5.*a9 -a2.*a6.*a10;
    D24 = a1.*a7.*a10 + a3.*a5t.*a8 + a6.*a6t.*a9 - a6t.*a7.*a8 - a1.*a3.*a9 - a5t.*a6.*a10;
    D34 = a1.*a2.*a10 + a5t.*a7t.*a8 + a5.*a6t.*a9 - a2.*a6t.*a8 - a1.*a7t.*a9 - a5.*a5t.*a10;
    D44 = a1.*a2.*a3 + a5t.*a6.*a7t + a5.*a6t.*a7 - a2.*a6.*a6t - a1.*a7.*a7t - a3.*a5.*a5t;

    det_A = a1.*D11 - a5.*D21 + a6.*D31 - a8.*D41;
    det_Ainv = 1 ./ (eps+det_A);

    % Allocate memory for first order gradient
    s1_1 = zeros([N Ne],'single'); z1_1 = zeros([N Ne],'single');
    s1_2 = zeros([N Ne],'single'); z1_2 = zeros([N Ne],'single');
    s1_3 = zeros([N Ne],'single'); z1_3 = zeros([N Ne],'single');
    
    % Allocate memory for symmetrized gradient
    s0_1 = zeros([N Ne],'single'); z0_1 = zeros([N Ne],'single'); 
    s0_2 = zeros([N Ne],'single'); z0_2 = zeros([N Ne],'single');
    s0_3 = zeros([N Ne],'single'); z0_3 = zeros([N Ne],'single'); 
    s0_4 = zeros([N Ne],'single'); z0_4 = zeros([N Ne],'single');
    s0_5 = zeros([N Ne],'single'); z0_5 = zeros([N Ne],'single');
    s0_6 = zeros([N Ne],'single'); z0_6 = zeros([N Ne],'single');

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
    
        
    
    if t < num_iter
            for echo = 1:Ne   

    %tx = E1t .* fftn(z_dx(:,:,:,echo) - s_dx(:,:,:,echo));
    %ty = E2t .* fftn(z_dy(:,:,:,echo) - s_dy(:,:,:,echo));
    %tz = E3t .* fftn(z_dz(:,:,:,echo) - s_dz(:,:,:,echo));
    
    
        F_z0_minus_s0_1 = fftn(z0_1(:,:,:,echo) - s0_1(:,:,:,echo));
        F_z0_minus_s0_2 = fftn(z0_2(:,:,:,echo) - s0_2(:,:,:,echo));
        F_z0_minus_s0_3 = fftn(z0_3(:,:,:,echo) - s0_3(:,:,:,echo));
        F_z0_minus_s0_4 = fftn(z0_4(:,:,:,echo) - s0_4(:,:,:,echo));
        F_z0_minus_s0_5 = fftn(z0_5(:,:,:,echo) - s0_5(:,:,:,echo));
        F_z0_minus_s0_6 = fftn(z0_6(:,:,:,echo) - s0_6(:,:,:,echo));
        

        F_z1_minus_s1_1 = fftn(z1_1(:,:,:,echo) - s1_1(:,:,:,echo));
        F_z1_minus_s1_2 = fftn(z1_2(:,:,:,echo) - s1_2(:,:,:,echo));
        F_z1_minus_s1_3 = fftn(z1_3(:,:,:,echo) - s1_3(:,:,:,echo));
        
        rhs0 = mu2*(fftn( z2(:,:,:,echo)-s2(:,:,:,echo) ));
        rhs1    =  rhs0                + mu1*(Et1.*F_z1_minus_s1_1 + Et2.*F_z1_minus_s1_2 + Et3.*F_z1_minus_s1_3);
        rhs2    = -mu1*F_z1_minus_s1_1 + mu0*(Et1.*F_z0_minus_s0_1 + Et2.*F_z0_minus_s0_4 + Et3.*F_z0_minus_s0_5);
        rhs3    = -mu1*F_z1_minus_s1_2 + mu0*(Et2.*F_z0_minus_s0_2 + Et1.*F_z0_minus_s0_4 + Et3.*F_z0_minus_s0_6);
        rhs4    = -mu1*F_z1_minus_s1_3 + mu0*(Et3.*F_z0_minus_s0_3 + Et1.*F_z0_minus_s0_5 + Et2.*F_z0_minus_s0_6);
        
        % Cramer's rule
        Fx = (rhs1.*D11 - rhs2.*D21 + rhs3.*D31 - rhs4.*D41) .* det_Ainv;
        Fv1 = (-rhs1.*D12 + rhs2.*D22 - rhs3.*D32 + rhs4.*D42) .* det_Ainv;
        Fv2 = (rhs1.*D13 - rhs2.*D23 + rhs3.*D33 - rhs4.*D43) .* det_Ainv;
        Fv3 = (-rhs1.*D14 +rhs2.*D24 - rhs3.*D34 + rhs4.*D44) .* det_Ainv;   
    
    
    
        z2(:,:,:,echo) = min(max(real( ifftn( Fx ) ),0.0),1.0);
        v1 = real(ifftn(Fv1));
        v2 = real(ifftn(Fv2));
        v3 = real(ifftn(Fv3));
    
    
        % update z 
        Fx = fftn(z2(:,:,:,echo));
        
        
        Dx1 = real(ifftn(E1.*Fx));
        Dx2 = real(ifftn(E2.*Fx));
        Dx3 = real(ifftn(E3.*Fx));

        E_v1 = real(ifftn(E1.*Fv1));
        E_v2 = real(ifftn(E2.*Fv2));
        E_v3 = real(ifftn(E3.*Fv3));
        E_v4 = real(ifftn(E1.*Fv2 + E2.*Fv1))/2;
        E_v5 = real(ifftn(E1.*Fv3 + E3.*Fv1))/2;
        E_v6 = real(ifftn(E2.*Fv3 + E3.*Fv2))/2;
        
        % Update z0: Symm grad
        z0_1(:,:,:,echo) = max(abs(E_v1 + s0_1(:,:,:,echo))-alpha0/mu0,0).*sign(E_v1 + s0_1(:,:,:,echo));
        z0_2(:,:,:,echo) = max(abs(E_v2 + s0_2(:,:,:,echo))-alpha0/mu0,0).*sign(E_v2 + s0_2(:,:,:,echo));
        z0_3(:,:,:,echo) = max(abs(E_v3 + s0_3(:,:,:,echo))-alpha0/mu0,0).*sign(E_v3 + s0_3(:,:,:,echo));
        z0_4(:,:,:,echo) = max(abs(E_v4 + s0_4(:,:,:,echo))-alpha0/mu0,0).*sign(E_v4 + s0_4(:,:,:,echo));
        z0_5(:,:,:,echo) = max(abs(E_v5 + s0_5(:,:,:,echo))-alpha0/mu0,0).*sign(E_v5 + s0_5(:,:,:,echo));
        z0_6(:,:,:,echo) = max(abs(E_v6 + s0_6(:,:,:,echo))-alpha0/mu0,0).*sign(E_v6 + s0_6(:,:,:,echo));
        
        % Update z1: Grad
        z1_1(:,:,:,echo) = max(abs(Dx1-v1+s1_1(:,:,:,echo))-regweight(:,:,:,1)*alpha1/mu1,0).*sign(Dx1-v1+s1_1(:,:,:,echo));
        z1_2(:,:,:,echo) = max(abs(Dx2-v2+s1_2(:,:,:,echo))-regweight(:,:,:,2)*alpha1/mu1,0).*sign(Dx2-v2+s1_2(:,:,:,echo));
        z1_3(:,:,:,echo) = max(abs(Dx3-v3+s1_3(:,:,:,echo))-regweight(:,:,:,3)*alpha1/mu1,0).*sign(Dx3-v3+s1_3(:,:,:,echo));
        
        
        
        % update s : Lagrange multiplier  
        s0_1(:,:,:,echo) = s0_1(:,:,:,echo) + E_v1-z0_1(:,:,:,echo);
        s0_2(:,:,:,echo) = s0_2(:,:,:,echo) + E_v2-z0_2(:,:,:,echo);
        s0_3(:,:,:,echo) = s0_3(:,:,:,echo) + E_v3-z0_3(:,:,:,echo);
        s0_4(:,:,:,echo) = s0_4(:,:,:,echo) + E_v4-z0_4(:,:,:,echo);
        s0_5(:,:,:,echo) = s0_5(:,:,:,echo) + E_v5-z0_5(:,:,:,echo);
        s0_6(:,:,:,echo) = s0_6(:,:,:,echo) + E_v6-z0_6(:,:,:,echo);
        
        s1_1(:,:,:,echo) = s1_1(:,:,:,echo) + Dx1-v1-z1_1(:,:,:,echo);
        s1_2(:,:,:,echo) = s1_2(:,:,:,echo) + Dx2-v2-z1_2(:,:,:,echo);
        s1_3(:,:,:,echo) = s1_3(:,:,:,echo) + Dx3-v3-z1_3(:,:,:,echo);
        
        
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

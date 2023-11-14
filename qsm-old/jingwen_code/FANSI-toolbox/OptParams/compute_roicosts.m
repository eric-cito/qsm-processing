function [ data_cost, reg_cost ] = compute_roicosts( chi, phase, kernel, roi )
% Calculates the costs associated with a TV regularization
% Useful to build a L-curve plot, with costs calculated in a region of
% interest to avoid unstabilities due to artifacts or noise.
%%
% Parameters:
% chi: calculated susceptibility map
% phase: acquired phase or local phase map.
% kernel: dipole kernel in the frequency space
% roi: binary mask selecting a Region Of Interest (ROI).
%
% Output:
% data_cost: Cost of the L2 (linear) data fidelity term.
% reg_cost: Cost of the TV functional.
%
% Last modified by Carlos Milovic in 2017.03.30
%
    delta_phi = phase-real(ifftn(kernel.*(fftn(chi))));
    data_cost = sum(delta_phi(roi>0).*delta_phi(roi>0));

    N = size(chi);
    [k1, k2, k3] = ndgrid(0:N(1)-1, 0:N(2)-1, 0:N(3)-1);
    E1 = 1 - exp(2i .* pi .* k1 / N(1));
    E2 = 1 - exp(2i .* pi .* k2 / N(2));
    E3 = 1 - exp(2i .* pi .* k3 / N(3)); 
    
    G1 = real(ifftn(E1.*fftn(chi)));
    G2 = real(ifftn(E2.*fftn(chi)));
    G3 = real(ifftn(E3.*fftn(chi)));
    GX = abs(G1) + abs(G2) + abs(G3);
    
    reg_cost = sum(GX(roi>0));

end


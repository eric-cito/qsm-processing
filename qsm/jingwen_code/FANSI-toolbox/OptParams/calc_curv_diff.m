function [ Kappa ] = calc_curv_diff( Lambda, regularization, consistency, isLinear )
% Calculate the curvature of the L-curve by using a finite differences approach
% Von Neumann boundary conditions are imposed.
% This function outputs the curvature Kappa.
%
% Parameters:
% Lambda: vector with regularization weight values
% regularization: vector with regularization costs (use the appropiate function to evaluate such costs)
% consistency: vector with data fidelity costs (use the appropiate function to evaluate such costs)
% isLinear: true for a calculation in the linear domain and false to use the log domain.
%
% Last modified by Carlos Milovic, 08.07.2020 



% finite differences differentiation to find Kappa (curvature) 

if isLinear == true
eta = regularization; 
rho = consistency;
else
eta = log(regularization.^2); 
rho = log(consistency.^2);
end

eta_del = ([eta(2:end) eta(end)] - [eta(1) eta(1:(end-1))]);
rho_del = ([rho(2:end) rho(end)] - [rho(1) rho(1:(end-1))]);

eta_del2 = ([eta(3:end) eta(end) eta(end)] + [eta(1) eta(1) eta(1:(end-2))])-2*eta;
rho_del2 = ([rho(3:end) rho(end) rho(end)] + [rho(1) rho(1) rho(1:(end-2))])-2*rho;


Kappa = 2 * (rho_del2 .* eta_del - eta_del2 .* rho_del) ./ (rho_del.^2 + eta_del.^2+eps).^1.5;

end


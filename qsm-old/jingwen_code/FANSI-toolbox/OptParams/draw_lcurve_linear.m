function [ Kappa ] = draw_lcurve_linear( Lambda, regularization, consistency, fig )
% Draw the L-curve and calculate the curvature in the linear domain.
% This function outputs the curvature Kappa.
%
% Parameters:
% Lambda: vector with regularization weight values
% regularization: vector with regularization costs (use the appropiate function to evaluate such costs)
% consistency: vector with data fidelity costs (use the appropiate function to evaluate such costs)
% fig: figure number to display the L-curve and the curvature.
%
% Last modified by Carlos Milovic, 08.07.2020 


% First plot the L-curve in the linear domain
figure(fig), subplot(1,2,1), plot((consistency), (regularization), 'marker', '*')

% Calculate the curvature in the linear domain, using splines
[ Kappa ] = calc_curv_splines( Lambda, regularization, consistency, true )

index_opt = find(Kappa == max(Kappa));
disp(['Maximum lambda, consistency, regularization: ', num2str([Lambda(index_opt), consistency(index_opt), regularization(index_opt)])])

% Plot the curvature as function of the regularization weight
figure(fig), subplot(1,2,2), semilogx(Lambda, Kappa, 'marker', '*')
end


function [ Kappa, etaP, rhoP ] = calc_curv_parametric( Lambda, regularization, consistency )
% Calculate the curvature of the L-curve by using a parametric approach
% Cost functions are approximated by a cuadratic polynomial as a function of the 
% logarithm of Lambda.
% This may be used if few points were calculated, as an approximation.
% This function outputs the curvature Kappa and the parametric functions.
%
% Parameters:
% Lambda: vector with regularization weight values
% regularization: vector with regularization costs (use the appropiate function to evaluate such costs)
% consistency: vector with data fidelity costs (use the appropiate function to evaluate such costs)
%
% Last modified by Carlos Milovic, 08.07.2020 

eta = (regularization);
rho = (consistency);

A = zeros([3 3]);
A(1,1) = length(Lambda);
A(1,2) = sum( log(Lambda(:)) );
A(2,1) = A(1,2);
A(2,2) = sum( log(Lambda(:)).^2 );
A(3,1) = A(2,2);
A(1,3) = A(2,2);
A(3,2) = sum( log(Lambda(:)).^3 );
A(2,3) = A(3,2);
A(3,3) = sum( log(Lambda(:)).^4 );

AtAiAt = inv(A'*A)*A';

be(1) = sum( eta );
be(2) = sum( eta.*log(Lambda) );
be(3) = sum( eta.*log(Lambda).*log(Lambda) );

br(1) = sum( rho );
br(2) = sum( rho.*log(Lambda) );
br(3) = sum( rho.*log(Lambda).*log(Lambda) );

Pe = AtAiAt*be';
Pr = AtAiAt*br';


etaP = Pe(1)+Pe(2)*log(Lambda)+Pe(3)*log(Lambda).*log(Lambda);
rhoP = Pr(1)+Pr(2)*log(Lambda)+Pr(3)*log(Lambda).*log(Lambda);

eta_del = 2*Pe(3)*log(Lambda)+Pe(2);
rho_del = 2*Pr(3)*log(Lambda)+Pr(2);

eta_del2 = 2*Pe(3);
rho_del2 = 2*Pr(3);


Kappa = 2 * (rho_del2 .* eta_del - eta_del2 .* rho_del) ./ (rho_del.^2 + eta_del.^2+eps).^1.5;

end


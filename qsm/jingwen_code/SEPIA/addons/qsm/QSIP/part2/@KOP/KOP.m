function res = KOP(k, b0z, chi0, delta)
% 
% kernel - is the second z derivative of 1/r kernel which is the same for all segmentations,
%     but differs in overall size. extract the kernel by running
%     ~/matlab/b0calc with the -k option set
%
% b0z - the value of main b0 field (ie 3 Tesla) 
%
% delta - is the perturbation parameter in the equation: chi = chi0 + delta*chi1
%
% chi0 - the constant in the equation: chi = chi0 + delta*chi1
%
%-------------------------------------------------------------------
% If chi1 labels exist on [0,1] set the following:
%
% chi0 = 0.4e-6
% delta = -9.45e-6
%-------------------------------------------------------------------
  

res.kernel = k;
res.b0z = b0z;
res.chi0 = chi0;
res.delta = delta;


res = class(res,'KOP');  %creates an object of class KOP from structure res


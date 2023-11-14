function [ out ] = unwrap( phase, voxel_size )
% Unwrap phase using Laplacian operation (integer solution)
% Schofield and Zhu, 2003 M.A. Schofield and Y. Zhu, Fast phase unwrapping
% algorithm for interferometric applications, Opt. Lett. 28 (2003), pp. 1194?196
%
% This function uses the approximate solution to find 2pi integers additive factors 
% that better approximate the unwrapped data.
%
%   Last modified by Carlos Milovic on 2020.07.14

% Approximate solution
puw = unwrapLaplacian(phase, size(phase) ,voxel_size);

out = phase;
% Iterate to find the propper integer factors
for i = 1:150
    out_old = out;
out = out + 2*pi*round( (puw - out)/(2*pi) );
% Stop if the residual is low enough
if sum(abs(out_old(:)-out(:))) < 1
    break;
end

end

end


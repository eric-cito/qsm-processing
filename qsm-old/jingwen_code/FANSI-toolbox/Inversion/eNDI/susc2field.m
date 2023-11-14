function [phi] = susc2field(D,X)
%Susceptibility to Field calculation

phi = real(ifftn( D.* fftn( X)));
end


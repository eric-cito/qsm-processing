function [ res ] = ifft3c( x )

%res = sqrt(length(x(:)))*ifftshift(ifftn(fftshift(x)));

res = sqrt(length(x(:)))*fftshift(ifftn(ifftshift(x)));

end


function [ res ] = fft3c( x )

res = 1/sqrt(length(x(:)))*fftshift(fftn(ifftshift(x)));

end


function [x, residual] = CG_dipole2(b, x, D, msk, num_iter)
    % note: this is the same as CG_dipole.m except it returns the residual,
    %       rsold
    % b : field map
    % x : initial guess
    % kernel : susceptibility kernel
    % msk : brain ROI mask
    % num_iter : number of CG iterations
    
    residual = zeros(num_iter,1);
    n = size(b);
    msk_out = 1-msk;
    
    Ah_b = msk_out .* ifft3c(conj(D).*fft3c(msk.*b));
    
    if length(find(x(:)~=0)) > 0 
        A_x = msk .* ifft3c(D .* fft3c(msk_out.*x));         
        A_x = msk_out .* ifft3c(conj(D) .* fft3c(A_x));
    else
        A_x = 0;
    end
    
    r = Ah_b - A_x;
    p = r;
    rsold = dot(r(:),r(:));
 
    for t = 1:num_iter
        disp(['Dipole fitting iter: ', num2str(t), ' / ', num2str(num_iter), '  Residual: ', num2str(rsold)])
        
        residual(t) = rsold;
        
        %Ap = A*p;
        A_p = msk .* ifft3c(D .* fft3c(msk_out.*p));         
        A_p = msk_out .* ifft3c(conj(D) .* fft3c(A_p));
        
        
        alpha = rsold / dot(p(:),A_p(:));
        x = x + alpha * p;
        r = r - alpha * A_p;
        rsnew = dot(r(:),r(:));
        
        if sqrt(rsnew) < 1e-10
              break;
        end
        
        p = r + rsnew/rsold*p;
        rsold = rsnew;
        
    end
    
end

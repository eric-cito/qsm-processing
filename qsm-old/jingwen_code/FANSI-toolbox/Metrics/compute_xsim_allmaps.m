function [mssim, ssim_map, lmap, cmap, smap] = compute_xsim_allmaps(img1, img2, mask, sw)
% This function reimplements the SSIM metric for QSM comparisons.
% Please use this function instead of "compute_ssim".
% This function returns all cost functions (luminance, contrast and structural)
% for (visual) evaluation.
%
% Changes made to the 2016 QSM Reconstruction Challenge implementation:
% - QSM maps are no longer rescaled. Keeping the reference and target in the 
% native [ppm] is recommended.
% - L parameter is set to 1 [ppm] accordingly
% - K = [0.01 0.001] to promote detection of streaking artifacts.
% - mask parameter allows to set a custom ROI for evaluation.
% - sw parameter sets the size of the gaussian kernel (default=[3 3 3])
%
% See README.txt for more details.
%
% Last modified by Carlos Milovic in 2019.05.26
% See below part of the original header of this function:
%========================================================================
%SSIM Index, Version 1.0
%Copyright(c) 2003 Zhou Wang
%All Rights Reserved.
%
%The author was with Howard Hughes Medical Institute, and Laboratory
%for Computational Vision at Center for Neural Science and Courant
%Institute of Mathematical Sciences, New York University, USA. He is
%currently with Department of Electrical and Computer Engineering,
%University of Waterloo, Canada.
%
%----------------------------------------------------------------------
%Permission to use, copy, or modify this software and its documentation
%for educational and research purposes only and without fee is hereby
%granted, provided that this copyright notice and the original authors'
%names appear on all copies and supporting documentation. This program
%shall not be used, rewritten, or adapted as the basis of a commercial
%software or hardware product without first obtaining permission of the
%authors. The authors make no representations about the suitability of
%this software for any purpose. It is provided "as is" without express
%or implied warranty.
%----------------------------------------------------------------------
%
%This is an implementation of the algorithm for calculating the
%Structural SIMilarity (SSIM) index between two images. Please refer
%to the following paper:
%
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%quality assessment: From error measurement to structural similarity"
%IEEE Transactios on Image Processing, vol. 13, no. 4, Apr. 2004.
%
%Kindly report any suggestions or corrections to zhouwang@ieee.org
%
%----------------------------------------------------------------------
%
%Output: (1) mssim: the mean SSIM index value between 2 images.
%            If one of the images being compared is regarded as 
%            perfect quality, then mssim can be considered as the
%            quality measure of the other image.
%            If img1 = img2, then mssim = 1.
%        (2) ssim_map: the SSIM index map of the test image. The map
%            has a smaller size than the input images. The actual size:
%            size(img1) - size(window) + 1.
%
%========================================================================


if (nargin < 2 | nargin > 4)
   mssim = -Inf;
   ssim_map = -Inf;
   return;
end

if (size(img1) ~= size(img2))
   mssim = -Inf;
   ssim_map = -Inf;
   return;
end

s = size(img1);



%--------------------------------------------------------------------------




if (nargin == 2)
    
   sw = [3 3 3];

   if ((s(1) < sw(1)) | (s(2) < sw(2)) | (s(3) < sw(3)))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
    ind = find(img2 ~=0);
    mask = zeros(size(img1));
    mask(ind) = 1;
%    window = fspecial('gaussian', 11, 1.5);	%
   window = gkernel(1.5,sw);	%
   K(1) = 0.01;								      % default settings
   K(2) = 0.001;								      %
   L = 1;                                  %
end

if (nargin == 3)
   sw = [3 3 3];

   if ((s(1) < sw(1)) | (s(2) < sw(2)) | (s(3) < sw(3)))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
%    ind = find(img1 ~=0);
%    window = fspecial('gaussian', 11, 1.5);	%
   window = gkernel(1.5,sw);	%
   K(1) = 0.01;								      % default settings
   K(2) = 0.001;								      %
   L = 1;                                  %
end

if (nargin == 4)
   if ((s(1) < sw(1)) | (s(2) < sw(2)) | (s(3) < sw(3)))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
%    window = fspecial('gaussian', 11, 1.5);
     window = gkernel(1.5,sw);	%
%     ind = find(img1 ~=0);
   K(1) = 0.01;								      % default settings
   K(2) = 0.001;	  
   L = 1;
   if (length(K) == 2)
      if (K(1) < 0 | K(2) < 0)
		   mssim = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   mssim = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end



C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;

window = window/sum(window(:));

img1 = double(img1);
img2 = double(img2);

mu1   = convn( img1,window, 'same');
mu2   = convn( img2,window, 'same');

mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;

sigma1_sq = convn( img1.*img1,window, 'same') - mu1_sq;
sigma2_sq = convn( img2.*img2,window, 'same') - mu2_sq;
sigma12 = convn( img1.*img2,window, 'same') - mu1_mu2;


if (C1 > 0 && C2 > 0)
    lmap = (2*mu1_mu2 + C1)./(mu1_sq + mu2_sq + C1);
    cmap = (2*sqrt(sigma1_sq).*sqrt(sigma2_sq) + C2)./(sigma1_sq + sigma2_sq + C2);
    smap = (2*sigma12 + C2)./(2*sqrt(sigma1_sq).*sqrt(sigma2_sq) + C2);
   ssim_map = (lmap.^1).*(cmap.^1).*(smap.^1);   % You may try changing the exponent of each cost function.
   
else
   numerator1 = 2*mu1_mu2 + C1;
   numerator2 = 2*sigma12 + C2;

   denominator1 = mu1_sq + mu2_sq + C1;
   denominator2 = sigma1_sq + sigma2_sq + C2;
   
   ssim_map = ones(size(mu1));
   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
end

%temp = zeros(size(img1));
%temp(ind) = 1;
%iind = find(temp ==0);
%ssim_map(mask==0)=0;
mssim = mean(ssim_map(mask>0));

return

function [gaussKernel] = gkernel(sigma,sk)

% Pierrick Coupe - pierrick.coupe@gmail.com                                                                         
% Brain Imaging Center, Montreal Neurological Institute.                     
% Mc Gill University                                                         
%                                                                            
% Copyright (C) 2008 Pierrick Coupe             

for x = 1:(2*sk(1)+1)
    for y=1:(2*sk(2)+1)
        for z=1:(2*sk(3)+1)
            radiusSquared = (x-(sk(1)+1))^2 + (y-(sk(2)+1))^2 + (z-(sk(3)+1))^2;
            gaussKernel(x, y, z) = exp(-radiusSquared/(2*sigma^2));
        end
    end
end

